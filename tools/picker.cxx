/*ckwg +29
 * Copyright 2012 by Kitware, Inc.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 *  * Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 *  * Neither name of Kitware, Inc. nor the names of any contributors may be used
 *    to endorse or promote products derived from this software without specific
 *    prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS ``AS IS''
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHORS OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include "picker.h"

#include <vtkCellPicker.h>
#include <vtkIdTypeArray.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRendererCollection.h>
#include <vtkAxis.h>
#include <vtkMath.h>
#include <vtkTable.h>
#include <vtkPlot.h>
#include <vtkPlotPoints.h>
#include <vtkFloatArray.h>
#include <vtkPen.h>
#include <vtkTextProperty.h>
#include <vtkObjectFactory.h>

#include <super3d/depth/cost_volume.h>
#include <super3d/depth/depth_map.h>

//*****************************************************************************

Picker::Picker()
{
  selectedMapper = vtkSmartPointer<vtkDataSetMapper>::New();
  selectedActor = vtkSmartPointer<vtkActor>::New();
  num_clicks = 0;
  costvol = NULL;
}

//*****************************************************************************

Picker::~Picker()
{
  if (costvol)
    fclose(costvol);
}

//*****************************************************************************

void Picker::SetUpChart(vtkRenderer *chart_renderer,
                        const char *costvol_filename,
                        const char *depth_filename)
{
  costvol = fopen(costvol_filename, "rb");
  fread(&dims[0], sizeof(unsigned int), 1, costvol);
  fread(&dims[1], sizeof(unsigned int), 1, costvol);
  fread(&dims[2], sizeof(unsigned int), 1, costvol);

  vcl_cout << dims[0] << " " << dims[1] << " " << dims[2] << "\n";
  // Now the chart
  chart = vtkSmartPointer<vtkChartXY>::New();
  chartScene = vtkSmartPointer<vtkContextScene>::New();
  chartActor = vtkSmartPointer<vtkContextActor>::New();

  chart->SetAutoSize(true);
  chart->SetSize(vtkRectf(0.0, 0.0, 1024, 200));
  chart->GetAxis(0)->SetTitle("COST");
  chart->GetAxis(1)->SetTitle("NORMALIZED DEPTH");

  vtkSmartPointer<vtkPlotPoints> plot = vtkSmartPointer<vtkPlotPoints>::New();
  chart->AddPlot(plot);

  for (unsigned int i = 0; i <= 1; i++)
  {
    chart->GetAxis(i)->GetPen()->SetColor(255, 0, 0, 255);
    chart->GetAxis(i)->GetGridPen()->SetColor(50, 0, 0, 255);
    chart->GetAxis(i)->GetTitleProperties()->SetColor(1.0, 0.0, 0.0);
    chart->GetAxis(i)->GetLabelProperties()->SetColor(1.0, 0.0, 0.0);
  }

  chartScene->AddItem(chart);
  chartActor->SetScene(chartScene);

  //both needed
  chart_renderer->AddActor(chartActor.GetPointer());
  chartScene->SetRenderer(chart_renderer);

  super3d::load_depth(depth, depth_filename);
}

//*****************************************************************************

void Picker::UpdateChart(unsigned int picked)
{
  int numplots = chart->GetNumberOfPlots();
  for (int i = numplots-1; i >= 0; i--)
    chart->RemovePlot(i);

  unsigned int i = picked % dims[0];
  unsigned int j = picked / dims[0];

  vcl_cout << i << " " << j << "\n";

  vnl_vector<double> values;
  super3d::read_cost_volume_at(costvol, dims, i, j, values);

  vtkSmartPointer<vtkTable> table = vtkSmartPointer<vtkTable>::New();
  vtkSmartPointer<vtkFloatArray> arrD = vtkSmartPointer<vtkFloatArray>::New();
  arrD->SetName("Depth Label");
  table->AddColumn(arrD);

  vtkSmartPointer<vtkFloatArray> arrC = vtkSmartPointer<vtkFloatArray>::New();
  arrC->SetName("Cost");
  table->AddColumn(arrC);

  table->SetNumberOfRows(values.size());
  for (unsigned int i = 0; i < values.size(); i++)
  {
    table->SetValue(i, 0, ((double)i+0.5)/(double)values.size());
    table->SetValue(i, 1, values(i));
  }

  vtkSmartPointer<vtkTable> observed = vtkSmartPointer<vtkTable>::New();
  vtkSmartPointer<vtkFloatArray> d = vtkSmartPointer<vtkFloatArray>::New();
  d->SetName("Depth Label");
  observed->AddColumn(d);

  vtkSmartPointer<vtkFloatArray> c = vtkSmartPointer<vtkFloatArray>::New();
  c->SetName("Cost");
  observed->AddColumn(c);
  observed->SetNumberOfRows(2);

  observed->SetValue(0, 0, depth(i,j));
  observed->SetValue(0, 1, 0);
  observed->SetValue(1, 0, depth(i,j));
  observed->SetValue(1, 1, values.max_value());

  vtkPlot *plot = chart->AddPlot(vtkChart::LINE);
  plot->SetInputData(observed, 0, 1);
  plot->SetColor(255, 0, 0, 255);
  plot->SetWidth(2.0);
  plot->Update();

  // Add multiple line plots, setting the colors etc
  plot = chart->AddPlot(vtkChart::LINE);
  plot->SetInputData(table, 0, 1);
  plot->SetColor(255, 255, 0, 255);
  plot->SetWidth(2.0);
  plot->Update();
}

//*****************************************************************************

void Picker::OnLeftButtonDown()
{
  int* pos = this->GetInteractor()->GetEventPosition();

  if (num_clicks == 1)
  {
    int xdist = pos[0] - this->last_pos[0];
    int ydist = pos[1] - this->last_pos[1];
    int moveDistance = xdist*xdist + ydist*ydist;
    if(moveDistance > 25)
      num_clicks = 0;
  }

  if (num_clicks == 0)
  {
    last_pos[0] = pos[0];
    last_pos[1] = pos[1];
    num_clicks = 1;
    vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
    return;
  }

  num_clicks = 0;

  vtkSmartPointer<vtkCellPicker> picker = vtkSmartPointer<vtkCellPicker>::New();
  picker->SetTolerance(0.0005);

  // Pick from this location.
  picker->Pick(pos[0], pos[1], 0, this->GetDefaultRenderer());

  double* worldPosition = picker->GetPickPosition();
  vnl_double_3 pick_loc;
  pick_loc.set(worldPosition);

  if (picker->GetCellId() != -1)
  {
    vtkSmartPointer<vtkIdTypeArray> ids = vtkSmartPointer<vtkIdTypeArray>::New();
    ids->SetNumberOfComponents(1);
    ids->InsertNextValue(picker->GetPointId());

    UpdateChart(picker->GetPointId());
  }
  // Forward events
  vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
}

//*****************************************************************************

void Picker::PlaceArrow(double *picked)
{
  // Generate a random start and end point
  double startPoint[3], endPoint[3];
  endPoint[0] = picked[0];
  endPoint[1] = picked[1];
  endPoint[2] = picked[2];
  startPoint[0] = endPoint[0];
  startPoint[1] = endPoint[1];
  startPoint[2] = endPoint[2] + 0.002;

  // Compute a basis
  double normalizedX[3];
  double normalizedY[3];
  double normalizedZ[3];

  // The X axis is a vector from start to end
  vtkMath::Subtract(endPoint, startPoint, normalizedX);
  double length = vtkMath::Norm(normalizedX);
  vtkMath::Normalize(normalizedX);

  // The Z axis is an arbitrary vector cross X
  double arbitrary[3];
  arbitrary[0] = vtkMath::Random(-10,10);
  arbitrary[1] = vtkMath::Random(-10,10);
  arbitrary[2] = vtkMath::Random(-10,10);
  vtkMath::Cross(normalizedX, arbitrary, normalizedZ);
  vtkMath::Normalize(normalizedZ);

  // The Y axis is Z cross X
  vtkMath::Cross(normalizedZ, normalizedX, normalizedY);
  vtkSmartPointer<vtkMatrix4x4> matrix = vtkSmartPointer<vtkMatrix4x4>::New();

  // Create the direction cosine matrix
  matrix->Identity();
  for (unsigned int i = 0; i < 3; i++)
  {
    matrix->SetElement(i, 0, normalizedX[i]);
    matrix->SetElement(i, 1, normalizedY[i]);
    matrix->SetElement(i, 2, normalizedZ[i]);
  }

  // Apply the transforms
  vtkSmartPointer<vtkTransform> xform = vtkSmartPointer<vtkTransform>::New();
  xform->Translate(startPoint);
  xform->Concatenate(matrix);
  xform->Scale(length, length, length);

  vtkSmartPointer<vtkArrowSource> arrow = vtkSmartPointer<vtkArrowSource>::New();
  arrow->SetShaftRadius(.1);
  arrow->SetTipLength(.2);
  arrow->SetTipRadius(.2);
  arrow->Update();

  selectedMapper->SetInputConnection(arrow->GetOutputPort());

  selectedActor->SetUserMatrix(xform->GetMatrix());
  selectedActor->SetMapper(selectedMapper);
  selectedActor->GetProperty()->SetColor(1,0,0);
  this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->AddActor(selectedActor);
}

//*****************************************************************************

vtkStandardNewMacro(Picker);

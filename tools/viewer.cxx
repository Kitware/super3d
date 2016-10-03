/*ckwg +29
 * Copyright 2012-2016 by Kitware, Inc.
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

#include <vtkPolyDataMapper.h>
#include <vtkRenderWindow.h>
#include <vtkXMLPolyDataReader.h>

#include <vul/vul_arg.h>
#include <string>

#include "picker.h"

int main (int argc, char *argv[])
{
  //Required
  vul_arg<std::string> vtp_file(0, "VTP file", std::string("world_vol.vtp"));
  vul_arg<std::string> cost_vol_file( 0, "Cost volume file", std::string("cost_volume.dat"));
  vul_arg<std::string> depth_map_file( 0, "depth image save file", std::string("depth_map.dat"));

  vul_arg_parse( argc, argv );

  vtkSmartPointer<vtkXMLPolyDataReader> polyreader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
  polyreader->SetFileName(vtp_file().c_str());
  polyreader->Update();
  vtkPolyData *polydata = polyreader->GetOutput();

  vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInputData(polydata);

  vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);

  vtkSmartPointer<vtkRenderer> data_renderer = vtkSmartPointer<vtkRenderer>::New();
  vtkSmartPointer<vtkRenderer> chart_renderer = vtkSmartPointer<vtkRenderer>::New();
  vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();

  chart_renderer->SetViewport(0.0, 0.0, 1.0, 0.25);
  data_renderer->SetViewport(0.0, 0.25, 1.0, 1.0);
  renderWindow->AddRenderer(data_renderer);
  renderWindow->AddRenderer(chart_renderer);
  renderWindow->SetSize(1024, 768);
  renderWindow->SetWindowName("Cost Viewer");

  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow(renderWindow);
  renderWindowInteractor->Initialize();

  // Set the custom stype to use for interaction.
  vtkSmartPointer<Picker> picker = vtkSmartPointer<Picker>::New();
  picker->SetDefaultRenderer(data_renderer);
  picker->Data = polydata;
  picker->SetUpChart(chart_renderer, cost_vol_file().c_str(), depth_map_file().c_str());

  renderWindowInteractor->SetInteractorStyle(picker);
  data_renderer->AddActor(actor);
  data_renderer->ResetCamera();
  data_renderer->SetBackground(0, 0, 0);
  chart_renderer->SetBackground(0.0, 0.0, 0.0);

  renderWindow->Render();
  renderWindowInteractor->Start();

  return EXIT_SUCCESS;
}

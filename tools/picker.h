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

#ifndef PICKER_H_
#define PICKER_H_

#include <vtkSmartPointer.h>
#include <vtkInteractorStyleTrackballCamera.h>

#include <vtkDataSetMapper.h>
#include <vtkArrowSource.h>
#include <vtkMatrix4x4.h>
#include <vtkTransform.h>
#include <vtkContextScene.h>
#include <vtkContextActor.h>
#include <vtkChartXY.h>
#include <vtkRenderer.h>

#include <vtkRenderWindowInteractor.h>

#include <stdio.h>
#include <vil/vil_image_view.h>


#include <vnl/vnl_double_3.h>

// Catch mouse events
class Picker : public vtkInteractorStyleTrackballCamera
{
  public:
  static Picker* New();

  Picker();
  ~Picker();

  void SetUpChart(vtkRenderer *chart_renderer,
                  const char *costvol_filename,
                  const char *depth_filename);

  void UpdateChart(unsigned int picked);

  //Pick only on double click
  virtual void OnLeftButtonDown();

  void PlaceArrow(double *picked);

  vtkSmartPointer<vtkChartXY> chart;
  vtkSmartPointer<vtkContextScene> chartScene;
  vtkSmartPointer<vtkContextActor> chartActor;
  vtkSmartPointer<vtkPolyData> Data;
  vtkSmartPointer<vtkDataSetMapper> selectedMapper;
  vtkSmartPointer<vtkActor> selectedActor;

  FILE *costvol;
  unsigned int dims[3];
  int num_clicks;
  int last_pos[2];
  vil_image_view<double> depth;
};

#endif

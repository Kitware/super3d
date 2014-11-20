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

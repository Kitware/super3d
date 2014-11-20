#include <vtkPolyDataMapper.h>
#include <vtkRenderWindow.h>
#include <vtkXMLPolyDataReader.h>

#include <vul/vul_arg.h>
#include <vcl_string.h>

#include "picker.h"

int main (int argc, char *argv[])
{
  //Required
  vul_arg<vcl_string> vtp_file(0, "VTP file", vcl_string("world_vol.vtp"));
  vul_arg<vcl_string> cost_vol_file( 0, "Cost volume file", vcl_string("cost_volume.dat"));
  vul_arg<vcl_string> depth_map_file( 0, "depth image save file", vcl_string("depth_map.dat"));

  vul_arg_parse( argc, argv );

  vtkSmartPointer<vtkXMLPolyDataReader> polyreader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
  polyreader->SetFileName(vtp_file().c_str());
  polyreader->Update();
  vtkPolyData *polydata = polyreader->GetOutput();

  vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInput(polydata);

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

/*
 * Copyright 2011 Kitware, Inc.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of this project nor the names of its contributors
 *       may be used to endorse or promote products derived from this software
 *       without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#include <vcl_iostream.h>
#include <vcl_fstream.h>
#include <vcl_string.h>
#include <vcl_iomanip.h>

#ifdef HAVE_VTK
#include <vtkXMLPolyDataWriter.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkFloatArray.h>
#include <vtkIntArray.h>
#include <vtkUnsignedCharArray.h>
#include <vtkSmartPointer.h>
#include <vtkPointData.h>
#include <vtkDelaunay2D.h>
#include <vtkCellArray.h>
#endif

#include <vpgl/vpgl_perspective_camera.h>
#include <vgl/vgl_vector_3d.h>
#include <vgl/algo/vgl_rotation_3d.h>
#include <vnl/vnl_vector_fixed.h>
#include <vnl/vnl_matrix_fixed.h>
#include <vil/vil_load.h>
#include <vil/vil_image_view.h>
#include <vcl_sstream.h>

#include <vpgl/algo/vpgl_bundle_adjust.h>




int main(int argc, char *argv[])
{
  vcl_string dir(argv[1]);
  vcl_string bundlerfilename = dir + "/bundle/bundle.out";
  vcl_string listfilename = dir + "/list.txt";

  vcl_ifstream bout(bundlerfilename.c_str());
  vcl_ifstream list(listfilename.c_str());
  vcl_ofstream outfile(argv[2]);

  vcl_string line;
  vcl_getline(bout, line);

  unsigned int cams, pts;
  bout >> cams >> pts;

  vcl_ofstream vrml("cameras.wrl");

#ifdef HAVE_VTK
  vtkSmartPointer<vtkPoints> vtkpts = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> verts = vtkSmartPointer<vtkCellArray>::New();
  vtkIdType vert[1];
#endif

  vcl_ofstream framelist("framelist.txt");

  vcl_vector<vpgl_perspective_camera<double> > cameras;
  for (unsigned int i = 0; i < cams; i++)
  {
    double fl, k1, k2;
    bout >> fl >> k1 >> k2;

    if (fabs(fl) < 1e-9)
      continue;

    vcl_getline(list, line);
    vcl_istringstream lss(line);
    lss >> line;
    line.erase(line.begin());

    line = dir + line;
    vcl_cout << "loading image: " << line << "\n";
    vil_image_view<vxl_byte> img = vil_load(line.c_str());

    vpgl_calibration_matrix<double> K;
    K.set_focal_length(fl);
    //K.set_focal_length(1000);
    K.set_principal_point(vgl_point_2d<double>(img.ni()/2.0, img.nj()/2.0));


  //  K.set_principal_point(vgl_point_2d<double>(0.0, 0.0));
    vnl_matrix_fixed<double, 3, 3> R;
    vnl_vector_fixed<double, 3> R1, R2, R3;
    vnl_vector_fixed<double, 3> t;
    bout >> R1 >> R2 >> R3 >> t;
    R.set_row(0, R1);
    R.set_row(1, R2);
    R.set_row(2, R3);

    vnl_vector_fixed<double, 3> camcenter = -(R.transpose() * t);

    vnl_matrix_fixed<double, 3, 3> Rx;
    Rx.set_identity();
    Rx(1,1) = -1;
    Rx(2,2) = -1;

#ifdef HAVE_VTK
    vtkpts->InsertNextPoint(camcenter.data_block());
    vert[0] = i;
    verts->InsertNextCell(1, vert);
#endif

    vcl_cout << camcenter << "\n";
    vpgl_perspective_camera<double> cam(K, vgl_point_3d<double>(camcenter.data_block()), vgl_rotation_3d<double>(Rx * R));
    vcl_cout << "Cam Center: " << cam.get_camera_center() << "\n";

    vcl_cout << "Principle axis bundle: " << R.transpose() * vnl_vector_fixed<double, 3>(0, 0, -1)  << "\n";

    outfile << i << "\n" << vcl_setprecision(11) << cam << "\n";
    framelist << line << "\n";
    cameras.push_back(cam);
    vrml_write(vrml, cam, 0.01);
  }

  vcl_vector<vgl_point_3d<double> > worldpts;
  vpgl_bundle_adjust::write_vrml(vcl_string("camerasnew.wrl"), cameras, worldpts);

#ifdef HAVE_VTK
  vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
  polydata->SetPoints(vtkpts);
  polydata->SetVerts(verts);

  vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  writer->SetInputData(polydata);
  writer->SetFileName("cameras.vtp");
  writer->Update();
#endif

  vrml.close();
  outfile.close();
  list.close();
  bout.close();

  return 0;
}

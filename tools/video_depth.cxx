/*ckwg +29
* Copyright 2016 by Kitware, Inc.
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

#include <super3d/depth/super_config.h>

#include <super3d/depth/tv_refine_search.h>
#include <super3d/depth/cost_volume.h>

#include <super3d/depth/file_io.h>
#include <super3d/depth/depth_map.h>
#include <super3d/depth/multiscale.h>
#include <super3d/depth/tv_refine_plane.h>
#include <super3d/depth/world_rectilinear.h>
#include <super3d/depth/world_frustum.h>
#include <super3d/depth/exposure.h>


// VXL includes
#include <vil/vil_convert.h>
#include <vil/vil_save.h>
#include <vil/vil_crop.h>
#include <vil/vil_load.h>
#include <vil/vil_copy.h>
#include <vil/vil_decimate.h>
#include <vul/vul_timer.h>
#include <super3d/imesh/imesh_mesh.h>
#include <super3d/imesh/imesh_fileio.h>

#include <vpgl/vpgl_perspective_camera.h>

#include <vtkXMLImageDataWriter.h>
#include <vtkImageData.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>

#include <sstream>

int main(int argc, char* argv[])
{
  try
  {
    std::unique_ptr<super3d::config> cfg(new super3d::config);
    cfg->read_config(argv[1]);

    std::string frame_file = cfg->get_value<std::string>("frame_list");
    std::string dir("");
    if (cfg->is_set("directory"))
      dir = cfg->get_value<std::string>("directory");
    std::string camera_dir = cfg->get_value<std::string>("camera_dir");
    std::cout << "Using frame file: " << frame_file << " to find images and " << camera_dir << " to find cameras.\n";

    std::ifstream infile(frame_file.c_str());
    std::vector<std::string> filenames;
    std::string x;
    while (infile >> x) filenames.push_back(x);

    int numsupport = cfg->get_value<int>("support_frames");
    int stride = cfg->get_value<int>("stride");

    std::cout << "Read " << filenames.size() << " filenames.\n";

    int halfsupport = numsupport / 2;
    for (int i = halfsupport; i < static_cast<int>(filenames.size()) - halfsupport; i += stride)
    {
      std::cout << "Computing depth map on frame: " << i << "\n";
      std::vector<std::string> support_frames(filenames.begin() + (i - halfsupport), filenames.begin() + (i + halfsupport));
      std::cout << support_frames.size() << std::endl;

      std::vector<vil_image_view<double> > frames;
      std::vector<vpgl_perspective_camera<double> >  cameras;

      super3d::load_frames(support_frames, frames, cfg->get_value<bool>("use_color"), cfg->get_value<bool>("use_rgb12"));
      for (unsigned int f = 0; f < support_frames.size(); f++)
      {
        std::string camname = support_frames[f];
        unsigned int found = camname.find_last_of("/\\");
        camname = camname.substr(found + 1, camname.size() - 4 - found - 1);
        camname = cfg->get_value<std::string>("camera_dir") + "/" + camname + ".krtd";
        cameras.push_back(super3d::load_cam(camname));
      }

      unsigned int ref_frame = halfsupport;
      vpgl_perspective_camera<double> ref_cam = cameras[ref_frame];

      super3d::world_space *ws = NULL;
      int ni = frames[ref_frame].ni(), nj = frames[ref_frame].nj();
      double depth_min, depth_max;


      std::cout << "Computing depth range from " << cfg->get_value<std::string>("landmarks_path") << "\n";
      std::vector<vnl_double_3> landmarks;
      super3d::read_landmark_file(cfg->get_value<std::string>("landmarks_path"), landmarks);
      super3d::compute_depth_range(cameras[ref_frame], 0, ni, 0, nj, landmarks, depth_min, depth_max);
      std::cout << "Max estimated depth: " << depth_max << "\n";
      std::cout << "Min estimated depth: " << depth_min << "\n";

      ws = new super3d::world_frustum(cameras[ref_frame], depth_min, depth_max, ni, nj);

      std::cout << "Refining depth" << std::endl;
      unsigned int S = cfg->get_value<unsigned int>("num_slices");
      double theta0 = cfg->get_value<double>("theta_start");
      double theta_end = cfg->get_value<double>("theta_end");
      //double beta = cfg->get_value<double>("beta");
      double lambda = cfg->get_value<double>("lambda");
      double gw_alpha = cfg->get_value<double>("gw_alpha");
      double epsilon = cfg->get_value<double>("epsilon");

      vil_image_view<double> g;
      vil_image_view<double> cost_volume;

      double iw = cfg->get_value<double>("intensity_cost_weight");
      double gw = cfg->get_value<double>("gradient_cost_weight");
      double cw = cfg->get_value<double>("census_cost_weight");
      super3d::compute_world_cost_volume(frames, cameras, ws, ref_frame, S, cost_volume, iw, gw, cw);
      //compute_cost_volume_warp(frames, cameras, ref_frame, S, depth_min, depth_max, cost_volume);
      ws->compute_g(frames[ref_frame], g, gw_alpha, 1.0);


      std::cout << "Refining Depth. ..\n";
      vil_image_view<double> depth(cost_volume.ni(), cost_volume.nj(), 1);

      vul_timer timer;

      unsigned int iterations = 2000;
      if (cfg->is_set("iterations"))
        iterations = cfg->get_value<unsigned int>("iterations");
      super3d::refine_depth(cost_volume, g, depth, iterations, theta0, theta_end, lambda, epsilon);

      double sec = 1e-3 * timer.real();
      std::cout << "super3d took " << sec << " seconds.\n";

      std::string outdir = cfg->get_value<std::string>("outdir");
      std::ostringstream depth_name;
      depth_name << outdir << "/" << i;

      vil_image_view<vxl_byte> dmap;
      vil_convert_stretch_range_limited(depth, dmap, 0.0, 1.0);
      // depth map are drawn inverted (white == closest) for viewing
      vil_math_scale_and_offset_values(dmap, -1.0, 255);
      std::string depthmap_file = depth_name.str() + ".png";
      vil_save(dmap, depthmap_file.c_str());

      super3d::save_depth_to_vtp((depth_name.str() + ".vtp").c_str(), depth, frames[ref_frame], ref_cam, ws);

      // map depth from normalized range back into true depth
      double depth_scale = depth_max - depth_min;
      vil_math_scale_and_offset_values(depth, depth_scale, depth_min);

      double minv, maxv;
      vil_math_value_range(depth, minv, maxv);
      std::cout << "Depth range: " << minv << " - " << maxv << "\n";

      vtkNew<vtkDoubleArray> uniquenessRatios;
      uniquenessRatios->SetName("Uniqueness Ratios");
      uniquenessRatios->SetNumberOfValues(ni*nj);

      vtkNew<vtkDoubleArray> bestCost;
      bestCost->SetName("Best Cost Values");
      bestCost->SetNumberOfValues(ni*nj);

      vtkNew<vtkUnsignedCharArray> color;
      color->SetName("Color");
      color->SetNumberOfComponents(3);
      color->SetNumberOfTuples(ni*nj);

      vtkNew<vtkDoubleArray> depths;
      depths->SetName("Depths");
      depths->SetNumberOfComponents(1);
      depths->SetNumberOfTuples(ni*nj);

      vil_image_view<vxl_byte> ref_img_color = vil_load(support_frames[ref_frame].c_str());

      vtkIdType pt_id = 0;

      for (int y = nj - 1; y >= 0; y--)
      {
        for (int x = 0; x < ni; x++)
        {
          uniquenessRatios->SetValue(pt_id, 0);
          bestCost->SetValue(pt_id, 0);
          depths->SetValue(pt_id, depth(x, y));
          color->SetTuple3(pt_id, (int)ref_img_color(x, y, 0), (int)ref_img_color(x, y, 1), (int)ref_img_color(x, y, 2));
          pt_id++;
        }
      }

      vtkNew<vtkImageData> imageData;
      imageData->SetSpacing(1, 1, 1);
      imageData->SetOrigin(0, 0, 0);
      imageData->SetDimensions(ni, nj, 1);
      imageData->GetPointData()->AddArray(depths.Get());
      imageData->GetPointData()->AddArray(color.Get());
      imageData->GetPointData()->AddArray(uniquenessRatios.Get());
      imageData->GetPointData()->AddArray(bestCost.Get());

      vtkNew<vtkXMLImageDataWriter> writerI;
      std::string depthmapImageFileName = depth_name.str() + ".vti";

      writerI->SetFileName(depthmapImageFileName.c_str());
      writerI->AddInputDataObject(imageData.Get());
      writerI->SetDataModeToBinary();
      writerI->Write();
      std::cout << "Saved : " << depthmapImageFileName << std::endl;


      if (ws) delete ws;
    }

  }
  catch (const super3d::config::cfg_exception &e)
  {
    std::cout << "Error in config: " << e.what() << "\n";
  }

  return 0;
}

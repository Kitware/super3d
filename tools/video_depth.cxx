/*ckwg +29
 *
* Copyright 2017 by Kitware, Inc.
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
#include <super3d/depth/world_angled_frustum.h>
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
#include <fstream>
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

    std::ifstream frame_infile(frame_file.c_str());
    std::vector<std::string> filenames;
    std::string x;
    while (frame_infile >> x) filenames.push_back(x);
    frame_infile.close();

    int numsupport = cfg->get_value<int>("support_frames");
    int stride = cfg->get_value<int>("stride");

    std::cout << "Read " << filenames.size() << " filenames.\n";

    std::vector<std::string> masknames;
    if (cfg->is_set("mask_list"))
    {
      std::string mask_list = cfg->get_value<std::string>("mask_list");
      std::ifstream mask_file_stream(mask_list.c_str());
      while (mask_file_stream >> x) masknames.push_back(x);
      mask_file_stream.close();
    }

    std::string outdir = cfg->get_value<std::string>("outdir");
    std::ofstream vtiList(outdir + "/vtiList.txt");
    std::ofstream kList(outdir + "/kList.txt");
    if( !vtiList || !kList )
    {
      std::cerr << "unable to open vtiList.txt or kList.txt for writing" << std::endl;
      return EXIT_FAILURE;
    }

    bool use_world_planes = false;
    vnl_double_3 normal;
    if (cfg->is_set("world_plane_normal"))
    {
      // use world coordinate slices in this direction instead of depth
      std::istringstream ss(cfg->get_value<std::string>("world_plane_normal"));
      ss >> normal;
      normal.normalize();
      use_world_planes = true;
    }

    int halfsupport = numsupport / 2;
    for (int i = halfsupport; i < static_cast<int>(filenames.size()) - halfsupport; i += stride)
    {
      std::cout << "Computing depth map on frame: " << i << "\n";
      std::vector<std::string> support_frames(filenames.begin() + (i - halfsupport), filenames.begin() + (i + halfsupport));

      std::vector<std::string> support_masks;
      if (masknames.size() == 1)
      {
        support_masks = std::vector<std::string>(support_frames.size(), masknames[0]);
      }
      else if (!masknames.empty())
      {
        support_masks = std::vector<std::string>(masknames.begin() + (i - halfsupport), masknames.begin() + (i + halfsupport));
      }
      std::cout << support_frames.size() << std::endl;

      std::vector<vil_image_view<double> > frames;
      std::vector<vil_image_view<double> > masks;
      std::vector<vpgl_perspective_camera<double> >  cameras;

      super3d::load_frames(support_frames, frames, cfg->get_value<bool>("use_color"), cfg->get_value<bool>("use_rgb12"));
      std::string ref_cam_name;
      for (unsigned int f = 0; f < support_frames.size(); f++)
      {
        std::string camname = support_frames[f];
        unsigned int found = camname.find_last_of("/\\");
        camname = camname.substr(found + 1, camname.size() - 4 - found - 1);
        if (f==halfsupport)
        {
          ref_cam_name = camname;
        }
        camname = cfg->get_value<std::string>("camera_dir") + "/" + camname + ".krtd";
        cameras.push_back(super3d::load_cam(camname));

        vil_image_view<double> mask;
        if (support_masks.empty())
        {
          mask = vil_image_view<double>(frames[f].ni(), frames[f].nj(), 1);
          mask.fill(1.0);
          vil_crop(mask, 2, mask.ni()-4, 2, mask.nj()-4).fill(0.0);
        }
        else
        {
          std::cout << support_masks[f] << "\n";
          vil_image_view<bool> maskb = vil_load(support_masks[f].c_str());
          vil_convert_cast(maskb, mask);
        }
        masks.push_back(mask);
      }

      unsigned int ref_frame = halfsupport;
      vpgl_perspective_camera<double> ref_cam = cameras[ref_frame];

      super3d::world_space *ws = NULL;
      int ni = frames[ref_frame].ni(), nj = frames[ref_frame].nj();
      double depth_min, depth_max;


      std::cout << "Computing depth range from " << cfg->get_value<std::string>("landmarks_path") << "\n";
      std::vector<vnl_double_3> landmarks;
      super3d::read_landmark_file(cfg->get_value<std::string>("landmarks_path"), landmarks);
      std::vector<vnl_double_3> visible_landmarks =
        super3d::filter_visible_landmarks(cameras[ref_frame], 0, ni, 0, nj, landmarks);
      if (use_world_planes)
      {
        super3d::compute_offset_range(visible_landmarks, normal, depth_min, depth_max, 0.0, 0.5);
        std::cout << "Max estimated offset: " << depth_max << "\n";
        std::cout << "Min estimated offset: " << depth_min << "\n";
        ws = new super3d::world_angled_frustum(cameras[ref_frame], normal, depth_min, depth_max, ni, nj);
      }
      else
      {
        super3d::compute_depth_range(visible_landmarks, cameras[ref_frame], depth_min, depth_max);
        std::cout << "Max estimated depth: " << depth_max << "\n";
        std::cout << "Min estimated depth: " << depth_min << "\n";
        ws = new super3d::world_frustum(cameras[ref_frame], depth_min, depth_max, ni, nj);
      }

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
      super3d::compute_world_cost_volume(frames, cameras, ws, ref_frame, S, cost_volume, iw, gw, cw,&masks);
      //compute_cost_volume_warp(frames, cameras, ref_frame, S, depth_min, depth_max, cost_volume);
      ws->compute_g(frames[ref_frame], g, gw_alpha, 1.0, &masks[ref_frame]);


      std::cout << "Refining Depth. ..\n";
      vil_image_view<double> depth(cost_volume.ni(), cost_volume.nj(), 1);

      vul_timer timer;

      unsigned int iterations = 2000;
      if (cfg->is_set("iterations"))
        iterations = cfg->get_value<unsigned int>("iterations");
      super3d::refine_depth(cost_volume, g, depth, iterations, theta0, theta_end, lambda, epsilon);

      double sec = 1e-3 * timer.real();
      std::cout << "super3d took " << sec << " seconds.\n";

      std::string depth_name = outdir + "/" + ref_cam_name;

      super3d::save_depth_to_vtp((depth_name + ".vtp").c_str(), depth, frames[ref_frame], ref_cam, ws);

      // map depth from normalized range back into true depth
      double depth_scale = depth_max - depth_min;
      vil_math_scale_and_offset_values(depth, depth_scale, depth_min);

      vil_image_view<double> height_map;
      if (use_world_planes)
      {
        height_map = depth;
        depth = vil_image_view<double>();
        super3d::height_map_to_depth_map(cameras[ref_frame], height_map, depth);
      }
      else
      {
        super3d::depth_map_to_height_map(cameras[ref_frame], depth, height_map);
      }

      // save byte depth map
      vil_image_view<vxl_byte> bmap;
      vil_convert_stretch_range(depth, bmap);
      // depth map are drawn inverted (white == closest) for viewing
      vil_math_scale_and_offset_values(bmap, -1.0, 255);
      std::string depthmap_file = depth_name + "_depth.png";
      vil_save(bmap, depthmap_file.c_str());
      // save byte height map
      vil_convert_stretch_range(height_map, bmap);
      std::string heightmap_file = depth_name + "_height.png";
      vil_save(bmap, heightmap_file.c_str());

      double minv, maxv;
      vil_math_value_range(depth, minv, maxv);
      std::cout << "Depth range: " << minv << " - " << maxv << "\n";
      vil_math_value_range(height_map, minv, maxv);
      std::cout << "Height range: " << minv << " - " << maxv << "\n";

      vil_image_view<vxl_byte> ref_img_color = vil_load(support_frames[ref_frame].c_str());
      super3d::save_depth_to_vti((depth_name + ".vti").c_str(), depth, ref_img_color);
      std::cout << "Saved : " << depth_name + ".vti" << std::endl;

      vtiList << i << " " << ref_cam_name + ".vti\n";
      kList << i << " ../krtd/" << ref_cam_name << ".krtd\n";


      if (ws) delete ws;
    }
    vtiList.close();
    kList.close();

  }
  catch (const super3d::config::cfg_exception &e)
  {
    std::cout << "Error in config: " << e.what() << "\n";
  }

  return 0;
}

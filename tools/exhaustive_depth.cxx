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
#include <super3d/imesh/imesh_mesh.h>
#include <super3d/imesh/imesh_fileio.h>

// VXL includes
#include <vul/vul_timer.h>
#include <vil/vil_convert.h>
#include <vil/vil_save.h>
#include <vil/vil_crop.h>
#include <vil/vil_load.h>
#include <vil/vil_copy.h>
#include <vil/vil_decimate.h>
#include <vpgl/vpgl_perspective_camera.h>

#ifdef HAVE_VISCL
#include <super3d/depth_cl/refine_depth.h>
#endif



int main(int argc, char* argv[])
{
  try
  {
    std::unique_ptr<super3d::config> cfg(new super3d::config);
    cfg->read_config(argv[1]);

#ifndef HAVE_VISCL
    if (cfg->is_set("use_gpu") && cfg->get_value<bool>("use_gpu"))
    {
      std::cerr << "use_gpu is true but not built with viscl.\n";
      return 1;
    }
#endif

  std::vector<vil_image_view<double> > frames;
  std::vector<vpgl_perspective_camera<double> >  cameras;
  std::vector<std::string> filenames;

  std::string frame_file = cfg->get_value<std::string>("frame_list");
  std::string dir("");
  if (cfg->is_set("directory"))
    dir = cfg->get_value<std::string>("directory");
  std::vector<int> frameindex;

  if (cfg->is_set("camera_file"))
  {
    std::string camera_file = cfg->get_value<std::string>("camera_file");
    std::cout << "Using frame file: " << frame_file << " to find images and " << camera_file  << " to find cameras.\n";
    super3d::load_from_frame_file(frame_file.c_str(), dir, filenames, frameindex, frames,
                         cfg->get_value<bool>("use_color"), cfg->get_value<bool>("use_rgb12"));
    super3d::load_cams(camera_file.c_str(), frameindex, cameras);
  }
  else if (cfg->is_set("camera_dir"))
  {
    std::string camera_dir = cfg->get_value<std::string>("camera_dir");
    std::cout << "Using frame file: " << frame_file << " to find images and " << camera_dir  << " to find cameras.\n";

    super3d::load_from_frame_file(frame_file.c_str(), dir, filenames, frameindex, frames,
                          cfg->get_value<bool>("use_color"), cfg->get_value<bool>("use_rgb12"));
    for (unsigned int i = 0; i < filenames.size(); i++)
    {
      std::string camname = filenames[i];
      unsigned int found = camname.find_last_of("/\\");
      camname = camname.substr(found+1, camname.size() - 4 - found - 1);
      camname = cfg->get_value<std::string>("camera_dir") + "/" + camname + ".krtd";
      cameras.push_back(super3d::load_cam(camname));
    }
  }
  else
  {
    std::cerr << "Error: must use camera_file or camera_dir.\n";
    return -1;
  }

  if (cfg->is_set("exposure_file"))
  {
    std::vector<std::pair<double, double> > exposures;
    exposures = super3d::load_exposure(cfg->get_value<std::string>("exposure_file"), frameindex);

    for (unsigned int i = 0; i < frames.size(); i++)
      super3d::apply_exposure_correction(frames[i], exposures[i]);
  }
  else
  {
    for (unsigned int i = 0; i < frames.size(); i++)
      vil_math_scale_and_offset_values(frames[i], 1.0/255.0, 0.0);
  }

  unsigned int ref_frame = cfg->get_value<unsigned int>("ref_frame");
  vpgl_perspective_camera<double> ref_cam = cameras[ref_frame];

  super3d::world_space *ws = NULL;
  int i0, ni, j0, nj;
  double depth_min, depth_max;

  if (cfg->is_set("world_volume"))
  {
    std::istringstream wvstream(cfg->get_value<std::string>("world_volume"));
    vnl_double_3 origin, dimensions;
    unsigned int ni, nj;
    for (unsigned int i = 0; i < 3; i++)  wvstream >> origin[i];
    for (unsigned int i = 0; i < 3; i++)  wvstream >> dimensions[i];
    wvstream >> ni >> nj;
    ws = new super3d::world_rectilinear(origin, dimensions, ni, nj);
  }
  else
  {
    double camera_scale = 1.0;
    if (cfg->is_set("camera_scale"))
    {
      camera_scale = cfg->get_value<double>("camera_scale");
      for (unsigned int i = 0; i < cameras.size(); i++)
        cameras[i] = super3d::scale_camera(cameras[i], camera_scale);
    }
    //Compute the window cropping, scale the cropping by the specified scale so that we do not
    //need to recompute cropping input for super resolution.
    if (cfg->is_set("crop_window"))
    {
      std::istringstream cwstream(cfg->get_value<std::string>("crop_window"));
      cwstream >> i0 >> ni >> j0 >> nj;
      i0 = (int)(i0*camera_scale);
      j0 = (int)(j0*camera_scale);
      ni = (int)(ni*camera_scale);
      nj = (int)(nj*camera_scale);
      std::cout << "Crop window: " << i0 << " " << ni << " " << j0 << " " << nj << "\n";
    }
    else
    {
      i0 = j0 = 0;
      ni = frames[ref_frame].ni();
      nj = frames[ref_frame].nj();
    }

    frames[ref_frame] = vil_crop(frames[ref_frame], i0, ni, j0, nj);
    cameras[ref_frame] = super3d::crop_camera(cameras[ref_frame], i0, j0);

    if (cfg->is_set("landmarks_path"))
    {
      std::cout << "Computing depth range from " << cfg->get_value<std::string>("landmarks_path") << "\n";
      super3d::compute_depth_range(cameras[ref_frame], 0, ni, 0, nj, cfg->get_value<std::string>("landmarks_path"), depth_min, depth_max);
      std::cout << "Max estimated depth: " << depth_max << "\n";
      std::cout << "Min estimated depth: " << depth_min << "\n";
    }
    else
    {
      depth_min = cfg->get_value<double>("depth_min");
      depth_max = cfg->get_value<double>("depth_max");
    }



    ws = new super3d::world_frustum(cameras[ref_frame], depth_min, depth_max, ni, nj);
  }

  std::cout << "Refining depth"<<std::endl;
  unsigned int S = cfg->get_value<unsigned int>("num_slices");
  double theta0 = cfg->get_value<double>("theta_start");
  double theta_end = cfg->get_value<double>("theta_end");
#ifdef HAVE_VISCL
  double beta = cfg->get_value<double>("beta");
#endif
  double lambda = cfg->get_value<double>("lambda");
  double gw_alpha = cfg->get_value<double>("gw_alpha");
  double epsilon = cfg->get_value<double>("epsilon");

  vil_image_view<double> g;
  vil_image_view<double> cost_volume;

  if (!cfg->is_set("compute_cost_volume") ||
      cfg->get_value<bool>("compute_cost_volume") )
  {
    double iw = cfg->get_value<double>("intensity_cost_weight");
    double gw = cfg->get_value<double>("gradient_cost_weight");
    double cw = cfg->get_value<double>("census_cost_weight");
    super3d::compute_world_cost_volume(frames, cameras, ws, ref_frame, S, cost_volume, iw, gw, cw);
    //compute_cost_volume_warp(frames, cameras, ref_frame, S, depth_min, depth_max, cost_volume);
    ws->compute_g(frames[ref_frame], g, gw_alpha, 1.0);

    if(cfg->is_set("cost_volume_file"))
    {
      std::string cost_file = cfg->get_value<std::string>("cost_volume_file");
      super3d::save_cost_volume(cost_volume, g, cost_file.c_str());
    }
  }
  else if (cfg->is_set("cost_volume_file"))
  {
    std::string cost_file = cfg->get_value<std::string>("cost_volume_file");
    if (cfg->is_set("mix_cost_volumes") &&
        cfg->get_value<bool>("mix_cost_volumes"))
    {
      double iw = cfg->get_value<double>("intensity_cost_weight");
      double gw = cfg->get_value<double>("gradient_cost_weight");
      double cw = cfg->get_value<double>("census_cost_weight");
      int pos = cost_file.find_last_of(".");
      std::string ext = cost_file.substr(pos);
      std::string basename = cost_file.substr(0, pos);
      vil_image_view<double> tmp;
      super3d::load_cost_volume(cost_volume, g, (basename+"_intensity"+ext).c_str());
      super3d::load_cost_volume(tmp, g, (basename+"_gradient"+ext).c_str());
      vil_math_add_image_fraction(cost_volume, iw, tmp, gw);
      super3d::load_cost_volume(tmp, g, (basename+"_census"+ext).c_str());
      vil_math_add_image_fraction(cost_volume, 1.0, tmp, cw);
    }
    else
    {
      super3d::load_cost_volume(cost_volume, g, cost_file.c_str());
    }
  }
  else
  {
    std::cerr << "Error: must either set compute_cost_volume"
             << " or specify a cost_volume_file to load." << std::endl;
    return -1;
  }


  std::cout << "Refining Depth. ..\n";
  vil_image_view<double> depth(cost_volume.ni(), cost_volume.nj(), 1);

#ifndef USE_BP

  if (!cfg->is_set("compute_init_depthmap") ||
       cfg->get_value<bool>("compute_init_depthmap") )
  {
#ifdef HAVE_VISCL
    if (cfg->is_set("use_gpu") &&
        cfg->get_value<bool>("use_gpu") )
    {
      super3d::cl::refine_depth_cl_t rd = NEW_VISCL_TASK(super3d::cl::refine_depth_cl);
      vil_image_view<float> depth_f(cost_volume.ni(), cost_volume.nj(), 1);
      vil_image_view<float> cv_f(cost_volume.ni(), cost_volume.nj(), 1, cost_volume.nplanes());
      vil_image_view<float> g_f;
      vil_convert_cast<double, float>(g, g_f);

      {
        vil_image_view<float> cv_float;
        vil_convert_cast<double, float>(cost_volume, cv_float);
        vil_copy_reformat<float>(cv_float, cv_f);
      }

      vul_timer timer;
      rd->refine(cv_f, depth_f, g_f, beta, theta0, theta_end, lambda);
      double sec = 1e-3 * timer.real();
      std::cout << "viscl took " << sec << " seconds.\n";
      vil_convert_cast<float, double>(depth_f, depth);
    }
    else
#endif
    {
      vul_timer timer;
      unsigned int iterations = 2000;
      if (cfg->is_set("iterations"))
        iterations = cfg->get_value<unsigned int>("iterations");
      super3d::refine_depth(cost_volume, g, depth, iterations, theta0, theta_end, lambda, epsilon);

      double sec = 1e-3 * timer.real();
      std::cout << "super3d took " << sec << " seconds.\n";
    }
    if(cfg->is_set("init_depthmap_file"))
    {
      std::string depthmap_file = cfg->get_value<std::string>("init_depthmap_file");
      super3d::save_depth(depth, depthmap_file.c_str());
    }
  }
  else if (cfg->is_set("init_depthmap_file"))
  {
    std::string depthmap_file = cfg->get_value<std::string>("init_depthmap_file");
    super3d::load_depth(depth, depthmap_file.c_str());
  }
  else
  {
    std::cerr << "Error: must either set compute_init_depthmap"
             << " or specify an init_depthmap_file to load." << std::endl;
    return -1;
  }

  super3d::save_depth(depth, "depth_map_normals.dat");

#else
  bp_refine(cost_volume, depth_min, depth_max, depth);    //TODO: need to check if this works after idepth->depth change
#endif

  if (cfg->is_set("output_depthmap"))
  {
    vil_image_view<vxl_byte> dmap;
    vil_convert_stretch_range_limited(depth, dmap, 0.0, 1.0);
    // depth map are drawn inverted (white == closest) for viewing
    vil_math_scale_and_offset_values(dmap, -1.0, 255);
    std::string depthmap_file = cfg->get_value<std::string>("output_depthmap");
    vil_save(dmap, depthmap_file.c_str());
  }

  std::cout << "writing mesh"<<std::endl;
#ifdef HAVE_VTK
  //vtp depth writer uses 0..1 depth scaling
  vil_image_view<double> d_texture;

  if (cfg->get_value<bool>("use_rgb12"))
  {
    vil_image_resource_sptr img_rsc = vil_load_image_resource((dir + filenames[ref_frame]).c_str());
    vil_image_view<unsigned short> us_texture = img_rsc->get_view();
    vil_convert_cast(us_texture, d_texture);
    vil_math_scale_values(d_texture, 255.0 / 4095.0);
  }
  else
  {
    vil_image_view<vxl_byte> b_texture;
    b_texture = vil_load((dir + filenames[ref_frame]).c_str());
    vil_convert_cast<vxl_byte, double>(b_texture, d_texture);
  }

  std::string output_file_name = cfg->get_value<std::string>("output_file");
  save_depth_to_vtp(output_file_name.c_str(), depth, d_texture, ref_cam, ws);
#endif

  // map depth from normalized range back into true depth
  double depth_scale = depth_max - depth_min;
  vil_math_scale_and_offset_values(depth, depth_scale, depth_min);

  if (cfg->is_set("obj_file"))
  {
    double output_decimate = 1.0;
    imesh_mesh nm = super3d::depth_map_to_mesh(super3d::scale_camera(cameras[ref_frame], 1.0/output_decimate),
                                               vil_decimate(depth,output_decimate));
    imesh_write_obj(cfg->get_value<std::string>("obj_file"), nm);
  }

  if (cfg->is_set("output_float_depthmap"))
  {
    std::string depthmap_file = cfg->get_value<std::string>("output_float_depthmap");
    vil_save(depth, depthmap_file.c_str());
  }

  vil_image_view<double> gt;
  if (cfg->is_set("ground_truth"))
  {
    std::string gtfile = cfg->get_value<std::string>("ground_truth");
    gt = vil_load(gtfile.c_str());
    gt = vil_crop(gt, i0, ni, j0, nj);
    vil_image_view<double> diff;
    vil_math_image_difference(gt, depth, diff);
    vil_image_view<vxl_byte> diff_byte;
    vil_convert_stretch_range_limited(diff, diff_byte, -depth_scale, depth_scale);
    vil_save(diff_byte, "diff.png");

    double threshold = depth_scale / S;
    super3d::score_vs_gt(depth, gt, threshold);
    std::cout << "cost of estimated depth: " << super3d::eval_hessian_frob(depth, cost_volume, lambda) << "\n";
    std::cout << "cost of gt depth: " << super3d::eval_hessian_frob(gt, cost_volume, lambda) << "\n";
  }

  if (ws) delete ws;

  } catch (const super3d::config::cfg_exception &e)
  {
    std::cout << "Error in config: " << e.what() << "\n";
  }

  return 0;
}

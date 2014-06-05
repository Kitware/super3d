/*
 * Copyright 2012 Kitware, Inc.
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

#include "super_config.h"

#include "tv_refine_search.h"
#include "cost_volume.h"

#include "file_io.h"
#include "depth_map.h"
#include "multiscale.h"
#include "tv_refine_plane.h"
#include "world_rectilinear.h"
#include "world_frustum.h"

#include <boost/scoped_ptr.hpp>

#include "exposure.h"

// VXL includes
#include <vil/vil_convert.h>
#include <vil/vil_save.h>
#include <vil/vil_crop.h>
#include <vil/vil_load.h>
#include <vil/vil_copy.h>
#include <vil/vil_decimate.h>
#include <imesh/imesh_mesh.h>
#include <imesh/imesh_fileio.h>

#ifdef HAVE_VISCL
#include "depth_cl/refine_depth.h"
#endif

#include <boost/chrono.hpp>

#include <vpgl/vpgl_perspective_camera.h>


int main(int argc, char* argv[])
{
  try
  {
    boost::scoped_ptr<config> cfg(new config);
    cfg->read_config(argv[1]);

#ifndef HAVE_VISCL
    if (cfg->is_set("use_gpu") && cfg->get_value<bool>("use_gpu"))
    {
      vcl_cerr << "use_gpu is true but not built with viscl.\n";
      return 1;
    }
#endif

  vcl_vector<vil_image_view<double> > frames;
  vcl_vector<vpgl_perspective_camera<double> >  cameras;
  vcl_vector<vcl_string> filenames;
  vcl_string camera_file = cfg->get_value<vcl_string>("camera_file");

  if (cfg->is_set("frame_format"))
  {
    vcl_cout << "Using formated string to find images/cameras.\n";
    vul_sequence_filename_map frame_seq(cfg->get_value<vcl_string>("frame_format"));

    //Read Cameras
    cameras = load_cams(camera_file, frame_seq);

    //Read Images
    frames = load_frames(frame_seq, filenames);
    if (frames.empty())
    {
      vcl_cerr << "No frames found"<<vcl_endl;
      return -1;
    }
    else if (frames.size() < 2)
    {
      vcl_cerr << "At least 2 frames are required"<<vcl_endl;
      return -1;
    }

    //Apply exposure correction
    if (cfg->is_set("exposure_file"))
    {
      vcl_vector<vcl_pair<double, double> > exposures;
      exposures = load_exposure(cfg->get_value<vcl_string>("exposure_file"), frame_seq);

      for (unsigned int i = 0; i < frames.size(); i++)
        apply_exposure_correction(frames[i], exposures[i]);
    }
  }
  else if (cfg->is_set("frame_list"))
  {
    vcl_string frame_file = cfg->get_value<vcl_string>("frame_list");
    vcl_string dir = cfg->get_value<vcl_string>("directory");
    vcl_cout << "Using frame file: " << frame_file << " to find images and " << camera_file  << " to find cameras.\n";
    vcl_vector<int> frameindex;
    load_from_frame_file(frame_file.c_str(), camera_file.c_str(), dir, filenames, frameindex, frames, cameras);
    if (cfg->is_set("exposure_file"))
    {
      vcl_vector<vcl_pair<double, double> > exposures;
      exposures = load_exposure(cfg->get_value<vcl_string>("exposure_file"), frameindex);

      for (unsigned int i = 0; i < frames.size(); i++)
        apply_exposure_correction(frames[i], exposures[i]);
    }
  }
  else
  {
    vcl_cerr << "Error: must use either frame list (-fl) or frame format string (-ff).\n";
    return -1;
  }

  if (!cfg->is_set("exposure_file"))
  {
    for (unsigned int i = 0; i < frames.size(); i++)
      vil_math_scale_and_offset_values(frames[i], 1.0/255.0, 0.0);
  }

  unsigned int ref_frame = cfg->get_value<unsigned int>("ref_frame");
  vpgl_perspective_camera<double> ref_cam = cameras[ref_frame];

  world_space *ws = NULL;
  int i0, ni, j0, nj;
  double depth_min, depth_max;

  if (cfg->is_set("landmarks_path"))
  {
    vcl_cout << "Computing depth range from " << cfg->get_value<vcl_string>("landmarks_path") << "\n";
    compute_depth_range(cameras[ref_frame], cfg->get_value<vcl_string>("landmarks_path"), depth_min, depth_max);
    vcl_cout << "Max estimated depth: " << depth_max << "\n";
    vcl_cout << "Min estimated depth: " << depth_min << "\n";
  }
  else
  {
    depth_min = cfg->get_value<double>("depth_min");
    depth_max = cfg->get_value<double>("depth_max");

  }

  if (cfg->is_set("world_volume"))
  {
    vcl_istringstream wvstream(cfg->get_value<vcl_string>("world_volume"));
    vnl_double_3 origin, dimensions;
    unsigned int ni, nj;
    for (unsigned int i = 0; i < 3; i++)  wvstream >> origin[i];
    for (unsigned int i = 0; i < 3; i++)  wvstream >> dimensions[i];
    wvstream >> ni >> nj;
    ws = new world_rectilinear(origin, dimensions, ni, nj);
  }
  else
  {
    double camera_scale = 1.0;
    if (cfg->is_set("camera_scale"))
    {
      camera_scale = cfg->get_value<double>("camera_scale");
      for (unsigned int i = 0; i < cameras.size(); i++)
        cameras[i] = scale_camera(cameras[i], camera_scale);
    }
    //Compute the window cropping, scale the cropping by the specified scale so that we do not
    //need to recompute cropping input for super resolution.
    if (cfg->is_set("crop_window"))
    {
      vcl_istringstream cwstream(cfg->get_value<vcl_string>("crop_window"));
      cwstream >> i0 >> ni >> j0 >> nj;
      i0 = (int)(i0*camera_scale);
      j0 = (int)(j0*camera_scale);
      ni = (int)(ni*camera_scale);
      nj = (int)(nj*camera_scale);
      vcl_cout << "Crop window: " << i0 << " " << ni << " " << j0 << " " << nj << "\n";
      frames[ref_frame] = vil_crop(frames[ref_frame], i0, ni, j0, nj);
      cameras[ref_frame] = crop_camera(cameras[ref_frame], i0, j0);
    }
    else
    {
      i0 = j0 = 0;
      ni = frames[ref_frame].ni();
      nj = frames[ref_frame].nj();
    }

    ws = new world_frustum(cameras[ref_frame], depth_min, depth_max, ni, nj);
  }

  vcl_cout << "Refining depth"<<vcl_endl;
  unsigned int S = cfg->get_value<unsigned int>("num_slices");
  double theta0 = cfg->get_value<double>("theta_start");
  double theta_end = cfg->get_value<double>("theta_end");
  double beta = cfg->get_value<double>("beta");
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
    compute_world_cost_volume(frames, cameras, ws, ref_frame, S, cost_volume, iw, gw, cw);
    //compute_cost_volume_warp(frames, cameras, ref_frame, S, depth_min, depth_max, cost_volume);
    ws->compute_g(frames[ref_frame], g, gw_alpha, 1.0);

    if(cfg->is_set("cost_volume_file"))
    {
      vcl_string cost_file = cfg->get_value<vcl_string>("cost_volume_file");
      save_cost_volume(cost_volume, g, cost_file.c_str());
    }
  }
  else if (cfg->is_set("cost_volume_file"))
  {
    vcl_string cost_file = cfg->get_value<vcl_string>("cost_volume_file");
    if (cfg->is_set("mix_cost_volumes") &&
        cfg->get_value<bool>("mix_cost_volumes"))
    {
      double iw = cfg->get_value<double>("intensity_cost_weight");
      double gw = cfg->get_value<double>("gradient_cost_weight");
      double cw = cfg->get_value<double>("census_cost_weight");
      int pos = cost_file.find_last_of(".");
      vcl_string ext = cost_file.substr(pos);
      vcl_string basename = cost_file.substr(0, pos);
      vil_image_view<double> tmp;
      load_cost_volume(cost_volume, g, (basename+"_intensity"+ext).c_str());
      load_cost_volume(tmp, g, (basename+"_gradient"+ext).c_str());
      vil_math_add_image_fraction(cost_volume, iw, tmp, gw);
      load_cost_volume(tmp, g, (basename+"_census"+ext).c_str());
      vil_math_add_image_fraction(cost_volume, 1.0, tmp, cw);
    }
    else
    {
      load_cost_volume(cost_volume, g, cost_file.c_str());
    }
  }
  else
  {
    vcl_cerr << "Error: must either set compute_cost_volume"
             << " or specify a cost_volume_file to load." << vcl_endl;
    return -1;
  }


  vcl_cout << "Refining Depth. ..\n";
  vil_image_view<double> depth(cost_volume.ni(), cost_volume.nj(), 1);

#ifndef USE_BP

  if (!cfg->is_set("compute_init_depthmap") ||
       cfg->get_value<bool>("compute_init_depthmap") )
  {
#ifdef HAVE_VISCL
    if (cfg->is_set("use_gpu") &&
        cfg->get_value<bool>("use_gpu") )
    {
      refine_depth_cl_t rd = NEW_VISCL_TASK(refine_depth_cl);
      vil_image_view<float> depth_f(cost_volume.ni(), cost_volume.nj(), 1);
      vil_image_view<float> cv_f(cost_volume.ni(), cost_volume.nj(), 1, cost_volume.nplanes());
      vil_image_view<float> g_f;
      vil_convert_cast<double, float>(g, g_f);

      {
        vil_image_view<float> cv_float;
        vil_convert_cast<double, float>(cost_volume, cv_float);
        vil_copy_reformat<float>(cv_float, cv_f);
      }

      boost::chrono::system_clock::time_point start = boost::chrono::system_clock::now();
      rd->refine(cv_f, depth_f, g_f, beta, theta0, theta_end, lambda);
      boost::chrono::duration<double> sec = boost::chrono::system_clock::now() - start;
      vcl_cout << "viscl took " << sec.count() << " seconds.\n";
      vil_convert_cast<float, double>(depth_f, depth);
    }
    else
#endif
    {
      boost::chrono::system_clock::time_point start = boost::chrono::system_clock::now();
      refine_depth(cost_volume, g, depth, 2000, theta0, theta_end, lambda, epsilon);
      boost::chrono::duration<double> sec = boost::chrono::system_clock::now() - start;
      vcl_cout << "super3d took " << sec.count() << " seconds.\n";
    }
    if(cfg->is_set("init_depthmap_file"))
    {
      vcl_string depthmap_file = cfg->get_value<vcl_string>("init_depthmap_file");
      save_depth(depth, depthmap_file.c_str());
    }
  }
  else if (cfg->is_set("init_depthmap_file"))
  {
    vcl_string depthmap_file = cfg->get_value<vcl_string>("init_depthmap_file");
    load_depth(depth, depthmap_file.c_str());
  }
  else
  {
    vcl_cerr << "Error: must either set compute_init_depthmap"
             << " or specify an init_depthmap_file to load." << vcl_endl;
    return -1;
  }

  //refine_depth_planar(ws, depth, g, frames[ref_frame], 1, 0.25);
  //refine_depth_planar(cost_volume, ws, depth, g, frames[ref_frame], beta, theta0, theta_end, lambda);
  save_depth(depth, "depth_map_normals.dat");

#else
  bp_refine(cost_volume, depth_min, depth_max, depth);    //TODO: need to check if this works after idepth->depth change
#endif

  if (cfg->is_set("output_depthmap"))
  {
    vil_image_view<vxl_byte> dmap;
    vil_convert_stretch_range_limited(depth, dmap, 0.0, 1.0);
    // depth map are drawn inverted (white == closest) for viewing
    vil_math_scale_and_offset_values(dmap, -1.0, 255);
    vcl_string depthmap_file = cfg->get_value<vcl_string>("output_depthmap");
    vil_save(dmap, depthmap_file.c_str());
  }

  vcl_cout << "writing mesh"<<vcl_endl;
#ifdef HAVE_VTK
  //vtp depth writer uses 0..1 depth scaling
  vil_image_view<vxl_byte> b_texture;
  if (cfg->is_set("directory"))
    b_texture = vil_load((cfg->get_value<vcl_string>("directory") + filenames[ref_frame]).c_str());
  else
    b_texture = vil_load(filenames[ref_frame].c_str());

  vil_image_view<double> d_texture;
  vil_convert_cast<vxl_byte, double>(b_texture, d_texture);
  vcl_string output_file_name = cfg->get_value<vcl_string>("output_file");
  save_depth_to_vtp(output_file_name.c_str(), depth, d_texture, ref_cam, ws);
#endif

  // map depth from normalized range back into true depth
  double depth_scale = depth_max - depth_min;
  vil_math_scale_and_offset_values(depth, depth_scale, depth_min);

  if (cfg->is_set("obj_file"))
  {
    double output_decimate = 1.0;
    imesh_mesh nm = depth_map_to_mesh(scale_camera(cameras[ref_frame], 1.0/output_decimate),
                                      vil_decimate(depth,output_decimate));
    imesh_write_obj(cfg->get_value<vcl_string>("obj_file"), nm);
  }

  if (cfg->is_set("output_float_depthmap"))
  {
    vcl_string depthmap_file = cfg->get_value<vcl_string>("output_float_depthmap");
    vil_save(depth, depthmap_file.c_str());
  }

  vil_image_view<double> gt;
  if (cfg->is_set("ground_truth"))
  {
    vcl_string gtfile = cfg->get_value<vcl_string>("ground_truth");
    gt = vil_load(gtfile.c_str());
    gt = vil_crop(gt, i0, ni, j0, nj);
    vil_image_view<double> diff;
    vil_math_image_difference(gt, depth, diff);
    vil_image_view<vxl_byte> diff_byte;
    vil_convert_stretch_range_limited(diff, diff_byte, -depth_scale, depth_scale);
    vil_save(diff_byte, "diff.png");

    double threshold = depth_scale / S;
    score_vs_gt(depth, gt, threshold);
    vcl_cout << "cost of estimated depth: " << eval_hessian_frob(depth, cost_volume, lambda) << "\n";
    vcl_cout << "cost of gt depth: " << eval_hessian_frob(gt, cost_volume, lambda) << "\n";
  }

  if (ws) delete ws;

  } catch (const config::cfg_exception &e)
  {
    vcl_cout << "Error in config: " << e.what() << "\n";
  }

  return 0;
}

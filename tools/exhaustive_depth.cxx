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

#include "depth_cl/refine_depth.h"

#include <boost/chrono.hpp>

#include <vpgl/vpgl_perspective_camera.h>

//C:/Data/gisr/cameras.txt c:/Data/meshes/crop.obj -fl C:/Data/gisr/frames_gisr.txt -dir C:/Data/gisr/ -f 0 -wv "0.026 0.030895 0.002363 0.123 0.09583 0.02121 1000 779"

//C:/Data/epfl_dataset/fountain-P11/cams.txt c:/Data/meshes/crop.obj -fl C:/Data/epfl_dataset/fountain-P11/frames.txt -dir C:/Data/epfl_dataset/fountain-P11/urd/ -f 3 -cw "563 1843 416 987"
//C:/Data/synthetic_data/simple2/cams.txt c:/Data/meshes/crop.obj -fl C:/Data/synthetic_data/simple2/frames.txt -dir C:/Data/synthetic_data/simple2/images/ -f 1

//C:/Data/ground_medians/cameras.txt c:/Data/meshes/crop.obj -fl frames.txt -dir C:/Data/ground_medians/full_window15_reset50/ -f 0 -e C:/Data/exposure/exposure.txt -wv "-0.0060024489 0.0049344006 0.00042798251 0.028599165 0.021470546 0.0051992279 1320 954"

//C:/Data/ground_medians/cameras.txt c:/Data/meshes/crop.obj -ff C:/Data/ground_medians/full_window15_reset50/frame###.png,000:100:950 -f 0 -e C:/Data/exposure/exposure.txt -cw "767 1320 757 954"
//C:/Data/ground_medians/cameras.txt c:/Data/meshes/crop.obj -ff C:/Data/ground_medians/full_window15_reset50/frame###.png,000:100:950 -f 0 -e C:/Data/exposure/exposure.txt -wv "-0.0060024489 0.0049344006 0.00042798251 0.028599165 0.021470546 0.0051992279 1320 954"

//C:/Data/ground_medians/cameras.txt c:/Data/meshes/crop.obj -ff C:/Data/ground_medians/full_window15_reset50/frame###.png,000:100:950 -f 0 -e C:/Data/exposure/exposure.txt -cw "962 334 1504 205"
//C:/Data/ground_medians/cameras.txt c:/Data/meshes/crop.obj -ff C:/Data/ground_medians/full_window15_reset50/frame###.png,000:100:950 -f 0 -e C:/Data/exposure/exposure.txt -wv "0.00068379926 0.012289556 0.00042798251 0.0071305364 0.0047269018 0.0051992279 334 205"
int main(int argc, char* argv[])
{
  try
  {
    config::inst()->read_config(argv[1]);


  vcl_vector<vil_image_view<double> > frames;
  vcl_vector<vpgl_perspective_camera<double> >  cameras;
  vcl_vector<vcl_string> filenames;
  vcl_string camera_file = config::inst()->get_value<vcl_string>("camera_file");

  if (config::inst()->is_set("frame_format"))
  {
    vcl_cout << "Using formated string to find images/cameras.\n";
    vul_sequence_filename_map frame_seq(config::inst()->get_value<vcl_string>("frame_format"));

    //Read Cameras
    cameras = load_cams(camera_file, frame_seq);

    load_from_frame_file(frame_file.c_str(), dir, filenames, frameindex, frames);
    for (unsigned int i = 0; i < filenames.size(); i++)
    {
      vcl_string camname = filenames[i];
      unsigned int found = camname.find_last_of("/\\");
      camname = camname.substr(found+1, camname.size() - 4 - found - 1);
      camname = config::inst()->get_value<vcl_string>("camera_dir") + "/" + camname + ".krtd";
      cameras.push_back(load_cam(camname));
    }

    }
    else if (frames.size() < 2)
    {
      vcl_cerr << "At least 2 frames are required"<<vcl_endl;
      return -1;
    }

    //Apply exposure correction
    if (config::inst()->is_set("exposure_file"))
    {
      vcl_vector<vcl_pair<double, double> > exposures;
      exposures = load_exposure(config::inst()->get_value<vcl_string>("exposure_file"), frame_seq);

      for (unsigned int i = 0; i < frames.size(); i++)
        apply_exposure_correction(frames[i], exposures[i]);
    }
  }
  else if (config::inst()->is_set("frame_list"))
  {
    vcl_string frame_file = config::inst()->get_value<vcl_string>("frame_list");
    vcl_string dir = config::inst()->get_value<vcl_string>("directory");
    vcl_cout << "Using frame file: " << frame_file << " to find images and " << camera_file  << " to find cameras.\n";
    vcl_vector<int> frameindex;
    load_from_frame_file(frame_file.c_str(), camera_file.c_str(), dir, filenames, frameindex, frames, cameras);
    if (config::inst()->is_set("exposure_file"))
    {
      vcl_vector<vcl_pair<double, double> > exposures;
      exposures = load_exposure(config::inst()->get_value<vcl_string>("exposure_file"), frameindex);

      for (unsigned int i = 0; i < frames.size(); i++)
        apply_exposure_correction(frames[i], exposures[i]);
    }
  }
  else
  {
    vcl_cerr << "Error: must use either frame list (-fl) or frame format string (-ff).\n";
    return -1;
  }

  if (!config::inst()->is_set("exposure_file"))
  {
    for (unsigned int i = 0; i < frames.size(); i++)
      vil_math_scale_and_offset_values(frames[i], 1.0/255.0, 0.0);
  }

  unsigned int ref_frame = config::inst()->get_value<unsigned int>("ref_frame");
  vpgl_perspective_camera<double> ref_cam = cameras[ref_frame];

  world_space *ws = NULL;
  int i0, ni, j0, nj;

  if (config::inst()->is_set("world_volume"))
  {
    vcl_istringstream wvstream(config::inst()->get_value<vcl_string>("world_volume"));
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
    if (config::inst()->is_set("camera_scale"))
    {
      camera_scale = config::inst()->get_value<double>("camera_scale");
      for (unsigned int i = 0; i < cameras.size(); i++)
        cameras[i] = scale_camera(cameras[i], camera_scale);
    }
    //Compute the window cropping, scale the cropping by the specified scale so that we do not
    //need to recompute cropping input for super resolution.
    if (config::inst()->is_set("crop_window"))
    {
      vcl_istringstream cwstream(config::inst()->get_value<vcl_string>("crop_window"));
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

    if (config::inst()->get_value<bool>("compute_depth_range") && config::inst()->is_set("landmarks_path"))
    {
      compute_depth_range(cameras[ref_frame], i0, ni, j0, nj, config::inst()->get_value<vcl_string>("landmarks_path"), depth_min, depth_max);
      vcl_cout << "depth range: " << depth_min << " - " << depth_max << "\n";
    }
    else
    {
      depth_min = config::inst()->get_value<double>("depth_min");
      depth_max = config::inst()->get_value<double>("depth_max");
    }

    ws = new world_frustum(cameras[ref_frame], depth_min, depth_max, ni, nj);
  }

  vcl_cout << "Refining depth"<<vcl_endl;
  unsigned int S = config::inst()->get_value<unsigned int>("num_slices");
  double theta0 = config::inst()->get_value<double>("theta_start");
  double theta_end = config::inst()->get_value<double>("theta_end");
  double beta = config::inst()->get_value<double>("beta");
  double lambda = config::inst()->get_value<double>("lambda");
  double gw_alpha = config::inst()->get_value<double>("gw_alpha");

  vil_image_view<double> g;
  vil_image_view<double> cost_volume;

  if (!config::inst()->is_set("compute_cost_volume") ||
      config::inst()->get_value<bool>("compute_cost_volume") )
  {
    double iw = config::inst()->get_value<double>("intensity_cost_weight");
    double gw = config::inst()->get_value<double>("gradient_cost_weight");
    double cw = config::inst()->get_value<double>("census_cost_weight");
    compute_world_cost_volume(frames, cameras, ws, ref_frame, S, cost_volume, iw, gw, cw);
    ws->compute_g(frames[ref_frame], g, gw_alpha, 1.0);

    if(config::inst()->is_set("cost_volume_file"))
    {
      vcl_string cost_file = config::inst()->get_value<vcl_string>("cost_volume_file");
      save_cost_volume(cost_volume, g, cost_file.c_str());
    }
  }
  else if (config::inst()->is_set("cost_volume_file"))
  {
    vcl_string cost_file = config::inst()->get_value<vcl_string>("cost_volume_file");
    if (config::inst()->is_set("mix_cost_volumes") &&
        config::inst()->get_value<bool>("mix_cost_volumes"))
    {
      double iw = config::inst()->get_value<double>("intensity_cost_weight");
      double gw = config::inst()->get_value<double>("gradient_cost_weight");
      double cw = config::inst()->get_value<double>("census_cost_weight");
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

  if (!config::inst()->is_set("compute_init_depthmap") ||
       config::inst()->get_value<bool>("compute_init_depthmap") )
  {
    if (config::inst()->is_set("use_gpu") &&
        config::inst()->get_value<bool>("use_gpu") )
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
    {
      boost::chrono::system_clock::time_point start = boost::chrono::system_clock::now();
      refine_depth(cost_volume, g, depth, beta, theta0, theta_end, lambda);
      boost::chrono::duration<double> sec = boost::chrono::system_clock::now() - start;
      vcl_cout << "super3d took " << sec.count() << " seconds.\n";
    }
    if(config::inst()->is_set("init_depthmap_file"))
    {
      vcl_string depthmap_file = config::inst()->get_value<vcl_string>("init_depthmap_file");
      save_depth(depth, depthmap_file.c_str());
    }
  }
  else if (config::inst()->is_set("init_depthmap_file"))
  {
    vcl_string depthmap_file = config::inst()->get_value<vcl_string>("init_depthmap_file");
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

  if (config::inst()->is_set("output_depthmap"))
  {
    vil_image_view<vxl_byte> dmap;
    vil_convert_stretch_range_limited(depth, dmap, 0.0, 1.0);
    // depth map are drawn inverted (white == closest) for viewing
    vil_math_scale_and_offset_values(dmap, -1.0, 255);
    vcl_string depthmap_file = config::inst()->get_value<vcl_string>("output_depthmap");
    vil_save(dmap, depthmap_file.c_str());
  }

  vcl_cout << "writing mesh"<<vcl_endl;
#ifdef HAVE_VTK
  //vtp depth writer uses 0..1 depth scaling
  vil_image_view<vxl_byte> b_texture;
  if (config::inst()->is_set("directory"))
    b_texture = vil_load((config::inst()->get_value<vcl_string>("directory") + filenames[ref_frame]).c_str());
  else
    b_texture = vil_load(filenames[ref_frame].c_str());

  vil_image_view<double> d_texture;
  vil_convert_cast<vxl_byte, double>(b_texture, d_texture);
  vcl_string output_file_name = config::inst()->get_value<vcl_string>("output_file");
  save_depth_to_vtp(output_file_name.c_str(), depth, d_texture, ref_cam, ws);
#endif

  // map depth from normalized range back into true depth
  double depth_scale = depth_max - depth_min;
  vil_math_scale_and_offset_values(depth, depth_scale, depth_min);

  if (config::inst()->is_set("obj_file"))
  {
    double output_decimate = 1.0;
    imesh_mesh nm = depth_map_to_mesh(scale_camera(cameras[ref_frame], 1.0/output_decimate),
                                      vil_decimate(depth,output_decimate));
    imesh_write_obj(config::inst()->get_value<vcl_string>("obj_file"), nm);
  }

  if (config::inst()->is_set("output_float_depthmap"))
  {
    vcl_string depthmap_file = config::inst()->get_value<vcl_string>("output_float_depthmap");
    vil_save(depth, depthmap_file.c_str());
  }

  vil_image_view<double> gt;
  if (config::inst()->is_set("ground_truth"))
  {
    vcl_string gtfile = config::inst()->get_value<vcl_string>("ground_truth");
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

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
#include <super3d/depth/file_io.h>
#include <super3d/depth/multiscale.h>
#include <super3d/depth/super_res.h>
#include <super3d/depth/flow_manip.h>
#include <super3d/depth/depth_map.h>
#include <super3d/depth/normal_map.h>
#include <super3d/depth/weighted_dbw.h>

#include <super3d/image/adjoint_flow_warp.h>
#include <super3d/image/adjoint_resample.h>

#include <iostream>
#include <boost/bind.hpp>
#include <boost/scoped_ptr.hpp>

#include <vil/vil_crop.h>
#include <vil/vil_bicub_interp.h>
#include <vil/vil_load.h>
#include <vil/vil_save.h>
#include <vil/vil_convert.h>
#include <vil/vil_resample_bilin.h>
#include <vil/vil_resample_bicub.h>
#include <vil/vil_crop.h>
#include <vil/vil_clamp.h>

#include <vpgl/vpgl_perspective_camera.h>
#include <vgl/vgl_intersection.h>


/// Use the valid region of the flow destinations to crop the frames and translate the flow.
void crop_frames_and_flows(std::vector<vil_image_view<double> > &flows,
                             std::vector<vil_image_view<double> > &frames,
                             int scale_factor,
                             int margin=0);

void create_warps_from_flows(const std::vector<vil_image_view<double> > &flows,
                             const std::vector<vil_image_view<double> > &weights,
                             const std::vector<vil_image_view<double> > &frames,
                             std::vector<super3d::adjoint_image_ops_func<double> > &warps,
                             int scale_factor,
                             super3d::config *cfg);

void difference_from_flow(const vil_image_view<double> &I0,
                          const vil_image_view<double> &I1,
                          const vil_image_view<double> &flow,
                          vil_image_view<double> &diff,
                          double scale_factor,
                          super3d::config *cfg);

void create_low_res(std::vector<vil_image_view<double> > &frames,
                    int scale,
                    super3d::config *cfg);

void load_flow(const char *flow_list, const std::string &dir, std::vector<vil_image_view<double> > &flows);

//*****************************************************************************

int main(int argc, char* argv[])
{
  try
  {
    boost::scoped_ptr<super3d::config> cfg(new super3d::config);
    cfg->read_config(argv[1]);
    std::vector<vil_image_view<double> > frames;
    std::vector<std::string> filenames;

    super3d::super_res_params srp;
    read_super_res_params( cfg, srp );

    std::string frame_file = cfg->get_value<std::string>("frame_list");
    std::string dir("");
    if (cfg->is_set("directory"))
      dir = cfg->get_value<std::string>("directory");
    std::vector<int> frameindex;

    std::vector<vpgl_perspective_camera<double> >  cameras;

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

    const unsigned int ref_frame = cfg->get_value<unsigned int>("ref_frame");
    const double scale_factor = cfg->get_value<double>("scale_factor");

    vil_image_view<double> gt;

    vil_image_view<double> ref_image;
    int i0 = 0;
    int j0 = 0;
    int ni, nj;

    std::vector<vil_image_view<double> > weights;
    std::vector<super3d::adjoint_image_ops_func<double> > warps;

    const double normalizer = 1.0/255.0;


    double camera_scale;

    //This executable can create its own low-res images to compare against the high res
    //images defined in the config or read from specific images
    if (cfg->get_value<bool>("create_low_res"))
    {
      gt.deep_copy(frames[ref_frame]);
      create_low_res(frames, scale_factor, cfg.get());
      camera_scale = scale_factor;
    }
    else
    {
      camera_scale = cfg->get_value<double>("camera_scale");
      if (cfg->is_set("ground_truth"))
      {
        std::string ground_truth = cfg->get_value<std::string>("ground_truth");
        vil_image_view<vxl_byte> gt_byte;
        gt_byte = vil_load(ground_truth.c_str());
        vil_convert_planes_to_grey(gt_byte, gt);
      }
    }

    std::vector<vpgl_perspective_camera<double> > scaled_cams;
    for (unsigned int i = 0; i < frames.size(); i++)
    {
      vil_math_scale_values(frames[i], normalizer);
      scaled_cams.push_back(super3d::scale_camera(cameras[i], scale_factor / camera_scale));
    }

    vil_image_view<double> depth;
    std::string depthfile = cfg->get_value<std::string>("high_res_depth");
    depth = vil_load(depthfile.c_str());
    int dni = depth.ni(), dnj = depth.nj();
    if (cfg->get_value<bool>("upsample_depth"))
    {
      dni *= scale_factor;
      dnj *= scale_factor;
      vil_image_view<double> temp;
      vil_resample_bilin(depth, temp, dni, dnj);
      depth = temp;
    }

    std::cout << depth.ni() << " " << depth.nj() << "\n";
    double min, max;
    vil_math_value_range(depth, min, max);
    std::cout << min << " " << max << "\n";

    //Crop regions are defined in the coordinates of the LOW RES reference frame
    vpgl_perspective_camera<double> ref_cam = scaled_cams[ref_frame];
    if (cfg->is_set("crop_window"))
    {
      std::istringstream cwstream(cfg->get_value<std::string>("crop_window"));
      cwstream >> i0 >> ni >> j0 >> nj;
      std::cout << frames[ref_frame].ni() << " " << frames[ref_frame].nj() << "\n";
      ref_image.deep_copy(vil_crop(frames[ref_frame], i0, ni, j0, nj));
      i0 = (int)(i0*scale_factor/camera_scale);
      j0 = (int)(j0*scale_factor/camera_scale);
      ni = (int)(ni*scale_factor/camera_scale);
      nj = (int)(nj*scale_factor/camera_scale);
      std::cout << "Crop window: " << i0 << " " << ni << " " << j0 << " " << nj << "\n";
      if (cfg->is_set("ground_truth"))
        gt = vil_crop(gt, i0, ni, j0, nj);

      if (cfg->is_set("crop_depth") && cfg->get_value<bool>("crop_depth"))
        depth = vil_crop(depth, i0, ni, j0, nj);
      ref_cam = super3d::crop_camera(ref_cam, i0, j0);
    }
    else
    {
      ni = frames[ref_frame].ni();
      nj = frames[ref_frame].nj();
      ref_image.deep_copy(frames[ref_frame]);
    }

    weights.resize(cameras.size());

    if (cfg->get_value<bool>("use_normal_weighting"))
    {
      vil_image_view<double> location_map, normal_map, ref_angle_map;
      super3d::depth_map_to_normal_map_inv_len(ref_cam, depth, normal_map);
      super3d::depth_map_to_location_map(ref_cam, depth, location_map);
      super3d::viewing_angle_map(ref_cam.camera_center(), location_map, normal_map, ref_angle_map);

      if( srp.debug )
      {
        vil_image_view<double> rotated;
        vil_image_view<vxl_byte> byte_normals;
        super3d::rotate_normal_map(normal_map, ref_cam.get_rotation(), rotated);
        super3d::byte_normal_map(rotated, byte_normals);
        vil_save(byte_normals, "images/normal_map.png");
      }

      for (unsigned int i = 0; i < cameras.size(); i++)
      {
        super3d::viewing_angle_map(cameras[i].camera_center(), location_map, normal_map, weights[i]);
        vil_math_image_ratio(weights[i], ref_angle_map, weights[i]);
        vil_clamp_below(weights[i], 0.0);
      }
    }
    else
    {
      for (unsigned int i = 0; i < cameras.size(); i++)
      {
        weights[i].set_size(dni, dnj);
        weights[i].fill(1.0);
      }
    }

    std::cout << "Computing flow from depth\n";
    std::vector<vil_image_view<double> > flows;
    super3d::compute_occluded_flows_from_depth(scaled_cams, ref_cam, depth, flows);
    crop_frames_and_flows(flows, frames, scale_factor, 0);
    create_warps_from_flows(flows, frames, weights, warps, scale_factor, cfg.get());


    if ( srp.debug )
    {

      for (unsigned int i = 0; i < warps.size(); i++)
      {
        // test if the DBW operator for image i is adjoint
        std::cout << "is adjoint "<<i<<" "<<warps[i].is_adjoint()<<std::endl;

        char buf[50];
        vil_image_view<vxl_byte> output;

        // cropped frames
        sprintf(buf, "images/frames-%03d.png", i);
        vil_convert_stretch_range_limited(frames[i], output, 0.0, 1.0);
        vil_save(output, buf);

        // apply DBW to ground truth super res image to predict frame i
        vil_image_view<double> pred;
        if (cfg->is_set("ground_truth"))
        {
          vil_image_view<double> gts;
          vil_convert_stretch_range_limited(gt, gts, 0.0, 255.0, 0.0, 1.0);
          warps[i].apply_A(gts, pred);
          vil_convert_stretch_range_limited(pred, output, 0.0, 1.0);
          sprintf(buf, "images/predicted-%03d.png", i);
          vil_save(output, buf);
        }

        // write out the alpha weight map
        vil_image_view<double> map = warps[i].weight_map();
        vil_convert_stretch_range_limited(map, output, 0.0, 2.0);
        sprintf(buf, "images/weight-%03d.png", i);
        vil_save(output, buf);

        // apply DBW alpha map weighting to frame i
        vil_image_view<double> wframe;
        vil_math_image_product(frames[i], map, wframe);
        vil_convert_stretch_range_limited(wframe, output, 0.0, 1.0);
        sprintf(buf, "images/weighted-frame-%03d.png", i);
        vil_save(output, buf);

        // write weighted difference between prediction and data
        if (cfg->is_set("gronund_truth"))
        {
          vil_math_image_difference(wframe, pred, pred);
          vil_convert_stretch_range_limited(pred, output, -1.0, 1.0);
          sprintf(buf, "images/predition-error%03d.png", i);
          vil_save(output, buf);
        }

        // write out the viewing angle based weights
        vil_convert_stretch_range_limited(weights[i], output, 0.0, 2.0);
        sprintf(buf, "images/angle-weight-%03d.png", i);
        vil_save(output, buf);

        // apply DBW inverse to frames
        vil_image_view<double> register_frames;
        warps[i].apply_At(frames[i], register_frames);
        vil_convert_stretch_range_limited(register_frames, output, 0.0, 1.0);
        sprintf(buf, "images/register-frames-%03d.png", i);
        vil_save(output, buf);
      }
    }

    std::cout << "Computing super resolution\n";

    //Initilize super resolution parameters
    vil_image_view<double> super_u;
    srp.scale_factor = scale_factor;
    srp.ref_frame = ref_frame;
    srp.s_ni = warps[ref_frame].src_ni();
    srp.s_nj = warps[ref_frame].src_nj();
    srp.l_ni = warps[ref_frame].dst_ni();
    srp.l_nj = warps[ref_frame].dst_nj();

    const unsigned int iterations = cfg->get_value<unsigned int>("iterations");
    std::string output_image = cfg->get_value<std::string>("output_image");

    std::vector< vil_image_view<double> > As;
    super3d::super_resolve_robust(frames, warps, super_u, srp, iterations, As, output_image);

    if (cfg->is_set("ground_truth"))
    {
      vil_math_scale_values(gt, normalizer);
      super3d::compare_to_original(ref_image, super_u, gt, scale_factor);
    }
    else
    {
      vil_image_view<double> upsamp;
      super3d::upsample(ref_image, upsamp, scale_factor, super3d::warp_image_parameters::CUBIC);

      vil_image_view<vxl_byte> output;
      vil_convert_stretch_range_limited(upsamp, output, 0.0, 1.0);
      vil_save(output, "bicub.png");
    }

    vil_image_view<vxl_byte> output;
    vil_convert_stretch_range_limited(super_u, output, 0.0, 1.0);
    vil_save(output, output_image.c_str());

    if( srp.illumination_prior )
    {
      for (unsigned int i = 0; i < As.size(); i++)
      {
        char buf[50];
        vil_image_view<vxl_byte> output;
        sprintf(buf, "A%d-%03d.png", i%2, i/2 );
        vil_convert_stretch_range_limited(As[i], output, -1.0, 1.0);
        vil_save(output, buf);
      }
    }
  }
  catch (const super3d::config::cfg_exception &e)  {
    std::cout << "Error in config: " << e.what() << "\n";
  }

  return 0;
}

//*****************************************************************************

/// Use the valid region of the flow destinations to crop the frames and translate the flow.
void crop_frames_and_flows(std::vector<vil_image_view<double> > &flows,
                           std::vector<vil_image_view<double> > &frames,
                           int scale_factor, int margin)
{
  assert(flows.size() == frames.size());
  for (unsigned int i=0; i<flows.size(); ++i)
  {
    // get a bounding box around the flow and expand by a margin
    vgl_box_2d<double> bounds = super3d::flow_destination_bounds(flows[i]);
    vgl_box_2d<int> bbox(std::floor(bounds.min_x()), std::ceil(bounds.max_x()),
                         std::floor(bounds.min_y()), std::ceil(bounds.max_y()));
    bbox.expand_about_centroid(2*margin);

    // get a bounding box around the (up sampled) image and intersect
    vgl_box_2d<int> img_box(0, frames[i].ni()-1, 0, frames[i].nj()-1);
    img_box.scale_about_origin(scale_factor);
    bbox = vgl_intersection(bbox, img_box);
    bbox.scale_about_origin(1.0 / scale_factor);

    // crop the image and translated the flow accordingly
    frames[i] = vil_crop(frames[i], bbox.min_x(), bbox.width(), bbox.min_y(), bbox.height());
    bbox.scale_about_origin(scale_factor);
    super3d::translate_flow(flows[i], -bbox.min_x(), -bbox.min_y());
  }
}

//*****************************************************************************

void create_warps_from_flows(const std::vector<vil_image_view<double> > &flows,
                             const std::vector<vil_image_view<double> > &frames,
                             const std::vector<vil_image_view<double> > &weights,
                             std::vector<super3d::adjoint_image_ops_func<double> > &warps,
                             int scale_factor,
                             super3d::config *cfg)
{
  assert(flows.size() == frames.size());
  bool down_sample_averaging = cfg->get_value<bool>("down_sample_averaging");
  bool bicubic_warping = cfg->get_value<bool>("bicubic_warping");
  double sensor_sigma = cfg->get_value<double>("sensor_sigma");

  warps.clear();
  warps.reserve(flows.size());
  for (unsigned int i=0; i<flows.size(); ++i)
  {
    warps.push_back(super3d::create_dbw_from_flow(flows[i], weights[i],
                                                  frames[i].ni(), frames[i].nj(),
                                                  frames[i].nplanes(), scale_factor,
                                                  sensor_sigma,
                                                  down_sample_averaging,
                                                  bicubic_warping));
  }
}


//*****************************************************************************

void difference_from_flow(const vil_image_view<double> &I0,
                          const vil_image_view<double> &I1,
                          const vil_image_view<double> &flow,
                          vil_image_view<double> &diff,
                          double scale_factor,
                          super3d::config *cfg)
{
  bool bicubic_warping = cfg->get_value<bool>("bicubic_warping");
  vil_image_view<double> I0_x, I1_x;
  vil_resample_bicub(I0, I0_x, scale_factor * I0.ni(), scale_factor * I0.nj());
  vil_resample_bicub(I1, I1_x, scale_factor * I1.ni(), scale_factor * I1.nj());
  diff.set_size(I0_x.ni(), I0_x.nj(), 1);

  vil_image_view<double> temp;
  if( bicubic_warping )
  {
    super3d::warp_back_with_flow_bicub(I1_x, flow, temp);
  }
  else
  {
    super3d::warp_back_with_flow_bilin(I1_x, flow, temp);
  }
  vil_math_image_difference(temp, I0_x, diff);
}

//*****************************************************************************

void create_low_res(std::vector<vil_image_view<double> > &frames,
                    int scale,
                    super3d::config *cfg)
{
  bool down_sample_averaging = cfg->get_value<bool>("down_sample_averaging");
  for (unsigned int i = 0; i < frames.size(); i++)
  {
    vil_image_view<double> temp;
    double sensor_sigma = 0.25 * sqrt(scale * scale - 1.0);
    vil_gauss_filter_2d(frames[i], temp, sensor_sigma, 3.0*sensor_sigma);
    if( down_sample_averaging )
    {
      super3d::down_scale(temp, frames[i], scale);
    }
    else
    {
      super3d::down_sample(temp, frames[i], scale);
    }
  }
}

//*****************************************************************************

void load_flow(const char *flow_list, const std::string &dir, std::vector<vil_image_view<double> > &flows)
{
  std::ifstream infile(flow_list);
  std::string flowname;

  while (infile >> flowname)
  {
    vil_image_view<double> flowimg;
    super3d::read_flow_file(flowimg, (dir + flowname).c_str());
    flows.push_back(flowimg);
  }
  std::cout << "Read " << flows.size() << " flows.\n";
}

//*****************************************************************************

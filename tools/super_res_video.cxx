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

#include <super3d/depth/super_config.h>
#include <super3d/depth/file_io.h>
#include <super3d/depth/multiscale.h>
#include <super3d/depth/super_res.h>
#include <super3d/depth/flow_manip.h>
#include <super3d/depth/depth_map.h>
#include <super3d/depth/normal_map.h>
#include <super3d/depth/weighted_dbw.h>

#include <video_transforms/adjoint_flow_warp.h>
#include <video_transforms/adjoint_resample.h>

#include <vcl_iostream.h>
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
#include <video_transforms/warp_image.h>

//#define DEBUG


/// Use the valid region of the flow destinations to crop the frames and translate the flow.
void crop_frames_and_flows(vcl_vector<vil_image_view<double> > &flows,
                             vcl_vector<vil_image_view<double> > &frames,
                             int scale_factor,
                             int margin=0);

void create_warps_from_flows(const vcl_vector<vil_image_view<double> > &flows,
                             const vcl_vector<vil_image_view<double> > &weights,
                             const vcl_vector<vil_image_view<double> > &frames,
                             vcl_vector<vidtk::adjoint_image_ops_func<double> > &warps,
                             int scale_factor,
                             super3d::config *cfg);

void difference_from_flow(const vil_image_view<double> &I0,
                          const vil_image_view<double> &I1,
                          const vil_image_view<double> &flow,
                          vil_image_view<double> &diff,
                          double scale_factor,
                          super3d::config *cfg);

void create_low_res(vcl_vector<vil_image_view<double> > &frames,
                    int scale,
                    super3d::config *cfg);

void upsample(const vil_image_view<double> &src, vil_image_view<double> &dest,
              double scale_factor, vidtk::warp_image_parameters::interp_type interp)
{
  vidtk::warp_image_parameters wip;
  wip.set_fill_unmapped(true);
  wip.set_unmapped_value(-1.0);
  wip.set_interpolator(interp);

  vnl_double_3x3 Sinv;
  Sinv.set_identity();
  Sinv(0,0) = 1.0/scale_factor;
  Sinv(1,1) = 1.0/scale_factor;

  dest.set_size(src.ni() * scale_factor, src.nj() * scale_factor, src.nplanes());
  vidtk::warp_image(src, dest, Sinv, wip);
}

//*****************************************************************************

int main(int argc, char* argv[])
{
  try  {
    boost::scoped_ptr<super3d::config> cfg(new super3d::config);
    cfg->read_config(argv[1]);
    vcl_vector<vil_image_view<double> > frames;
    vcl_vector<vpgl_perspective_camera<double> >  cameras;
    vcl_vector<vcl_string> filenames;

    vcl_string camera_file = cfg->get_value<vcl_string>("camera_file");
    vcl_string frame_file = cfg->get_value<vcl_string>("frame_list");
    vcl_string dir = cfg->get_value<vcl_string>("directory");
    vcl_cout << "Using frame file: " << frame_file << " to find images and "
             << camera_file  << " to find cameras.\n";
    vcl_vector<int> frameindex;
    super3d::load_from_frame_file(frame_file.c_str(), camera_file.c_str(), dir,
                                  filenames, frameindex, frames, cameras, cfg->get_value<bool>("use_color"));

    const unsigned int ref_frame = cfg->get_value<unsigned int>("ref_frame");
    const double scale_factor = cfg->get_value<double>("scale_factor");
    double camera_scale;

    vil_image_view<double> gt;

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
        vcl_string ground_truth = cfg->get_value<vcl_string>("ground_truth");
        vil_image_view<vxl_byte> gt_byte;
        gt_byte = vil_load(ground_truth.c_str());
        vil_convert_planes_to_grey(gt_byte, gt);
      }
    }

    const double normalizer = 1.0/255.0;
    vcl_vector<vpgl_perspective_camera<double> > scaled_cams;
    for (unsigned int i = 0; i < frames.size(); i++)
    {
      vil_math_scale_values(frames[i], normalizer);
      scaled_cams.push_back(super3d::scale_camera(cameras[i], scale_factor / camera_scale));
    }

    //Should check that this depth map is the correct dimensions
    vil_image_view<double> depth;
    vcl_string depthfile = cfg->get_value<vcl_string>("high_res_depth");
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

    vcl_cout << depth.ni() << " " << depth.nj() << "\n";
    double min, max;
    vil_math_value_range(depth, min, max);
    vcl_cout << min << " " << max << "\n";

    int i0, ni, j0, nj;
    vil_image_view<double> ref_image;

    //Crop regions are defined in the coordinates of the LOW RES reference frame
    vpgl_perspective_camera<double> ref_cam = scaled_cams[ref_frame];
    if (cfg->is_set("crop_window"))
    {
      vcl_istringstream cwstream(cfg->get_value<vcl_string>("crop_window"));
      cwstream >> i0 >> ni >> j0 >> nj;
      vcl_cout << frames[ref_frame].ni() << " " << frames[ref_frame].nj() << "\n";
      ref_image.deep_copy(vil_crop(frames[ref_frame], i0, ni, j0, nj));
      i0 = (int)(i0*scale_factor/camera_scale);
      j0 = (int)(j0*scale_factor/camera_scale);
      ni = (int)(ni*scale_factor/camera_scale);
      nj = (int)(nj*scale_factor/camera_scale);
      vcl_cout << "Crop window: " << i0 << " " << ni << " " << j0 << " " << nj << "\n";
      if (cfg->is_set("ground_truth"))
        gt = vil_crop(gt, i0, ni, j0, nj);

      if (cfg->is_set("crop_depth") && cfg->get_value<bool>("crop_depth"))
        depth = vil_crop(depth, i0, ni, j0, nj);
      ref_cam = super3d::crop_camera(ref_cam, i0, j0);
    }
    else
    {
      i0 = j0 = 0;
      ni = frames[ref_frame].ni();
      nj = frames[ref_frame].nj();
      ref_image.deep_copy(frames[ref_frame]);
    }

    vcl_vector<vil_image_view<double> > weights(cameras.size());

    if (cfg->get_value<bool>("use_normal_weighting"))
    {
      vil_image_view<double> location_map, normal_map, ref_angle_map;
      super3d::depth_map_to_normal_map_inv_len(ref_cam, depth, normal_map);
      super3d::depth_map_to_location_map(ref_cam, depth, location_map);
      super3d::viewing_angle_map(ref_cam.camera_center(), location_map, normal_map, ref_angle_map);

#ifdef DEBUG
      vil_image_view<double> rotated;
      vil_image_view<vxl_byte> byte_normals;
      rotate_normal_map(normal_map, ref_cam.get_rotation(), rotated);
      byte_normal_map(rotated, byte_normals);
      vil_save(byte_normals, "images/normal_map.png");
#endif

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

    vcl_cout << "Computing flow from depth\n";
    vcl_vector<vil_image_view<double> > flows;
    super3d::compute_occluded_flows_from_depth(scaled_cams, ref_cam, depth, flows);
    crop_frames_and_flows(flows, frames, scale_factor, 3);
    vcl_vector<vidtk::adjoint_image_ops_func<double> > warps;
    create_warps_from_flows(flows, frames, weights, warps, scale_factor, cfg.get());

 #ifdef DEBUG
    for (unsigned int i = 0; i < warps.size(); i++)
    {
      // test if the DBW operator for image i is adjoint
      vcl_cout << "is adjoint "<<i<<" "<<warps[i].is_adjoint()<<vcl_endl;

      char buf[50];
      vil_image_view<vxl_byte> output;

      // apply DBW to ground truth super res image to predict frame i
      vil_image_view<double> gts, pred;
      vil_convert_stretch_range_limited(gt, gts, 0.0, 255.0, 0.0, 1.0);
      warps[i].apply_A(gts, pred);
      vil_convert_stretch_range_limited(pred, output, 0.0, 1.0);
      sprintf(buf, "images/predicted-%03d.png", i);
      vil_save(output, buf);

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
      vil_math_image_difference(wframe, pred, pred);
      vil_convert_stretch_range_limited(pred, output, -1.0, 1.0);
      sprintf(buf, "images/predition-error%03d.png", i);
      vil_save(output, buf);

      // write out the viewing angle based weights
      vil_convert_stretch_range_limited(weights[i], output, 0.0, 2.0);
      sprintf(buf, "images/angle-weight-%03d.png", i);
      vil_save(output, buf);
    }
#endif

    vcl_cout << "Computing super resolution\n";

    //Initilize super resolution parameters
    vil_image_view<double> super_u;
    super3d::super_res_params srp;
    srp.scale_factor = scale_factor;
    srp.ref_frame = ref_frame;
    srp.s_ni = warps[ref_frame].src_ni();
    srp.s_nj = warps[ref_frame].src_nj();
    srp.l_ni = warps[ref_frame].dst_ni();
    srp.l_nj = warps[ref_frame].dst_nj();

    srp.lambda = cfg->get_value<double>("lambda");
    srp.epsilon_data = cfg->get_value<double>("epsilon_data");
    srp.epsilon_reg = cfg->get_value<double>("epsilon_reg");
    srp.sigma = cfg->get_value<double>("sigma");
    srp.tau = cfg->get_value<double>("tau");
    const unsigned int iterations = cfg->get_value<unsigned int>("iterations");
    vcl_string output_image = cfg->get_value<vcl_string>("output_image");
//    super3d::super_resolve(frames, warps, super_u, srp, iterations, output_image);

    // additional parameters
    vcl_string tv_method_str = cfg::inst()->get_value<vcl_string>("tv_method");
    if( tv_method_str.compare("SUPER3D_BASELINE") == 0 )
    {
      srp.tv_method = super_res_params::SUPER3D_BASELINE;
    }
    else if ( tv_method_str.compare("IMAGEDATA_IMAGEPRIOR") == 0 )
    {
      srp.tv_method = super_res_params::IMAGEDATA_IMAGEPRIOR;
    }
    else if ( tv_method_str.compare("GRADIENTDATA_IMAGEPRIOR") == 0 )
    {
      srp.tv_method = super_res_params::GRADIENTDATA_IMAGEPRIOR;
    }
    else if ( tv_method_str.compare("IMAGEDATA_IMAGEPRIOR_ILLUMINATIONPRIOR") == 0 )
    {
      srp.tv_method = super_res_params::IMAGEDATA_IMAGEPRIOR_ILLUMINATIONPRIOR;
    }
    else if ( tv_method_str.compare("IMAGEDATA_GRADIENTDATA_IMAGEPRIOR_ILLUMINATIONPRIOR") == 0 )
    {
      srp.tv_method = super_res_params::IMAGEDATA_GRADIENTDATA_IMAGEPRIOR_ILLUMINATIONPRIOR;
    }
    else
    {
      vcl_cerr << "unknown tv method\n";
      return(-1);
    }
    vcl_cout << "tv method : " << tv_method_str << vcl_endl;

    vcl_string str = cfg::inst()->get_value<vcl_string>("cost_function");
    if( str.compare("HUMBER_NORM") == 0 )
    {
      srp.cost_function = super_res_params::HUBER_NORM;
    }
    else if( str.compare("TRUNCATED_QUADRATIC") == 0 )
    {
      srp.cost_function = super_res_params::TRUNCATED_QUADRATIC;
    }
    else if( str.compare("GENERALIZED_HUBER") == 0 )
    {
      srp.cost_function = super_res_params::GENERALIZED_HUBER;
    }
    else
    {
      vcl_cerr << "unknown cost function\n";
      return(-1);
    }
    vcl_cout << "cost function : " << str << vcl_endl;

    srp.alpha_a = cfg::inst()->get_value<double>("alpha_a");
    srp.gamma_a = cfg::inst()->get_value<double>("gamma_a");
    srp.beta_a = cfg::inst()->get_value<double>("beta_a");
    srp.sigma_a = cfg::inst()->get_value<double>("sigma_a");

    srp.lambda_g = cfg::inst()->get_value<double>("lambda_g");
    srp.alpha_g = cfg::inst()->get_value<double>("alpha_g");
    srp.gamma_g = cfg::inst()->get_value<double>("gamma_g");
    srp.beta_g = cfg::inst()->get_value<double>("beta_g");
    srp.sigma_g = cfg::inst()->get_value<double>("sigma_g");

    srp.lambda_r = cfg::inst()->get_value<double>("lambda_r");
    srp.alpha_r = cfg::inst()->get_value<double>("alpha_r");
    srp.gamma_r = cfg::inst()->get_value<double>("gamma_r");
    srp.beta_r = cfg::inst()->get_value<double>("beta_r");
    srp.sigma_r = cfg::inst()->get_value<double>("sigma_r");

    srp.lambda_l = cfg::inst()->get_value<double>("lambda_l");
    srp.alpha_l = cfg::inst()->get_value<double>("alpha_l");
    srp.gamma_l = cfg::inst()->get_value<double>("gamma_l");
    srp.beta_l = cfg::inst()->get_value<double>("beta_l");
    srp.sigma_l = cfg::inst()->get_value<double>("sigma_l");

    srp.sigma_pr = cfg::inst()->get_value<double>("sigma_pr");
    srp.sigma_pl = cfg::inst()->get_value<double>("sigma_pl");
    srp.sigma_qa = cfg::inst()->get_value<double>("sigma_qa");
    srp.sigma_qg = cfg::inst()->get_value<double>("sigma_qg");
    srp.sigma_A = cfg::inst()->get_value<double>("sigma_A");
    srp.sigma_Y = cfg::inst()->get_value<double>("sigma_Y");

    vcl_vector< vil_image_view<double> > As;
    super3d::super_resolve_robust(frames, warps, super_u, srp, iterations, As);

    if (cfg::inst()->is_set("ground_truth"))
    {
      vil_math_scale_values(gt, normalizer);
      super3d::compare_to_original(ref_image, super_u, gt, scale_factor);
    }
    else
    {
      vil_image_view<double> upsamp;
      upsample(ref_image, upsamp, scale_factor, vidtk::warp_image_parameters::CUBIC);

      vil_image_view<vxl_byte> output;
      vil_convert_stretch_range_limited(upsamp, output, 0.0, 1.0);
      vil_save(output, "bicub.png");
    }

    vil_image_view<vxl_byte> output;
    vil_convert_stretch_range_limited(super_u, output, 0.0, 1.0);
    vil_save(output, output_image.c_str());
  }
  catch (const super3d::cfg::cfg_exception &e)  {
    vcl_cout << "Error in config: " << e.what() << "\n";
  }

  return 0;
}

//*****************************************************************************

/// Use the valid region of the flow destinations to crop the frames and translate the flow.
void crop_frames_and_flows(vcl_vector<vil_image_view<double> > &flows,
                           vcl_vector<vil_image_view<double> > &frames,
                           int scale_factor, int margin)
{
  assert(flows.size() == frames.size());
  for (unsigned int i=0; i<flows.size(); ++i)
  {
    // get a bounding box around the flow and expand by a margin
    vgl_box_2d<double> bounds = super3d::flow_destination_bounds(flows[i]);
    vgl_box_2d<int> bbox(vcl_floor(bounds.min_x()), vcl_ceil(bounds.max_x()),
                         vcl_floor(bounds.min_y()), vcl_ceil(bounds.max_y()));
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

void create_warps_from_flows(const vcl_vector<vil_image_view<double> > &flows,
                             const vcl_vector<vil_image_view<double> > &frames,
                             const vcl_vector<vil_image_view<double> > &weights,
                             vcl_vector<vidtk::adjoint_image_ops_func<double> > &warps,
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
    vidtk::warp_back_with_flow_bicub(I1_x, flow, temp);
  }
  else
  {
    vidtk::warp_back_with_flow_bilin(I1_x, flow, temp);
  }
  vil_math_image_difference(temp, I0_x, diff);
}

//*****************************************************************************

void create_low_res(vcl_vector<vil_image_view<double> > &frames,
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
      vidtk::down_scale(temp, frames[i], scale);
    }
    else
    {
      vidtk::down_sample(temp, frames[i], scale);
    }
  }
}

//*****************************************************************************

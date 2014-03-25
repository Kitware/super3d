/*
 * Copyright 2013 Kitware, Inc.
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
#include "file_io.h"
#include "multiscale.h"
#include "super_res.h"
#include "flow_manip.h"
#include "depth_map.h"

#include "depth_cl/super_res.h"

#include <video_transforms/adjoint_flow_warp.h>
#include <video_transforms/adjoint_resample.h>
#include <video_transforms/adjoint_dbw.h>

#include <vcl_iostream.h>
#include <boost/bind.hpp>

#include <vil/vil_crop.h>
#include <vil/vil_bicub_interp.h>
#include <vil/vil_load.h>
#include <vil/vil_save.h>
#include <vil/vil_convert.h>
#include <vil/vil_math.h>
#include <vil/vil_resample_bilin.h>
#include <vil/vil_resample_bicub.h>
#include <vnl/vnl_inverse.h>

#include <vgl/algo/vgl_h_matrix_2d.h>
#include <vgl/vgl_homg_point_2d.h>
#include <vgl/vgl_intersection.h>

#include <video_transforms/warp_image.h>
#include <video_transforms/warp_and_average.h>
#include <tracking/refine_homography.h>

//#define GPU

void load_flow(const char *flow_list, const vcl_string &dir, vcl_vector<vil_image_view<double> > &flows);
void create_warps_from_flows(const vcl_vector<vil_image_view<double> > &flows,
                             const vcl_vector<vil_image_view<double> > &frames,
                             vcl_vector<vidtk::adjoint_image_ops_func<double> > &warps,
                             int scale_factor);
void load_homogs(const char *homog_list, vcl_vector<vgl_h_matrix_2d<double> > &homogs,
                 const vcl_vector<int> &frameindex);
void refine_homogs(vcl_vector<vgl_h_matrix_2d<double> > &homogs,
                   const vcl_vector<vil_image_view<double> > &frames,
                   unsigned int ref_frame,
                   int i0, int ni, int j0, int nj,
                   unsigned int margin);
void homogs_to_flows(const vcl_vector<vgl_h_matrix_2d<double> > &homogs,
                     int ref_frame, const vcl_vector<vil_image_view<double> > &frames,
                     int i0, int j0, int ni, int nj,
                     int scale_factor, vcl_vector<vil_image_view<double> > &flows);
void crop_frames_and_flows(vcl_vector<vil_image_view<double> > &flows,
                           vcl_vector<vil_image_view<double> > &frames,
                           int scale_factor, int margin);
void upsample(const vil_image_view<double> &src, vil_image_view<double> &dest,
              double scale_factor, vidtk::warp_image_parameters::interp_type interp);
void write_warped_frames(const vcl_vector<vil_image_view<double> > &frames,
                         const vcl_vector<vgl_h_matrix_2d<double> > &homogs,
                         unsigned int i0, unsigned int j0, unsigned int ni, unsigned int nj,
                         unsigned int ref_frame);
///Converts optical flow 2 plane image to an RGB image
void flow_to_colormap(const vil_image_view<double> &flow,
                      const vil_image_view<double> &ref_flow,
                      vil_image_view<vil_rgb<vxl_byte> > &colormap_flow,
                      double max_flow = 0.0);

void create_low_res(vcl_vector<vil_image_view<double> > &frames,
                    int scale);

//*****************************************************************************

int main(int argc, char* argv[])
{
  try  {
    config::inst()->read_config(argv[1]);
    config::inst()->read_argument_updates(argc, argv);

    vcl_vector<vil_image_view<double> > frames;
    vcl_vector<vpgl_perspective_camera<double> >  cameras;
    vcl_vector<vcl_string> filenames;
    vcl_vector<int> frameindex;
    vcl_string frame_file = config::inst()->get_value<vcl_string>("frame_list");
    vcl_string dir = config::inst()->get_value<vcl_string>("directory");
    vcl_cout << "Using frame file: " << frame_file << " to find images and flows.\n";

    const unsigned int ref_frame = config::inst()->get_value<unsigned int>("ref_frame");
    const double scale_factor = config::inst()->get_value<double>("scale_factor");

    load_from_frame_file(frame_file.c_str(), dir, filenames, frameindex, frames,
                         config::inst()->get_value<bool>("use_color"), config::inst()->get_value<bool>("use_rgb12"));

    int i0, ni, j0, nj;
    vil_image_view<double> ref_image;

    if (config::inst()->is_set("crop_window"))
    {
      vcl_istringstream cwstream(config::inst()->get_value<vcl_string>("crop_window"));
      cwstream >> i0 >> ni >> j0 >> nj;
    }
    else
    {
      i0 = j0 = 0;
      ni = frames[ref_frame].ni();
      nj = frames[ref_frame].nj();
    }

    vil_image_view<double> original, waa;
    original.deep_copy(frames[ref_frame]);
    if (config::inst()->is_set("create_low_res") && config::inst()->get_value<bool>("create_low_res"))
    {
      vcl_cout << "Creating low resolution data\n";
      original = vil_crop(original, i0, ni, j0, nj);
      vil_image_view<unsigned short> out;
      vil_convert_stretch_range_limited(original, out, 0.0, 255.0, 0, 65535);
      vil_save(out, "original.png");
      create_low_res(frames, scale_factor);
    }

    ref_image.deep_copy(frames[ref_frame]);

#ifndef GPU
    vcl_vector<vidtk::adjoint_image_ops_func<double> > warps;
#endif

    vcl_vector<vil_image_view<double> > flows;

    if (config::inst()->is_set("flow_file"))
    {
      load_flow( config::inst()->get_value<vcl_string>("flow_file").c_str(), dir, flows);
#ifndef GPU
      create_warps_from_flows(flows, frames, warps, scale_factor);
#endif
    }
    else if (config::inst()->is_set("homog_file"))
    {
      vcl_vector<vgl_h_matrix_2d<double> > homogs;
      if (config::inst()->is_set("crop_window"))
      {
        vcl_cout << frames[ref_frame].ni() << " " << frames[ref_frame].nj() << "\n";
        if (config::inst()->get_value<bool>("create_low_res"))
        {
          i0 /= scale_factor;
          ni /= scale_factor;
          j0 /= scale_factor;
          nj /= scale_factor;
        }
        load_homogs(config::inst()->get_value<vcl_string>("homog_file").c_str(), homogs, frameindex);
        refine_homogs(homogs, frames, ref_frame, i0, ni, j0, nj, 50);
        ref_image.deep_copy(vil_crop(frames[ref_frame], i0, ni, j0, nj));
      }
      else
      {
        load_homogs(config::inst()->get_value<vcl_string>("homog_file").c_str(), homogs, frameindex);
        //refine_homogs(homogs, frames, ref_frame, i0, ni, j0, nj, 100);
      }

      vidtk::warp_image_parameters wip;
      wip.set_fill_unmapped(true);
      wip.set_unmapped_value(-1.0);
      wip.set_interpolator(vidtk::warp_image_parameters::CUBIC);
      vidtk::warp_and_average<double>(frames, waa, homogs, ref_frame, i0, j0, ni, nj, wip, scale_factor);
      vil_image_view<unsigned short> out;
      vil_convert_stretch_range_limited(waa, out, 0.0, 255.0, 0, 65535);
      vil_save(out, "mfavg.png");

      //Should be called on the full images to avoid black borders
      write_warped_frames(frames, homogs, i0, j0, ni, nj, ref_frame);

      homogs_to_flows(homogs, ref_frame, frames, i0, j0, ni, nj, scale_factor, flows);
      crop_frames_and_flows(flows, frames, scale_factor, 0);

      //for (unsigned int i = 0; i < flows.size(); i++)
      //{
      //  vcl_cout << flows[i].ni() << " " << flows[i].nj() << "\n";
      //  vil_image_view<vil_rgb<vxl_byte> > colormap_flow;
      //  flow_to_colormap(flows[i], flows[ref_frame], colormap_flow, 3);
      //  char buf[40];
      //  sprintf(buf, "color%d.png", i);
      //  vil_save(colormap_flow, buf);
      //  sprintf(buf, "img%d.png", i);
      //  vil_image_view<vxl_byte> img_b;
      //  vil_convert_cast(frames[i], img_b);
      //  vil_save(img_b, buf);
      //}

#ifndef GPU
      create_warps_from_flows(flows, frames, warps, scale_factor);
#endif
    }
    else
    {
      vcl_cerr << "Must specify flow file or homog file.\n";
      return 1;
    }

    if (!config::inst()->get_value<bool>("use_rgb12"))
    {
      double normalizer = 1.0/255.0;
      for (unsigned int i = 0; i < frames.size(); i++)
      {
        vil_math_scale_values(frames[i], normalizer);
      }
    }

    //Initilize super resolution parameters
    vil_image_view<double> super_u;
#ifdef GPU
    super_res_cl::params srp;
    srp.sdim.s[0] = scale_factor * frames[ref_frame].ni();
    srp.sdim.s[1] = scale_factor * frames[ref_frame].nj();
    srp.ldim.s[0] = frames[ref_frame].ni();
    srp.ldim.s[1] = frames[ref_frame].nj();
#else
    super_res_params srp;
    srp.s_ni = warps[ref_frame].src_ni();
    srp.s_nj = warps[ref_frame].src_nj();
    srp.l_ni = warps[ref_frame].dst_ni();
    srp.l_nj = warps[ref_frame].dst_nj();
    srp.ref_frame = ref_frame;
#endif

    srp.scale_factor = scale_factor;
    srp.lambda = config::inst()->get_value<double>("lambda");
    srp.epsilon_data = config::inst()->get_value<double>("epsilon_data");
    srp.epsilon_reg = config::inst()->get_value<double>("epsilon_reg");
    srp.sigma = config::inst()->get_value<double>("sigma");
    srp.tau = config::inst()->get_value<double>("tau");
    const unsigned int iterations = config::inst()->get_value<unsigned int>("iterations");
#ifdef GPU
    {
    super_res_cl srcl;
    vil_image_view<float> super_u_flt;
    vcl_vector<vil_image_view<float> > frames_flt, flows_flt;
    frames_flt.resize(frames.size());
    flows_flt.resize(frames.size());
    for (unsigned int i = 0; i < frames.size(); i++)
    {
      vil_convert_cast(frames[i], frames_flt[i]);
      vil_convert_cast(flows[i], flows_flt[i]);
    }

    vcl_cout << "Starting GPU Super Res.\n";
    srcl.super_resolve(frames_flt, flows_flt, super_u_flt, srp, iterations);
    vil_convert_cast(super_u_flt, super_u);
    vcl_cout << "finished GPU!\n";
    return 0;
    }
#else
    super_resolve(frames, warps, super_u, srp, iterations);
#endif
    vil_image_view<double> upsamp;
    upsample(ref_image, upsamp, scale_factor, vidtk::warp_image_parameters::CUBIC);
    vil_image_view<unsigned short> output;
    vil_convert_stretch_range_limited(upsamp, output, 0.0, 255.0, 0, 65535);
    vil_save(output, "bicub.png");

    if (config::inst()->get_value<bool>("create_low_res"))
    {
      vil_math_scale_values(original, 1.0/255.0);
      vil_math_scale_values(upsamp, 1.0/255.0);
      vil_math_scale_values(waa, 1.0/255.0);
      vcl_cout << "ssd bicub: " << vil_math_ssd(upsamp, original, double()) << "\n";
      vcl_cout << "ssd mfavg: " << vil_math_ssd(waa, original, double()) << "\n";
      vcl_cout << "ssd super: " << vil_math_ssd(super_u, original, double()) << "\n";
    }

    vil_convert_stretch_range_limited(super_u, output, 0.0, 1.0, 0, 65535);
    vil_save(output, config::inst()->get_value<vcl_string>("output_image").c_str());
  }
  catch (const config::cfg_exception &e)  {
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
    vgl_box_2d<double> bounds = flow_destination_bounds(flows[i]);
    vgl_box_2d<int> bbox(vcl_floor(bounds.min_x()), vcl_ceil(bounds.max_x()),
                         vcl_floor(bounds.min_y()), vcl_ceil(bounds.max_y()));
    bbox.expand_about_centroid(2*margin);

    // get a bounding box around the (up sampled) image and intersect
    vgl_box_2d<int> img_box(0, frames[i].ni()-1, 0, frames[i].nj()-1);
    img_box.scale_about_origin(scale_factor);
    bbox = vgl_intersection(bbox, img_box);
    bbox.scale_about_origin(1.0 / scale_factor);

    // crop the image and translated the flow accordingly
    frames[i].deep_copy(vil_crop(frames[i], bbox.min_x(), bbox.width(), bbox.min_y(), bbox.height()));
    bbox.scale_about_origin(scale_factor);
    translate_flow(flows[i], -bbox.min_x(), -bbox.min_y());
  }
}

//*****************************************************************************

void create_warps_from_flows(const vcl_vector<vil_image_view<double> > &flows,
                             const vcl_vector<vil_image_view<double> > &frames,
                             vcl_vector<vidtk::adjoint_image_ops_func<double> > &warps,
                             int scale_factor)
{
  assert(flows.size() == frames.size());
  bool down_sample_averaging = config::inst()->get_value<bool>("down_sample_averaging");
  bool bicubic_warping = config::inst()->get_value<bool>("bicubic_warping");
  double smoothing_sigma = config::inst()->get_value<double>("sensor_sigma");

  warps.clear();
  warps.reserve(flows.size());
  for (unsigned int i = 0; i < flows.size(); ++i)
  {
    warps.push_back(vidtk::create_dbw_from_flow(flows[i],
                                         frames[i].ni(),
                                         frames[i].nj(),
                                         frames[i].nplanes(),
                                         scale_factor,
                                         smoothing_sigma,
                                         down_sample_averaging,
                                         bicubic_warping));
  }
}

//*****************************************************************************

void load_flow(const char *flow_list, const vcl_string &dir, vcl_vector<vil_image_view<double> > &flows)
{
  vcl_ifstream infile(flow_list);
  vcl_string flowname;

  while (infile >> flowname)
  {
    vil_image_view<double> flowimg;
    read_flow_file(flowimg, (dir + flowname).c_str());
    flows.push_back(flowimg);
  }
  vcl_cout << "Read " << flows.size() << " flows.\n";
}

//*****************************************************************************

void load_homogs(const char *homog_list, vcl_vector<vgl_h_matrix_2d<double> > &homogs,
                 const vcl_vector<int> &frameindex)
{
  vcl_ifstream infile(homog_list);
  vgl_h_matrix_2d<double> H;
  int index = 0; int framenum = 0;
  while (infile >> H)
  {
    if (frameindex[index] == framenum)
    {
      if (config::inst()->get_value<bool>("create_low_res"))
      {
        double scale_factor = config::inst()->get_value<double>("scale_factor");
        vnl_double_3x3 S, Sinv;
        S.set_identity();
        S(0,0) = scale_factor;
        S(1,1) = scale_factor;
        Sinv.set_identity();
        Sinv(0,0) = 1.0/scale_factor;
        Sinv(1,1) = 1.0/scale_factor;
        H = Sinv * H.get_matrix() * S;
      }

      homogs.push_back(H);
      index++;
      if (index == frameindex.size())
        break;
    }

    framenum++;
  }
}

//*****************************************************************************

void refine_homogs(vcl_vector<vgl_h_matrix_2d<double> > &homogs,
                   const vcl_vector<vil_image_view<double> > &frames,
                   unsigned int ref_frame,
                   int i0, int ni, int j0, int nj,
                   unsigned int margin)
{
  vcl_cout << "Refining homog ";
  vil_image_view<double> refimg;

  if (frames[ref_frame].nplanes() > 1)
    vil_convert_planes_to_grey(frames[ref_frame], refimg);
  else
    refimg = frames[ref_frame];

  refimg = vil_crop(refimg, i0 - margin, ni + 2*margin, j0 - margin, nj + 2*margin);
  vil_image_view<vxl_byte> output;
  vil_convert_cast(refimg, output);
  vil_save(output, "refhomg/test.png");

  vnl_double_3x3 ref_H = homogs[ref_frame].get_matrix();
  vnl_double_3x3 ref_H_inv = homogs[ref_frame].get_inverse().get_matrix();
  for (unsigned int i = 0; i < homogs.size(); i++)
  {
    vcl_cout << i << " ";
    vil_image_view<double> img;

    if (frames[i].nplanes() > 1)
      vil_convert_planes_to_grey(frames[i], img);
    else
      img = frames[i];

    vnl_double_3x3 H = ref_H_inv * homogs[i].get_matrix();
    vidtk::crop_homography_source(H, i0 - margin, j0 - margin);
    vidtk::refine_homography(refimg, img, H, 4, 100, 15);
    vidtk::crop_homography_source(H, -i0 + margin, -j0 + margin);

    vil_image_view<double> warped(ni+2*margin, nj+2*margin);
    vnl_double_3x3 Hinv = vnl_inverse<double>(H);
    vidtk::warp_image(img, warped, Hinv);
    vil_convert_cast(warped, output);
    char buf[50];
    sprintf(buf, "refhomgwand/ref%d.png", i);
    vil_save(output, buf);

    //vidtk::crop_homography(H, margin, margin);
    homogs[i].set(H);
  }
  //exit(1);
  vcl_cout << "\n";
}

//*****************************************************************************

void homogs_to_flows(const vcl_vector<vgl_h_matrix_2d<double> > &homogs,
                     int ref_frame, const vcl_vector<vil_image_view<double> > &frames,
                     int i0, int j0, int ni, int nj,
                     int scale_factor, vcl_vector<vil_image_view<double> > &flows)
{
  vcl_string flowname;

  const vgl_h_matrix_2d<double> &ref_H = homogs[ref_frame];

  i0 *= scale_factor;  j0 *= scale_factor;
  ni *= scale_factor;  nj *= scale_factor;

  vnl_double_3x3 S, Sinv;
  S.set_identity();
  S(0,0) = scale_factor;
  S(1,1) = scale_factor;
  Sinv.set_identity();
  Sinv(0,0) = 1.0/scale_factor;
  Sinv(1,1) = 1.0/scale_factor;

  for (unsigned int f = 0; f < homogs.size(); f++)
  {
    vgl_h_matrix_2d<double> ref_to_f = S * homogs[f].get_inverse().get_matrix() * ref_H.get_matrix() * Sinv;

    vil_image_view<double> flow;
    flow.set_size(ni, nj, 2);
    flow.fill(vcl_numeric_limits<double>::quiet_NaN());
    for (int i = 0; i < ni; i++)
    {
      for (int j = 0; j < nj; j++)
      {
        vgl_point_2d<double> p = ref_to_f * vgl_homg_point_2d<double>(i+i0, j+j0);
        flow(i, j, 0) = p.x() - i;
        flow(i, j, 1) = p.y() - j;
      }
    }

    flows.push_back(flow);
  }

  vcl_cout << "Converted " << flows.size() << " Homogs to flows.\n";
}

//*****************************************************************************

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

void write_warped_frames(const vcl_vector<vil_image_view<double> > &frames,
                         const vcl_vector<vgl_h_matrix_2d<double> > &homogs,
                         unsigned int i0, unsigned int j0, unsigned int ni, unsigned int nj,
                         unsigned int ref_frame)
{
  vnl_double_3x3 T;
  T.set_identity();
  T(0,2) = i0;
  T(1,2) = j0;

  //pre computed post multiply matrix
  const vnl_double_3x3 M = homogs[ref_frame].get_matrix() * T;

  vidtk::warp_image_parameters wip;
  wip.set_fill_unmapped(true);
  wip.set_unmapped_value(-1.0);
  wip.set_interpolator(vidtk::warp_image_parameters::LINEAR);

  char buf[64];
  for (unsigned int i = 0; i < frames.size(); i++)
  {
    vgl_h_matrix_2d<double> H = homogs[i].get_inverse().get_matrix() * M;
    vil_image_view<double> warped(ni, nj, frames[i].nplanes());

    vidtk::warp_image(frames[i], warped, H, wip);
    vil_image_view<vxl_byte> output;
    vil_convert_stretch_range_limited(warped, output, 0.0, 1.0);
    sprintf(buf, "images/warped_frame%02d.png", i);
    vil_save(output, buf);
  }
}


///Converts optical flow 2 plane image to an RGB image
void flow_to_colormap(const vil_image_view<double> &flow,
                      const vil_image_view<double> &ref_flow,
                      vil_image_view<vil_rgb<vxl_byte> > &colormap_flow,
                      double max_flow)
{
  const unsigned ni = flow.ni();
  const unsigned nj = flow.nj();
  //If no max flow distance was specified, then use the max for the flow.
  if (max_flow <= 0.0)
  {
    for (unsigned int j = 0; j < flow.nj(); j++)
    {
      for (unsigned int i = 0; i < flow.ni(); i++)
      {
        double len = vnl_double_2(flow(i, j, 0), flow(i, j, 1)).two_norm();
        if (len > max_flow)
          max_flow = len;
      }
    }
  }

  vcl_cout << "Using max flow: " << max_flow << ".\n";

  colormap_flow.set_size(ni,nj);
  //Convert flow to HSV to RGB
  for (unsigned int j = 0; j < flow.nj(); j++)
  {
    for (unsigned int i = 0; i < flow.ni(); i++)
    {
       vnl_double_2 uv(flow(i, j, 0) - ref_flow(i, j, 0), flow(i, j, 1) - ref_flow(i, j, 1));
       double sat = uv.two_norm() / max_flow;
       if (sat > 1.0)
       {
         sat = 1.0;
       }
       if (sat < 1e-6)
       {
         colormap_flow(i, j) = vil_rgb<vxl_byte>(255, 255, 255);
         continue;
       }
       if (vnl_math::isnan(sat))
       {
         colormap_flow(i, j) = vil_rgb<vxl_byte>(0, 0, 0);
         continue;
       }

       // compute hue as flow angle in degrees
       double hue = (atan2(uv(1), uv(0)) * vnl_math::deg_per_rad);
       // rotate negative angles by 360 degrees to positive angles
       if (hue < 0.0)
         hue += 360.0;
       double hueprime = hue / 60.0;
       double C = sat;
       double X = C * (1.0 - fabs(fmod(hueprime,2.0) - 1.0));
       double m = 1.0 - C;  //Use white for center
       vxl_byte Cmb = static_cast<vxl_byte>((C+m)*255.0);
       vxl_byte Xmb = static_cast<vxl_byte>((X+m)*255.0);
       vxl_byte mb = static_cast<vxl_byte>(m*255.0);
       if (hueprime < 1.0)
         colormap_flow(i, j) = vil_rgb<vxl_byte>(Cmb, Xmb, mb);
       else if (hueprime < 2.0)
         colormap_flow(i, j) = vil_rgb<vxl_byte>(Xmb, Cmb, mb);
       else if (hueprime < 3.0)
         colormap_flow(i, j) = vil_rgb<vxl_byte>(mb, Cmb, Xmb);
       else if (hueprime < 4.0)
         colormap_flow(i, j) = vil_rgb<vxl_byte>(mb, Xmb, Cmb);
       else if (hueprime < 5.0)
         colormap_flow(i, j) = vil_rgb<vxl_byte>(Xmb, mb, Cmb);
       else if (hueprime < 6.0)
         colormap_flow(i, j) = vil_rgb<vxl_byte>(Cmb, mb, Xmb);
    }
  }
}

void create_low_res(vcl_vector<vil_image_view<double> > &frames,
                    int scale)
{
  bool down_sample_averaging = config::inst()->get_value<bool>("down_sample_averaging");

  for (unsigned int i = 0; i < frames.size(); i++)
  {
    vil_image_view<double> temp;
    double sensor_sigma = config::inst()->get_value<double>("sensor_sigma");
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

/*ckwg +29
 * Copyright 2012 by Kitware, Inc.
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
#include <super3d/depth/super_res.h>
#include <super3d/depth/resample.h>

#include <super3d/image/adjoint_dbw.h>
#include <super3d/image/adjoint_flow_warp.h>
#include <cstdio>

#include <boost/bind.hpp>
#include <boost/scoped_ptr.hpp>

#include <vil/vil_load.h>
#include <vil/vil_save.h>
#include <vil/vil_math.h>
#include <vil/vil_convert.h>


#define DEBUG

vcl_vector<super3d::adjoint_image_ops_func<double> >
create_warps_simple(int sni, int snj, int scale_factor, super3d::config *cfg);

vcl_vector<super3d::adjoint_image_ops_func<double> >
create_warps(int sni, int snj, int scale_factor, super3d::config *cfg);

void create_downsampled_frames(const vil_image_view<double> &high_res,
                               const vcl_vector<super3d::adjoint_image_ops_func<double> > &warps,
                               vcl_vector<vil_image_view<double> > &frames);

//*****************************************************************************

int main(int argc, char* argv[])
{
  try  {
    boost::scoped_ptr<super3d::config> cfg(new super3d::config);
    cfg->read_config(argv[1]);
    double scale_factor = cfg->get_value<double>("scale_factor");
    vcl_string img_name = cfg->get_value<vcl_string>("single_frame");
        vcl_cout << "Reading image " << img_name << "\n";
    vil_image_view<double> img;
    vil_image_resource_sptr img_rsc = vil_load_image_resource(img_name.c_str());
    if (img_rsc != NULL)
    {
      vil_image_view<vxl_byte> img_ = img_rsc->get_view();

      if (img_.nplanes() == 3)
        vil_convert_planes_to_grey(img_, img);
      else
        vil_convert_cast(img_, img);
    }
    else
      return 0;

    unsigned int ref_frame = 0;

    vil_math_scale_values(img, 1.0/255.0);

    vcl_cout << "Creating downsampled frames\n";
    vcl_vector<vil_image_view<double> > frames;

    vcl_vector<super3d::adjoint_image_ops_func<double> > warps;
    warps = create_warps_simple(img.ni(), img.nj(), scale_factor, cfg.get());
    create_downsampled_frames(img, warps, frames);

#ifdef DEBUG
    for (unsigned int i = 0; i < frames.size(); i++)
    {
      vil_image_view<vxl_byte> output;
      vil_convert_stretch_range_limited(frames[i], output, 0.0, 1.0);
      char buf[50];
      sprintf(buf, "images/ds%d.png", i);
      vil_save(output, buf);
    }
#endif

    vcl_cout << "Computing super resolution\n";
    vil_image_view<double> super_u;
    super3d::super_res_params srp;
    srp.scale_factor = scale_factor;
    srp.ref_frame = ref_frame;
    srp.s_ni = img.ni();
    srp.s_nj = img.nj();
    srp.l_ni = static_cast<unsigned int>(srp.s_ni/scale_factor);
    srp.l_nj = static_cast<unsigned int>(srp.s_nj/scale_factor);
    srp.lambda = cfg->get_value<double>("lambda");
    srp.epsilon_data = cfg->get_value<double>("epsilon_data");
    srp.epsilon_reg = cfg->get_value<double>("epsilon_reg");
    srp.sigma = cfg->get_value<double>("sigma");
    srp.tau = cfg->get_value<double>("tau");
    const unsigned int iterations = cfg->get_value<unsigned int>("iterations");
    vcl_string output_image = cfg->get_value<vcl_string>("output_image");
    super_resolve(frames, warps, super_u, srp, iterations, output_image);

    super3d::compare_to_original(frames[ref_frame], super_u, img, scale_factor);

    vil_image_view<vxl_byte> output;
    vil_convert_stretch_range_limited(super_u, output, 0.0, 1.0);
    vil_save(output, output_image.c_str());
  }
  catch (const super3d::config::cfg_exception &e)  {
    vcl_cout << "Error in config: " << e.what() << "\n";
  }

  return 0;
}

//*****************************************************************************

vcl_vector<super3d::adjoint_image_ops_func<double> >
create_warps_simple(int sni, int snj, int scale_factor, super3d::config *cfg)
{
  //int num_imgs = scale_factor * scale_factor;
  int dni = sni / scale_factor;
  int dnj = snj / scale_factor;
  bool down_sample_averaging = cfg->get_value<bool>("down_sample_averaging");
  bool bicubic_warping = cfg->get_value<bool>("bicubic_warping");
  double sensor_sigma = cfg->get_value<double>("sensor_sigma");

  typedef super3d::adjoint_image_ops_func<double>::func_t func_t;
  vcl_vector<super3d::adjoint_image_ops_func<double> > warps;

  for (unsigned int i = 0; i < static_cast<unsigned int>(scale_factor); i++)
  {
    for (unsigned int j = 0; j < static_cast<unsigned int>(scale_factor); j++)
    {
      vil_image_view<double> flow(sni, snj, 2);
      vil_plane<double>(flow, 0).fill(static_cast<double>(i));
      vil_plane<double>(flow, 1).fill(static_cast<double>(j));

      warps.push_back(super3d::create_dbw_from_flow(flow, dni, dnj, 1, scale_factor, sensor_sigma,
                                           down_sample_averaging,
                                           bicubic_warping));
    }
  }
  return warps;
}

//*****************************************************************************

vcl_vector<super3d::adjoint_image_ops_func<double> >
create_warps(int sni, int snj, int scale_factor, super3d::config *cfg)
{
  int dni = sni / scale_factor;
  int dnj = snj / scale_factor;
  bool down_sample_averaging = cfg->get_value<bool>("down_sample_averaging");
  bool bicubic_warping = cfg->get_value<bool>("bicubic_warping");
  double sensor_sigma = cfg->get_value<double>("sensor_sigma");

  typedef super3d::adjoint_image_ops_func<double>::func_t func_t;
  vcl_vector<super3d::adjoint_image_ops_func<double> > warps;
  for (double i = 0.0; i < scale_factor; i+=0.5)
  {
    for (double j = 0.0; j < scale_factor; j+=0.5)
    {
      vil_image_view<double> flow(sni, snj, 2);
      vil_plane<double>(flow, 0).fill(i);
      vil_plane<double>(flow, 1).fill(j);

      warps.push_back(super3d::create_dbw_from_flow(flow, dni, dnj, 1, scale_factor, sensor_sigma,
                                           down_sample_averaging,
                                           bicubic_warping));
    }
  }
  return warps;
}

//*****************************************************************************


void create_downsampled_frames(const vil_image_view<double> &high_res,
                               const vcl_vector<super3d::adjoint_image_ops_func<double> > &warps,
                               vcl_vector<vil_image_view<double> > &frames)
{
  frames.clear();
  for (unsigned i=0; i<warps.size(); ++i)
  {
    vil_image_view<double> frame(warps[i].dst_ni(),
                                 warps[i].dst_nj(),
                                 warps[i].dst_nplanes());
    warps[i].apply_A(high_res, frame);
    frames.push_back(frame);
  }
}

//*****************************************************************************

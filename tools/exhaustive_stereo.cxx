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

#include <super3d/depth/tv_refine_search.h>
#include <super3d/depth/cost_volume.h>

// VXL includes
#include <vil/vil_image_view.h>
#include <vil/vil_convert.h>
#include <vil/vil_save.h>
#include <vil/vil_bilin_interp.h>
#include <vil/vil_crop.h>
#include <vil/algo/vil_sobel_3x3.h>
#include <vil/vil_save.h>
#include <vil/vil_load.h>
#include <vul/vul_arg.h>
#include <vnl/vnl_double_3.h>
#include <vnl/vnl_inverse.h>

#include <vil/vil_decimate.h>
#include <vil/vil_resample_bilin.h>

// vcl includes
#include <vcl_iostream.h>
#include <vcl_string.h>


//C:/Data/ground_medians/cameras.txt C:/Data/ground_medians/full_window15_reset50/frame###.png,000:100:950 c:/Data/meshes/crop.obj -f 0  -e C:/Data/exposure/exposure.txt -cw "767 1320 757 954"
int main(int argc, char* argv[])
{
  vul_arg<vcl_string> image1_file( 0, "image 1 name", "" );
  vul_arg<vcl_string> image2_file( 0, "image 2 name", "" );
  vul_arg<vcl_string> output_file( 0, "output disparity map (image file)", "");

  vul_arg_parse( argc, argv );

  vil_image_view<double> img1, img2;

  //Read Images
  {
    vil_image_resource_sptr img1_rsc = vil_load_image_resource(image1_file().c_str());
    if (img1_rsc == NULL)
    {
      vcl_cerr << "cannot load image 1\n";
      return 1;
    }

    vil_image_resource_sptr img2_rsc = vil_load_image_resource(image2_file().c_str());
    if (img2_rsc == NULL)
    {
      vcl_cerr << "cannot load image 2\n";
      return 1;
    }

    vil_convert_planes_to_grey<vxl_byte, double>(img1_rsc->get_view(), img1);
    vil_convert_planes_to_grey<vxl_byte, double>(img2_rsc->get_view(), img2);
    vil_math_scale_and_offset_values(img1, 1.0/255.0, 0);
    vil_math_scale_and_offset_values(img2, 1.0/255.0, 0);
  }

  double idepth_min = 0.25;
  double idepth_max = 16;

  vcl_cout << "Initializing depth map"<<vcl_endl;
  vil_image_view<double> depth(img1.ni(), img1.nj(), 1);

  vcl_cout << "Refining depth"<<vcl_endl;
  unsigned int S = 64;
  double theta0 = 100.0;
  double beta_end = 1e-7;
  double beta = .1;
  double lambda = 100;
  double epsilon = 0.01;

  vil_image_view<double> g;
  vil_image_view<double> cost_volume;

#if 1
  vcl_vector<vil_image_view<double> > frames;
  frames.push_back(img1);
  frames.push_back(img2);
  super3d::compute_cost_volume_rectified(frames, 0, S, idepth_min, idepth_max, cost_volume);

  g.set_size(depth.ni(), depth.nj(), 1);

  vil_image_view<double> img1_g;
  vil_sobel_3x3(img1, img1_g);

  const double alpha = 5.0;
  vcl_cout << "Computing g weighting.\n";
  for (unsigned int i = 0; i < img1_g.ni(); i++)
  {
    for (unsigned int j = 0; j < img1_g.nj(); j++)
    {
      double dx = img1_g(i,j,0);
      double dy = img1_g(i,j,1);
      double mag = sqrt(dx*dx + dy*dy);
      if (mag > .5)
        mag = .5;
      g(i,j) = exp(-alpha * mag);
    }
  }

  super3d::save_cost_volume(cost_volume, g, "cost_volume_stereo.dat");
#else
  super3d::load_cost_volume(cost_volume, g, "cost_volume_stereo.dat");
#endif

  super3d::refine_depth(cost_volume, g, depth, beta, theta0, beta_end, lambda, epsilon);

  // convert inverse depths back to depths
  for (unsigned int j = 0; j < depth.nj(); j++) {
    for (unsigned int i = 0; i < depth.ni(); i++) {
      double idepth = depth(i,j) * (idepth_max - idepth_min) + idepth_min;
      depth(i,j) = 1.0 / idepth;
    }
  }

  vil_image_view<vxl_byte> output(depth.ni(), depth.nj(), 1);
  for (unsigned int i = 0; i < depth.ni(); i++)
  {
    for (unsigned int j = 0; j < depth.nj(); j++)
    {
      output(i,j) = (vxl_byte)((1.0 / depth(i,j)) * 16.0);
    }
  }

  vil_save(output, output_file().c_str());

  return 0;
}


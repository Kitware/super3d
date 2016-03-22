/*ckwg +29
 * Copyright 2011-2014 by Kitware, Inc.
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

#include <viscl/core/manager.h>

#include <vcl_iostream.h>

#include <vil/vil_load.h>
#include <vil/vil_image_view.h>
#include <vil/vil_convert.h>
#include <vil/vil_save.h>
#include <vil/algo/vil_gauss_filter.h>
#include <vil/vil_math.h>

#include <vul/vul_arg.h>

#include <super3d/depth_cl/dual_rof.h>

const char *print_cl_errstring(cl_int err);

int main(int argc, char *argv[])
{
  vul_arg<vcl_string> input_image( 0, "input image", "" );

  vul_arg<vcl_string> output_image( "-o", "output image", "denoised.png");


  vul_arg_parse( argc, argv );

  vil_image_view<unsigned short> img_byte = vil_load(input_image().c_str());
  vil_image_view<float> img, out_img;
  if (img_byte.nplanes() > 1)
    vil_convert_planes_to_grey(img_byte, img);
  else
    vil_convert_cast<unsigned short, float>(img_byte, img);

  vil_math_scale_and_offset_values(img, 1.0/65535.0, 0.0);

  out_img.set_size(img.ni(), img.nj());
  vil_image_view<float> g(img.ni(), img.nj());
  g.fill(1.0f);

  super3d::cl::dual_rof_t rof = NEW_VISCL_TASK(super3d::cl::dual_rof);
  rof->denoise(img, out_img, g, 20000, 200.0f, .25, 0.01f);

  float min, max;
  vil_math_value_range(out_img, min, max);
  float scale = 65535.0f;
  vcl_cout << min << " " << max << "\n";
  vil_math_scale_and_offset_values(out_img, scale, 0);

  vil_image_view<unsigned short> out_byte;
  vil_convert_cast(out_img, out_byte);

  vil_save(out_byte, output_image().c_str());

  return 0;
}

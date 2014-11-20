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


// VXL includes
#include <vil/vil_image_view.h>
#include <vil/vil_crop.h>
#include <vil/vil_convert.h>
#include <vil/vil_math.h>
#include <vil/vil_load.h>
#include <vul/vul_arg.h>
#include <vnl/vnl_math.h>


/// Compute the minimum and maximum value ignoring infinite and NaN values.
void finite_value_range(const vil_image_view<double>& img,
                        double& min_value, double& max_value)
{
  min_value = 0;
  max_value = 0;

  if (img.size()==0)
  {
    return;
  }

  const unsigned ni = img.ni();
  const unsigned nj = img.nj();
  assert(img.nplanes() == 1);

  vcl_ptrdiff_t istep=img.istep(),  jstep=img.jstep();
  bool first_finite = true;

  const double* row = img.top_left_ptr();
  for (unsigned j=0; j<nj; ++j, row += jstep)
  {
    const double* pixel = row;
    for (unsigned i=0; i<ni; ++i, pixel+=istep)
    {
      if (vnl_math::isfinite(*pixel))
      {
        if (first_finite)
        {
          min_value = *pixel;
          max_value = *pixel;
          first_finite = false;
        }
        else if (*pixel<min_value)
        {
          min_value=*pixel;
        }
        else if (*pixel>max_value)
        {
          max_value=*pixel;
        }
      }
    }
  }
}


int main(int argc, char* argv[])
{
  vul_arg<std::string> depth_image( 0, "Input depth image", "" );
  vul_arg<std::string> crop_window("-cw", "Crop window i0 ni j0 nj (e.g., \"0 100 0 200\"", "");

  vul_arg_parse( argc, argv );

  vil_image_view<double> depth_map = vil_convert_cast(double(), vil_load(depth_image().c_str()));

  if (!crop_window().empty())
  {
    std::istringstream cwstream(crop_window());
    int i0, ni, j0, nj;
    cwstream >> i0 >> ni >> j0 >> nj;
    std::cout << "Crop window: " << i0 << " " << ni << " " << j0 << " " << nj
              << std::endl;
    depth_map = vil_crop(depth_map, i0, ni, j0, nj);
  }

  double min_d, max_d;
  finite_value_range(depth_map, min_d, max_d);
  std::cout << "min depth: "<<min_d<<" max depth: "<<max_d<<std::endl;

  return 0;
}

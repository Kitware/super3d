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

#ifndef image_utils_h_
#define image_utils_h_

#include <vil/vil_image_view.h>
#include <vil/vil_math.h>
#include <vcl_cstdlib.h>


namespace super3d
{

/// fill an image with random values between min_v and max_v
template <typename T>
void fill_random(vil_image_view<T>& img, T min_v, T max_v)
{
  const T scale = max_v - min_v;
  for (unsigned p=0; p < img.nplanes(); ++p)
  {
    for (unsigned j=0; j < img.nj(); ++j)
    {
      for (unsigned i=0; i < img.ni(); ++i)
      {
        img(i,j,p) = static_cast<T>(vcl_rand()) / RAND_MAX * scale + min_v;
      }
    }
  }
}


/// compute the 2-norm of an image
template <typename T>
double image_norm2(const vil_image_view<T>& img)
{
  double sum = 0.0;
  for (unsigned p=0; p < img.nplanes(); ++p)
  {
    for (unsigned j=0; j < img.nj(); ++j)
    {
      for (unsigned i=0; i < img.ni(); ++i)
      {
        const T& val = img(i,j,p);
        sum += val * val;
      }
    }
  }
  return vcl_sqrt(sum);
}


/// compute the dot product between two vil_image_views
template <typename T>
double dot_product(const vil_image_view<T>& img1,
                   const vil_image_view<T>& img2)
{
  assert(img1.ni() == img2.ni());
  assert(img1.nj() == img2.nj());
  assert(img1.nplanes() == img2.nplanes());
  const unsigned np = img1.nplanes();

  vil_image_view<double> tmp;
  vil_math_image_product(img1, img2, tmp);
  double sum = 0.0;
  for (unsigned p=0; p<np; ++p)
  {
    double s = 0;
    vil_math_sum(s, tmp, p);
    sum += s;
  }
  return sum;
}

} // end namespace super3d

#endif //norm_estimation_h_

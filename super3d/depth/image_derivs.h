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

#ifndef IMAGE_DERIVS_H_
#define IMAGE_DERIVS_H_

#include <vil/vil_image_view.h>
#include <vil/vil_math.h>
#include <vcl_cstdlib.h>


namespace super3d
{

/// Compute the forward gradient of an image
/// produces an image with 2x the number of planes
template <typename T>
void forward_gradient(const vil_image_view<T>& src, vil_image_view<T>& dst)
{
  const unsigned ni = src.ni();
  const unsigned nj = src.nj();
  const unsigned np = src.nplanes();
  dst.set_size(ni, nj, 2*np);
  dst.fill(T(0));
  for (unsigned p=0; p < np; ++p)
  {
    const unsigned px = 2*p;
    const unsigned py = 2*p+1;
    for (unsigned j=0; j < nj; ++j)
    {
      for (unsigned i=0; i < ni; ++i)
      {
        const double& s = src(i,j,p);
        if( i < ni-1 )
        {
          dst(i,j,px) = src(i+1,j,p) - s;
        }
        if( j < nj-1 )
        {
          dst(i,j,py) = src(i,j+1,p) - s;
        }
      }
    }
  }
}


/// Compute the backward divergence of an image
/// produces an image with 0.5x the number of planes
template <typename T>
void backward_divergence(const vil_image_view<T>& src, vil_image_view<T>& dst)
{
  const unsigned ni = src.ni();
  const unsigned nj = src.nj();
  const unsigned np = src.nplanes() / 2;
  assert((src.nplanes() > 0) && (src.nplanes()%2 == 0));
  dst.set_size(ni, nj, np);
  for (unsigned p=0; p < np; ++p)
  {
    const unsigned px = 2*p;
    const unsigned py = 2*p+1;
    for (unsigned j=0; j < nj; ++j)
    {
      for (unsigned i=0; i < ni; ++i)
      {
        T& d = dst(i,j,p);
        d = T(0);
        if (i == 0)
          d += src(i,j,px);
        else if (i < ni - 1)
          d += src(i,j,px) - src(i-1,j,px);
        else
          d += -src(i-1,j,px);
        if (j == 0)
          d += src(i,j,py);
        else if (j < nj - 1)
          d += src(i,j,py) - src(i,j-1,py);
        else
          d += -src(i,j-1,py);
      }
    }
  }
}

} // end namespace super3d

#endif //IMAGE_DERIVS_H_

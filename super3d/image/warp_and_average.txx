/*ckwg +29
 * Copyright 2013-2015 by Kitware, Inc.
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

#include "warp_and_average.h"

#include <vil/vil_math.h>

namespace vidtk
{

template < typename PixType >
bool warp_and_average(const std::vector<vil_image_view<PixType> > &src_frames,
                      vil_image_view<PixType> &dest,
                      const std::vector<vgl_h_matrix_2d<double> > &homogs,
                      unsigned int ref_frame, unsigned int i0, unsigned int j0,
                      unsigned int ni, unsigned int nj,
                      const warp_image_parameters &wip, double scale_factor)
{
  const vil_image_view<double> &ref_img = src_frames[ref_frame];

  //Check that the crop region is inside the ref image
  if (i0 + ni > ref_img.ni() || j0 + nj > ref_img.nj())
  {
    return false;
  }

  const vgl_h_matrix_2d<double> &ref_H = homogs[ref_frame];

  vnl_double_3x3 Sinv;
  Sinv.set_identity();
  Sinv(0,0) = 1.0/scale_factor;
  Sinv(1,1) = 1.0/scale_factor;

  vnl_double_3x3 T;
  T.set_identity();
  T(0,2) = i0;
  T(1,2) = j0;

  //pre computed post multiply matrix
  const vnl_double_3x3 M = ref_H.get_matrix() * T * Sinv;

  const unsigned int dest_ni = static_cast<unsigned int>(ni * scale_factor);
  const unsigned int dest_nj = static_cast<unsigned int>(nj * scale_factor);
  const unsigned int nplanes = ref_img.nplanes();

  //Need to count the images that contribute to each pixel since warps may have unmapped pixels
  vil_image_view<double> avg(dest_ni, dest_nj, nplanes);
  vil_image_view<int> count(dest_ni, dest_nj, 1);
  avg.fill(0.0);
  count.fill(0);

  for (unsigned int f = 0; f < src_frames.size(); f++)
  {
    vil_image_view<PixType> warp;
    warp.set_size(dest_ni, dest_nj, nplanes);
    vgl_h_matrix_2d<double> H = homogs[f].get_inverse().get_matrix() * M;
    warp_image(src_frames[f], warp, H, wip);

    for (unsigned int j = 0; j < dest_nj; j++)
    {
      for (unsigned int i = 0; i < dest_ni; i++)
      {
        if (warp(i,j) != wip.unmapped_value_)
        {
          for (unsigned int k = 0; k < nplanes; k++)
          {
            avg(i,j,k) += static_cast<double>(warp(i,j,k));
          }

          count(i,j)++;
        }
      }
    }
  }

  vil_math_image_ratio(avg, count, dest);
  return true;
}


template < typename PixType >
void warp_and_average(const std::vector<vil_image_view<PixType> > &src_frames,
                      vil_image_view<PixType> &dest,
                      const std::vector<vgl_h_matrix_2d<double> > &homogs,
                      unsigned int ref_frame, const warp_image_parameters &wip,
                      double scale_factor)
{
  warp_and_average(src_frames, dest, homogs, ref_frame, 0, 0,
                   src_frames[ref_frame].ni(), src_frames[ref_frame].nj(), wip, scale_factor);
}

}

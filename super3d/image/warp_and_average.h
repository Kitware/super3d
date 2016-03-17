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

#ifndef super3d_warp_and_average_h_
#define super3d_warp_and_average_h_

#include <vil/vil_image_view.h>
#include <vgl/algo/vgl_h_matrix_2d.h>
#include <vector>

#include "warp_image.h"

//template < typename PixType >
//class vil_image_view;

namespace super3d
{

/// \brief Warps images to a reference frame and averages them with a crop region
/// returns false if the crop region is not fully contained in ref_frame
///
/// \param src_frames The input images.
/// \param dest the output averaged image
/// \param homogs the homographies for src_frames, these are assumed to warp to the same plane
/// \param ref_frame the index into src_frames of the image to warp to
/// \param i0 left crop pixel index
/// \param j0 top crop pixel index
/// \param ni width of crop region
/// \param nj height of crop region
/// \param wip parameters to warp image, set the interp method and the unwarped value
/// \param scale_factor the output scale of dest
template < typename PixType >
bool warp_and_average(const std::vector<vil_image_view<PixType> > &src_frames,
                      vil_image_view<PixType> &dest,
                      const std::vector<vgl_h_matrix_2d<double> > &homogs,
                      unsigned int ref_frame, unsigned int i0, unsigned int j0,
                      unsigned int ni, unsigned int nj, const warp_image_parameters &wip,
                      double scale_factor = 1.0);

/// \brief Warps images to a reference frame and averages them
///
/// \param src_frames The input images.
/// \param dest the output averaged image
/// \param homogs the homographies for src_frames, these are assumed to warp to the same plane
/// \param ref_frame the index into src_frames of the image to warp to
/// \param wip parameters to warp image, set the interp method and the unwarped value
/// \param scale_factor the output scale of dest
template < typename PixType >
void warp_and_average(const std::vector<vil_image_view<PixType> > &src_frames,
                      vil_image_view<PixType> &dest,
                      const std::vector<vgl_h_matrix_2d<double> > &homogs,
                      unsigned int ref_frame, const warp_image_parameters &wip,
                      double scale_factor = 1.0);
}

#endif // super3d_warp_and_average

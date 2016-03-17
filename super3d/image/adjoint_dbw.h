/*ckwg +29
 * Copyright 2012-2015 by Kitware, Inc.
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


#ifndef adjoint_dbw_h_
#define adjoint_dbw_h_

#include <vil/vil_image_view.h>

#include "adjoint_image_op.h"


namespace vidtk
{

/// Create DBW adjoint operator from a flow field and various options.
/// \param flow The optical flow field used in the warping step, also defined input resolution.
/// \param ni The width of the output image after applying DBW
/// \param nj The height of the output image after applying DBW
/// \param np The number of planes (e.g. color channels) in input and output
/// \param scale_factor The scale factor for the down sampling step
/// \param smoothing_sigma The amount of bluring for the blur step.
/// \param down_sample_averaging If true, average pixels in a box to down scale
///        rather than just picking one pixel to sample.
/// \param bicubic_warping toggle between bilinear (false) and bicubic (true)
adjoint_image_ops_func<double>
create_dbw_from_flow(const vil_image_view<double> &flow,
                     const unsigned ni, const unsigned nj, const unsigned np,
                     int scale_factor,
                     double smoothing_sigma,
                     bool down_sample_averaging = false,
                     bool bicubic_warping = false);

} // end namespace vidtk

#endif //adjoint_dbw_h_

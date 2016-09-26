/*ckwg +29
 * Copyright 2012-2016 by Kitware, Inc.
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

#ifndef DBW_H_
#define DBW_H_

#include "depth_config.h"

#include <vil/vil_image_view.h>
#include <vector>
#include <super3d/image/adjoint_dbw.h>


namespace super3d
{

/// Create DBW adjoint operator from a flow field, weights, and various options
/// \param flow The optical flow field used in the warping step, also defined input resolution.
/// \param weights The weights image multiplied by the input before warping
/// \param ni The width of the output image after applying DBW
/// \param nj The height of the output image after applying DBW
/// \param np The number of planes (e.g. color channels) in input and output
/// \param scale_factor The scale factor for the down sampling step
/// \param down_sample_averaging If true, average pixels in a box to down scale
///        rather than just picking one pixel to sample.
SUPER3D_DEPTH_EXPORT
super3d::adjoint_image_ops_func<double>
create_dbw_from_flow(const vil_image_view<double> &flow,
                     const vil_image_view<double> &weights,
                     const unsigned ni, const unsigned nj, const unsigned np,
                     int scale_factor, double sensor_sigma, bool down_sample_averaging,
                     bool bicubic_warping);

} // end namespace super3d

#endif

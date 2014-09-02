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

#ifndef multiscale_h_
#define multiscale_h_

#include "depth_config.h"

#include <vil/vil_image_view.h>
#include <vil/algo/vil_gauss_filter.h>

#include <vpgl/vpgl_perspective_camera.h>
#include <video_transforms/warp_image.h>

namespace super3d
{

/// Produce the camera corresponding to a downsampled image
/// \param camera The input camera
/// \param scale The scale of the downsampled image (default is 0.5)
/// \return A camera corresponding to the downsampled image
SUPER3D_DEPTH_EXPORT
vpgl_perspective_camera<double>
scale_camera(const vpgl_perspective_camera<double>& camera,
                  double scale=0.5);


/// Produce the camera corresponding to a cropped image
/// \param camera The input camera
/// \param left the left coordinate of the cropping
/// \param top the left coordinate of the cropping
/// \return A camera corresponding to the cropped image
SUPER3D_DEPTH_EXPORT
vpgl_perspective_camera<double>
crop_camera(const vpgl_perspective_camera<double>& camera, double left, double top);


/// Compute the scaling of a homography for downsampled images.
/// Assumes the same \a scale has been used to sample each image
/// \param H The homography for the original images
/// \param scale The scale of the downsampled images (default is 0.5)
/// \returns The homography mapping between the downsampled images
SUPER3D_DEPTH_EXPORT
vnl_matrix_fixed<double,3,3>
scale_homography(const vnl_matrix_fixed<double,3,3>& H,
                      double scale = 0.5);


/// Compute the scaling of a homogeneous point in a downsampled image.
/// \param p The homogeneous point in the original image
/// \param scale The scale of the downsampled image (default is 0.5)
/// \returns The homogeneous point in the downsampled image
SUPER3D_DEPTH_EXPORT
vnl_vector_fixed<double,3>
scale_point(const vnl_vector_fixed<double,3>& p,
                 double scale = 0.5);

/// Rescale image src to image dest by a factor of scale_factor
/// using interpolation method interp
/// \param src The input/source image
/// \param dest The rescale/output images
/// \param scale_factor The rescale factor, e.g. 0.5, 2.0
/// \param interp The interpolation method, e.g. bilinear, bicubic
SUPER3D_DEPTH_EXPORT
void upsample(const vil_image_view<double> &src, vil_image_view<double> &dest,
              double scale_factor, vidtk::warp_image_parameters::interp_type interp);

} // end namespace super3d

#endif // multiscale_h_

/*ckwg +29
 * Copyright 2010-2013 by Kitware, Inc.
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

#ifndef dual_rof_denoise_h_
#define dual_rof_denoise_h_

#include <vil/vil_image_view.h>

/// \file
/// Denoise images using total variation with a L1 norm.
///
/// This denoising code is also useful in dense optical flow
/// and dense depth estimation.  The implementation below is
/// based on the equations in the following paper:
///
/// A. Wedel, T. Pock, C. Zach, H. Bischof, and D. Cremers,
/// "An improved algorithm for TV-L1 optical flow",
/// Statistical and Geometrical Approaches to Visual Motion Analysis
/// pg 23-25, Springer, 2009.
/// http://cvpr.in.tum.de/old/pub/pub/DagstuhlOpticalFlowChapter.pdf
///
/// Templates are instantiated with float and double in the .cxx file
/// and are intended to be only instantiated with double and float.


namespace vidtk
{

/// Apply several iterations of the dual Rudin, Osher and Fatemi model.
/// Solves the total variation minimization with L1 norm over the image.
/// The resulting image is denoised, but sharp edges are preserved.
/// Solves min_u Int { |grad(u(i,j))| + 1/(2*theta) (u(i,j)-v(i,j))^2} di dj
/// where u is \a dest and v is \a src.
/// \param src The source image
/// \retval dest The destination image
/// \param iterations The number of iterations
/// \param theta A tuning parameter to control amount of smoothing.
///              Larger values produce more smoothing.
/// \param step The step size for each iteration
template <typename T>
void
dual_rof_denoise(const vil_image_view<T>& src,
                 vil_image_view<T>& dest,
                 unsigned iterations,
                 T theta, T step = 0.25);


/// Apply several iterations of the dual Rudin, Osher and Fatemi model.
/// Solves the total variation minimization with L1 norm over the image.
/// The resulting image is denoised, but sharp edges are preserved.
/// This version starts from and updates an existing dual variable
/// \param src The source image
/// \retval dest The destination image
/// \param dual the dual variable
/// \param iterations The number of iterations
/// \param theta A tuning parameter to control amount of smoothing.
///              Larger values produce more smoothing.
/// \param step The step size for each iteration
template <typename T>
void
dual_rof_denoise(const vil_image_view<T>& src,
                 vil_image_view<T>& dest,
                 vil_image_view<T>& dual,
                 unsigned iterations,
                 T theta, T step = 0.25);


/// Apply several iterations of the weighted dual Rudin, Osher and Fatemi model.
/// Solves the total variation minimization with L1 norm over the image.
/// The resulting image is denoised, but sharp edges are preserved.
/// A weight image allows control over the amount of denoising at each pixel.
/// Solves min_u Int { w(i,j)*|grad(u(i,j))| + 1/(2*theta) (u(i,j)-v(i,j))^2} di dj
/// where u is \a dest and v is \a src and w is \a weights.
/// \param src The source image
/// \param weights An image of weights in [0,1].  A value of 1 results in denoising
///                as usual.  Smaller values reduce the denoising amount
/// \retval dest The destination image
/// \param iterations The number of iterations
/// \param theta A tuning parameter to control amount of smoothing.
///              Larger values produce more smoothing.
/// \param step The step size for each iteration
template <typename T>
void
dual_rof_weighted_denoise(const vil_image_view<T>& src,
                          const vil_image_view<T>& weights,
                          vil_image_view<T>& dest,
                          unsigned iterations,
                          T theta, T step = 0.25);


/// Add the scaled divergence of a vector field to the source image.
/// Compute dest = src + scale*div(vec)
/// \param vec The vector field (x,y in planes 0,1)
/// \param src The source image
/// \param scale The scale factor applied to the divergence term
/// \retval dest The destination image
template <typename T>
void
add_scaled_divergence(const vil_image_view<T>& vec,
                      const vil_image_view<T>& src,
                      T scale,
                      vil_image_view<T>& dest);


/// Add (in place) the scaled gradient of the source image to the vector field.
/// Compute vec += scale*grad(src)
/// \param src The source image
/// \param vec The vector field (x,y in planes 0,1)
/// \param scale The scale factor applied to the gradient vectors
template <typename T>
void
add_scaled_gradient(const vil_image_view<T>& src,
                    vil_image_view<T>& vec,
                    T scale);


/// Truncate all vectors greater than unit length to unit length.
/// \retval vec The input vector field, modified in place (x,y in planes 0,1).
template <typename T>
void
truncate_vectors(vil_image_view<T>& vec);


/// Truncate all vectors greater in length than \a weight to a magnitude of \a weight.
/// \param weight is a image of weights in [0,1] and varies at each pixel location.
/// \retval vec The input vector field, modified in place (x,y in planes 0,1).
template <typename T>
void
truncate_vectors(vil_image_view<T>& vec,
                 const vil_image_view<T>& weights);

}

#endif // dual_rof_denoise_h_

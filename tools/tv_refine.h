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

#ifndef tv_refine_h_
#define tv_refine_h_

#include <vector>
#include <vil/vil_image_view.h>
#include <vnl/vnl_matrix_fixed.h>
#include <vnl/vnl_vector_fixed.h>


/// compute the depth refinement using multiple images.
/// \note The size of vectors \a I1, \a H, and \a e should be the same.
/// \param I0 The reference image.
/// \param I1 The vector of destination images.
/// \param H The vector of homographies induced by the planes at inifinity.
/// \param e The vector eipoles in the destination images.
/// \retval depth The current depth map image (relative to \a I0)
/// \param outer_iterations The number of outer (warping) iterations
/// \param inner_iterations The nummber of inner (propagating) iterations
void
refine_depths(const vil_image_view<double>& I0,
              const std::vector<vil_image_view<double> >& I1,
              const std::vector<vnl_matrix_fixed<double,3,3> >& H,
              const std::vector<vnl_vector_fixed<double,3> >& e,
              vil_image_view<double>& depth,
              unsigned outer_iterations = 10,
              unsigned inner_iterations = 10,
              double theta = 0.00125,
              double lambda = 1e-2);


/// compute the linear brightness constancy constraint at each pixel
/// \param I0 The reference image.
/// \param I1 The destination image.
/// \param I1xy The destination image gradients (x in plane 0, y in plane 1).
/// \param depth The current depth map image.
/// \param H The homography induced by the plane at inifinity.
/// \param e The eipole in the destination image.
/// \retval bcc The linear brighness constance constraints (2-plane image).
///             BCC is satisfied when bcc(i,j,0) + depth(i,j)*bcc(i,j,1) == 0
void
compute_linear_bcc(const vil_image_view<double>& I0,
                   const vil_image_view<double>& I1,
                   const vil_image_view<double>& I1xy,
                   const vil_image_view<double>& depth,
                   const vnl_matrix_fixed<double,3,3>& H,
                   const vnl_vector_fixed<double,3>& e,
                   vil_image_view<double>& mask,
                   vil_image_view<double>& bcc);


/// Apply multiple image BCCs for a least squares depth refinement
/// \param bcc The linear brighness constance constraints (2 planes per image).
///            BCC is satisfied for image k when
///            bcc(i,j,2*k) + depth(i,j)*bcc(i,j,2*k+1) == 0
/// \param lambda Determines the influence of the data (BCC) term.
/// \retval depth The current depth map image, modified in place.
void
apply_bcc_to_depth_lst_sqr_two(const vil_image_view<double>& bcc,
                               const vil_image_view<double>& mask,
                               double lambda,
                               vil_image_view<double>& depth);

void
apply_bcc_to_depth_all(const vil_image_view<double>& bcc,
                       const vil_image_view<double>& mask,
                       double lambda, double theta,
                       vil_image_view<double>& depth);

double
eval_vij(const vil_image_view<double>& bcc,
         double lambda, double theta,
         double depth,
         double v, unsigned i, unsigned j);

void output_function(const vil_image_view<double>& bcc,
                     double lambda, double theta,
                     double depth,
                     unsigned i, unsigned j);

#endif // tv_refine_h_

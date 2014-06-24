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

#ifndef SUPER_RES_H_
#define SUPER_RES_H_

#include "depth_config.h"

#include <vil/vil_image_view.h>
#include <vcl_vector.h>

#include <video_transforms/adjoint_image_op.h>


namespace super3d
{

struct super_res_params {
  double lambda, epsilon_data, epsilon_reg, tau, sigma;
  unsigned int s_ni, s_nj, l_ni, l_nj;
  double scale_factor;
  unsigned int ref_frame;
};


SUPER3D_DEPTH_EXPORT
void super_resolve(const vcl_vector<vil_image_view<double> > &frames,
                   const vcl_vector<vidtk::adjoint_image_ops_func<double> > &warps,
                   vil_image_view<double> &u,
                   const super_res_params &srp,
                   unsigned int iterations,
                   const vcl_string &output_image = "");


SUPER3D_DEPTH_EXPORT
void compare_to_original(const vil_image_view<double> &ref_img,
                         const vil_image_view<double> &super,
                         const vil_image_view<double> &original,
                         unsigned int scale_factor);

} // end namespace super3d

#endif

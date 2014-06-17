/*
 * Copyright 2013-2014 Kitware, Inc.
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

#ifndef CL_SUPER_RES_H_
#define CL_SUPER_RES_H_

#include <viscl/core/task.h>
#include <viscl/core/image.h>
#include <viscl/core/buffer.h>
#include <vcl_vector.h>

#include <vil/vil_image_view.h>

class super_res_cl : public viscl::task
{
public:

  super_res_cl();

  struct params {
    float lambda, epsilon_data, epsilon_reg, tau, sigma;
    cl_int2 sdim, ldim;
    float scale_factor;
    float sensor_sigma;
    float dual_p_denom, dual_q_denom, sf_2;
  };


void super_resolve(const vcl_vector<vil_image_view<float> > &frames,
                   const vcl_vector<vil_image_view<float> > &flows,
                   vil_image_view<float> &u,
                   params &srp,
                   unsigned int iterations);


private:

  viscl::cl_kernel_t dual_step_q_k, dual_step_p_k,
                     primal_step_u_k, down_sample_k, up_sample_k,
                     warp_foward_k, warp_backward_k, blur1D_horiz_k, blur1D_vert_k,
                     zero_k, zero2_k;

};

#endif

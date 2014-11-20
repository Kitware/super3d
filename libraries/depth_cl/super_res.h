/*ckwg +5
 * Copyright 2013 by Kitware, Inc. All Rights Reserved. Please refer to
 * KITWARE_LICENSE.TXT for licensing information, or contact General Counsel,
 * Kitware, Inc., 28 Corporate Drive, Clifton Park, NY 12065.
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

/*ckwg +5
 * Copyright 2012 by Kitware, Inc. All Rights Reserved. Please refer to
 * KITWARE_LICENSE.TXT for licensing information, or contact General Counsel,
 * Kitware, Inc., 28 Corporate Drive, Clifton Park, NY 12065.
 */

#ifndef CL_DUAL_ROF_H_
#define CL_DUAL_ROF_H_

#include <vcl_vector.h>
#include <vcl_string.h>

#include <viscl/core/task.h>
#include <viscl/core/image.h>
#include <viscl/core/buffer.h>

#include <vil/vil_image_view.h>

class dual_rof : public viscl::task
{
public:

  dual_rof();

  viscl::buffer create_dual(size_t ni, size_t nj);

  void denoise(const vil_image_view<float> &src,
               vil_image_view<float> &dest,
               vil_image_view<float> &g,
               int iterations,
               float lambda,
               float step,
               float epsilon);

  void denoise(const viscl::buffer &denoised,
               const viscl::buffer &dual,
               const viscl::buffer &src,
               const viscl::image &g,
               size_t ni, size_t nj,
               int iterations,
               float lambda,
               float step,
               float epsilon);

private:

  viscl::cl_kernel_t gradient_k, divergence_k, init_dual_k;
};

typedef boost::shared_ptr<dual_rof> dual_rof_t;

#endif

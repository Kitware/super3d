/*ckwg +29
 * Copyright 2012-2014 by Kitware, Inc.
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

#ifndef CL_DUAL_ROF_H_
#define CL_DUAL_ROF_H_

#include "depth_cl_config.h"

#include <vcl_vector.h>
#include <vcl_string.h>

#include <viscl/core/task.h>
#include <viscl/core/image.h>
#include <viscl/core/buffer.h>

#include <vil/vil_image_view.h>


namespace super3d
{

namespace cl
{

class SUPER3D_DEPTH_CL_EXPORT dual_rof : public viscl::task
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

} // end namespace cl

} // end namespace super3d

#endif

/*ckwg +5
 * Copyright 2012 by Kitware, Inc. All Rights Reserved. Please refer to
 * KITWARE_LICENSE.TXT for licensing information, or contact General Counsel,
 * Kitware, Inc., 28 Corporate Drive, Clifton Park, NY 12065.
 */

#ifndef CL_REFINE_DEPTH_H_
#define CL_REFINE_DEPTH_H_

#include <vcl_vector.h>
#include <vcl_string.h>

#include <viscl/core/task.h>
#include <viscl/core/image.h>
#include <viscl/core/buffer.h>

#include <vil/vil_image_view.h>

class refine_depth_cl : public viscl::task
{
public:

  refine_depth_cl();

  void refine(const vil_image_view<float> &cost_volume,
              vil_image_view<float> &d,
              const vil_image_view<float> &g,
              float beta,
              float theta0,
              float beta_end,
              float lambda);

private:

  viscl::cl_queue_t queue;
  viscl::cl_kernel_t init_depth_k, search_k;
};

typedef boost::shared_ptr<refine_depth_cl> refine_depth_cl_t;

#endif

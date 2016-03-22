/*ckwg +29
 * Copyright 2012 by Kitware, Inc.
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

#ifndef tv_refine_search_h_
#define tv_refine_search_h_

#include "depth_config.h"

#include <vcl_vector.h>
#include <vil/vil_image_view.h>
#include <vnl/vnl_matrix_fixed.h>
#include <vnl/vnl_vector_fixed.h>
#include <vpgl/vpgl_perspective_camera.h>

#include <boost/function.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/locks.hpp>


namespace super3d
{

class SUPER3D_DEPTH_EXPORT depth_refinement_monitor;


SUPER3D_DEPTH_EXPORT
void
refine_depth(vil_image_view<double> &cost_volume,
             const vil_image_view<double> &g,
             vil_image_view<double> &d,
             unsigned int iterations,
             double theta0,
             double beta_end,
             double lambda,
             double epsilon,
             depth_refinement_monitor *drm = NULL);


class SUPER3D_DEPTH_EXPORT depth_refinement_monitor
{
public:

  struct update_data
  {
    vil_image_view<double> current_result;
    unsigned int num_iterations;
  };

  depth_refinement_monitor(boost::function<void (update_data)> callback,
                           unsigned int interval,
                           boost::shared_ptr<bool> interrupted) : callback_(callback),
                                                                  interval_(interval),
                                                                  interrupted_(interrupted) {}

private:

  friend void refine_depth(vil_image_view<double> &cost_volume,
                           const vil_image_view<double> &g,
                           vil_image_view<double> &d,
                           unsigned int iterations,
                           double theta0,
                           double beta_end,
                           double lambda,
                           double epsilon,
                           depth_refinement_monitor *drm);

  boost::function<void (update_data)> callback_;
  unsigned int interval_;
  boost::shared_ptr<bool const> interrupted_;
};


//semi-implicit gradient ascent on q and descent on d
SUPER3D_DEPTH_EXPORT
void huber(vil_image_view<double> &q,
           vil_image_view<double> &d,
           const vil_image_view<double> &a,
           const vil_image_view<double> &g,
           double theta,
           double step,
           double epsilon);


SUPER3D_DEPTH_EXPORT
void huber_central(vil_image_view<double> &q,
           vil_image_view<double> &d,
           const vil_image_view<double> &a,
           const vil_image_view<double> &g,
           double theta,
           double step,
           double epsilon);


SUPER3D_DEPTH_EXPORT
void hessian_frob(vil_image_view<double> &q,
                  vil_image_view<double> &d,
                  const vil_image_view<double> &a,
                  const vil_image_view<double> &g,
                  double theta,
                  double step,
                  double epsilon);


SUPER3D_DEPTH_EXPORT
double eval_hessian_frob(const vil_image_view<double> &d,
                         const vil_image_view<double> &costvol,
                         double lambda);

} // end namespace super3d

#endif // tv_refine_search_h_

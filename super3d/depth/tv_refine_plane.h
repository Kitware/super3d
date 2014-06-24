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

#ifndef tv_refine_plane_h_
#define tv_refine_plane_h_

#include "depth_config.h"
#include "world_space.h"

#include <vcl_vector.h>
#include <vil/vil_image_view.h>
#include <vnl/vnl_matrix_fixed.h>
#include <vnl/vnl_vector_fixed.h>
#include <vpgl/vpgl_perspective_camera.h>


namespace super3d
{

SUPER3D_DEPTH_EXPORT
void refine_depth_planar(world_space *ws,
                         vil_image_view<double> &d,
                         vil_image_view<double> &g,
                         const vil_image_view<double> &ref,
                         double lambda,
                         double step,
                         double epsilon);


SUPER3D_DEPTH_EXPORT
void refine_depth_planar(const vil_image_view<double> &cost_volume,
                         world_space *ws,
                         vil_image_view<double> &d,
                         const vil_image_view<double> &g,
                         const vil_image_view<double> &ref,
                         double beta,
                         double theta0,
                         double theta_end,
                         double lambda,
                         double epsilon);


SUPER3D_DEPTH_EXPORT
void
huber_planar_rof(vil_image_view<double> &q,
           vil_image_view<double> &d,
           const vil_image_view<double> &g,
           const vil_image_view<double> &a,
           const vil_image_view<double> &n,
           double lambda,
           double step,
           double epsilon);


SUPER3D_DEPTH_EXPORT
void
huber_planar_coupled(vil_image_view<double> &q,
           vil_image_view<double> &d,
           const vil_image_view<double> &g,
           const vil_image_view<double> &a,
           const vil_image_view<double> &n,
           double theta,
           double step,
           double epsilon);

} // end namespace super3d

#endif // tv_refine_plane_h_

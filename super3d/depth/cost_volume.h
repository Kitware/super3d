/*ckwg +29
 * Copyright 2012-2016 by Kitware, Inc.
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


#ifndef COST_VOLUME_H_
#define COST_VOLUME_H_

#include "depth_config.h"

#include <vector>
#include <vil/vil_image_view.h>
#include <vpgl/vpgl_perspective_camera.h>

#include "world_space.h"

namespace super3d
{

SUPER3D_DEPTH_EXPORT
bool compute_depth_range(const vpgl_perspective_camera<double> &ref_cam,
                         int i0, int ni, int j0, int nj, const std::string &landmark_file,
                         double &min_depth, double &max_depth);

SUPER3D_DEPTH_EXPORT

void
compute_world_cost_volume(const std::vector<vil_image_view<double> > &frames,
                          const std::vector<vpgl_perspective_camera<double> > &cameras,
                          world_space *ws,
                          unsigned int ref_frame,
                          unsigned int S,
                          vil_image_view<double> &cost_volume,
                          double intesity_weight,
                          double gradient_weight,
                          double census_weight);

SUPER3D_DEPTH_EXPORT
void
compute_cost_volume_warp(const std::vector<vil_image_view<double> > &frames,
                         const std::vector<vpgl_perspective_camera<double> > &cameras,
                         unsigned int ref_frame,
                         unsigned int S,
                         double depth_min,
                         double depth_max,
                         vil_image_view<double> &cost_volume);

SUPER3D_DEPTH_EXPORT
void
compute_cost_volume_bp(const std::vector<vil_image_view<double> > &frames,
                       const std::vector<vpgl_perspective_camera<double> > &cameras,
                       unsigned int ref_frame,
                       unsigned int S,
                       double depth_min,
                       double depth_max,
                       vil_image_view<double> &cost_volume);


//disparity increases along negative x
SUPER3D_DEPTH_EXPORT
void compute_cost_volume_rectified(const std::vector<vil_image_view<double> > &frames,
                                  unsigned int ref_frame,
                                  unsigned int S,
                                  double idepth_min,
                                  double idepth_max,
                                  vil_image_view<double> &cost_volume);

SUPER3D_DEPTH_EXPORT
vxl_uint_64 compute_census(const vil_image_view<double> &ref, int u, int v);

SUPER3D_DEPTH_EXPORT
unsigned int hamming_distance(vxl_uint_64 l, vxl_uint_64 r);

struct g_census {
  vxl_uint_64 ori, mag;
};

SUPER3D_DEPTH_EXPORT
g_census compute_g_census(const vil_image_view<double> &grad, int u, int v);

SUPER3D_DEPTH_EXPORT
void save_cost_volume(const vil_image_view<double> &cost_volume,
                      const vil_image_view<double> &g_weight,
                      const char *file_name);

SUPER3D_DEPTH_EXPORT
void load_cost_volume(vil_image_view<double> &cost_volume,
                      vil_image_view<double> &g_weight,
                      const char *file_name);

SUPER3D_DEPTH_EXPORT
void
read_cost_volume_at(FILE *file,
                    unsigned int *dims,
                    unsigned int i,
                    unsigned int j,
                    vnl_vector<double> &values);

} // end namespace suepr3d


#endif

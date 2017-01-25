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
void
compute_world_cost_volume(const std::vector<vil_image_view<double> > &frames,
                          const std::vector<vpgl_perspective_camera<double> > &cameras,
                          world_space *ws,
                          unsigned int ref_frame,
                          unsigned int S,
                          vil_image_view<double> &cost_volume,
                          double intesity_weight,
                          double gradient_weight,
                          double census_weight,
                          std::vector<vil_image_view<double> > *masks = NULL);

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

/// Return a subset of landmark points that project into the given region of interest
/// \param camera is the camera used to project the points
/// \param i0 is the horizontal coordinate of upper left corner of the ROI
/// \param ni is the width of the ROI
/// \param j0 is the vertical coordinate of upper left corner of the ROI
/// \param nj is the height of the ROI
/// \param landmarks is the set of 3D landmark points to project
/// \return the subset of \p landmarks that project into the ROI
SUPER3D_DEPTH_EXPORT
std::vector<vnl_double_3>
filter_visible_landmarks(const vpgl_perspective_camera<double> &camera,
                         int i0, int ni, int j0, int nj,
                         const std::vector<vnl_double_3> &landmarks);

/// Robustly compute the bounding planes of the landmarks in a given direction
/// \param  landmarks is the set of 3D landmark points
/// \param  normal is the normal vector of the plane
/// \retval min_offset is the minimum plane offset
/// \retval max_offset is the maximum plane offset
/// \param  outlier_thresh is the threshold for fraction of outlier offsets to
///         reject at both the top and bottom
/// \param  safety_margin_factor is the fraction of total offset range to pad
///         both top and bottom to account for insufficient landmark samples
SUPER3D_DEPTH_EXPORT
void
compute_offset_range(const std::vector<vnl_double_3> &landmarks,
                     const vnl_vector_fixed<double,3> &normal,
                     double &min_offset, double &max_offset,
                     const double outlier_thresh = 0.05,
                     const double safety_margin_factor = 0.33);

/// Robustly compute the depth range of the landmarks with respect to a camera
/// \param  landmarks is the set of 3D landmark points
/// \param  camera is the camera used relative to which depths are computed
/// \retval min_depth is the computed minimum depth
/// \retval max_depth is the computed maximum depth
/// \param  outlier_thresh is the threshold for fraction of outlier offsets to
///         reject at both the near and far planes
/// \param  safety_margin_factor is the fraction of total depth range to pad
///         both near and far to account for insufficient landmark samples
SUPER3D_DEPTH_EXPORT
void
compute_depth_range(const std::vector<vnl_double_3> &landmarks,
                    const vpgl_perspective_camera<double> &camera,
                    double &min_depth, double &max_depth,
                    const double outlier_thresh = 0.05,
                    const double safety_margin_factor = 0.33);

} // end namespace suepr3d


#endif

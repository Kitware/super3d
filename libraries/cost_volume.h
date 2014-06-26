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


#ifndef COST_VOLUME_H_
#define COST_VOLUME_H_

#include <vcl_vector.h>
#include <vil/vil_image_view.h>
#include <vpgl/vpgl_perspective_camera.h>

#include "world_space.h"

void
compute_world_cost_volume(const vcl_vector<vil_image_view<double> > &frames,
                          const vcl_vector<vpgl_perspective_camera<double> > &cameras,
                          world_space *ws,
                          unsigned int ref_frame,
                          unsigned int S,
                          vil_image_view<double> &cost_volume,
                          double intesity_weight,
                          double gradient_weight,
                          double census_weight);

void
compute_cost_volume_warp(const vcl_vector<vil_image_view<double> > &frames,
                         const vcl_vector<vpgl_perspective_camera<double> > &cameras,
                         unsigned int ref_frame,
                         unsigned int S,
                         double depth_min,
                         double depth_max,
                         vil_image_view<double> &cost_volume);

void
compute_cost_volume_bp(const vcl_vector<vil_image_view<double> > &frames,
                       const vcl_vector<vpgl_perspective_camera<double> > &cameras,
                       unsigned int ref_frame,
                       unsigned int S,
                       double depth_min,
                       double depth_max,
                       vil_image_view<double> &cost_volume);


//disparity increases along negative x
void compute_cost_volume_rectified(const vcl_vector<vil_image_view<double> > &frames,
                                  unsigned int ref_frame,
                                  unsigned int S,
                                  double idepth_min,
                                  double idepth_max,
                                  vil_image_view<double> &cost_volume);

vxl_uint_64 compute_census(const vil_image_view<double> &ref, int u, int v);

unsigned int hamming_distance(vxl_uint_64 l, vxl_uint_64 r);

struct g_census {
  vxl_uint_64 ones, twos, fours;
};

g_census compute_g_census(const vil_image_view<double> &grad, int u, int v, double thresh);

unsigned int hamming_distance(const g_census &l, const g_census &r);

void save_cost_volume(const vil_image_view<double> &cost_volume,
                      const vil_image_view<double> &g_weight,
                      const char *file_name);

void load_cost_volume(vil_image_view<double> &cost_volume,
                      vil_image_view<double> &g_weight,
                      const char *file_name);

void
read_cost_volume_at(FILE *file,
                    unsigned int *dims,
                    unsigned int i,
                    unsigned int j,
                    vnl_vector<double> &values);

void compute_depth_range(const vpgl_perspective_camera<double> &ref_cam, unsigned int i0,
                         unsigned int ni, unsigned int j0, unsigned int nj,
                         const vcl_string &landmark_file, double &min_depth, double &max_depth);


#endif

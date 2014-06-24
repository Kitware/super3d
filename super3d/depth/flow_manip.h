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

#ifndef FLOW_MANIP_H_
#define FLOW_MANIP_H_

#include "depth_config.h"

#include <vcl_vector.h>

#include <vil/vil_image_view.h>
#include <vpgl/vpgl_perspective_camera.h>
#include <vgl/vgl_box_2d.h>


namespace super3d
{

/// Compute optical flows induced by camera pair and a depth map.
/// One flow field is produced for each camera in \a cameras.
/// \param cameras The cameras viewing destination images of the flow fields
/// \param ref_cam The reference camera viewing the source of all flow fields.
/// \param depth  The depth map relative to \a ref_cam.
/// \retval flows The computed optical flow fields from \a ref_cam to all \a cameras.
SUPER3D_DEPTH_EXPORT
void compute_flows_from_depth(const vcl_vector<vpgl_perspective_camera<double> > &cameras,
                              const vpgl_perspective_camera<double> &ref_cam,
                              const vil_image_view<double> &depth,
                              vcl_vector<vil_image_view<double> > &flows);


/// Compute optical flows induced by camera pair and a depth map.
/// One flow field is produced for each camera in \a cameras.
/// This is like the \c compute_flows_from_depth except it tries to account
/// for occulsions by invalidating flow that is occluded by other flow.
/// \param cameras The cameras viewing destination images of the flow fields
/// \param ref_cam The reference camera viewing the source of all flow fields.
/// \param depth  The depth map relative to \a ref_cam.
/// \retval flows The computed optical flow fields from \a ref_cam to all \a cameras.
SUPER3D_DEPTH_EXPORT
void compute_occluded_flows_from_depth(const vcl_vector<vpgl_perspective_camera<double> > &cameras,
                                       const vpgl_perspective_camera<double> &ref_cam,
                                       const vil_image_view<double> &depth,
                                       vcl_vector<vil_image_view<double> > &flows);


/// Compute the axis-aligned bounding box contain all flow vectors destinations.
/// \note this ignores all pixels with invalid flow
SUPER3D_DEPTH_EXPORT
vgl_box_2d<double>
flow_destination_bounds(const vil_image_view<double> &flow);


/// Translate all flow vectors in the field by (dx, dy) in place.
/// \param flow The flow field to translate in place
/// \param dx The value added to the first (X) channel
/// \param dy The value added to the second (Y) channel
SUPER3D_DEPTH_EXPORT
void translate_flow(vil_image_view<double> &flow, double dx, double dy);

} // end namespace super3d

#endif

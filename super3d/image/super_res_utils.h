/*ckwg +29
 * Copyright 2013-2015 by Kitware, Inc.
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

#ifndef super_res_utils_h_
#define super_res_utils_h_

#include <vgl/vgl_intersection.h>
#include <vil/vil_image_view.h>
#include <vgl/algo/vgl_h_matrix_2d.h>
#include <super3d/image/adjoint_image_op.h>

namespace super3d
{

/// Computes an optical flow field from a homography at a crop region
/// \param homogs the homographies to be converted to flow
/// \param ref_frame the index of the target frame for the flow
/// \param frames the vector of images parallel with homogs
/// \param i0 \param j0 \param ni \param nj specify the crop region in the ref frame
/// \param scale_factor how much to scale the images before computing the flow
/// \param flows the output optical flow, a 2 plane image for x and y components
void homogs_to_flows(const std::vector<vgl_h_matrix_2d<double> > &homogs,
                     const int ref_frame,
                     int i0, int j0, int ni, int nj,
                     const int scale_factor,
                     std::vector<vil_image_view<double> > &flows);


/// Compute the axis-aligned bounding box contain all flow vectors destinations.
/// \note this ignores all pixels with invalid flow
vgl_box_2d<double> flow_destination_bounds(const vil_image_view<double> &flow);


/// Crops the images and flows to what is needed from the crop region
/// (deep copies the crop to remove the full frames from being stored in memory)
/// \param flows the optical flow for each frame
/// \param frames the images
/// \param scale_factor the amount of scaling used in the flows
/// \param margin an expansion tolerance on the cropping region
void crop_frames_and_flows(std::vector<vil_image_view<double> > &flows,
                           std::vector<vil_image_view<double> > &frames,
                           const int scale_factor,
                           const int margin);


/// Crops the image boxes and flows to what is needed from the crop region
/// This function is the same as \a crop_frames_and_flows except that instead
/// of croping the frames it produces the cropped bounding boxes that can
/// be used to crop the frames later.
/// \param flows the optical flow for each frame
/// \param boxes the bounding boxes corresponding to image sizes
/// \param scale_factor the amount of scaling used in the flows
/// \param margin an expansion tolerance on the cropping region
void crop_boxes_and_flows(std::vector<vil_image_view<double> > &flows,
                          std::vector<vgl_box_2d<int> > &boxes,
                          const int scale_factor,
                          const int margin);


/// Builds DBW ops from flows
/// \param flows the input optical flow for each image
/// \param frames the input original images
/// \param warps the output warp DBW ops
/// \param scale_factor the amount of scaling for the super res
/// \param sensor_sigma the camera blur estimate
/// \param down_scaling toggle between simple down sampling (false) and down scaling (true)
/// \param bicubic_warping toggle between bilinear (false) and bicubic (true)
void create_warps_from_flows(const std::vector<vil_image_view<double> > &flows,
                             const std::vector<vil_image_view<double> > &frames,
                             std::vector<adjoint_image_ops_func<double> > &warps,
                             const int scale_factor,
                             const double sensor_sigma,
                             const bool down_scaling,
                             const bool bicubic_warping);

} // end namespace super3d

#endif // super_res_utils_h_

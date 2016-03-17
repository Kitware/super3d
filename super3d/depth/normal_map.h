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

#ifndef normal_map_h_
#define normal_map_h_

#include "depth_config.h"

#include <vil/vil_image_view.h>
#include <vpgl/vpgl_perspective_camera.h>
#include <vgl/vgl_homg_point_3d.h>


namespace super3d
{

/// Compute world normals at each pixel using a 4-connected neighborhood
/// This version uses weighting of normals by triangle area.  These weights
/// come out automatically from the cross product, as a result this algorithm
/// is very simple and fast to compute.
/// \param camera The camera to use for back projection
/// \param depth_map The depths at each pixel
/// \retval normal_map The computed normal map (a 3-plane image)
SUPER3D_DEPTH_EXPORT
void
depth_map_to_normal_map(const vpgl_perspective_camera<double>& camera,
                        const vil_image_view<double>& depth_map,
                        vil_image_view<double>& normal_map);


/// Compute world normals at each pixel using a 4-connected neighborhood
/// In this function, normals at adjacent triangles are weighted inversely
/// by the length of their sides.  See
/// Nelson Max, "Weights for computing vertex normals from facet normals",
/// Journal of Graphics Tools archive Vol. 4 Issue 2, March 1999 Pg. 1-6
/// \param camera The camera to use for back projection
/// \param depth_map The depths at each pixel
/// \retval normal_map The computed normal map (a 3-plane image)
/// \note currently assumes a zero skew and square pixel camera
SUPER3D_DEPTH_EXPORT
void
depth_map_to_normal_map_inv_len(const vpgl_perspective_camera<double>& camera,
                                const vil_image_view<double>& depth_map,
                                vil_image_view<double>& normal_map);


/// Compute world point location at each pixel by back projecting to depth
/// \param camera The camera to use for back projection
/// \param depth_map The depths at each pixel
/// \retval location_map The 3D world pixel locations (a 3-plane image)
SUPER3D_DEPTH_EXPORT
void
depth_map_to_location_map(const vpgl_perspective_camera<double>& camera,
                          const vil_image_view<double>& depth_map,
                          vil_image_view<double>& location_map);


/// Compute world normals at each pixel using a 4-connected neighborhood.
/// In this function, normals at adjacent triangles are weighted inversely
/// by the length of their sides.  See
/// Nelson Max, "Weights for computing vertex normals from facet normals",
/// Journal of Graphics Tools archive Vol. 4 Issue 2, March 1999 Pg. 1-6
/// \param location_map The 3D world pixel locations (a 3-plane image)
/// \retval normal_map The computed normal map (a 3-plane image)
SUPER3D_DEPTH_EXPORT
void
location_map_to_normal_map(const vil_image_view<double>& location_map,
                           vil_image_view<double>& normal_map);


/// Compute dot product (cos of angle) between rays to center and normals.
/// At each pixel the normalized ray to \a center from world location (i,j)
/// is dotted with the normal at that location.
/// \param center The camera or light center at which rays point
/// \param location_map The 3D world pixel locations (a 3-plane image)
/// \param normal_map The normal map (a 3-plane image)
/// \retval angle_map The map of dot products between normals and rays
SUPER3D_DEPTH_EXPORT
void
viewing_angle_map(const vgl_homg_point_3d<double>& center,
                  const vil_image_view<double>& location_map,
                  const vil_image_view<double>& normal_map,
                  vil_image_view<double>& angle_map);


/// Compute dot product (cos of angle) between two normal maps.
/// \param normal_map1 The first normal map (a 3-plane image)
/// \param normal_map2 The second normal map (a 3-plane image)
/// \retval angle_map The map of dot products between normals and rays
SUPER3D_DEPTH_EXPORT
void
dot_product_map(const vil_image_view<double>& normal_map1,
                const vil_image_view<double>& normal_map2,
                vil_image_view<double>& angle_map);


/// Rotate all normals in the normal map by rotatation matrix \a R
/// \param normal_map is the normal map to rotate
/// \param R is the rotation matrix
/// \retval rotated is the rotated normal map
/// \note \a rotated can be set to \a normal_map to rotate in place
SUPER3D_DEPTH_EXPORT
void rotate_normal_map(const vil_image_view<double>& normal_map,
                       const vgl_rotation_3d<double>& R,
                       vil_image_view<double>& rotated);


/// Scale the normal map to byte range using the typical conventions
/// X maps (-1.0, 1.0) to (0, 255)
/// Y maps (1.0, -1.0) to (0, 255)
/// Z maps (0.0, 1.0) to (0, 255)
/// Normals with negative Z map to (0,0,0)
/// \param normal_map is the input normal map
/// \retval dest is the byte mapped destination image
SUPER3D_DEPTH_EXPORT
void byte_normal_map(const vil_image_view<double>& normal_map,
                     vil_image_view<vxl_byte>& dest);

} // end namespace super3d

#endif // normal_map_h_

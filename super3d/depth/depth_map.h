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

#ifndef depth_map_h_
#define depth_map_h_

#include "depth_config.h"

#include <super3d/imesh/imesh_mesh.h>

#include <vil/vil_image_view.h>
#include <vpgl/vpgl_perspective_camera.h>


namespace super3d
{

class world_space;


/// Compute the average of the finite values in the image.
/// \param img The image. It may contain some infinite values.
/// \return The average of finite depths
SUPER3D_DEPTH_EXPORT
double average_finite(const vil_image_view<double>& img);


/// Fill in infinite values in the image with constant \a value.
/// \retval img The image containing some infinite values
///             (modified in place)
/// \param value The value used to replace infinite values.
SUPER3D_DEPTH_EXPORT
void fill_infinite(vil_image_view<double>& img, double value);


/// Copy finite values from \a src into \a dest
/// For infinite pixels in \a src, the matching pixels in \a dst are unchanaged
/// \retval src The source image containing finite and infinite values
/// \param dest The destination image.
SUPER3D_DEPTH_EXPORT
void copy_finite(const vil_image_view<double>& src,
                 vil_image_view<double>& dest);


/// Fill in infinite values in \a src by diffusing from adjacent finite values.
/// The output \a dest is equal to \a src at all finite values of \a src
/// and has diffused finite values where ever \a src in infinite.
/// If \a dest is not initialized it will be filled with the average finite
/// depth computed from \a src.
/// \param src The source image containing finite and infinite values.
/// \retval dest The destination image (modified in place).
/// \param iterations The number of diffusion iterations
SUPER3D_DEPTH_EXPORT
void diffuse_into_infinite(const vil_image_view<double>& src,
                           vil_image_view<double>& dest,
                           unsigned iterations = 10);


/// Fill in missing (infinite) depths in the depth map.
/// Initially missing depths are filled with the average of finite depths.
/// Then the know depths are defused into the unknown depths using
/// a multiscale application of the heat equation.
/// \retval depth_map The depth map containing some infinite values
///                   (modified in place)
/// \param levels The number of 2x pyramid levels to use
/// \param iterations The number of iterations to use at each level
SUPER3D_DEPTH_EXPORT
void fill_missing_depths(vil_image_view<double>& depth_map,
                         unsigned levels = 7,
                         unsigned iterations = 10);


/// Create dense mesh vertices by back projecting pixels to their given depth.
/// \param camera The camera to use for back projection
/// \param depth_map The depths at which to project each pixel
/// \retval index_image Image of indices into the array of mesh vertices
///                     The index is unsigned(-1) for infinite depths.
/// \return An array of mesh vertices covering all finite depth pixels
SUPER3D_DEPTH_EXPORT
std::unique_ptr<imesh_vertex_array<3> >
depth_map_to_vertices(const vpgl_perspective_camera<double>& camera,
                      const vil_image_view<double>& depth_map,
                      vil_image_view<unsigned>& index_image);


/// Create dense mesh vertices from a height map image (e.g. GeoTIFF LIDAR)
/// \param height_map The heights at which to place each pixel
/// \retval index_image Image of indices into the array of mesh vertices
///                     The index is unsigned(-1) for non-finite heights.
/// \param z_scale Amount to scale the height values
/// \param x_scale Amount to scale the unit horizontal distance between pixels
/// \param y_scale Amount to scale the unit vertical distance between pixels
/// \return An array of mesh vertices covering all finite height pixels
SUPER3D_DEPTH_EXPORT
std::unique_ptr<imesh_vertex_array<3> >
height_map_to_vertices(const vil_image_view<double>& height_map,
                       vil_image_view<unsigned>& index_image,
                       double z_scale = 1.0,
                       double x_scale = 1.0,
                       double y_scale = 1.0);


/// Create a triangulation from an image of vertex indices
/// \param index_image Image of indices into an array of mesh vertices
///                    The index is unsigned(-1) for missing vertices
/// \return An array of mesh faces (triangles) using the indices
SUPER3D_DEPTH_EXPORT
std::unique_ptr<imesh_regular_face_array<3> >
triangulate_index_image(const vil_image_view<unsigned>& index_image);


/// Create a dense mesh by back projecting pixels to their given depth.
/// \param camera The camera to use for back projection
/// \param depth_map The depths at which to project each pixel
/// \return A dense mesh with one vertex per finite depth pixel
SUPER3D_DEPTH_EXPORT
imesh_mesh
depth_map_to_mesh(const vpgl_perspective_camera<double>& camera,
                  const vil_image_view<double>& depth_map);


/// Create a dense mesh from a height map image (e.g. GeoTIFF LIDAR)
/// \param height_map The heights at which to place each pixel
/// \param z_scale Amount to scale the height values
/// \param x_scale Amount to scale the unit horizontal distance between pixels
/// \param y_scale Amount to scale the unit vertical distance between pixels
/// \return A dense mesh with one vertex per finite height pixel
SUPER3D_DEPTH_EXPORT
imesh_mesh
height_map_to_mesh(const vil_image_view<double>& height_map,
                   double z_scale = 1.0,
                   double x_scale = 1.0,
                   double y_scale = 1.0);


/// Convert a depth map into a height map
SUPER3D_DEPTH_EXPORT
void depth_map_to_height_map(const vpgl_perspective_camera<double>& camera,
                             const vil_image_view<double>& depth_map,
                                   vil_image_view<double>& height_map);


/// Compute the minimum and maximum value ignoring infinite and NaN values.
SUPER3D_DEPTH_EXPORT
void finite_value_range(const vil_image_view<double>& img,
                        double& min_value, double& max_value);


SUPER3D_DEPTH_EXPORT
void save_depth(const vil_image_view<double> &depth, const char *filename);


SUPER3D_DEPTH_EXPORT
void load_depth(vil_image_view<double> &depth, const char *filename);


/// Compares an image (depth) with its ground truth then reports the average
/// error between them.  It also reports the percentage of pixels
/// that by more than a threshold.
SUPER3D_DEPTH_EXPORT
void score_vs_gt(const vil_image_view<double> &depth,
                 const vil_image_view<double> &gt,
                 double threshold);

#ifdef HAVE_VTK

/// Writes a vtp (vtk) from a depth image
/// \param filename .vtp file to write to
/// \param depth the depth image
/// \param ref reference image that the depth was computed by, used for texturing
/// \param cam perspective camera associated with the depth map
/// \param world_space the world space class used in the depth estimation
SUPER3D_DEPTH_EXPORT
void save_depth_to_vtp(const char *filename,
                       const vil_image_view<double> &depth,
                       const vil_image_view<double> &ref,
                       const vpgl_perspective_camera<double> &cam,
                       world_space *ws);

#endif

SUPER3D_DEPTH_EXPORT
void write_points_to_vtp(std::vector<vnl_double_3> &points, const char *filename);

} // end namespace super3d

#endif // depth_map_h_

/*ckwg +29
 * Copyright 2010 by Kitware, Inc.
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

#include "render.h"


#include <vnl/vnl_inverse.h>
#include <imesh/algo/imesh_project.h>


/// Compute the texture coordinates for \a mesh by projected with \a camera
/// \retval mesh The mesh for computing texture coordinates (modified in place)
/// \param camera The camera projecting the vertices to generate texture coords.
/// \param width The width of the texture image in pixels
/// \param height The height of the texture image in pixels
/// \note \a width and \a height are used to normalize the coordinates to [0,1]
void compute_tex_coords(imesh_mesh& mesh,
                        const vpgl_perspective_camera<double>& camera,
                        unsigned width,
                        unsigned height)
{
  // compute the texture coordinates
  std::vector<vgl_point_2d<double> > tex_coords;
  imesh_project_verts(mesh.vertices<3>(), camera, tex_coords);

  // map texture coordinates to the unit square and flip vertically
  for (int i=0; i<tex_coords.size(); ++i)
  {
    tex_coords[i].x() /= width;
    tex_coords[i].y() = 1.0 - tex_coords[i].y() / height;
  }

  // update the mesh texture coordinates
  mesh.set_tex_coords(tex_coords);
}


/// Offset the principal point in \a camera to account for a cropped output
/// \param camera The camera to offset
/// \param offset A vector specifying the offset of the crop
/// \returns A camera with modified principal point
vpgl_perspective_camera<double>
offset_camera(const vpgl_perspective_camera<double>& camera,
              const vgl_vector_2d<double>& offset)
{
  vpgl_perspective_camera<double> cam(camera);
  vpgl_calibration_matrix<double> K = camera.get_calibration();
  K.set_principal_point(K.principal_point()-offset);
  cam.set_calibration(K);
  return cam;
}


/// Apply \c offset_camera to each camera in a camera sequence
/// \param cameras The input camera sequence
/// \param offset A vector specifying the offset of the crop
/// \returns A modified sequence of cameras
std::map<int, vpgl_perspective_camera<double> >
offset_cameras(const std::map<int, vpgl_perspective_camera<double> >& cameras,
               const vgl_vector_2d<double>& offset)
{
  typedef std::map<int, vpgl_perspective_camera<double> >::const_iterator Cam_itr;
  std::map<int, vpgl_perspective_camera<double> > crop_cameras;
  for (Cam_itr citr = cameras.begin(); citr != cameras.end(); ++citr)
  {
    crop_cameras[citr->first] = offset_camera(citr->second, offset);
  }
  return crop_cameras;
}


/// Compute the plane plus parallax model for two cameras.
/// Produces homography \a H and epipole \a e such that a point p = w*(u,v,1),
/// in the image of \a src_camera maps to H*p + e in the image of \a dst_camera.
/// Here w is the projective depth relative to \a src_camera and \a H is the
/// homography induced by the plane at inifinity.
void compute_plane_parallax(const vpgl_proj_camera<double>& src_camera,
                            const vpgl_proj_camera<double>& dst_camera,
                            vnl_matrix_fixed<double,3,3>& H,
                            vnl_vector_fixed<double,3>& e)
{
  // if src_camera matrix is P1 = [M1 | t1] and dst_camera is P2 = [M2 | t2]
  // the H = M2*(M1^-1) and e = -H*t1 + t2
  const vnl_matrix_fixed<double,3,4>& P1 = src_camera.get_matrix();
  const vnl_matrix_fixed<double,3,4>& P2 = dst_camera.get_matrix();
  H = P2.extract(3,3);
  H *= vnl_inverse(P1.extract(3,3));

  e = P2.get_column(3);
  e -= H*P1.get_column(3);
}

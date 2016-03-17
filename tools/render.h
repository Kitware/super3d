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

#ifndef render_h_
#define render_h_


#include <vnl/vnl_matrix_fixed.h>
#include <vnl/vnl_vector_fixed.h>
#include <vnl/vnl_double_3.h>
#include <vil/vil_image_view.h>
#include <vil/vil_bilin_interp.h>
#include <imesh/imesh_mesh.h>
#include <vpgl/vpgl_perspective_camera.h>



/// Compute the texture coordinates for \a mesh by projected with \a camera
/// \retval mesh The mesh for computing texture coordinates (modified in place)
/// \param camera The camera projecting the vertices to generate texture coords.
/// \param width The width of the texture image in pixels
/// \param height The height of the texture image in pixels
/// \note \a width and \a height are used to normalize the coordinates to [0,1]
void compute_tex_coords(imesh_mesh& mesh,
                        const vpgl_perspective_camera<double>& camera,
                        unsigned width,
                        unsigned height);


/// Offset the principal point in \a camera to account for a cropped output
/// \param camera The camera to offset
/// \param offset A vector specifying the offset of the crop
/// \returns A camera with modified principal point
vpgl_perspective_camera<double>
offset_camera(const vpgl_perspective_camera<double>& camera,
              const vgl_vector_2d<double>& offset);


/// Apply \c offset_camera to each camera in a camera sequence
/// \param cameras The input camera sequence
/// \param offset A vector specifying the offset of the crop
/// \returns A modified sequence of cameras
std::map<int, vpgl_perspective_camera<double> >
offset_cameras(const std::map<int, vpgl_perspective_camera<double> >& cameras,
               const vgl_vector_2d<double>& offset);


/// Compute the plane plus parallax model for two cameras.
/// Produces homography \a H and epipole \a e such that a point p = w*(u,v,1),
/// in the image of \a src_camera maps to H*p + e in the image of \a dst_camera.
/// Here w is the projective depth relative to \a src_camera and \a H is the
/// homography induced by the plane at inifinity.
void compute_plane_parallax(const vpgl_proj_camera<double>& src_camera,
                            const vpgl_proj_camera<double>& dst_camera,
                            vnl_matrix_fixed<double,3,3>& H,
                            vnl_vector_fixed<double,3>& e);


/// Warp \a other_view to the current view given \a depth_map in the current
/// view, homography induced by the plane at infinity \a H, and epipole \a e
template <typename T>
void
parallax_warp(const vil_image_view<double>& depth_map,
              const vnl_matrix_fixed<double,3,3>& H,
              const vnl_vector_fixed<double,3>& e,
              const vil_image_view<T>& other_view,
              vil_image_view<T>& rendered)
{
  const unsigned ni = depth_map.ni();
  const unsigned nj = depth_map.nj();
  const unsigned np = other_view.nplanes();
  rendered.set_size(ni,nj,np);

  vcl_ptrdiff_t istepD=depth_map.istep(), jstepD=depth_map.jstep();
  vcl_ptrdiff_t istepR=rendered.istep(),  jstepR=rendered.jstep();
  vcl_ptrdiff_t pstepR=rendered.planestep();

  const double* rowD = depth_map.top_left_ptr();
  T*            rowR = rendered.top_left_ptr();
  for (unsigned j=0; j<nj; ++j, rowD+=jstepD, rowR+=jstepR)
  {
    const double* pixelD = rowD;
    T*            pixelR = rowR;
    for (unsigned i=0; i<ni; ++i, pixelD+=istepD, pixelR+=istepR)
    {
      const double& depth = *pixelD;
      vnl_double_3 p(i,j,1.0);
      vnl_double_3 q = H*p;
      if (vnl_math_isfinite(depth))
      {
        q *= depth;
        q += e;
      }
      q /= q[2];

      if (q[0] < 0.0 || q[0] > other_view.ni()-1 ||
          q[1] < 0.0 || q[1] > other_view.nj()-1)
      {
        T* planeR = pixelR;
        for (unsigned p=0; p<np; ++p, planeR+=pstepR)
        {
          *planeR = T(0);
        }
      }
      else
      {
        T* planeR = pixelR;
        for (unsigned p=0; p<np; ++p, planeR+=pstepR)
        {
          *planeR = vil_bilin_interp(other_view,q[0],q[1],p);
        }
      }
    }
  }
}


#endif // render_h_

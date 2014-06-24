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

#include "normal_map.h"

#include <vnl/vnl_double_3.h>
#include <vnl/vnl_inverse.h>


namespace super3d
{

/// Compute world normals at each pixel using a 4-connected neighborhood
/// This version uses weighting of normals by triangle area.  These weights
/// come out automatically from the cross product, as a result this algorithm
/// is very simple and fast to compute.
/// \param camera The camera to use for back projection
/// \param depth_map The depths at each pixel
/// \retval normal_map The computed normal map (a 3-plane image)
void
depth_map_to_normal_map(const vpgl_perspective_camera<double>& camera,
                        const vil_image_view<double>& depth_map,
                        vil_image_view<double>& normal_map)
{
  const unsigned ni = depth_map.ni();
  const unsigned nj = depth_map.nj();
  assert(depth_map.nplanes() == 1);

  normal_map.set_size(ni,nj,3);
  normal_map.fill(0.0);

  vnl_matrix_fixed<double,3,3> M_t = camera.get_matrix().extract(3,3).transpose();

  const vcl_ptrdiff_t istepD=depth_map.istep(),   jstepD=depth_map.jstep();
  const vcl_ptrdiff_t istepN=normal_map.istep(),  jstepN=normal_map.jstep();
  const vcl_ptrdiff_t pstepN=normal_map.planestep();
  const vcl_ptrdiff_t pstep2N=2*pstepN;

  // currently ignores the boundary cases
  const double*   rowD = depth_map.top_left_ptr() + jstepD;
  double*         rowN = normal_map.top_left_ptr() + jstepN;
  for (unsigned j=1; j<nj-1; ++j, rowD += jstepD, rowN += jstepN)
  {
    const double* pixelD = rowD + istepD;
    double*       pixelN = rowN + istepN;
    for (unsigned i=1; i<ni-1; ++i, pixelD+=istepD, pixelN+=istepN)
    {
      const double& d1 = *(pixelD - jstepD);
      const double& d2 = *(pixelD + istepD);
      const double& d3 = *(pixelD + jstepD);
      const double& d4 = *(pixelD - istepD);
      const double d12 = d1*d2;
      const double d14 = d1*d4;
      const double d23 = d2*d3;
      const double d34 = d3*d4;
      const double ax =  d12 - d14 + d23 - d34;
      const double ay = -d12 - d14 + d23 + d34;
      const double a0 = -d12 - d14 - d23 - d34;
      vnl_double_3 normal = M_t * vnl_double_3(ax, ay, a0 - ax*i - ay*j);
      normal.normalize();
      *pixelN = normal[0];
      *(pixelN + pstepN) = normal[1];
      *(pixelN + pstep2N) = normal[2];
    }
  }
}


/// Compute world normals at each pixel using a 4-connected neighborhood
/// In this function, normals at adjacent triangles are weighted inversely
/// by the length of their sides.  See
/// Nelson Max, "Weights for computing vertex normals from facet normals",
/// Journal of Graphics Tools archive Vol. 4 Issue 2, March 1999 Pg. 1-6
/// \param camera The camera to use for back projection
/// \param depth_map The depths at each pixel
/// \retval normal_map The computed normal map (a 3-plane image)
/// \note currently assumes a zero skew and square pixel camera
void
depth_map_to_normal_map_inv_len(const vpgl_perspective_camera<double>& camera,
                                const vil_image_view<double>& depth_map,
                                vil_image_view<double>& normal_map)
{
  const unsigned ni = depth_map.ni();
  const unsigned nj = depth_map.nj();
  assert(depth_map.nplanes() == 1);

  normal_map.set_size(ni,nj,3);
  normal_map.fill(0.0);

  const vpgl_calibration_matrix<double>& K = camera.get_calibration();
  const double pp_x =  K.principal_point().x();
  const double pp_y =  K.principal_point().y();
  const double f = K.focal_length() * K.x_scale();
  assert(K.skew() == 0);
  assert(K.x_scale() / K.y_scale() == 1.0);
  vnl_matrix_fixed<double,3,3> M_t = camera.get_matrix().extract(3,3).transpose();

  const vcl_ptrdiff_t istepD=depth_map.istep(),   jstepD=depth_map.jstep();
  const vcl_ptrdiff_t istepN=normal_map.istep(),  jstepN=normal_map.jstep();
  const vcl_ptrdiff_t pstepN=normal_map.planestep();
  const vcl_ptrdiff_t pstep2N=2*pstepN;

  // currently ignores the boundary cases
  const double*   rowD = depth_map.top_left_ptr() + jstepD;
  double*         rowN = normal_map.top_left_ptr() + jstepN;
  for (unsigned j=1; j<nj-1; ++j, rowD += jstepD, rowN += jstepN)
  {
    const double* pixelD = rowD + istepD;
    double*       pixelN = rowN + istepN;
    const double y = (j - pp_y) / f;
    for (unsigned i=1; i<ni-1; ++i, pixelD+=istepD, pixelN+=istepN)
    {
      const double x = (i - pp_x) / f;
      const double& d0 = *pixelD;
      const double& d1 = *(pixelD - jstepD);
      const double& d2 = *(pixelD + istepD);
      const double& d3 = *(pixelD + jstepD);
      const double& d4 = *(pixelD - istepD);
      const double d1_0 = (d1-d0);
      const double d2_0 = (d2-d0);
      const double d3_0 = (d3-d0);
      const double d4_0 = (d4-d0);
      const double x2_y2_1 = x*x + y*y + 1;
      const double len1 = d1_0*d1_0*x2_y2_1 - (2*d1*d1_0*y - d1*d1/f)/f;
      const double len2 = d2_0*d2_0*x2_y2_1 + (2*d2*d2_0*x + d2*d2/f)/f;
      const double len3 = d3_0*d3_0*x2_y2_1 + (2*d3*d3_0*y + d3*d3/f)/f;
      const double len4 = d4_0*d4_0*x2_y2_1 - (2*d4*d4_0*x - d4*d4/f)/f;
      const double w12 = 1.0 / (len1*len2);
      const double w23 = 1.0 / (len2*len3);
      const double w34 = 1.0 / (len3*len4);
      const double w41 = 1.0 / (len4*len1);
      const double ax =  d2_0*d1*w12 + d2_0*d3*w23 - d4_0*d3*w34 - d4_0*d1*w41;
      const double ay = -d1_0*d2*w12 + d3_0*d2*w23 + d3_0*d4*w34 - d1_0*d4*w41;
      const double a0 = -d1*d2*w12   - d2*d3*w23   - d3*d4*w34   - d4*d1*w41;
      vnl_double_3 normal = M_t * vnl_double_3(ax, ay, a0 - ax*i - ay*j);
      normal.normalize();
      *pixelN = normal[0];
      *(pixelN + pstepN) = normal[1];
      *(pixelN + pstep2N) = normal[2];
    }
  }
}


/// Compute world point location at each pixel by back projecting to depth
/// \param camera The camera to use for back projection
/// \param depth_map The depths at each pixel
/// \retval location_map The 3D world pixel locations (a 3-plane image)
void
depth_map_to_location_map(const vpgl_perspective_camera<double>& camera,
                          const vil_image_view<double>& depth_map,
                          vil_image_view<double>& location_map)
{
  const unsigned ni = depth_map.ni();
  const unsigned nj = depth_map.nj();
  assert(depth_map.nplanes() == 1);

  location_map.set_size(ni,nj,3);

  vnl_matrix_fixed<double,3,3> M_inv = vnl_inverse(camera.get_matrix().extract(3,3));
  vnl_double_3 p4 = camera.get_matrix().get_column(3);

  const vcl_ptrdiff_t istepD=depth_map.istep(),    jstepD=depth_map.jstep();
  const vcl_ptrdiff_t istepL=location_map.istep(), jstepL=location_map.jstep();
  const vcl_ptrdiff_t pstepL=location_map.planestep();
  const vcl_ptrdiff_t pstep2L=2*pstepL;

  const double*   rowD = depth_map.top_left_ptr();
  double*         rowL = location_map.top_left_ptr();
  for (unsigned j=0; j<nj; ++j, rowD += jstepD, rowL += jstepL)
  {
    const double* pixelD = rowD;
    double*       pixelL = rowL;
    for (unsigned i=0; i<ni; ++i, pixelD+=istepD, pixelL+=istepL)
    {
      vnl_double_3 pt = M_inv * ((*pixelD) * vnl_double_3(i,j,1) - p4);
      *pixelL = pt[0];
      *(pixelL + pstepL) = pt[1];
      *(pixelL + pstep2L) = pt[2];
    }
  }
}


/// Compute world normals at each pixel using a 4-connected neighborhood.
/// In this function, normals at adjacent triangles are weighted inversely
/// by the length of their sides.  See
/// Nelson Max, "Weights for computing vertex normals from facet normals",
/// Journal of Graphics Tools archive Vol. 4 Issue 2, March 1999 Pg. 1-6
/// \param location_map The 3D world pixel locations (a 3-plane image)
/// \retval normal_map The computed normal map (a 3-plane image)
void
location_map_to_normal_map(const vil_image_view<double>& location_map,
                           vil_image_view<double>& normal_map)
{
  const unsigned ni = location_map.ni();
  const unsigned nj = location_map.nj();
  assert(location_map.nplanes() == 3);

  normal_map.set_size(ni,nj,3);
  normal_map.fill(0.0);

  const vcl_ptrdiff_t istepL=location_map.istep(), jstepL=location_map.jstep();
  const vcl_ptrdiff_t pstepL=location_map.planestep(), pstep2L=2*pstepL;
  const vcl_ptrdiff_t istepN=normal_map.istep(),  jstepN=normal_map.jstep();
  const vcl_ptrdiff_t pstepN=normal_map.planestep(), pstep2N=2*pstepN;
  typedef vgl_vector_3d<double> vec_t;
  typedef vgl_point_3d<double> pnt_t;

  // currently ignores the boundary cases
  const double*   rowL = location_map.top_left_ptr() + jstepL;
  double*         rowN = normal_map.top_left_ptr() + jstepN;
  for (unsigned j=1; j<nj-1; ++j, rowL += jstepL, rowN += jstepN)
  {
    const double* pixelL = rowL + istepL;
    double*       pixelN = rowN + istepN;
    for (unsigned i=1; i<ni-1; ++i, pixelL+=istepL, pixelN+=istepN)
    {
      const pnt_t pt0(*pixelL, *(pixelL+pstepL), *(pixelL+pstep2L));
      const double* L = pixelL - jstepL;
      const vec_t v1 = pnt_t(*L, *(L+pstepL), *(L+pstep2L)) - pt0;
      L = pixelL + istepL;
      const vec_t v2 = pnt_t(*L, *(L+pstepL), *(L+pstep2L)) - pt0;
      L = pixelL + jstepL;
      const vec_t v3 = pnt_t(*L, *(L+pstepL), *(L+pstep2L)) - pt0;
      L = pixelL - istepL;
      const vec_t v4 = pnt_t(*L, *(L+pstepL), *(L+pstep2L)) - pt0;

      vec_t normal(0,0,0);
      normal += cross_product(v2, v1) / (v2.sqr_length() * v1.sqr_length());
      normal += cross_product(v3, v2) / (v3.sqr_length() * v2.sqr_length());
      normal += cross_product(v4, v3) / (v4.sqr_length() * v3.sqr_length());
      normal += cross_product(v1, v4) / (v1.sqr_length() * v4.sqr_length());
      normalize(normal);
      *pixelN = normal.x();
      *(pixelN + pstepN) = normal.y();
      *(pixelN + pstep2N) = normal.z();
    }
  }
}


/// Compute dot product (cos of angle) between rays to center and normals.
/// At each pixel the normalized ray to \a center from world location (i,j)
/// is dotted with the normal at that location.
/// \param center The camera or light center at which rays point
/// \param location_map The 3D world pixel locations (a 3-plane image)
/// \param normal_map The computed normal map (a 3-plane image)
/// \retval angle_map The map of dot products between normals and rays
void
viewing_angle_map(const vgl_homg_point_3d<double>& center,
                  const vil_image_view<double>& location_map,
                  const vil_image_view<double>& normal_map,
                  vil_image_view<double>& angle_map)
{
  const unsigned ni = location_map.ni();
  const unsigned nj = location_map.nj();
  assert(location_map.nplanes() == 3);
  assert(normal_map.ni() == ni);
  assert(normal_map.nj() == nj);
  assert(normal_map.nplanes() == 3);

  angle_map.set_size(ni,nj);

  const vcl_ptrdiff_t istepL=location_map.istep(), jstepL=location_map.jstep();
  const vcl_ptrdiff_t pstepL=location_map.planestep();
  const vcl_ptrdiff_t pstep2L=2*pstepL;
  const vcl_ptrdiff_t istepN=normal_map.istep(),  jstepN=normal_map.jstep();
  const vcl_ptrdiff_t pstepN=normal_map.planestep();
  const vcl_ptrdiff_t pstep2N=2*pstepN;
  const vcl_ptrdiff_t istepA=angle_map.istep(), jstepA=angle_map.jstep();

  const double*   rowL = location_map.top_left_ptr();
  const double*   rowN = normal_map.top_left_ptr();
  double*         rowA = angle_map.top_left_ptr();
  for (unsigned j=0; j<nj; ++j, rowL += jstepL, rowN += jstepN, rowA += jstepA)
  {
    const double* pixelL = rowL;
    const double* pixelN = rowN;
    double*       pixelA = rowA;
    for (unsigned i=0; i<ni; ++i, pixelL+=istepL, pixelN+=istepN, pixelA+=istepA)
    {
      vnl_double_3 n(*pixelN, *(pixelN+pstepN), *(pixelN+pstep2N));
      vnl_double_3 v(*pixelL, *(pixelL+pstepL), *(pixelL+pstep2L));
      v *= center.w();
      v -= vnl_double_3(center.x(), center.y(), center.z());
      v.normalize();
      *pixelA = dot_product(-v,n);
    }
  }
}


/// Compute dot product (cos of angle) between two normal maps.
/// \param normal_map1 The first normal map (a 3-plane image)
/// \param normal_map2 The second normal map (a 3-plane image)
/// \retval angle_map The map of dot products between normals and rays
void
dot_product_map(const vil_image_view<double>& normal_map1,
                const vil_image_view<double>& normal_map2,
                vil_image_view<double>& angle_map)
{
  const unsigned ni = normal_map1.ni();
  const unsigned nj = normal_map1.nj();
  assert(normal_map1.nplanes() == 3);
  assert(normal_map2.ni() == ni);
  assert(normal_map2.nj() == nj);
  assert(normal_map2.nplanes() == 3);

  angle_map.set_size(ni,nj);

  const vcl_ptrdiff_t istepB=normal_map1.istep(), jstepB=normal_map1.jstep();
  const vcl_ptrdiff_t pstepB=normal_map1.planestep();
  const vcl_ptrdiff_t pstep2B=2*pstepB;
  const vcl_ptrdiff_t istepC=normal_map2.istep(),  jstepC=normal_map2.jstep();
  const vcl_ptrdiff_t pstepC=normal_map2.planestep();
  const vcl_ptrdiff_t pstep2C=2*pstepC;
  const vcl_ptrdiff_t istepA=angle_map.istep(), jstepA=angle_map.jstep();

  const double*   rowB = normal_map1.top_left_ptr();
  const double*   rowC = normal_map2.top_left_ptr();
  double*         rowA = angle_map.top_left_ptr();
  for (unsigned j=0; j<nj; ++j, rowB += jstepB, rowC += jstepC, rowA += jstepA)
  {
    const double* pixelB = rowB;
    const double* pixelC = rowC;
    double*       pixelA = rowA;
    for (unsigned i=0; i<ni; ++i, pixelB+=istepB, pixelC+=istepC, pixelA+=istepA)
    {
      vnl_double_3 n1(*pixelB, *(pixelB+pstepB), *(pixelB+pstep2B));
      vnl_double_3 n2(*pixelC, *(pixelC+pstepC), *(pixelC+pstep2C));
      *pixelA = dot_product(n1,n2);
    }
  }
}


/// Rotate all normals in the normal map by rotatation matrix \a R
/// \param normal_map is the normal map to rotate
/// \param R is the rotation matrix
/// \retval rotated is the rotated normal map
/// \note \a rotated can be set to \a normal_map to rotate in place
void rotate_normal_map(const vil_image_view<double>& normal_map,
                       const vgl_rotation_3d<double>& R,
                       vil_image_view<double>& rotated)
{
  const unsigned ni = normal_map.ni();
  const unsigned nj = normal_map.nj();
  assert(normal_map.nplanes() == 3);
  rotated.set_size(ni, nj, 3);
  for(unsigned j=0; j<nj; ++j)
  {
    for(unsigned i=0; i<ni; ++i)
    {
      vnl_double_3 n(normal_map(i,j,0), normal_map(i,j,1), normal_map(i,j,2));
      n = R * n;
      rotated(i,j,0) = n[0];
      rotated(i,j,1) = n[1];
      rotated(i,j,2) = n[2];
    }
  }
}


/// Scale the normal map to byte range using the typical conventions
/// X maps (-1.0, 1.0) to (0, 255)
/// Y maps (1.0, -1.0) to (0, 255)
/// Z maps (0.0, 1.0) to (0, 255)
/// Normals with negative Z map to (0,0,0)
/// \param normal_map is the input normal map
/// \retval dest is the byte mapped destination image
void byte_normal_map(const vil_image_view<double>& normal_map,
                     vil_image_view<vxl_byte>& dest)
{
  const unsigned ni = normal_map.ni();
  const unsigned nj = normal_map.nj();
  assert(normal_map.nplanes() == 3);
  dest.set_size(ni, nj, 3);
  dest.fill(0);
  for(unsigned j=0; j<nj; ++j)
  {
    for(unsigned i=0; i<ni; ++i)
    {
      vnl_double_3 n(normal_map(i,j,0), normal_map(i,j,1), normal_map(i,j,2));
      n.normalize();
      if (n[2] <= 0.0)
      {
        dest(i,j,0) = static_cast<vxl_byte>(127.5 * (n[0]+1));
        dest(i,j,1) = static_cast<vxl_byte>(127.5 * (1-n[1]));
        dest(i,j,2) = static_cast<vxl_byte>(255 * -n[2]);
      }
    }
  }
}

} // end namespace super3d

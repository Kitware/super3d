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

#include "world_rectilinear.h"


namespace super3d
{

world_rectilinear::world_rectilinear(const vnl_double_3 &origin,
                                     const vnl_double_3 &dimensions,
                                     unsigned int pixel_width,
                                     unsigned int pixel_height) :
  world_space(pixel_width, pixel_height)
{
  double halfx = 0.5 * dimensions[0];
  double halfy = 0.5 * dimensions[1];
  double halfz = 0.5 * dimensions[2];
  b_corners[0] = origin + vnl_double_3(-halfx, -halfy, -halfz);
  b_corners[1] = origin + vnl_double_3( halfx, -halfy, -halfz);
  b_corners[2] = origin + vnl_double_3( halfx,  halfy, -halfz);
  b_corners[3] = origin + vnl_double_3(-halfx,  halfy, -halfz);
  height = dimensions[2];

  wip.set_fill_unmapped(true);
  wip.set_unmapped_value(-1.0);

  xscale = dimensions[0] / pixel_width;
  yscale = dimensions[1] / pixel_height;

  vcl_cout << "x scale: " << xscale << "\n";
  vcl_cout << "y scale: " << yscale << "\n";
  vcl_cout << "z scale: " << height << "\n";
}


/// returns the corner points of an image slice at depth slice.
/// depth slice is a value between 0 and 1 over the depth range
vcl_vector<vnl_double_3> world_rectilinear::get_slice(double depth_slice) const
{
  vcl_vector<vnl_double_3> slice_corners;
  for (unsigned int i = 0; i < 4; i++)
  {
    vnl_double_3 corner = b_corners[i];
    corner[2] += depth_slice*height;
    slice_corners.push_back(corner);
  }

  return slice_corners;
}


//gradient weighting doesn't work with rectilinear world because there is no
void world_rectilinear::compute_g(const vil_image_view<double> &ref_img,
                                  vil_image_view<double> &g,
                                  double alpha,
                                  double beta)
{
  g.set_size(ref_img.ni(), ref_img.nj(), 1);
  g.fill(1.0);
}


vnl_double_3 world_rectilinear::point_at_depth(unsigned int i, unsigned int j, double d) const
{
  double x = (double)i/(double)ni_;
  double y = (double)j/(double)nj_;

  return vnl_double_3(x*b_corners[2][0] + (1.0-x)*b_corners[0][0],
                      y*b_corners[2][1] + (1.0-y)*b_corners[0][1],
                      b_corners[0][2]   + d*height);
}


vnl_double_3 world_rectilinear::map_normal_w2n(const vnl_double_3 &vec, const vnl_double_3 &loc) const
{
  vnl_double_3 new_vec = vec;
  new_vec[0] /= yscale*height;
  new_vec[1] /= xscale*height;
  new_vec[2] /= xscale*yscale;
  new_vec.normalize();
  return new_vec;
}


vnl_double_3 world_rectilinear::map_normal_n2w(const vnl_double_3 &vec, const vnl_double_3 &loc) const
{
  vnl_double_3 new_vec = vec;
  new_vec[0] *= yscale*height;
  new_vec[1] *= xscale*height;
  new_vec[2] *= xscale*yscale;
  new_vec.normalize();
  return new_vec;
}

} // end namespace super3d

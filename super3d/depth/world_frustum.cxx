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

#include "world_frustum.h"

world_frustum::world_frustum(const vpgl_perspective_camera<double> &cam,
                             double min_depth,
                             double max_depth,
                             unsigned int pixel_width,
                             unsigned int pixel_height)
: world_space(pixel_width, pixel_height),
  n(min_depth), f(max_depth)
{
  const vpgl_calibration_matrix<double> &K = cam.get_calibration();
  px = K.principal_point().x();
  py = K.principal_point().y();
  focal_length = K.focal_length();
  R = cam.get_rotation().as_matrix();
  vgl_point_3d<double> c = cam.get_camera_center();
  cam_center = vnl_double_3(c.x(), c.y(), c.z());
  refcam = cam;

  vnl_double_3 b1 = point_at_depth_on_axis(0.0, 0.0, 0.0);
  vnl_double_3 b2 = point_at_depth_on_axis(ni_, 0.0, 0.0);
  vnl_double_3 t1 = point_at_depth_on_axis(0.0, 0.0, 1.0);
  vnl_double_3 t2 = point_at_depth_on_axis(ni_, 0.0, 1.0);

  double w1 = 0.5 * ((b2 - b1).two_norm() + (t2 - t1).two_norm());
  dscale = fabs(f - n) * (double)ni_ / w1;
  vcl_cout << dscale << "\n";
}

/// returns the corner points of an image slice at depth slice.
/// depth slice is a value between 0 and 1 over the depth range
vcl_vector<vnl_double_3> world_frustum::get_slice(double depth_slice) const
{
  vcl_vector<vnl_double_3> slice;
  slice.push_back(point_at_depth_on_axis(0.0, 0.0, depth_slice));
  slice.push_back(point_at_depth_on_axis((double)ni_, 0.0, depth_slice));
  slice.push_back(point_at_depth_on_axis((double)ni_, (double)nj_, depth_slice));
  slice.push_back(point_at_depth_on_axis(0.0, (double)nj_, depth_slice));
  return slice;
}

vcl_vector<vpgl_perspective_camera<double> >
world_frustum::warp_cams(const vcl_vector<vpgl_perspective_camera<double> > &cameras, int ref_frame) const
{
  vcl_vector<vpgl_perspective_camera<double> > newcams(cameras.size());
  vgl_rotation_3d<double> R_ref = cameras[ref_frame].get_rotation().inverse();
  for (unsigned int i = 0; i < newcams.size(); i++)
  {
    newcams[i] = cameras[i];
    if (i == ref_frame)
    {
      vgl_rotation_3d<double> identity;
      identity.set_identity();
      newcams[i].set_rotation(identity);
      newcams[i].set_translation(vgl_vector_3d<double>(0.0, 0.0, 0.0));
    }
    else
    {
      vgl_rotation_3d<double> R_relative = cameras[i].get_rotation() * R_ref;
      vgl_vector_3d<double> t_relative = R_relative * -cameras[ref_frame].get_translation() + cameras[i].get_translation();
      newcams[i].set_rotation(R_relative);
      newcams[i].set_translation(t_relative);
    }
  }

  return newcams;
}

vnl_double_3 world_frustum::point_at_depth_on_axis(double i, double j, double depth) const
{
  double denomij = (focal_length * (depth * f - depth * n - f));
  vnl_double_3 pt = vnl_double_3( (f * n * (-i + px))/denomij,
                       (f * n * (-j + py))/denomij,
                       (f * n) / (-depth * f + depth * n + f) );
  return pt;
}

vnl_double_3 world_frustum::point_at_depth(unsigned int i, unsigned int j, double depth) const
{
  vnl_double_3 pt3d = point_at_depth_on_axis((double)i, (double)j, depth);
  return cam_center + R.transpose() * pt3d;
}

vnl_double_3 world_frustum::map_normal_w2n(const vnl_double_3 &vec, const vnl_double_3 &loc) const
{
  vnl_double_3 wpt = point_at_depth_on_axis(loc(0), loc(1), loc(2));
  double z_over_focal = wpt(2) / focal_length;
  double term = (wpt(2) * (f - n))/(f*n);
  vnl_double_3 normal(vec(0) * z_over_focal,
                      vec(1) * z_over_focal,
                      vec(0) * wpt(0) * term + vec(1) * wpt(1) * term + vec(2) * wpt(2) * term);
  return normal.normalize();
}

vnl_double_3 world_frustum::map_normal_n2w(const vnl_double_3 &vec, const vnl_double_3 &loc) const
{
  vnl_double_3 wpt = point_at_depth_on_axis(loc(0), loc(1), loc(2));
  double focal_over_z = focal_length / wpt(2);
  vnl_double_3 normal;
  normal(0) = vec(0) * focal_over_z;
  normal(1) = vec(1) * focal_over_z,
  normal(2) = -(normal(0) * wpt(0)) / wpt(2) - (normal(1) * wpt(1)) / wpt(2) + (vec[2]*f*n)/(wpt[2]*wpt[2]*(f-n));
  return normal.normalize();
}

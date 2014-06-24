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

#ifndef WORLD_FRUSTUM
#define WORLD_FRUSTUM

#include "depth_config.h"
#include "world_space.h"


namespace super3d
{

class SUPER3D_DEPTH_EXPORT world_frustum : public world_space
{
public:

  world_frustum(const vpgl_perspective_camera<double> &cam,
                double min_depth,
                double max_depth,
                unsigned int pixel_width,
                unsigned int pixel_height);

  /// returns the corner points of an image slice at depth slice.
  /// depth slice is a value between 0 and 1 over the depth range
  vcl_vector<vnl_double_3> get_slice(double depth_slice) const;

  vcl_vector<vpgl_perspective_camera<double> >
  warp_cams(const vcl_vector<vpgl_perspective_camera<double> > &cameras, int ref_frame) const;

  vnl_double_3 point_at_depth(unsigned int i, unsigned int j, double depth) const;

  vnl_double_3 map_normal_w2n(const vnl_double_3 &vec, const vnl_double_3 &loc) const;
  vnl_double_3 map_normal_n2w(const vnl_double_3 &vec, const vnl_double_3 &loc) const;

  vnl_double_3 point_at_depth_on_axis(double i, double j, double depth) const;

private:

  double dscale;
  double n, f, focal_length, px, py;
  vnl_double_3x3 R;
  vnl_double_3 cam_center;
  vpgl_perspective_camera<double> refcam;
};

} // end namespace super3d

#endif

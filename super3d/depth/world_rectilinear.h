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

#ifndef WORLD_RECTILINEAR_H_
#define WORLD_RECTILINEAR_H_

#include "depth_config.h"
#include "world_space.h"


namespace super3d
{

class SUPER3D_DEPTH_EXPORT world_rectilinear : public world_space
{
public:

  world_rectilinear(const vnl_double_3 &origin,
                    const vnl_double_3 &dimensions,
                    unsigned int pixel_width,
                    unsigned int pixel_height);

  world_rectilinear(const std::string &landmarks_file,
          unsigned int pixel_width,
          unsigned int pixel_height);

  /// returns the corner points of an image slice at depth slice.
  /// depth slice is a value between 0 and 1 over the depth range
  std::vector<vnl_double_3> get_slice(double depth_slice) const;

  void compute_g(const vil_image_view<double> &ref_img,
                 vil_image_view<double> &g,
                 double alpha,
                 double beta,
                 vil_image_view<double> *mask = NULL);

  vnl_double_3 point_at_depth(unsigned int i, unsigned int j, double depth) const;
  vnl_double_3 point_at_depth_on_axis(double i, double j, double depth) const { return point_at_depth((unsigned int)i, (unsigned int)j, depth); }

  double get_height() const { return height; }
  double get_xscale() const { return xscale; }
  double get_yscale() const { return yscale; }

  vnl_double_3 map_normal_w2n(const vnl_double_3 &vec, const vnl_double_3 &loc) const;
  vnl_double_3 map_normal_n2w(const vnl_double_3 &vec, const vnl_double_3 &loc) const;

private:

  vnl_double_3 b_corners[4]; //corners for bottom of world volume
  double height, xscale, yscale;

};

} // end namespace super3d

#endif

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

#include "world_space.h"

#include <vgl/algo/vgl_h_matrix_2d_compute_4point.h>
#include <vil/vil_save.h>
#include <vil/vil_convert.h>
#include <vil/vil_math.h>
#include <vil/algo/vil_sobel_3x3.h>

world_space::world_space(unsigned int pixel_width, unsigned int pixel_height)
  : ni_(pixel_width), nj_(pixel_height)
{
  wip.set_fill_unmapped(true);
  wip.set_unmapped_value(-1.0);
  wip.set_interpolator(vidtk::warp_image_parameters::LINEAR);
}

vcl_vector<vpgl_perspective_camera<double> >
world_space::warp_cams(const vcl_vector<vpgl_perspective_camera<double> > &cameras, int ref_frame) const
{
  return cameras;
}

/// warps image \in to the world volume at depth_slice,
/// uses ni and nj as out's dimensions
void world_space::warp_image_to_depth(const vil_image_view<double> &in,
                                      vil_image_view<double> &out,
                                      const vpgl_perspective_camera<double> &cam,
                                      double depth_slice, int f) const
{
  vcl_vector<vnl_double_3> wpts = this->get_slice(depth_slice);

  vcl_vector<vgl_homg_point_2d<double> > warp_pts;
  vcl_vector<vgl_homg_point_2d<double> > proj_pts;

  warp_pts.push_back(vgl_homg_point_2d<double>(0.0, 0.0, 1.0));
  warp_pts.push_back(vgl_homg_point_2d<double>(ni_, 0.0, 1.0));
  warp_pts.push_back(vgl_homg_point_2d<double>(ni_, nj_, 1.0));
  warp_pts.push_back(vgl_homg_point_2d<double>(0.0, nj_, 1.0));

  for (unsigned int i = 0; i < wpts.size(); i++)
  {
    double u, v;
    cam.project(wpts[i][0], wpts[i][1], wpts[i][2], u, v);
    proj_pts.push_back(vgl_homg_point_2d<double>(u, v, 1));
  }

  vgl_h_matrix_2d<double> H;
  vgl_h_matrix_2d_compute_4point dlt;
  dlt.compute(warp_pts, proj_pts, H);

  out.set_size(ni_, nj_, 1);
  vidtk::warp_image(in, out, vgl_h_matrix_2d<double>(H), wip);

#if 1
  vil_image_view<double> outwrite;
  outwrite.deep_copy(out);
  vil_math_scale_and_offset_values(outwrite, 255.0, 0.0);
  vil_image_view<vxl_byte> to_save;
  vil_convert_cast<double, vxl_byte>(outwrite, to_save);
  char buf[60];
  sprintf(buf, "images/slice%2f_frame%d_%d.png", wpts[0][2], f, wip.interpolator_);
  vil_save(to_save, buf);
#endif
}

void world_space::compute_g(const vil_image_view<double> &ref_img, vil_image_view<double> &g, double alpha, double beta)
{
  g.set_size(ref_img.ni(), ref_img.nj(), 1);

  vil_image_view<double> ref_img_g;
  vil_sobel_3x3(ref_img, ref_img_g);

  vcl_cout << "Computing g weighting.\n";
  for (unsigned int i = 0; i < ref_img_g.ni(); i++)
  {
    for (unsigned int j = 0; j < ref_img_g.nj(); j++)
    {
      double dx = ref_img_g(i,j,0);
      double dy = ref_img_g(i,j,1);
      double mag = sqrt(dx*dx + dy*dy);
      //if (mag > .5) //clamp gradients because of over exposure
      //  mag = .5;
      g(i,j) = exp(-alpha * mag);
    }
  }
}

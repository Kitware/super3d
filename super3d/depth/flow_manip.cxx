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

#include "flow_manip.h"
#include "depth_map.h"
#include "multiscale.h"

#include <vcl_limits.h>

#include <vnl/vnl_inverse.h>
#include <vnl/vnl_double_3x3.h>
#include <vnl/vnl_double_3.h>
#include <vnl/vnl_double_4.h>

#include <super3d/imesh/algo/imesh_project.h>

namespace super3d
{

void compute_flows_from_depth(const vcl_vector<vpgl_perspective_camera<double> > &cameras,
                              const vpgl_perspective_camera<double> &ref_cam,
                              const vil_image_view<double> &depth,
                              vcl_vector<vil_image_view<double> > &flows)
{
  flows.resize(cameras.size());
  const vnl_matrix_fixed<double,3,4> &P1 = ref_cam.get_matrix();
  vnl_matrix_fixed<double, 3, 3> inv1 = vnl_inverse(P1.extract(3,3));
  for (unsigned int c = 0; c < cameras.size(); c++)
  {
    const vnl_matrix_fixed<double,3,4> &P2 = cameras[c].get_matrix();
    vnl_double_3x3 H = P2.extract(3,3);
    H *= inv1;
    vnl_double_3 e = P2.get_column(3);
    e -= H*P1.get_column(3);

    vil_image_view<double> &flow = flows[c];
    flow.set_size(depth.ni(), depth.nj(), 2);
    flow.fill(vcl_numeric_limits<double>::quiet_NaN());

    //Use plane+parallax to find the flow between current image and the ref image
    for (unsigned int j = 0; j < depth.nj(); j++)
    {
      for (unsigned int i = 0; i < depth.ni(); i++)
      {
        vnl_double_3 src_p = vnl_double_3(i,j,1.0);
        vnl_double_3 p = depth(i,j) * src_p;
        vnl_double_3 q = H * p + e;
        q /= q(2);
        flow(i,j,0) = q(0) - src_p(0);
        flow(i,j,1) = q(1) - src_p(1);
      }
    }
  }
}


void compute_occluded_flows_from_depth(const vcl_vector<vpgl_perspective_camera<double> > &cameras,
                                       const vpgl_perspective_camera<double> &ref_cam,
                                       const vil_image_view<double> &depth,
                                       vcl_vector<vil_image_view<double> > &flows)
{
  const unsigned int ni = depth.ni();
  const unsigned int nj = depth.nj();

  flows.resize(cameras.size());
  const vnl_matrix_fixed<double,3,4> &P1 = ref_cam.get_matrix();
  imesh_mesh mesh = depth_map_to_mesh(ref_cam, depth);
  vnl_matrix_fixed<double, 3, 3> inv1 = vnl_inverse(P1.extract(3,3));
  vil_image_view<double> fdepth(ni, nj, 1);
  for (unsigned int c = 0; c < cameras.size(); c++)
  {
    const vnl_matrix_fixed<double,3,4> &P2 = cameras[c].get_matrix();
    vnl_double_3x3 H = P2.extract(3,3);
    H *= inv1;
    vnl_double_3 e = P2.get_column(3);
    e -= H*P1.get_column(3);

    vil_image_view<double> &flow = flows[c];
    flow.set_size(ni, nj, 2);
    flow.fill(vcl_numeric_limits<double>::quiet_NaN());

    //Use plane+parallax to find the flow between current image and the ref image
    for (unsigned int j = 0; j < nj; j++)
    {
      for (unsigned int i = 0; i < ni; i++)
      {
        vnl_double_3 src_p = vnl_double_3(i,j,1.0);
        vnl_double_3 p = depth(i,j) * src_p;
        vnl_double_3 q = H * p + e;
        fdepth(i,j) = q(2);
        q /= q(2);
        flow(i,j,0) = q(0) - src_p(0);
        flow(i,j,1) = q(1) - src_p(1);
      }
    }

    vgl_box_2d<double> bounds = flow_destination_bounds(flow);
    vgl_box_2d<int> bbox(vcl_floor(bounds.min_x()), vcl_ceil(bounds.max_x()),
                         vcl_floor(bounds.min_y()), vcl_ceil(bounds.max_y()));
    const int ox = bbox.min_x();
    const int oy = bbox.min_y();

    vil_image_view<double> warped_depth(bbox.width()+2, bbox.height()+2);
    warped_depth.fill(vcl_numeric_limits<double>::infinity());
    imesh_project_depth(mesh, crop_camera(cameras[c], ox, oy), warped_depth);

    for (unsigned int j = 0; j < nj; j++)
    {
      for (unsigned int i = 0; i < ni; i++)
      {
        int x = i + vcl_floor(flow(i,j,0)) - ox;
        int y = j + vcl_floor(flow(i,j,1)) - oy;
        assert(x>=0);  assert(x<=bbox.width());
        assert(y>=0);  assert(y<=bbox.height());
        const double& fd = fdepth(i,j);
        vnl_double_4 wd( warped_depth(x,y), warped_depth(x+1,y),
                         warped_depth(x,y+1), warped_depth(x+1,y+1));
        double min_d = wd.min_value();
        double max_d = wd.max_value();
        if( fd > max_d + 0.5 * (max_d-min_d))
        {
          flow(i,j,0) = flow(i,j,1) = vcl_numeric_limits<double>::quiet_NaN();
        }
      }
    }
  }
}


vgl_box_2d<double>
flow_destination_bounds(const vil_image_view<double> &flow)
{
  const unsigned int ni=flow.ni();
  const unsigned int nj=flow.nj();

  vgl_box_2d<double> bbox;
  for( unsigned j=0; j<nj; ++j)
  {
    for( unsigned i=0; i<ni; ++i)
    {
      const double& di = flow(i,j,0);
      const double& dj = flow(i,j,1);
      if( vnl_math::isfinite(di) && vnl_math::isfinite(dj) )
      {
        bbox.add(vgl_point_2d<double>(i + di, j + dj));
      }
    }
  }
  return bbox;
}


void translate_flow(vil_image_view<double> &flow, double dx, double dy)
{
  const unsigned int ni=flow.ni();
  const unsigned int nj=flow.nj();

  for( unsigned j=0; j<nj; ++j)
  {
    for( unsigned i=0; i<ni; ++i)
    {
      flow(i,j,0) += dx;
      flow(i,j,1) += dy;
    }
  }
}

} // end namespace super3d

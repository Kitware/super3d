/*ckwg +29
 * Copyright 2012-2013 by Kitware, Inc.
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


#include <vcl_iostream.h>
#include <vcl_limits.h>

#include <testlib/testlib_test.h>

#include <boost/bind.hpp>

#include <super3d/image/adjoint_flow_warp.h>
#include <super3d/image/adjoint_image_utils.h>
#include <super3d/image/adjoint_image_op.h>

#include <vil/vil_math.h>
#include <vil/vil_crop.h>

// Put everything in an anonymous namespace so that different tests
// won't conflict.
namespace
{

using namespace vidtk;

void
test_flow_warp_bilin()
{
  const unsigned ni = 100, nj = 200;
  vil_image_view<float> flow(ni, nj, 2);
  fill_random(flow, -10.0f, 10.0f);

  typedef adjoint_image_ops_func<float>::func_t func_t;
  func_t forward = boost::bind(warp_forward_with_flow_bilin<float,float,float>, _1, flow, _2);
  func_t backward = boost::bind(warp_back_with_flow_bilin<float,float,float>, _1, flow, _2);

  adjoint_image_ops_func<float> adj(forward, backward, ni, nj, 3);
  TEST("Bilin: forward and backward warping adjoint (same size images)", true, adj.is_adjoint());

  adj.set_dst_size(ni+20, nj+20, 3);
  TEST("Bilin: forward and backward warping adjoint (different size images)", true, adj.is_adjoint());

  vil_crop(flow, 10, 20, 10, 20).fill(vcl_numeric_limits<float>::quiet_NaN());
  TEST("Bilin: forward and backward warping adjoint (with invalid flow)", true, adj.is_adjoint());
}


void
test_flow_warp_bicub()
{
  const unsigned ni = 100, nj = 200;
  vil_image_view<float> flow(ni, nj, 2);
  fill_random(flow, -10.0f, 10.0f);

  typedef adjoint_image_ops_func<float>::func_t func_t;
  func_t forward = boost::bind(warp_forward_with_flow_bicub<float,float,float>, _1, flow, _2);
  func_t backward = boost::bind(warp_back_with_flow_bicub<float,float,float>, _1, flow, _2);

  adjoint_image_ops_func<float> adj(forward, backward, ni, nj, 3);
  TEST("Bicub: forward and backward warping adjoint (same size images)", true, adj.is_adjoint());

  adj.set_dst_size(ni+20, nj+20, 3);
  TEST("Bicub: forward and backward warping adjoint (different size images)", true, adj.is_adjoint());

  vil_crop(flow, 10, 20, 10, 20).fill(vcl_numeric_limits<float>::quiet_NaN());
  TEST("Bicub: forward and backward warping adjoint (with invalid flow)", true, adj.is_adjoint());
}


void
test_zero_flow_bilin()
{
  const unsigned ni = 100, nj = 200;
  vil_image_view<float> flow(ni, nj, 2);
  flow.fill(0.0f);

  typedef adjoint_image_ops_func<float>::func_t func_t;
  func_t forward = boost::bind(warp_forward_with_flow_bilin<float,float,float>, _1, flow, _2);
  func_t backward = boost::bind(warp_back_with_flow_bilin<float,float,float>, _1, flow, _2);

  adjoint_image_ops_func<float> adj(forward, backward, ni, nj, 3);
  TEST("Bilin: forward and backward warping with zero flow adjoint", true, adj.is_adjoint());

  vil_image_view<float> src(ni, nj, 1), dst(ni,nj,1);
  fill_random(src, 0.0f, 1.0f);
  warp_forward_with_flow_bilin(src, flow, dst);
  TEST("Bilin: forward warping with zero flow identity operation check",
       true, vil_image_view_deep_equality(src,dst));

  warp_back_with_flow_bilin(src, flow, dst);
  TEST("Bilin: backward warping with zero flow identity operation check",
       true, vil_image_view_deep_equality(src,dst));
}


void
test_zero_flow_bicub()
{
  const unsigned ni = 100, nj = 200;
  vil_image_view<float> flow(ni, nj, 2);
  flow.fill(0.0f);

  typedef adjoint_image_ops_func<float>::func_t func_t;
  func_t forward = boost::bind(warp_forward_with_flow_bicub<float,float,float>, _1, flow, _2);
  func_t backward = boost::bind(warp_back_with_flow_bicub<float,float,float>, _1, flow, _2);

  adjoint_image_ops_func<float> adj(forward, backward, ni, nj, 3);
  TEST("Bicub: forward and backward warping with zero flow adjoint", true, adj.is_adjoint());

  // Note: for bicubic interpolation the outer 1 pixel margin around the image
  // is not interpolated due to boundary conditions.  So zero flow is only the
  // identity warping if you ignore the outer 1 pixel boundary.
  vil_image_view<float> src(ni, nj, 1), dst(ni,nj,1);
  fill_random(src, 0.0f, 1.0f);
  warp_forward_with_flow_bicub(src, flow, dst);
  vil_image_view<float> src_crop = vil_crop(src, 1, ni-2, 1, nj-2);
  vil_image_view<float> dst_crop = vil_crop(dst, 1, ni-2, 1, nj-2);
  TEST("Bicub: forward warping with zero flow identity operation check",
       true, vil_image_view_deep_equality(src_crop, dst_crop));

  warp_back_with_flow_bicub(src, flow, dst);
  dst_crop = vil_crop(dst, 1, ni-2, 1, nj-2);
  TEST("Bicub: backward warping with zero flow identity operation check",
       true, vil_image_view_deep_equality(src_crop, dst_crop));
}

} // end anonymous namespace

int
test_flow_warp( int /*argc*/, char */*argv*/[] )
{
  testlib_test_start( "flow_warp_bilin" );

  test_flow_warp_bilin();
  test_flow_warp_bicub();
  test_zero_flow_bilin();
  test_zero_flow_bicub();

  return testlib_test_summary();
}

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

#include <iostream>

#include <vil/vil_math.h>
#include <boost/bind.hpp>

#include <super3d/image/adjoint_resample.h>
#include <super3d/image/adjoint_image_utils.h>
#include <super3d/image/adjoint_image_op.h>

#include <testlib/testlib_test.h>

// Put everything in an anonymous namespace so that different tests
// won't conflict.
namespace
{

using namespace super3d;

void
test_resampler()
{
  const unsigned ni = 100, nj = 200;
  const unsigned scale = 2;

  typedef adjoint_image_ops_func<float>::func_t func_t;
  func_t forward = boost::bind(down_sample<float>, _1, _2, scale, 0, 0);
  func_t backward = boost::bind(up_sample<float>, _1, _2, scale, 0, 0);

  adjoint_image_ops_func<float> adj(forward, backward, ni, nj, 3);

  TEST("Down sampling and up sampling adjoint", true, adj.is_adjoint());
}

void
test_rescaler()
{
  const unsigned ni = 100, nj = 200;
  vil_image_view<float> I0(ni, nj, 3), I1(ni/2, nj/2, 3);
  fill_random(I0, 0.0f, 255.0f);
  fill_random(I1, 0.0f, 255.0f);

  vil_image_view<float> d_I0, u_I1;
  down_scale(I0, d_I0, 2);
  up_scale(I1, u_I1, 2);

  double sum1 = dot_product(I0, u_I1);
  double sum2 = dot_product(I1, d_I0);
  double relative_error = std::abs(sum1 - sum2) / std::abs(sum1);

  TEST_NEAR("down scaling and up scaling adjoint", relative_error, 0.0, 1e-8);
}


} // end anonymous namespace

int
test_resample( int /*argc*/, char */*argv*/[] )
{
  testlib_test_start( "resample" );
  test_resampler();
  test_rescaler();

  return testlib_test_summary();
}

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
#include <testlib/testlib_test.h>

#include <super3d/image/adjoint_image_derivs.h>
#include <super3d/image/adjoint_image_utils.h>
#include <super3d/image/adjoint_image_op.h>

#include <vil/vil_math.h>

// Put everything in an anonymous namespace so that different tests
// won't conflict.
namespace {

using namespace super3d;

void
test_image_derivs_adjoint()
{
  const unsigned ni = 100, nj = 200;

  typedef adjoint_image_ops_func<double>::func_t func_t;
  func_t forward = forward_gradient<double>;
  func_t backward = backward_divergence<double>;

  adjoint_image_ops_func<double> adj(forward, backward, ni, nj, 1, ni, nj, 2);

  TEST("Forward gradient and backward divergence adjoint (1 plane)",
       true, adj.is_adjoint(1e-6));

  adj.set_src_size(ni, nj, 3);
  adj.set_src_size(ni, nj, 6);

  TEST("Forward gradient and backward divergence adjoint (3 planes)",
       true, adj.is_adjoint(1e-6));
}

void
test_image_derivs_norm()
{
  const unsigned ni = 100, nj = 100, np = 1;

  typedef adjoint_image_ops_func<double>::func_t func_t;
  func_t forward = forward_gradient<double>;
  func_t backward = backward_divergence<double>;

  adjoint_image_ops_func<double> adj(forward, backward, ni, nj, np);

  double norm = adj.norm_estimation(1e-6);
  TEST_NEAR("Gradient norm", norm, std::sqrt(8.0), 1e-2);
}


} // end anonymous namespace

int test_image_derivs( int /*argc*/, char */*argv*/[] )
{
  testlib_test_start( "image_derivs" );

  test_image_derivs_adjoint();
  test_image_derivs_norm();

  return testlib_test_summary();
}

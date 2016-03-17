/*ckwg +29
 * Copyright 2013-2015 by Kitware, Inc.
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

#include <vil/vil_load.h>
#include <vil/vil_convert.h>
#include <vil/vil_crop.h>

#include <vnl/vnl_inverse.h>

#include <testlib/testlib_test.h>

#include <super3d/image/refine_homography.h>
#include <super3d/image/warp_image.h>

// Put everything in an anonymous namespace so that different tests
// won't conflict.
namespace
{

using namespace super3d;

void test_refine_homog(const vil_image_view<double> &img,
                       vnl_double_3x3 &H)
{
  vil_image_view<double> warped(img.ni(), img.nj());
  vnl_double_3x3 Hinv = vnl_inverse<double>(H);
  warp_image(img, warped, Hinv);

  vnl_double_3x3 est_H;
  est_H.set_identity();

  refine_homography(warped, img, est_H, 4.0, 60, 4);
  est_H /= est_H.frobenius_norm();

  TEST_NEAR( "Homography difference (frobenius)", (est_H - H).frobenius_norm(), 0.0, 0.01);
}

} // end anonymous namespace


int test_refine_homography( int argc, char* argv[] )
{
  if( argc < 2 )
  {
    std::cerr << "Need the data directory as an argument!\n";
    return EXIT_FAILURE;
  }

  // load the input test image
  std::string dir = argv[1];
  std::string src_path = dir + "/ocean_city.png";
  vil_image_view<vxl_byte> input = vil_load(src_path.c_str());
  vil_image_view<double> img;
  vil_convert_cast<vxl_byte, double>(input, img);
  img.deep_copy(vil_crop(img, 270, 100, 150, 100));

  testlib_test_start( "refine_homography" );

  vnl_double_3x3 H;
  double angle = 0.1;
  H(0,0) = std::cos(angle); H(0,1) = -std::sin(angle); H(0,2) = -5.0;
  H(1,0) = std::sin(angle); H(1,1) = std::cos(angle); H(1,2) = -5.0;
  H(2,0) = 0.0; H(2,1) = 0.0; H(2,2) = 1.0;
  H /= H.frobenius_norm();

  test_refine_homog(img, H);

  return testlib_test_summary();
}

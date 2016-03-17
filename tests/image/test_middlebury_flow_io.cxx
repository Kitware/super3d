/*ckwg +29
 * Copyright 2012-2015 by Kitware, Inc.
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
#include <algorithm>
#include <exception>
#include <testlib/testlib_test.h>

#include <super3d/image/middlebury_flow_io.h>


// Put everything in an anonymous namespace so that different tests
// won't conflict.
namespace {

using namespace vidtk;

/// test the reading and writing Middlebury flow files
void
test_flow_io( const std::string& file_path )
{
  vil_image_view<float> flow;

  bool success = false;
  try
  {
    read_middlebury_flow(file_path, flow);
    success = true;
  }
  catch( std::exception& e )
  {
    std::cerr << e.what() << std::endl;
  }
  TEST( "Reading flow file successful", success, true );


  success = false;
  std::string out_file = "temp_flow.flo";
  try
  {
    write_middlebury_flow(out_file, flow);
    success = true;
  }
  catch( std::exception& e )
  {
    std::cerr << e.what() << std::endl;
  }
  TEST( "Writing flow file successful", success, true );


  success = false;
  vil_image_view<float> flow2;
  try
  {
    read_middlebury_flow(out_file, flow2);
    success = true;
  }
  catch( std::exception& e )
  {
    std::cerr << e.what() << std::endl;
  }
  TEST( "Reading written flow file successful", success, true );

  TEST( "Read flow equals written flow",
        vil_image_view_deep_equality(flow, flow2), true);

}

} // end anonymous namespace

int test_middlebury_flow_io( int argc, char* argv[] )
{
  if( argc < 2 )
  {
    std::cerr << "Need the data directory as an argument\n";
    return EXIT_FAILURE;
  }

  // load the input test image and the expected output images
  // This flow sample was take from Middlebury sequence RubberWhale
  // using ground truth flow file flow10.flo with crop string 20x15+100+100
  std::string dir = argv[1];
  std::string src_path = dir + '/' + "middlebury_sample.flo";

  testlib_test_start( "middlebury_flow_io" );

  test_flow_io(src_path);

  return testlib_test_summary();
}

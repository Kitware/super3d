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

#include <vcl_iostream.h>
#include <vcl_vector.h>

#include <vul/vul_arg.h>
#include <vil/algo/vil_gauss_filter.h>
#include <vil/vil_image_view.h>
#include <vil/vil_convert.h>
#include <vil/vil_save.h>

#include <super3d/image/adjoint_dbw.h>
#include <super3d/image/adjoint_resample.h>

#include <super3d/depth/file_io.h>

int main(int argc, char *argv[])
{
  vul_arg<vcl_string> directory( 0, "input directory prefix", "" );
  vul_arg<vcl_string> framefile( 0, "frame file", "" );
  vul_arg<vcl_string> output_dir(0, "output directory prefix", "");
  vul_arg<int> downsampling( "-d", "downsample factor", 2);

  vul_arg_parse( argc, argv );

  vcl_vector<vil_image_view<double> > frames;
  vcl_vector<vcl_string> filenames;
  vcl_vector<int> framelist;
  super3d::load_from_frame_file(framefile().c_str(), directory(), filenames, framelist, frames);

  double scale_factor = downsampling();
  double sensor_sigma = 0.25 * sqrt(scale_factor * scale_factor - 1.0);
  for (unsigned int i = 0; i < frames.size(); i++)
  {
    vil_image_view<double> temp;
    vil_gauss_filter_2d(frames[i], temp, sensor_sigma, 3.0*sensor_sigma);
    vidtk::down_sample(temp, frames[i], downsampling());
    vcl_string output_filename = output_dir() + filenames[i];
    vil_image_view<vxl_byte> output_image;
    vil_convert_cast(frames[i], output_image);
    vil_save(output_image, output_filename.c_str());
  }

  return 0;
}

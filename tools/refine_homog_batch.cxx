/*ckwg +29
 * Copyright 2013-2016 by Kitware, Inc.
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

#include "refine_homog.h"

#include <vul/vul_arg.h>
#include <vul/vul_file.h>

#include <vil/vil_image_view.h>
#include <vil/vil_convert.h>
#include <vil/vil_crop.h>
#include <vil/vil_load.h>
#include <vil/vil_save.h>

#include <vnl/vnl_double_3.h>
#include <vnl/vnl_inverse.h>

#include <vgl/algo/vgl_h_matrix_2d.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <utility>
#include <iomanip>

#include <super3d/image/warp_image.h>
#include <super3d/image/refine_homography.h>



void load_homogs_tagged(const std::string &dir, std::vector<std::pair<unsigned int, std::string> > &timestamps,
                        std::vector<std::string> &frame_names, const char *homog_list, std::vector<vnl_double_3x3> &homogs)
{
  std::ifstream infile(homog_list);
  std::string x;
  char buf[60];
  while (infile.good())
  {
    std::getline(infile, x);
    std::istringstream line(x);

    std::string c;
    line >> c;
    std::string garb;

    vgl_h_matrix_2d<double> H;
    if (c == std::string("TS:"))
    {
      std::pair<unsigned int, std::string> ts;
      line >> ts.first >> ts.second;
      std::getline(infile, x);
      std::getline(infile, x);
      std::istringstream h_line(x);
      h_line >> garb >> garb >> garb >> H;
      timestamps.push_back(ts);
      homogs.push_back(H.get_matrix());
      sprintf(buf, "im%04d.png", ts.first);
      frame_names.push_back(dir + std::string(buf));
    }
  }
  std::cout << homogs.size() << "\n";
  infile.close();
}

void load_homogs(const char *homog_list, unsigned int num_homogs, std::vector<std::pair<unsigned int, std::string> > &timestamps, std::vector<vnl_double_3x3> &homogs)
{
  std::ifstream infile(homog_list);
  std::string x;
  for (unsigned int i = 0; i < num_homogs; i++)
  {
    unsigned int index;
    std::string timestamp, nl;
    vgl_h_matrix_2d<double> H;
    infile >> index >> timestamp >> H;
    timestamps.push_back(std::make_pair(index, timestamp));
    homogs.push_back(H.get_matrix());
  }
  std::cout << homogs.size() << "\n";
  infile.close();
}

int main(int argc, char* argv[])
{
  vul_arg<std::string> dir(0, "dir", "");
  vul_arg<std::string> homog_file(0, "homog file", "");
  vul_arg<unsigned int> ref_frame(0, "ref_frame", 0);
  vul_arg<unsigned int> start_frame(0, "start frame", 0);
  vul_arg<unsigned int> end_frame(0, "end_frame", 0);
  vul_arg<std::string> frame_file("-f", "frame list", "");

  vul_arg_parse( argc, argv );

  std::vector<vnl_double_3x3> Hs;
  std::vector<std::string> frame_names;
  std::vector<std::pair<unsigned int, std::string> > timestamps;

  if (frame_file.set())
  {
    std::ifstream infile(frame_file().c_str());
    std::string x;
    while (infile >> x)
    {
      frame_names.push_back(x);
    }
    infile.close();
    load_homogs(homog_file().c_str(), frame_names.size(), timestamps, Hs);
  }
  else
  {
    load_homogs_tagged(dir(), timestamps, frame_names, homog_file().c_str(), Hs);
  }

  vil_image_view_base_sptr image_sptr_ref = vil_load(frame_names[ref_frame()].c_str());
  image_sptr_ref = vil_convert_to_grey_using_rgb_weighting(image_sptr_ref);
  vil_image_view<double> ref = vil_convert_cast<>(double(), image_sptr_ref);

  unsigned int end = std::min((unsigned int)frame_names.size(), (unsigned int)end_frame());
  char buf[100];
  sprintf(buf, "ref_homgs%d-%d.txt", start_frame(), end);
  std::ofstream outfile(buf);
  for (unsigned int i = start_frame(); i < end; i++)
  {
    std::cout << "Reading " << frame_names[i] << "\n";
    vil_image_view_base_sptr image_sptr = vil_load(frame_names[i].c_str());
    image_sptr = vil_convert_to_grey_using_rgb_weighting(image_sptr);
    vil_image_view<double> moving = vil_convert_cast<>(double(), image_sptr);

    vnl_double_3x3 H = Hs[i];
    super3d::refine_homography(ref, moving, H, 10, 100, 15);

    vil_image_view<vxl_byte> out;

    vil_image_view<double> warped(ref.ni(), ref.nj());
    vnl_double_3x3 Hinv = vnl_inverse<double>(H);
    super3d::warp_image(moving, warped, Hinv);
    sprintf(buf, "refined/img%03d.png", i);
    vil_convert_cast(warped, out);
    vil_save(out, buf);

    vil_image_view<double> original_warped(ref.ni(), ref.nj());
    vnl_double_3x3 original_Hinv = vnl_inverse<double>(Hs[i]);
    super3d::warp_image(moving, original_warped, original_Hinv);
    sprintf(buf, "klt/img%03d.png", i);
    vil_convert_cast(original_warped, out);
    vil_save(out, buf);

    outfile << timestamps[i].first << " " <<  timestamps[i].second << "\n" << std::setprecision(20) << H << "\n" << std::endl;
  }

  outfile.close();
}

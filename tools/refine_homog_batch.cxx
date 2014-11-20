/*ckwg +5
 * Copyright 2013 by Kitware, Inc. All Rights Reserved. Please refer to
 * KITWARE_LICENSE.TXT for licensing information, or contact General Counsel,
 * Kitware, Inc., 28 Corporate Drive, Clifton Park, NY 12065.
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

#include <vcl_iostream.h>
#include <vcl_fstream.h>
#include <vcl_sstream.h>
#include <vcl_utility.h>
#include <vcl_iomanip.h>

#include <video_transforms/warp_image.h>
#include <tracking/refine_homography.h>



void load_homogs_tagged(const vcl_string &dir, vcl_vector<vcl_pair<unsigned int, vcl_string> > &timestamps,
                        vcl_vector<vcl_string> &frame_names, const char *homog_list, vcl_vector<vnl_double_3x3> &homogs)
{
  vcl_ifstream infile(homog_list);
  vcl_string x;
  char buf[60];
  while (infile.good())
  {
    vcl_getline(infile, x);
    vcl_istringstream line(x);

    vcl_string c;
    line >> c;
    vcl_string garb;

    vgl_h_matrix_2d<double> H;
    if (c == vcl_string("TS:"))
    {
      vcl_pair<unsigned int, vcl_string> ts;
      line >> ts.first >> ts.second;
      vcl_getline(infile, x);
      vcl_getline(infile, x);
      vcl_istringstream h_line(x);
      h_line >> garb >> garb >> garb >> H;
      timestamps.push_back(ts);
      homogs.push_back(H.get_matrix());
      sprintf(buf, "im%04d.png", ts.first);
      frame_names.push_back(dir + vcl_string(buf));
    }
  }
  vcl_cout << homogs.size() << "\n";
  infile.close();
}

void load_homogs(const char *homog_list, unsigned int num_homogs, vcl_vector<vcl_pair<unsigned int, vcl_string> > &timestamps, vcl_vector<vnl_double_3x3> &homogs)
{
  vcl_ifstream infile(homog_list);
  vcl_string x;
  for (unsigned int i = 0; i < num_homogs; i++)
  {
    unsigned int index;
    vcl_string timestamp, nl;
    vgl_h_matrix_2d<double> H;
    infile >> index >> timestamp >> H;
    timestamps.push_back(std::make_pair(index, timestamp));
    homogs.push_back(H.get_matrix());
  }
  vcl_cout << homogs.size() << "\n";
  infile.close();
}

int main(int argc, char* argv[])
{
  vul_arg<vcl_string> dir(0, "dir", "");
  vul_arg<vcl_string> homog_file(0, "homog file", "");
  vul_arg<unsigned int> ref_frame(0, "ref_frame", 0);
  vul_arg<unsigned int> start_frame(0, "start frame", 0);
  vul_arg<unsigned int> end_frame(0, "end_frame", 0);
  vul_arg<vcl_string> frame_file("-f", "frame list", "");

  vul_arg_parse( argc, argv );

  vcl_vector<vnl_double_3x3> Hs;
  vcl_vector<vcl_string> frame_names;
  vcl_vector<vcl_pair<unsigned int, vcl_string> > timestamps;

  if (frame_file.set())
  {
    vcl_ifstream infile(frame_file().c_str());
    vcl_string x;
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

  unsigned int end = vcl_min((unsigned int)frame_names.size(), (unsigned int)end_frame());
  char buf[100];
  sprintf(buf, "ref_homgs%d-%d.txt", start_frame(), end);
  vcl_ofstream outfile(buf);
  for (unsigned int i = start_frame(); i < end; i++)
  {
    vcl_cout << "Reading " << frame_names[i] << "\n";
    vil_image_view_base_sptr image_sptr = vil_load(frame_names[i].c_str());
    image_sptr = vil_convert_to_grey_using_rgb_weighting(image_sptr);
    vil_image_view<double> moving = vil_convert_cast<>(double(), image_sptr);

    vnl_double_3x3 H = Hs[i];
    vidtk::refine_homography(ref, moving, H, 10, 100, 15);

    vil_image_view<vxl_byte> out;

    vil_image_view<double> warped(ref.ni(), ref.nj());
    vnl_double_3x3 Hinv = vnl_inverse<double>(H);
    vidtk::warp_image(moving, warped, Hinv);
    sprintf(buf, "refined/img%03d.png", i);
    vil_convert_cast(warped, out);
    vil_save(out, buf);

    vil_image_view<double> original_warped(ref.ni(), ref.nj());
    vnl_double_3x3 original_Hinv = vnl_inverse<double>(Hs[i]);
    vidtk::warp_image(moving, original_warped, original_Hinv);
    sprintf(buf, "klt/img%03d.png", i);
    vil_convert_cast(original_warped, out);
    vil_save(out, buf);

    outfile << timestamps[i].first << " " <<  timestamps[i].second << "\n" << vcl_setprecision(20) << H << "\n" << vcl_endl;
  }

  outfile.close();
}

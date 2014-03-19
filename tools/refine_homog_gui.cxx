/*
 * Copyright 2013 Kitware, Inc.
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

#include <vgui/vgui.h>
#include <vgui/vgui_style.h>
#include <vgui/vgui_shell_tableau.h>
#include <vgui/vgui_easy2D_tableau.h>
#include <vgui/vgui_viewer2D_tableau.h>
#include <vgui/vgui_selector_tableau.h>
#include <vgui/vgui_deck_tableau.h>

#include <vcl_iostream.h>
#include <vcl_fstream.h>
#include <vcl_sstream.h>

#include <video_transforms/warp_image.h>

#include <rrel/rrel_estimation_problem.h>
#include <rrel/rrel_tukey_obj.h>


class matches_tableau : public vgui_tableau
{
public:
  matches_tableau(const vcl_vector<match>& m, const vcl_vector<edgel> &e_fixed,
                  const vcl_vector<edgel> &e_moving, const vnl_double_3x3 &Homog)
  : matches(m), H(Homog), fixed_edgels(e_fixed), moving_edgels(e_moving) { }

  //: handle events
  bool handle(const vgui_event& e)
  {
    if(e.type == vgui_DRAW )
    {
      vgui_style::new_style(0.0,0.0,1.0,2.0,1.0)->apply_all();
#if 0
      for (unsigned int i = 0; i < fixed_edgels.size(); i++)
      {
        const edgel *f = &fixed_edgels[i];
        const double& xf = f->pos(0);
        const double& yf = f->pos(1);
        double dyf = f->n(0);
        double dxf = f->n(1);
        glBegin(GL_LINES);
            glColor3f(0.0f, 1.0f, 0.0f);
            glVertex2f(xf, yf);
            glVertex2f(xf+dxf, yf+dyf);
       glEnd();
       glBegin(GL_POINTS);
          glColor3f(0.0f, 0.0f, 1.0f);
          glVertex2f(xf, yf);
        glEnd();
      }
#endif

#if 0
      for (unsigned int i = 0; i < moving_edgels.size(); i++)
      {
        const edgel *m = &moving_edgels[i];
        const double& xm = m->pos(0);
        const double& ym = m->pos(1);
        double dym = m->n(0);
        double dxm = m->n(1);
        glBegin(GL_LINES);
            glColor3f(0.0f, 1.0f, 0.0f);
            glVertex2f(xm, ym);
            glVertex2f(xm+dxm, ym+dym);
       glEnd();
       glBegin(GL_POINTS);
          glColor3f(1.0f, 0.0f, 0.0f);
          glVertex2f(xm, ym);
        glEnd();
      }
#endif
#if 1
      vnl_double_3x3 Hinvt = vnl_inverse<double>(H).transpose();
      for(unsigned int i = 0; i < matches.size(); ++i)
      {
        const edgel *f = matches[i].fixed;
        const edgel *m = matches[i].moving;

        vnl_double_3 temp = H * vnl_double_3(m->pos(0), m->pos(1), 1.0);
        const double xm = temp(0) / temp(2);
        const double ym = temp(1) / temp(2);

        const double& xf = f->pos(0);
        const double& yf = f->pos(1);

        double dyf = -f->n(1);
        double dxf = f->n(0);



       //Warp as a directional tangent
      vnl_double_2 wt = mult_and_norm(H, m->pos + vnl_double_2(-m->n(1), m->n(0))) - mult_and_norm(H, m->pos);
      vnl_double_2 n = wt.normalize();

      //  vnl_double_2 n = warp_normal(H, m->pos, m->n);

        double dym = n(0);
        double dxm = n(1);

        glBegin(GL_LINES);

            glColor3f(1.0f, 1.0f, 0.0f);
            glVertex2f(xm, ym);
            glVertex2f(xf, yf);
            glColor3f(0.0f, 1.0f, 0.0f);
            glVertex2f(xf, yf);
            glVertex2f(xf+dxf, yf+dyf);
            glColor3f(1.0f, 0.0f, 1.0f);
            glVertex2f(xm, ym);
            glVertex2f(xm+dxm, ym+dym);
        glEnd();

        glBegin(GL_POINTS);
          glColor3f(1.0f, 0.0f, 0.0f);
          glVertex2f(xm, ym);
          glColor3f(0.0f, 0.0f, 1.0f);
          glVertex2f(xf, yf);
        glEnd();
      }
#endif

      return true;
    }
    return false;
  }

  const vcl_vector<match> &matches;
  const vnl_double_3x3 &H;
  const vcl_vector<edgel> &fixed_edgels;
  const vcl_vector<edgel> &moving_edgels;
  vgui_parent_child_link child;
};


typedef vgui_tableau_sptr_t<matches_tableau> matches_tableau_sptr;

struct matches_tableau_new : public matches_tableau_sptr
{
  typedef matches_tableau_sptr base;
  matches_tableau_new(const std::vector<match>& m, const vcl_vector<edgel> &e_fixed,
                  const vcl_vector<edgel> &e_moving, const vnl_double_3x3 &H)
  : base(new matches_tableau(m, e_fixed, e_moving, H)) {}
};

vnl_double_3x3 load_homog(const char *homog_list, int homog_index_f, int homog_index_m)
{
  vcl_ifstream infile(homog_list);
  vgl_h_matrix_2d<double> H_f, H_m, H;
  int index = 0; int framenum = 0;
  while (infile >> H)
  {
    if (homog_index_f == framenum)
    {
      H_f = H;
    }

    if (homog_index_m == framenum)
    {
      H_m = H;
    }

    framenum++;
  }

  return H_f.get_inverse().get_matrix() * H_m.get_matrix();
}

//"C:\Data\argus\seq8\107.png" "C:\Data\argus\seq8\092.png" -h  "C:\Data\argus\seq8\homogs.txt 107 92" -c "1366 400 918 380"

int main(int argc, char* argv[])
{
  vgui::init(argc, argv);

  vul_arg<vcl_string> input_image_f( 0, "input image fixed", "" );
  vul_arg<vcl_string> input_image_m( 0, "input image moving", "" );
  vul_arg<vcl_string> homog_info( "-h", "homog file fixed index moving index", "" );
  vul_arg<vcl_string> crop_string("-c", "crop string", "");

  vul_arg_parse( argc, argv );

  vil_image_view_base_sptr image_sptr1 = vil_load(input_image_f().c_str());
  image_sptr1 = vil_convert_to_grey_using_rgb_weighting(image_sptr1);
  vil_image_view<double> I_f = vil_convert_cast<>(double(), image_sptr1);

  vil_image_view_base_sptr image_sptr2 = vil_load(input_image_m().c_str());
  image_sptr2 = vil_convert_to_grey_using_rgb_weighting(image_sptr2);
  vil_image_view<double> I_m = vil_convert_cast<>(double(), image_sptr2);

  vnl_double_3x3 H;
  if (homog_info.set())
  {
    vcl_istringstream hstream(homog_info());
    vcl_string homog_file;
    unsigned int homog_index_f, homog_index_m;
    hstream >> homog_file >> homog_index_f >> homog_index_m;
    H = load_homog(homog_file.c_str(), homog_index_f, homog_index_m);
  }
  else
  {
    H.set_identity();
  }

  //H.set_identity();


  if (crop_string.set())
  {
    int i0f, nif, j0f, njf;
    vcl_istringstream cwstream(crop_string());
    cwstream >> i0f >> nif >> j0f >> njf;
    I_f.deep_copy(vil_crop(I_f, i0f, nif, j0f, njf));
    crop_homography(H, i0f, j0f);
  }

  vcl_vector<match> matches;
  vcl_vector<edgel> e_fixed, e_moving;
  vil_image_view<unsigned int> index_map;

  extract_driving_edgels(I_m, 5, 0.0, e_moving);
  extract_matchable_edgels(I_f, 5, 0.0, e_fixed, index_map);


  vnl_double_3x3 original_H = H;

  for (unsigned int i = 0; i < 50; i++)
  {
    match_edgels(H, e_fixed, e_moving, index_map, 15, matches);
    estimate_homog_lm(matches, H);
  }

  vil_image_view<double> warped(I_f.ni(), I_f.nj());
  vnl_double_3x3 Hinv = vnl_inverse<double>(H);
  vidtk::warp_image(I_m, warped, Hinv);

  vil_image_view<double> original_warped(I_f.ni(), I_f.nj());
  vnl_double_3x3 original_Hinv = vnl_inverse<double>(original_H);
  vidtk::warp_image(I_m, original_warped, original_Hinv);

  matches_tableau_new matches_tab(matches, e_fixed, e_moving, H);
  vgui_image_tableau_new imagef_tab(I_f);
  vgui_image_tableau_new imagew_tab(warped);
  vgui_image_tableau_new imagem_tab(original_warped);

  vgui_deck_tableau_sptr deck = vgui_deck_tableau_new();
  deck->add(imagem_tab);
  deck->add(imagew_tab);
  deck->add(imagef_tab);

  vgui_selector_tableau_new selector_tab;
  selector_tab->add(deck,"image");
  selector_tab->add(matches_tab,"matches");
  selector_tab->set_active("matches");
  selector_tab->set_active("image");

  vgui_viewer2D_tableau_new viewer_tab(selector_tab);
  vgui_shell_tableau_new shell_tab(viewer_tab);
  return vgui::run(shell_tab, I_f.ni(), I_f.nj());
}

/*
 * Copyright 2012 Kitware, Inc.
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


#include "cost_volume.h"

#include <vcl_cstdio.h>
#include <vcl_fstream.h>

#include <video_transforms/warp_image.h>
#include <vil/vil_bilin_interp.h>
#include <vil/algo/vil_sobel_3x3.h>
#include <vnl/vnl_double_3.h>
#include <vnl/vnl_double_2.h>
#include <vbl/vbl_array_2d.h>
#include <vil/algo/vil_gauss_filter.h>

#include <vil/vil_math.h>
#include <vil/vil_save.h>
#include <vil/vil_convert.h>

#include <vcl_limits.h>


const double two_pi = 2.0*vnl_math::pi;

double p(double c, double alpha)
{
  return 1.0 - exp(-c/alpha);
}

//*****************************************************************************

void
compute_world_cost_volume(const vcl_vector<vil_image_view<double> > &frames,
                          const vcl_vector<vpgl_perspective_camera<double> > &cameras,
                          world_space *ws,
                          unsigned int ref_frame,
                          unsigned int S,
                          vil_image_view<double> &cost_volume,
                          double intesity_weight,
                          double gradient_weight,
                          double census_weight)
{
  const vil_image_view<double> &ref = frames[ref_frame];
  cost_volume = vil_image_view<double>(ws->ni(), ws->nj(), 1, S);
  cost_volume.fill(0.0);

  vcl_vector<vpgl_perspective_camera<double> > warp_cams = ws->warp_cams(cameras, ref_frame);

  double s_step = 1.0/static_cast<double>(S);

  vcl_cout << "Computing cost volume of size (" << cost_volume.ni() << ", " << cost_volume.nj() << ", " << cost_volume.nplanes() << ").\n";

  double total_weight = intesity_weight + gradient_weight + census_weight;
  intesity_weight /= total_weight;
  gradient_weight /= total_weight;
  census_weight /= total_weight;

  int ni = ws->ni(), nj = ws->nj();
  vil_image_view<double> warp_ref(ni, nj, 1), warp(ni, nj, 1);
  vil_image_view<double> warp_ref_grad(ni, nj, 2), warp_grad(ni, nj, 2);
  double denom = 1.0 / (frames.size() - 1.0);

  vbl_array_2d<vxl_uint_64> warp_ref_census;
  //vbl_array_2d<g_census> warp_ref_g_census;

  const double g_census_thresh = 1e-4;
  vil_image_view<int> counts(ni, nj, 1);

  //Depths
  for (unsigned int k = 0; k < S; k++)
  {
    vcl_cout << k << " " << vcl_flush;
    double s = (k + 0.5) * s_step;

    //Warp ref image to world volume
    ws->warp_image_to_depth(ref, warp_ref, warp_cams[ref_frame], s, ref_frame);
    if (gradient_weight)
    {
      vil_sobel_3x3(warp_ref, warp_ref_grad);
    }
    if (census_weight)
    {
      warp_ref_census.resize(ref.ni(), ref.nj());
      //warp_ref_g_census.resize(ref.ni(), ref.nj());
      for (unsigned int i = 0; i < warp_ref.ni(); i++) {
        for (unsigned int j = 0; j < warp_ref.nj(); j++) {
        warp_ref_census(i, j) = compute_census(warp_ref, i, j);
        //warp_ref_g_census(i, j) = compute_g_census(warp_ref_grad, i, j, g_census_thresh);
        }
      }
    }

    counts.fill(0);
    for (unsigned int f = 0; f < frames.size(); f++)
    {
      if (f == ref_frame)
        continue;

      //Warp frame to world volume
      ws->warp_image_to_depth(frames[f], warp, warp_cams[f], s, f);
      if (gradient_weight)
      {
        vil_sobel_3x3(warp, warp_grad);
      }

      for (unsigned int j = 0; j < warp_ref.nj(); j++)
      {
        for (unsigned int i = 0; i < warp_ref.ni(); i++)
        {
          if (warp(i,j) == -1)
            continue;

          double Di = fabs(warp_ref(i,j) - warp(i,j));

          double cost = intesity_weight * Di;
          if (gradient_weight)
          {
            double Dgx = warp_ref_grad(i,j,0) - warp_grad(i,j,0);
            double Dgy = warp_ref_grad(i,j,1) - warp_grad(i,j,1);
            double Dg = sqrt(Dgx*Dgx + Dgy*Dgy);
            cost += gradient_weight * Dg;
          }
          if (census_weight)
          {
            double hamming = (double)hamming_distance(compute_census(warp, i, j), warp_ref_census(i,j));
            //double hamming = (double)hamming_distance(compute_g_census(warp_grad, i, j, g_census_thresh), warp_ref_g_census(i, j));
            cost += census_weight * hamming/63.0;
            //cost += p(hamming, 30) + p(Di, 10);
          }
          cost_volume(i, j, k) += cost;
          counts(i,j) += 1;
        }
      }
    }

    for (unsigned int j = 0; j < warp_ref.nj(); j++)
    {
      for (unsigned int i = 0; i < warp_ref.ni(); i++)
      {
        if (counts(i,j) == 0)
          cost_volume(i,j,k) = 1e6;
        else
          cost_volume(i, j, k) /= (double)counts(i,j);
      }
    }
  }
    //
  vcl_cout << "\n";
}

//*****************************************************************************

void
compute_cost_volume_warp(const vcl_vector<vil_image_view<double> > &frames,
                         const vcl_vector<vpgl_perspective_camera<double> > &cameras,
                         unsigned int ref_frame,
                         unsigned int S,
                         double depth_min,
                         double depth_max,
                         vil_image_view<double> &cost_volume)
{
  const vil_image_view<double> &ref = frames[ref_frame];
  cost_volume = vil_image_view<double>(ref.ni(), ref.nj(), 1, S);

  cost_volume.fill(0.0);

  double s_step = 1.0/static_cast<double>(S);

  vil_image_view<double> ref_grad(ref.ni(), ref.nj(), 2), grad(ref.ni(), ref.nj(), 2);
  vil_sobel_3x3(frames[ref_frame], ref_grad);

  vcl_cout << "Computing cost volume of size (" << ref.ni() << ", " << ref.nj() << ", " << cost_volume.nplanes() << ").\n";
  double bcc_mix = 0.3;

  const vpgl_perspective_camera<double>& camera_ref = cameras[ref_frame];
  vgl_rotation_3d<double> R_ref = camera_ref.get_rotation().inverse();
  vnl_svd<double> svd(camera_ref.get_calibration().get_matrix());
  vnl_matrix_fixed<double, 3, 3> Kref_v = svd.pinverse();

  vil_image_view<double> warped(frames[ref_frame].ni(), frames[ref_frame].nj(),1);
  vidtk::warp_image_parameters wip;
  wip.set_fill_unmapped(true);
  wip.set_unmapped_value(-1.0);

  double denom = 1.0 / (frames.size() - 1.0);

  for (unsigned int f = 0; f < frames.size(); f++)
  {
    if (f == ref_frame)
      continue;

    //Depths
    for (unsigned int k = 0; k < S; k++)
    {
      double s = (k + 0.5) * s_step;
      double idepth = 1.0/((1.0-s)*depth_min + s*depth_max);

      // Compute relative rotation and translation between cameras
      vgl_rotation_3d<double> R_relative = cameras[f].get_rotation() * R_ref;
      vgl_vector_3d<double> t_relative = R_relative * -camera_ref.get_translation() + cameras[f].get_translation();

      // Compute the homography between image planes induced by a plane parallel to the reference image plane
      // at an inverse depth of idepth.
      vnl_double_3x3 H = R_relative.as_matrix();
      H.set_column(2, H.get_column(2) + idepth * vnl_double_3(t_relative.x(), t_relative.y(), t_relative.z()));
      H = cameras[f].get_calibration().get_matrix() * H * Kref_v;
      vidtk::warp_image(frames[f], warped, vgl_h_matrix_2d<double>(H), wip );

      #if 0
      vil_image_view<double> outwrite;
      outwrite.deep_copy(warped);
      vil_math_scale_and_offset_values(outwrite, 255.0, 0.0);
      vil_image_view<vxl_byte> to_save;
      vil_convert_cast<double, vxl_byte>(outwrite, to_save);
      char buf[60];
      sprintf(buf, "images/slice%2f_frame%d_%d.png", s, f, wip.interpolator_);
      vil_save(to_save, buf);
      #endif

      vil_sobel_3x3(warped, grad);

      for (unsigned int j = 0; j < ref.nj(); j++)
      {
        for (unsigned int i = 0; i < ref.ni(); i++)
        {
          double cost = bcc_mix * fabs(ref(i,j) - warped(i,j)) + (1.0 - bcc_mix)*(fabs(ref_grad(i,j,0) - grad(i,j,0)) + fabs(ref_grad(i,j,1) - grad(i,j,1)));
          cost_volume(i, j, k) += cost * denom;
        }
      }
    }
  }
}

//*****************************************************************************

void
compute_cost_volume_bp(const vcl_vector<vil_image_view<double> > &frames,
                       const vcl_vector<vpgl_perspective_camera<double> > &cameras,
                       unsigned int ref_frame,
                       unsigned int S,
                       double depth_min,
                       double depth_max,
                       vil_image_view<double> &cost_volume)
{
  const vil_image_view<double> &ref = frames[ref_frame];
  cost_volume = vil_image_view<double>(ref.ni(), ref.nj(), 1, S);

  double s_step = 1.0/static_cast<double>(S);

  vcl_vector<vil_image_view<double> > grads(frames.size());
  for (unsigned int f = 0;  f < frames.size(); f++)
    vil_sobel_3x3(frames[f], grads[f]);

  vcl_cout << "Computing cost volume of size (" << ref.ni() << ", " << ref.nj() << ", " << cost_volume.nplanes() << ").\n";
  double bcc_mix = 1;

  for (unsigned int j = 0; j < ref.nj(); j++)
  {
    for (unsigned int i = 0; i < ref.ni(); i++)
    {
      vgl_line_3d_2_points<double> ray = cameras[ref_frame].backproject(vgl_point_2d<double>(i,j));
      vgl_vector_3d<double> dir = ray.direction();
      normalize(dir);

      //Depths
      for (unsigned int k = 0; k < S; k++)
      {
        double s = (k + 0.5) * s_step;
        double cost = 0.0, depth = (1.0-s)*depth_min + s*depth_max;
        double ref_val = ref(i,j), ref_gradx = grads[ref_frame](i,j,0), ref_grady = grads[ref_frame](i,j,1);
        double num_frames = (double)frames.size()-1;
        vgl_point_3d<double> pt = cameras[ref_frame].get_camera_center() + depth*dir;

        //Images
        for (unsigned int f = 0; f < frames.size(); f++)
        {
          if (f == ref_frame)
            continue;

          vgl_point_2d<double> uf = cameras[f].project(pt);

          //Intensity difference
          if (uf.x() > 0.0 && uf.x() < frames[f].ni()-1 && uf.y() > 0.0 && uf.y() < frames[f].nj()-1)
          {
            double val = vil_bilin_interp(frames[f], uf.x(), uf.y());
            double gradx = vil_bilin_interp(grads[f], uf.x(), uf.y(), 0);
            double grady = vil_bilin_interp(grads[f], uf.x(), uf.y(), 1);

            cost += fabs(ref_val - val) + bcc_mix*(fabs(ref_gradx - gradx) + fabs(ref_grady - grady));
          }
          else
            num_frames -= 1.0;
        }

        if (num_frames > 0)
          cost_volume(i, j, k) = cost/num_frames;
        else
          cost_volume(i, j, k) = 1e6;
      }
    }
  }
}

//*****************************************************************************

//disparity increases along negative x
void compute_cost_volume_rectified(const vcl_vector<vil_image_view<double> > &frames,
                                  unsigned int ref_frame,
                                  unsigned int S,
                                  double idepth_min,
                                  double idepth_max,
                                  vil_image_view<double> &cost_volume)
{
  const vil_image_view<double> &ref = frames[ref_frame];
  cost_volume = vil_image_view<double>(ref.ni(), ref.nj(), 1, S);
  vcl_cout << "Computing cost volume of size (" << ref.ni() << ", " << ref.nj() << ", " << cost_volume.nplanes() << ").\n";

  double s_step = 1.0/static_cast<double>(S);

  vcl_vector<vil_image_view<double> > grads(frames.size());
  for (unsigned int f = 0;  f < frames.size(); f++)
    vil_sobel_3x3(frames[f], grads[f]);

  double bcc_mix = .5;

  for (unsigned int i = 0; i < ref.ni(); i++)
  {
    for (unsigned int j = 0; j < ref.nj(); j++)
    {
      for (unsigned int k = 0; k < S; ++k)
      {
        double s = (k + 0.5) * s_step;
        double idepth = (1.0-s)*idepth_min + s*idepth_max;
        double ref_val = ref(i,j), ref_gradx = grads[ref_frame](i,j,0), ref_grady = grads[ref_frame](i,j,1);
        double numframes = (double)frames.size() - 1.0;
        double cost = 0.0;
        for (unsigned int f = 0; f < frames.size(); f++)
        {
          if (f == ref_frame)
            continue;

          double x = (double)i - idepth;
          //Intensity difference
          if (x >= 0.0)
          {
            double val = vil_bilin_interp(frames[f], x, j);
            double gradx = vil_bilin_interp(grads[f], x, j, 0);
            double grady = vil_bilin_interp(grads[f], x, j, 1);

            cost += bcc_mix*fabs(ref_val - val) + (1.0 - bcc_mix)*(fabs(ref_gradx - gradx) + fabs(ref_grady - grady));
          }
          else
            numframes -= 1.0;
        }

        if (numframes > 0.0)
          cost_volume(i, j, k) = cost/numframes;
        else
          cost_volume(i, j, k) = 1e6;
      }
    }
  }
}

//*****************************************************************************

vxl_uint_64 compute_census(const vil_image_view<double> &img, int u, int v)
{
  double val = img(u,v);
  vxl_uint_64 census = 0;
  int mini = vcl_max(0, u-4);
  int minj = vcl_max(0, v-3);
  int maxi = vcl_min((int)img.ni()-1, u+4);
  int maxj = vcl_min((int)img.nj()-1, v+3);
  for (int j = minj; j <= maxj; j++)
  {
    for (int i = mini; i <= maxi; i++)
    {
      if (img(i,j) <= val)
        census = (census << 1) | 1;
      else
        census = census << 1;
    }
  }

  return census;
}

//*****************************************************************************

unsigned int hamming_distance(vxl_uint_64 l, vxl_uint_64 r)
{
  vxl_uint_64 n = l ^ r;
  unsigned int count = 0;
  while (n) {
     count++;
     n &= (n - 1);
  }

  return count;
}

//*****************************************************************************

g_census compute_g_census(const vil_image_view<double> &grad, int u, int v, double thresh)
{
  const unsigned int n_bins = 16;
  vnl_double_2 center(grad(u,v,0), grad(u,v,1));
  g_census census;
  census.ones = census.twos = census.fours = 0;
  for (int j = v - 3; j <= v + 3; j++)
  {
    for (int i = u - 4; i <= u + 4; i++)
    {
      if (i < 0 || i >= grad.ni() || j < 0 || j >= grad.nj()) {
        census.ones  = census.ones  << 1;
        census.twos  = census.twos  << 1;
        census.fours = census.fours << 1;
        continue;
      }

      vnl_double_2 vec = vnl_double_2(grad(i,j,0), grad(i,j,1));
      vec = center - vec;

      unsigned char bin = 0;
      double angle = atan2(vec[1], vec[0]);
      if (angle < 0.0) angle += two_pi;
      bin = ((unsigned char)((angle/two_pi) * (double)n_bins));

      census.ones  = (census.ones  << 1) | (bin  & 0x1);
      census.twos  = (census.twos  << 1) | ((bin & 0x2) >> 1);
      census.fours = (census.fours << 1) | ((bin & 0x4) >> 2);
    }
  }

  return census;
}

//*****************************************************************************

unsigned int hamming_distance(const g_census &l, const g_census &r)
{
  vxl_uint_64 n = (l.ones ^ r.ones) | (l.twos ^ r.twos) | (l.fours ^ r.fours);
  unsigned int count = 0;
  while (n) {
     count++;
     n &= (n - 1);
  }

  return count;
}

//*****************************************************************************

void save_cost_volume(const vil_image_view<double> &cost_volume,
                      const vil_image_view<double> &g_weight,
                      const char *file_name)
{
  vcl_cout << "Saving cost volume to " << file_name << "\n";
  FILE *file = vcl_fopen(file_name, "wb");

  unsigned int ni = cost_volume.ni(), nj = cost_volume.nj();
  unsigned int np = cost_volume.nplanes();
  fwrite(&ni, sizeof(unsigned int), 1, file);
  fwrite(&nj, sizeof(unsigned int), 1, file);
  fwrite(&np, sizeof(unsigned int), 1, file);

  for (unsigned int i = 0; i < ni; i++)
  {
    for (unsigned int j = 0; j < nj; j++)
    {
      for (unsigned int s = 0; s < np; s++)
      {
        fwrite(&cost_volume(i,j,s), sizeof(double), 1, file);
      }
    }
  }

  for (unsigned int i = 0; i < ni; i++)
    for (unsigned int j = 0; j < nj; j++)
      fwrite(&g_weight(i,j), sizeof(double), 1, file);

  fclose(file);
}

//*****************************************************************************

void load_cost_volume(vil_image_view<double> &cost_volume,
                      vil_image_view<double> &g_weight,
                      const char *file_name)
{
  vcl_cout << "Loading cost volume from " << file_name << "\n";

  FILE *file = fopen(file_name, "rb");

  unsigned int ni, nj, np;
  fread(&ni, sizeof(unsigned int), 1, file);
  fread(&nj, sizeof(unsigned int), 1, file);
  fread(&np, sizeof(unsigned int), 1, file);

  g_weight.set_size(ni, nj, 1);
  cost_volume = vil_image_view<double>(ni, nj, 1, np);

  for (unsigned int i = 0; i < ni; i++)
  {
    for (unsigned int j = 0; j < nj; j++)
    {
      for (unsigned int s = 0; s < np; s++)
      {
        fread(&cost_volume(i,j,s), sizeof(double), 1, file);
      }
    }
  }

  for (unsigned int i = 0; i < ni; i++)
    for (unsigned int j = 0; j < nj; j++)
      fread(&g_weight(i,j), sizeof(double), 1, file);

  fclose(file);
}

//*****************************************************************************

void
read_cost_volume_at(FILE *file,
                    unsigned int *dims,
                    unsigned int i,
                    unsigned int j,
                    vnl_vector<double> &values)
{
  values.set_size(dims[2]);
  long offset = sizeof(unsigned int)*3 + (i*dims[1]*dims[2] + j*dims[2])*sizeof(double);
  fseek(file, offset, SEEK_SET);
  for (unsigned int s = 0; s < dims[2]; s++)
  {
    fread(&values(s), sizeof(double), 1, file);
  }
}

//*****************************************************************************

void compute_depth_range(const vpgl_perspective_camera<double> &ref_cam,
                         const vcl_string &landmark_file, double &min_depth, double &max_depth)
{
  vcl_ifstream infile(landmark_file.c_str());
  vcl_string x;
  unsigned int numverts;
  do
  {
    infile >> x;
    if (x == "element")
    {
      infile >> x;
      if (x == "vertex")
      {
        infile >> numverts;
      }
    }
  } while (x != vcl_string("end_header"));

  min_depth = vcl_numeric_limits<double>::infinity();
  max_depth = -vcl_numeric_limits<double>::infinity();

  vnl_double_3 c(ref_cam.get_camera_center().x(),
                 ref_cam.get_camera_center().y(),
                 ref_cam.get_camera_center().z());

  for (unsigned int i = 0; i < numverts; i++)
  {
    double x, y, z;
    unsigned int id;
    infile >> x >> y >> z >> id;
    vnl_vector_fixed<double, 4> pt(x, y, z, 1.0);
    vnl_double_3 res = ref_cam.get_matrix() * pt;
    double dist = res(2);
    if (min_depth > dist)
      min_depth = dist;
    if (max_depth < dist)
      max_depth = dist;
  }

  double diff = max_depth - min_depth;
  double offset = diff * 0.5;
  max_depth += offset;
  min_depth -= offset;

  infile.close();
}

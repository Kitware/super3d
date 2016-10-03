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

#include "cost_volume.h"
#include "depth_map.h"

#include <cstdio>
#include <fstream>

#include <super3d/image/warp_image.h>
#include <vil/vil_bilin_interp.h>
#include <vil/algo/vil_sobel_3x3.h>
#include <vnl/vnl_double_3.h>
#include <vnl/vnl_double_2.h>
#include <vbl/vbl_array_2d.h>
#include <vil/algo/vil_gauss_filter.h>
#include <vgl/vgl_box_2d.h>

#include <vil/vil_math.h>
#include <vil/vil_save.h>
#include <vil/vil_convert.h>

#include <limits>


namespace
{

const double two_pi = 2.0*vnl_math::pi;

double p(double c, double alpha)
{
  return 1.0 - exp(-c/alpha);
}

} // end anonymous namespace


namespace super3d
{

void
compute_world_cost_volume(const std::vector<vil_image_view<double> > &frames,
                          const std::vector<vpgl_perspective_camera<double> > &cameras,
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

  std::vector<vpgl_perspective_camera<double> > warp_cams = ws->warp_cams(cameras, ref_frame);

  double s_step = 1.0/static_cast<double>(S);

  std::cout << "Computing cost volume of size (" << cost_volume.ni() << ", " << cost_volume.nj() << ", " << cost_volume.nplanes() << ").\n";

  double total_weight = intesity_weight + gradient_weight + census_weight;
  intesity_weight /= total_weight;
  gradient_weight /= total_weight;
  census_weight /= total_weight;

  int ni = ws->ni(), nj = ws->nj();
  vil_image_view<double> warp_ref(ni, nj, 1), warp(ni, nj, 1);
  vil_image_view<double> warp_ref_grad(ni, nj, 2), warp_grad(ni, nj, 2);

  vbl_array_2d<vxl_uint_64> warp_ref_census;
  vbl_array_2d<g_census> warp_ref_g_census;

  vil_image_view<int> counts(ni, nj, 1);

  //Depths
  for (unsigned int k = 0; k < S; k++)
  {
    std::cout << k << " " << std::flush;
    double s = (k + 0.5) * s_step;

    //Warp ref image to world volume
    ws->warp_image_to_depth(ref, warp_ref, warp_cams[ref_frame], s, ref_frame);
    //if (gradient_weight)
    {
      vil_sobel_3x3(warp_ref, warp_ref_grad);
    }
    if (census_weight)
    {
      //warp_ref_census.resize(ref.ni(), ref.nj());
      warp_ref_g_census.resize(ref.ni(), ref.nj());
      for (unsigned int i = 0; i < warp_ref.ni(); i++) {
        for (unsigned int j = 0; j < warp_ref.nj(); j++) {
        //warp_ref_census(i, j) = compute_census(warp_ref, i, j);
        warp_ref_g_census(i, j) = compute_g_census(warp_ref_grad, i, j);
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
      //if (gradient_weight)
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
            vnl_double_2 warp_ref(warp_ref_grad(i,j,0), warp_ref_grad(i,j,1));
            vnl_double_2 warp(warp_grad(i,j,0), warp_grad(i,j,1));
            double Dg1 = (warp_ref - warp).two_norm();
            double Dg2 = (warp_ref + warp).two_norm();

            //double Dgx = warp_ref_grad(i,j,0) - warp_grad(i,j,0);
            //double Dgy = warp_ref_grad(i,j,1) - warp_grad(i,j,1);
            //double Dg = sqrt(Dgx*Dgx + Dgy*Dgy);
            cost += gradient_weight * (Dg1 < Dg2 ? Dg1 : Dg2);
          }
          if (census_weight)
          {
            //double hamming = (double)hamming_distance(compute_census(warp, i, j), warp_ref_census(i,j));
            //cost += census_weight * hamming/63.0;

            g_census gc = compute_g_census(warp_grad, i, j);
            //double hamming_mag = (double)hamming_distance(gc.mag, warp_ref_g_census(i, j).mag);
            double hamming_ori = (double)hamming_distance(gc.ori, warp_ref_g_census(i, j).ori);
            cost += census_weight * hamming_ori/63.0;

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
  std::cout << "\n";
}


void
compute_cost_volume_warp(const std::vector<vil_image_view<double> > &frames,
                         const std::vector<vpgl_perspective_camera<double> > &cameras,
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

  std::cout << "Computing cost volume of size (" << ref.ni() << ", " << ref.nj() << ", " << cost_volume.nplanes() << ").\n";
  double bcc_mix = 0.3;

  const vpgl_perspective_camera<double>& camera_ref = cameras[ref_frame];
  vgl_rotation_3d<double> R_ref = camera_ref.get_rotation().inverse();
  vnl_svd<double> svd(camera_ref.get_calibration().get_matrix());
  vnl_matrix_fixed<double, 3, 3> Kref_v = svd.pinverse();

  vil_image_view<double> warped(frames[ref_frame].ni(), frames[ref_frame].nj(),1);
  super3d::warp_image_parameters wip;
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
      super3d::warp_image(frames[f], warped, vgl_h_matrix_2d<double>(H), wip );

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


void
compute_cost_volume_bp(const std::vector<vil_image_view<double> > &frames,
                       const std::vector<vpgl_perspective_camera<double> > &cameras,
                       unsigned int ref_frame,
                       unsigned int S,
                       double depth_min,
                       double depth_max,
                       vil_image_view<double> &cost_volume)
{
  const vil_image_view<double> &ref = frames[ref_frame];
  cost_volume = vil_image_view<double>(ref.ni(), ref.nj(), 1, S);

  double s_step = 1.0/static_cast<double>(S);

  std::vector<vil_image_view<double> > grads(frames.size());
  for (unsigned int f = 0;  f < frames.size(); f++)
    vil_sobel_3x3(frames[f], grads[f]);

  std::cout << "Computing cost volume of size (" << ref.ni() << ", " << ref.nj() << ", " << cost_volume.nplanes() << ").\n";
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


//disparity increases along negative x
void compute_cost_volume_rectified(const std::vector<vil_image_view<double> > &frames,
                                  unsigned int ref_frame,
                                  unsigned int S,
                                  double idepth_min,
                                  double idepth_max,
                                  vil_image_view<double> &cost_volume)
{
  const vil_image_view<double> &ref = frames[ref_frame];
  cost_volume = vil_image_view<double>(ref.ni(), ref.nj(), 1, S);
  std::cout << "Computing cost volume of size (" << ref.ni() << ", " << ref.nj() << ", " << cost_volume.nplanes() << ").\n";

  double s_step = 1.0/static_cast<double>(S);

  std::vector<vil_image_view<double> > grads(frames.size());
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


vxl_uint_64 compute_census(const vil_image_view<double> &img, int u, int v)
{
  double val = img(u,v);
  vxl_uint_64 census = 0;
  for (int j = v - 3; j <= v + 3; j++)
  {
    for (int i = u - 4; i <= u + 4; i++)
    {
      if (i < 0 || i >= static_cast<int>(img.ni()) ||
          j < 0 || j >= static_cast<int>(img.nj())) {
        census = census << 1;
        continue;
      }

      if (img(i,j) <= val)
        census = (census << 1) | 1;
      else
        census = census << 1;
    }
  }

  return census;
}


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


g_census compute_g_census(const vil_image_view<double> &grad, int u, int v)
{
  vnl_double_2 center(grad(u,v,0), grad(u,v,1));
  g_census census;
  census.ori = 0;
  census.mag = 0;
  for (int j = v - 3; j <= v + 3; j++)
  {
    for (int i = u - 4; i <= u + 4; i++)
    {
      if (i < 0 || i >= static_cast<int>(grad.ni()) ||
          j < 0 || j >= static_cast<int>(grad.nj())) {
        census.ori = census.ori << 1;
        census.mag = census.mag << 1;
        continue;
      }

      vnl_double_2 neighbor(grad(i,j,0), grad(i,j,1));
      double angle = fabs(atan2(center[1], center[0]) - atan2(neighbor[1], neighbor[0]));
      if (angle > vnl_math::pi)
        angle = fabs(angle - two_pi);

      if (center.squared_magnitude() <= neighbor.squared_magnitude())
        census.mag = (census.mag << 1) | 1;
      else
        census.mag = census.mag << 1;

      if (angle <= vnl_math::pi_over_2)
        census.ori = (census.ori << 1) | 1;
      else
        census.ori = census.ori << 1;
    }
  }

  return census;
}

void save_cost_volume(const vil_image_view<double> &cost_volume,
                      const vil_image_view<double> &g_weight,
                      const char *file_name)
{
  std::cout << "Saving cost volume to " << file_name << "\n";
  FILE *file = std::fopen(file_name, "wb");

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


void load_cost_volume(vil_image_view<double> &cost_volume,
                      vil_image_view<double> &g_weight,
                      const char *file_name)
{
  std::cout << "Loading cost volume from " << file_name << "\n";

  FILE *file = fopen(file_name, "rb");

  unsigned int ni, nj, np;
  if (fread(&ni, sizeof(unsigned int), 1, file) != 1 ||
      fread(&nj, sizeof(unsigned int), 1, file) != 1 ||
      fread(&np, sizeof(unsigned int), 1, file) != 1 )
  {
    std::cerr << "Error loading cost volume" << std::endl;
    return;
  }

  g_weight.set_size(ni, nj, 1);
  cost_volume = vil_image_view<double>(ni, nj, 1, np);

  for (unsigned int i = 0; i < ni; i++)
  {
    for (unsigned int j = 0; j < nj; j++)
    {
      for (unsigned int s = 0; s < np; s++)
      {
        if (fread(&cost_volume(i,j,s), sizeof(double), 1, file) != 1)
        {
          std::cerr << "Error loading cost volume" << std::endl;
          return;
        }
      }
    }
  }

  for (unsigned int i = 0; i < ni; i++)
    for (unsigned int j = 0; j < nj; j++)
      if (fread(&g_weight(i,j), sizeof(double), 1, file) != 1)
      {
        std::cerr << "Error loading cost volume" << std::endl;
        return;
      }

  fclose(file);
}


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
    if (fread(&values(s), sizeof(double), 1, file) != 1)
    {
      std::cerr << "Error loading cost volume at location" << std::endl;
      return;
    }
  }
}

bool compute_depth_range(const vpgl_perspective_camera<double> &ref_cam,
                         int i0, int ni, int j0, int nj,
                         const std::vector<vnl_double_3> &landmarks, double &min_depth, double &max_depth)
{
  min_depth = std::numeric_limits<double>::infinity();
  max_depth = -std::numeric_limits<double>::infinity();

  std::vector<double> depths;

  vgl_box_2d<double> box(i0, i0+ni, j0, j0+nj);
  std::cout << box << "\n";
  for (unsigned int i = 0; i < landmarks.size(); i++)
  {
    const vnl_double_3 &p = landmarks[i];
    vnl_vector_fixed<double, 4> pt(p[0], p[1], p[2], 1.0);
    vnl_double_3 res = ref_cam.get_matrix() * pt;

    double u, v;
    ref_cam.project(p[0], p[1], p[2], u, v);

    if (box.contains(u, v))
    {
      depths.push_back(res(2));
    }
  }

  if (depths.size() < 3)
    return false;

  std::sort(depths.begin(), depths.end());

  // threshold for fraction of outlier depths to reject at the
  // near and far of the depth range.
  const double outlier_thresh = 0.05;
  // fraction of total depth range to pad both near and far to
  // account for depths not sampled by the point cloud
  const double safety_margin_factor = 0.33;

  const unsigned int min_index =
    static_cast<unsigned int>((depths.size()-1) * outlier_thresh);
  const unsigned int max_index =
    static_cast<unsigned int>((depths.size()-1) * (1.0 - outlier_thresh));
  min_depth = depths[min_index];
  max_depth = depths[max_index];

  const double safety_margin = safety_margin_factor * (max_depth - min_depth);
  max_depth += safety_margin;
  min_depth -= safety_margin;
  min_depth = std::max(0.0, min_depth);

  return true;
}

} // end namespace super3d

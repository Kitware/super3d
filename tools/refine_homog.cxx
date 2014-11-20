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

// VXL includes
#include <vil/algo/vil_sobel_3x3.h>
#include <vil/algo/vil_suppress_non_max_edges.h>

#include <vnl/vnl_double_3.h>
#include <vnl/vnl_math.h>
#include <vnl/algo/vnl_levenberg_marquardt.h>
#include <vnl/vnl_least_squares_function.h>
#include <vnl/vnl_inverse.h>

#include <vil/algo/vil_gauss_filter.h>

#include <vcl_iostream.h>
#include <vcl_limits.h>


void refine_homography(const vil_image_view<double> &fixed, const vil_image_view<double> &moving,
                       vnl_double_3x3 &H, double grad_thresh, unsigned int num_iters,
                       unsigned int search_radius, double normal_cutoff, double sigma)
{
  vcl_vector<match> matches;
  vcl_vector<edgel> e_fixed, e_moving;
  vil_image_view<unsigned int> index_map;

  extract_driving_edgels(moving, grad_thresh, sigma, e_moving);
  extract_matchable_edgels(fixed, grad_thresh, sigma, e_fixed, index_map);

  for (unsigned int i = 0; i < num_iters; i++)
  {
    match_edgels(H, e_fixed, e_moving, index_map, search_radius, matches, normal_cutoff);
    estimate_homog_lm(matches, H);
  }
}


void crop_homography(vnl_double_3x3 &H, unsigned int i0, unsigned int j0)
{
  vnl_double_3x3 Tinv;
  Tinv.set_identity();
  Tinv(0,2) = -static_cast<double>(i0);
  Tinv(1,2) = -static_cast<double>(j0);

  H = Tinv * H;
}


edgel::edgel(double x, double y, double angle, double m) : pos(x, y), mag(m)
{
  //Angle is normal to the gradient
  n = vnl_double_2(sin(angle), cos(angle)).normalize();
}


void extract_driving_edgels(const vil_image_view<double> &img,
                            double grad_thresh,
                            double sigma,
                            vcl_vector<edgel> &edgels)
{
  vil_image_view<double> grad_x, grad_y;
  if (sigma > 0.0)
  {
    vil_image_view<double> smoothed;
    vil_gauss_filter_2d(img, smoothed, sigma, static_cast<unsigned int>(2 * sigma));
    vil_sobel_3x3(smoothed, grad_x, grad_y);
  }
  else
  {
    vil_sobel_3x3(img, grad_x, grad_y);
  }

  vil_image_view<double> grad;
  vil_suppress_non_max_edges_subpixel(grad_x, grad_y, grad_thresh, grad);

  const unsigned int ni = grad.ni();
  const unsigned int nj = grad.nj();

  const vcl_ptrdiff_t istep = grad.istep();
  const vcl_ptrdiff_t jstep = grad.jstep();
  const vcl_ptrdiff_t pstep = grad.planestep();

  const double *row = grad.top_left_ptr();
  for (unsigned int j = 0; j < nj; ++j, row += jstep)
  {
    const double* pixel = row;
    for (unsigned int i = 0; i < ni; ++i, pixel += istep)
    {
      if (*pixel == 0.0)
      {
        continue;
      }

      const double &mag = *pixel;
      const double &angle = *(pixel + pstep);
      const double &offset = *(pixel + 2*pstep);
      double x = i + vcl_cos(angle)*offset;
      double y = j + vcl_sin(angle)*offset;
      edgels.push_back(edgel(x, y, angle, mag));
    }
  }
}


void extract_matchable_edgels(const vil_image_view<double> &img,
                              double grad_thresh,
                              double sigma,
                              vcl_vector<edgel> &edgels,
                              vil_image_view<unsigned int> &index)
{
  vil_image_view<double> grad_x, grad_y;
  if (sigma > 0.0)
  {
    vil_image_view<double> smoothed;
    vil_gauss_filter_2d(img, smoothed, sigma, static_cast<unsigned int>(2 * sigma));
    vil_sobel_3x3(smoothed, grad_x, grad_y);
  }
  else
  {
    vil_sobel_3x3(img, grad_x, grad_y);
  }

  vil_image_view<double> grad;
  vil_suppress_non_max_edges_subpixel(grad_x, grad_y, grad_thresh, grad);

  const unsigned int ni = grad.ni();
  const unsigned int nj = grad.nj();
  index.set_size(ni, nj, 1);
  index.fill(0);

  const vcl_ptrdiff_t istep_E = grad.istep();
  const vcl_ptrdiff_t istep_M = index.istep();
  const vcl_ptrdiff_t jstep_E = grad.jstep();
  const vcl_ptrdiff_t jstep_M = index.jstep();
  const vcl_ptrdiff_t pstep_E = grad.planestep();

  const double *row_E = grad.top_left_ptr();
  unsigned int *row_M = index.top_left_ptr();
  for (unsigned int j = 0; j < nj; ++j, row_E += jstep_E, row_M += jstep_M)
  {
    const double *pixel_E = row_E;
    unsigned int *pixel_M = row_M;
    for (unsigned int i = 0; i < ni; ++i, pixel_E += istep_E, pixel_M += istep_M)
    {
      if (*pixel_E == 0.0)
      {
        continue;
      }

      const double &mag = *pixel_E;
      const double &angle = *(pixel_E + pstep_E);
      const double &offset = *(pixel_E + 2*pstep_E);

      double x = i + vcl_cos(angle) * offset;
      double y = j + vcl_sin(angle) * offset;
      edgels.push_back(edgel(x, y, angle, mag));
      *pixel_M = static_cast<unsigned int>(edgels.size()); //index + 1 so that 0 is reserved for no edge
    }
  }
}

///Computes T * [x' 1]' then returns the dehomogenized result
vnl_double_2 mult_and_norm(const vnl_double_3x3 &T, const vnl_double_2 &x)
{
  vnl_double_3 xp_homg = T * vnl_double_3(x(0), x(1), 1.0);
  return vnl_double_2(xp_homg(0) / xp_homg(2), xp_homg(1) / xp_homg(2));
}


vnl_double_2 warp_normal(const vnl_double_3x3 &H, const vnl_double_2 &pt, const vnl_double_2 &normal)
{
  //Warp as a directional tangent
  vnl_double_2 wt = mult_and_norm(H, pt + vnl_double_2(-normal(1), normal(0))) - mult_and_norm(H, pt);
  wt.normalize();
  return vnl_double_2(wt(1), -wt(0));
}


void match_edgels(const vnl_double_3x3 &H,
                  const vcl_vector<edgel> &e_fixed,
                  const vcl_vector<edgel> &e_moving,
                  const vil_image_view<unsigned int> &index,
                  const double search_rad,
                  vcl_vector<match> &matches,
                  double normal_cutoff)
{
  double x_bound = static_cast<double>(index.ni()-1);
  double y_bound = static_cast<double>(index.nj()-1);

  matches.clear();

  for (unsigned int i = 0; i < e_moving.size(); i++)
  {
    const edgel *e = &e_moving[i];
    vnl_double_3 temp = H * vnl_double_3(e->pos(0), e->pos(1), 1.0);
    vnl_double_2 pp;
    pp(0) = temp(0) / temp(2);
    pp(1) = temp(1) / temp(2);

    //Remove matches from outside of the ROI
    if (pp(0) < 0.0 || pp(0) > x_bound || pp(1) < 0.0 || pp(1) > y_bound)
      continue;

    vnl_double_2 np = warp_normal(H, e->pos, e->n);

    int left = static_cast<int>(vcl_max(pp(0) - search_rad, 0.0));
    int right = static_cast<int>(vcl_min(pp(0) + search_rad, x_bound));
    int top = static_cast<int>(vcl_max(pp(1) - search_rad, 0.0));
    int bot = static_cast<int>(vcl_min(pp(1) + search_rad, y_bound));

    double closest_dist = vcl_numeric_limits<double>::max();
    unsigned int closest_index;
    for (int m = left; m <= right; m++)
    {
      for (int n = top; n <= bot; n++)
      {
        unsigned int ind = index(m,n);

        if (ind == 0)
          continue;

        --ind;
        const edgel &e = e_fixed[ind];

        //Check orientations are similar
        if (dot_product(e.n, np) < normal_cutoff)
          continue;

        //Euclidean distance to find matches
        double dist = (e.pos - pp).squared_magnitude();
        if (dist < closest_dist)
        {
          closest_dist = dist;
          closest_index = ind;
        }
      }
    }

    if (closest_dist != vcl_numeric_limits<double>::max())
    {
      matches.push_back(match(e, &e_fixed[closest_index]));
    }
  }
  vcl_cout << matches.size() << "\n";
}

///Class required by vnl_levenberg_marquardt to define the energy function
class energy_func : public vnl_least_squares_function
{
public:
  energy_func(unsigned int number_of_unknowns,
              unsigned int number_of_residuals,
              const vcl_vector<match> &correspondences,
              const vcl_vector<vnl_double_2> &normalized_fixed,
              const vcl_vector<vnl_double_2> &normalized_moving,
              const vnl_vector<double> &weights)
             : vnl_least_squares_function(number_of_unknowns, number_of_residuals),
               corresp(correspondences), norm_fixed(normalized_fixed),
               norm_moving(normalized_moving), wgts(weights) {}

  void f(const vnl_vector<double> &x, vnl_vector<double> &fx);
  void gradf(const vnl_vector<double> &x, vnl_matrix<double> &jacobian);

private:

  const vcl_vector<match> &corresp;
  const vcl_vector<vnl_double_2> &norm_fixed;
  const vcl_vector<vnl_double_2> &norm_moving;
  const vnl_vector<double> &wgts;
};

/// Converts parameter vector x to a 3x3 homography matrix
void x_to_H(const vnl_vector<double> &x, vnl_double_3x3 &H)
{
  H(0,0) = x[0];  H(0,1) = x[1];  H(0,2) = x[2];
  H(1,0) = x[3];  H(1,1) = x[4];  H(1,2) = x[5];
  H(2,0) = x[6];  H(2,1) = x[7];  H(2,2) = x[8];
}

///Computes the residuals using parameter vector x
void energy_func::f(const vnl_vector<double> &x, vnl_vector<double> &fx)
{
  vnl_double_3x3 H;
  x_to_H(x, H);

  for (unsigned int i = 0; i < corresp.size(); i++)
  {
    vnl_double_2 pp = mult_and_norm(H, norm_moving[i]);
    fx[i] = sqrt(wgts[i]) * dot_product(pp - norm_fixed[i], corresp[i].fixed->n);
  }
}

///Computes the jacobian of the energy function at x
void energy_func::gradf(const vnl_vector<double> &x, vnl_matrix<double> &jacobian)
{
  vnl_double_3x3 H;
  x_to_H(x, H);

  for (unsigned int i = 0; i < corresp.size(); i++)
  {
    const vnl_double_2 &n = corresp[i].fixed->n;
    const vnl_double_2 &p = norm_moving[i];

    vnl_double_3 p_homog = H * vnl_double_3(p(0), p(1), 1.0);
    double w = p_homog(2);
    double w_sqr = w * w;
    double sqrtwgt = sqrt(wgts[i]);

    jacobian(i, 0) = sqrtwgt * dot_product(vnl_double_2(p(0)/w, 0.0), n);
    jacobian(i, 1) = sqrtwgt * dot_product(vnl_double_2(p(1)/w, 0.0), n);
    jacobian(i, 2) = sqrtwgt * dot_product(vnl_double_2(1.0/w, 0.0), n);
    jacobian(i, 3) = sqrtwgt * dot_product(vnl_double_2(0.0, p(0)/w), n);
    jacobian(i, 4) = sqrtwgt * dot_product(vnl_double_2(0.0, p(1)/w), n);
    jacobian(i, 5) = sqrtwgt * dot_product(vnl_double_2(0.0, 1.0/w), n);
    jacobian(i, 6) = sqrtwgt * dot_product(vnl_double_2(-p(0) * p_homog(0) / w_sqr, -p(0) * p_homog(1) / w_sqr), n);
    jacobian(i, 7) = sqrtwgt * dot_product(vnl_double_2(-p(1) * p_homog(0) / w_sqr, -p(1) * p_homog(1) / w_sqr), n);
    jacobian(i, 8) = sqrtwgt * dot_product(vnl_double_2(-p_homog(0) / w_sqr, -p_homog(1) / w_sqr), n);
  }
}

///Computes normalization matrices for numerical stability in LM
void compute_normalization(const vcl_vector<match> &matches,
                           vnl_double_3x3 &T_fixed,
                           vnl_double_3x3 &T_moving)
{
  vnl_double_2 c_fixed(0.0, 0.0), c_moving(0.0, 0.0);
  double s_fixed = 0.0, s_moving = 0.0;

  double num = static_cast<double>(matches.size());

  //Compute centers
  for (unsigned int i = 0; i < matches.size(); i++)
  {
    c_fixed += matches[i].fixed->pos;
    c_moving += matches[i].moving->pos;
  }

  c_fixed /= num;
  c_moving /= num;

  //Compute avg scalings
  for (unsigned int i = 0; i < matches.size(); i++)
  {
    s_fixed += (matches[i].fixed->pos - c_fixed).two_norm();
    s_moving += (matches[i].moving->pos - c_moving).two_norm();
  }

  s_fixed /= num;
  s_moving /= num;

  //Build normalization matrices
  T_fixed.fill(0.0);
  T_fixed.put(0, 0, 1.0/s_fixed);  T_fixed.put(0, 2, -c_fixed(0)/s_fixed);
  T_fixed.put(1, 1, 1.0/s_fixed);  T_fixed.put(1, 2, -c_fixed(1)/s_fixed);
  T_fixed.put(2, 2, 1.0);

  T_moving.fill(0.0);
  T_moving.put(0, 0, 1.0/s_moving);  T_moving.put(0, 2, -c_moving(0)/s_moving);
  T_moving.put(1, 1, 1.0/s_moving);  T_moving.put(1, 2, -c_moving(1)/s_moving);
  T_moving.put(2, 2, 1.0);
}


void estimate_homog_lm(const vcl_vector<match> &corresp, vnl_double_3x3 &H)
{
  const unsigned int num_params = 9;
  const unsigned int num_resid = static_cast<unsigned int>(corresp.size());

  //Compute normalizations
  vnl_double_3x3 T_fixed, T_fixed_inv, T_moving, T_moving_inv;
  compute_normalization(corresp, T_fixed, T_moving);
  T_fixed_inv = vnl_inverse<double>(T_fixed);
  T_moving_inv = vnl_inverse<double>(T_moving);

  vcl_vector<vnl_double_2> norm_fixed, norm_moving;
  norm_fixed.reserve(corresp.size()), norm_moving.reserve(corresp.size());
  for (unsigned int i = 0; i < corresp.size(); i++)
  {
    norm_fixed.push_back(mult_and_norm(T_fixed, corresp[i].fixed->pos));
    norm_moving.push_back(mult_and_norm(T_moving, corresp[i].moving->pos));
  }

  H = T_fixed * H * T_moving_inv;

  vnl_vector<double> x(num_params);
  x[0] = H(0,0);  x[1] = H(0,1);  x[2] = H(0,2);
  x[3] = H(1,0);  x[4] = H(1,1);  x[5] = H(1,2);
  x[6] = H(2,0);  x[7] = H(2,1);  x[8] = H(2,2);

  vnl_vector<double> weights(static_cast<unsigned int>(corresp.size()), 1.0);
  energy_func ef(num_params, num_resid, corresp, norm_fixed, norm_moving, weights);
  vnl_levenberg_marquardt lm(ef);
  lm.minimize(x);
  x /= x.two_norm(); //normalize by frobenius

  x_to_H(x, H);
  H = T_fixed_inv * H * T_moving;
}

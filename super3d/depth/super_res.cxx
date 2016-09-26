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

#include "super_res.h"
#include "super_config.h"

#include <super3d/image/adjoint_image_derivs.h>
#include <super3d/image/warp_image.h>

#include <vil/algo/vil_gauss_filter.h>
#include <vil/vil_save.h>
#include <vil/vil_math.h>
#include <vil/vil_resample_bilin.h>
#include <vil/vil_resample_bicub.h>
#include <vil/vil_convert.h>

#include <vil/vil_resample_bicub.hxx>
VIL_RESAMPLE_BICUB_INSTANTIATE( double , double );

#define DEBUG


namespace super3d
{

namespace // anonymous
{

//step of dual variable for lambda * || grad(primal) ||_epsilon
void dual_step_grad_huber(const vil_image_view<double> &primal,
                          vil_image_view<double> &dual,
                          double lambda,
                          double sigma,
                          double epsilon)
{
  const double denom = 1.0 + (sigma * epsilon) / lambda;
  vil_image_view<double> work;
  super3d::forward_gradient(primal, work);
  vil_math_add_image_fraction(dual, 1.0/denom, work, sigma/denom);

  for (unsigned int j = 0; j < dual.nj(); j++)
  {
    for (unsigned int i = 0; i < dual.ni(); i++)
    {
      for (unsigned int k = 0; k < primal.nplanes(); k++)
      {
        double &x = dual(i,j,2*k), &y = dual(i,j,2*k+1);

        //truncate vectors
        const double mag = sqrt(x*x + y*y)/lambda;
        if (mag > 1.0)
        {
          x /= mag;
          y /= mag;
        }
      }
    }
  }
}


void dual_step_b(const vcl_vector<vil_image_view<double> > &a0,
                 vcl_vector<vil_image_view<double> > &b0,
                 const super_res_params &srp)
{
  for (unsigned int i = 0; i < a0.size(); i++)
  {
    if (i == srp.ref_frame)
      continue;

    dual_step_grad_huber(a0[i], b0[i], srp.lambda_l, srp.sigma, srp.epsilon_reg);
  }
}

//*****************************************************************************

void dual_step_q(const vcl_vector<vil_image_view<double> > &frames,
                 const vcl_vector<super3d::adjoint_image_ops_func<double> > &warps,
                 const vcl_vector<vil_image_view<double> > &weights,
                 const vil_image_view<double> &u_bar,
                 const vcl_vector<vil_image_view<double> > &a0,
                 vcl_vector<vil_image_view<double> > &q,
                 const super_res_params &srp)
{
  vil_image_view<double> l_u, l_a0;
  const double sf_2 = 1.0/(srp.scale_factor * srp.scale_factor);
  const double denom = 1.0 + (srp.sigma * srp.epsilon_data) / sf_2;

  for (unsigned int f = 0; f < frames.size(); f++)
  {
    const unsigned int ni = weights[f].ni();
    const unsigned int nj = weights[f].nj();

    // apply the linear operator to warp, blur, and downsample
    warps[f].apply_A(u_bar, l_u);
    warps[f].apply_A(a0[f], l_a0);

    for (unsigned int j = 0; j < nj; j++)
    {
      for (unsigned int i = 0; i < ni; i++)
      {
        for (unsigned int k = 0; k < q[f].nplanes(); k++)
        {
          double &qfijk = q[f](i, j, k);
          qfijk = (qfijk + srp.sigma * sf_2 * (l_u(i, j, k) + l_a0(i, j) - frames[f](i, j, k) * weights[f](i,j)))/denom;
          qfijk = vcl_max(qfijk, -sf_2);
          qfijk = vcl_min(qfijk, sf_2);
        }
      }
    }
  }
}

//*****************************************************************************

void dual_step_q(const vcl_vector<vil_image_view<double> > &frames,
                 const vcl_vector<super3d::adjoint_image_ops_func<double> > &warps,
                 const vcl_vector<vil_image_view<double> > &weights,
                 const vil_image_view<double> &u_bar,
                 vcl_vector<vil_image_view<double> > &q,
                 const super_res_params &srp)
{
  vil_image_view<double> l_u;
  const double sf_2 = 1.0/(srp.scale_factor * srp.scale_factor);
  const double denom = 1.0 + (srp.sigma * srp.epsilon_data) / sf_2;

  for (unsigned int f = 0; f < frames.size(); f++)
  {
    const unsigned int ni = weights[f].ni();
    const unsigned int nj = weights[f].nj();

    // apply the linear operator to warp, blur, and downsample
    warps[f].apply_A(u_bar, l_u);

    for (unsigned int j = 0; j < nj; j++)
    {
      for (unsigned int i = 0; i < ni; i++)
      {
        for (unsigned int k = 0; k < q[f].nplanes(); k++)
        {
          double &qfijk = q[f](i, j, k);
          qfijk = (qfijk + srp.sigma * sf_2 * (l_u(i, j, k) - frames[f](i, j, k) * weights[f](i,j)))/denom;
          qfijk = vcl_max(qfijk, -sf_2);
          qfijk = vcl_min(qfijk, sf_2);
        }
      }
    }
  }
}

void primal_step_u(const vcl_vector<vil_image_view<double> > &q,
                   const vcl_vector<super3d::adjoint_image_ops_func<double> > &warps,
                   const vil_image_view<double> &p,
                   vil_image_view<double> &u,
                   vil_image_view<double> &u_bar,
                   const super_res_params &srp)
{
  const double sf_2 = 1.0/(srp.scale_factor * srp.scale_factor);
  vil_image_view<double> sum_super_q(srp.s_ni, srp.s_nj, u.nplanes());
  sum_super_q.fill(0.0);
  vil_image_view<double> super_q(srp.s_ni, srp.s_nj, u.nplanes()), temp;
  for (unsigned int i = 0; i < q.size(); i++)
  {
    // apply transpose linear operator to upsample, blur, and warp
    warps[i].apply_At(q[i], super_q);
    vil_math_image_sum(sum_super_q, super_q, sum_super_q);
  }

  vil_image_view<double> work;
  super3d::backward_divergence(p, work);
  vil_math_add_image_fraction(work, srp.tau, sum_super_q, -srp.tau * sf_2);
  vil_math_image_sum(u, work, work);
  u_bar.deep_copy(work);
  vil_math_add_image_fraction(u_bar, 2.0, u, -1.0);
  vcl_swap(u, work);
}

//*****************************************************************************

void primal_step_a(const vcl_vector<vil_image_view<double> > &q,
                   const vcl_vector<super3d::adjoint_image_ops_func<double> > &warps,
                   const vil_image_view<double> &u,
                   const vcl_vector<vil_image_view<double> > &b0,
                   vcl_vector<vil_image_view<double> > &a0,
                   const super_res_params &srp)
{
  const double sf_2 = 1.0/(srp.scale_factor * srp.scale_factor);

  vil_image_view<double> super_q(srp.s_ni, srp.s_nj, u.nplanes());
  for (unsigned int i = 0; i < q.size(); i++)
  {
    if (i == srp.ref_frame)
      continue;

    vil_image_view<double> work;
    // apply transpose linear operator to upsample, blur, and warp
    warps[i].apply_At(q[i], super_q);

    super3d::backward_divergence(b0[i], work);
    vil_math_add_image_fraction(work, srp.tau, super_q, -srp.tau * sf_2);
    vil_math_image_sum(a0[i], work, a0[i]);

    //elementwise product
    //vil_math_image_product(super_q, u, super_q);
    //super3d::backward_divergence(b1[i], work);
    //vil_math_add_image_fraction(work, srp.tau, super_q, -srp.tau * sf_2);
    //vil_math_image_sum(a1[i], work, a1[i]);
  }
}

//*****************************************************************************

void save_image(const vil_image_view<double> &img, const char *filename)
{
  double minv, maxv;
  vil_image_view<vxl_byte> output;
  vil_math_value_range(img, minv, maxv);
  vil_image_view<double> outd;
  vil_convert_stretch_range_limited(img, outd, vcl_max(0.0,minv), vcl_min(1.0,maxv), 0.0, 255.0);
  vil_convert_cast(outd, output);
  vil_save(output, filename);
}

//*****************************************************************************

void upsamplep(const vil_image_view<double> &src, vil_image_view<double> &dest,
              double scale_factor, super3d::warp_image_parameters::interp_type interp)
{
  super3d::warp_image_parameters wip;
  wip.set_fill_unmapped(true);
  wip.set_unmapped_value(-1.0);
  wip.set_interpolator(interp);

  vnl_double_3x3 Sinv;
  Sinv.set_identity();
  Sinv(0,0) = 1.0/scale_factor;
  Sinv(1,1) = 1.0/scale_factor;

  dest.set_size(src.ni() * scale_factor, src.nj() * scale_factor, src.nplanes());
  super3d::warp_image(src, dest, Sinv, wip);
}

} // end anonymous namespace

//*****************************************************************************

void super_resolve(const vcl_vector<vil_image_view<double> > &frames,
                   const vcl_vector<super3d::adjoint_image_ops_func<double> > &warps,
                   vil_image_view<double> &u,
                   const super_res_params &srp,
                   unsigned int iterations,
                   const vcl_string &output_image,
                   super_res_monitor *srm)
{
  if (frames.empty())
    return;

  unsigned int np = frames[0].nplanes();
  for (unsigned int i = 1; i < frames.size(); i++)
  {
    if (frames[i].nplanes() != np)
    {
      vcl_cerr << "All frames must have the same number of planes.\n";
      return;
    }
  }

  //If the size was correct assume that it was pre-initialized
  if (u.ni() != srp.s_ni || u.nj() != srp.s_nj)
  {
    u.set_size(srp.s_ni, srp.s_nj, np);
    u.fill(0.0);
  }

  vil_image_view<double> p(srp.s_ni, srp.s_nj, 2*np);
  p.fill(0.0);

  vcl_vector<vil_image_view<double> > q(frames.size());
  vcl_vector<vil_image_view<double> > a0(frames.size());
  vcl_vector<vil_image_view<double> > b0(frames.size());
  vcl_vector<vil_image_view<double> > weights;
  for (unsigned int i = 0; i < frames.size(); i++)
  {
    vcl_cout << warps[i].dst_ni() << " " << warps[i].dst_nj() << "\n";
    q[i].set_size(warps[i].dst_ni(), warps[i].dst_nj(), np);
    q[i].fill(0.0);
    weights.push_back(warps[i].weight_map());
    a0[i].set_size(u.ni(), u.nj(), 1);
    a0[i].fill(0.0);
    b0[i].set_size(u.ni(), u.nj(), 2);
    b0[i].fill(0.0);
  }

  double ssd;
  vil_image_view<double> u_bar, last;
  u_bar.deep_copy(u);
  last.deep_copy(u);
  unsigned int i = 0;

  do
  {
    vcl_cout << "Iteration: " << i;
    if (srp.tv_method == super3d::super_res_params::IMAGEDATA_IMAGEPRIOR_ILLUMINATIONPRIOR)
    {
      dual_step_grad_huber(u, p, srp.lambda, srp.sigma, srp.epsilon_reg);
      dual_step_b(a0, b0, srp);
      dual_step_q(frames, warps, weights, u, a0, q, srp);
      primal_step_u(q, warps, p, u, u_bar, srp);
      primal_step_a(q, warps, u, b0, a0, srp);
    }
    else
    {
      dual_step_grad_huber(u, p, srp.lambda, srp.sigma, srp.epsilon_reg);
      dual_step_q(frames, warps, weights, u, q, srp);
      primal_step_u(q, warps, p, u, u_bar, srp);
    }

    if (srm)
    {
      if (*srm->interrupted_)
      {
        return;
      }

      if (srm->callback_ && !(i % srm->interval_))
      {
        super_res_monitor::update_data data;
        data.current_result.deep_copy(u);
        data.num_iterations = i;
        srm->callback_(data);
      }
    }
    else
    {
      double minv, maxv;
      vil_math_value_range(u, minv, maxv);
      if (!output_image.empty() && !(i % 15))
      {
        save_image(u, output_image.c_str());

        if (srp.tv_method == super3d::super_res_params::IMAGEDATA_IMAGEPRIOR_ILLUMINATIONPRIOR)
        {
          char buf[40];
          for (unsigned int i = 0; i < frames.size(); i++)
          {
            sprintf(buf, "a0-%d.png", i);
            save_image(a0[i], buf);
          }
        }
      }

      ssd = vil_math_ssd(last, u, double());
      vcl_cout << " SSD: " << ssd << " " << minv << " " << maxv;
      vcl_swap(last, u);
    }

    vcl_cout << "\n";
  } while (++i < iterations);
}


void compare_to_original(const vil_image_view<double> &ref_img,
                         const vil_image_view<double> &super,
                         const vil_image_view<double> &original,
                         unsigned int scale_factor)
{
  //Change to bicubic once warping becomes bicubic
  vil_image_view<double> upsample;
  vil_resample_bicub(ref_img, upsample, scale_factor * ref_img.ni(), scale_factor * ref_img.nj());

  vcl_cout << "SSD Super-res: " << vil_math_ssd(original, super, double()) << "\n";
  vcl_cout << "SSD Upsample: " << vil_math_ssd(original, upsample, double()) << "\n";

  vil_image_view<vxl_byte> output;
  vil_convert_stretch_range_limited(upsample, output, 0.0, 1.0);
  vil_save(output, "bilin.png");
  vil_convert_stretch_range_limited(original, output, 0.0, 1.0);
  vil_save(output, "original.png");
}


} // end namespace super3d

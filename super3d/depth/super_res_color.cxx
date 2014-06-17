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

#include "super_res.h"
#include "super_config.h"
#include "adjoint_image_op.h"
#include "image_derivs.h"


#include <vil/algo/vil_gauss_filter.h>
#include <vil/vil_save.h>
#include <vil/vil_math.h>
#include <vil/vil_resample_bilin.h>
#include <vil/vil_resample_bicub.h>
#include <vil/vil_convert.h>

#include <vil/vil_resample_bicub.txx>
VIL_RESAMPLE_BICUB_INSTANTIATE( double , double );

#define DEBUG


//*****************************************************************************

void dual_step_p(const vil_image_view<double> &u_bar,
                 vil_image_view<double> &p,
                 const super_res_params &srp)
{
  const double denom = 1.0 + (srp.sigma * srp.epsilon_reg) / srp.lambda;
  vil_image_view<double> work(srp.s_ni, srp.s_nj, p.nplanes());
  forward_gradient(u_bar, work);
  vil_math_add_image_fraction(p, 1.0/denom, work, srp.sigma/denom);

  for (unsigned int j = 0; j < srp.s_nj; j++)
  {
    for (unsigned int i = 0; i < srp.s_ni; i++)
    {
      for (unsigned int k = 0; k < u_bar.nplanes(); k++)
      {
        double &x = p(i,j,2*k), &y = p(i,j,2*k+1);

        //truncate vectors
        const double mag = sqrt(x*x + y*y)/srp.lambda;
        if (mag > 1.0f)
        {
          x /= mag;
          y /= mag;
        }
      }
    }
  }
}

//*****************************************************************************

void dual_step_q(const vcl_vector<vil_image_view<double> > &frames,
                 const vcl_vector<adjoint_image_ops_func<double> > &warps,
                 const vcl_vector<vil_image_view<double> > &weights,
                 const vil_image_view<double> &u_bar,
                 vcl_vector<vil_image_view<double> > &q,
                 const super_res_params &srp)
{
  vil_image_view<double> l_u, temp;
  const double sf_2 = 1.0/(srp.scale_factor * srp.scale_factor);
  const double denom = 1.0 + (srp.sigma * srp.epsilon_data) / sf_2;

  for (unsigned int f = 0; f < frames.size(); f++)
  {
    const unsigned int ni = weights[f].ni();
    const unsigned int nj = weights[f].nj();
    l_u.set_size(ni*srp.scale_factor, nj*srp.scale_factor, 1);
    // apply the linear operator to warp, blur, and downsample
    warps[f].apply_A(u_bar, l_u);

    for (unsigned int j = 0; j < nj; j++)
    {
      for (unsigned int i = 0; i < ni; i++)
      {
        for (unsigned int k = 0; k < u_bar.nplanes(); k++)
        {
          double &qfij = q[f](i, j, k);
          qfij = (qfij + srp.sigma * sf_2 * (l_u(i, j, k) - frames[f](i, j, k) * weights[f](i,j)))/denom;
          qfij = vcl_max(qfij, -sf_2);
          qfij = vcl_min(qfij, sf_2);
        }
      }
    }
  }
}

//*****************************************************************************

void primal_step_u(const vcl_vector<vil_image_view<double> > &q,
                   const vcl_vector<adjoint_image_ops_func<double> > &warps,
                   const vil_image_view<double> &p,
                   vil_image_view<double> &u,
                   vil_image_view<double> &u_bar,
                   const super_res_params &srp)
{
  const double sf_2 = 1.0/(srp.scale_factor * srp.scale_factor);
  vil_image_view<double> sum_super_q(srp.s_ni, srp.s_nj, 3);
  sum_super_q.fill(0.0);
  vil_image_view<double> super_q(srp.s_ni, srp.s_nj), temp;
  for (unsigned int i = 0; i < q.size(); i++)
  {
    // apply transpose linear operator to upsample, blur, and warp
    warps[i].apply_At(q[i], super_q);
    vil_math_image_sum(sum_super_q, super_q, sum_super_q);
  }

  vil_image_view<double> work(srp.s_ni, srp.s_nj, 3);
  backward_divergence(p, work);
  vil_math_add_image_fraction(work, srp.tau, sum_super_q, -srp.tau * sf_2);
  vil_math_image_sum(u, work, work);
  u_bar.deep_copy(work);
  vil_math_add_image_fraction(u_bar, 2.0, u, -1.0);
  vcl_swap(u, work);
}

//*****************************************************************************

void super_resolve_color(const vcl_vector<vil_image_view<double> > &frames,
                         const vcl_vector<adjoint_image_ops_func<double> > &warps,
                         vil_image_view<double> &u,
                         const super_res_params &srp,
                         unsigned int iterations)
{
  //If the size was correct assume that it was pre-initialized
  if (u.ni() != srp.s_ni || u.nj() != srp.s_nj)
  {
    u.set_size(srp.s_ni, srp.s_nj, 1);
    u.fill(0.0);
  }

  vil_image_view<double> p(srp.s_ni, srp.s_nj, 6);
  p.fill(0.0);

  vcl_vector<vil_image_view<double> > q(frames.size());
  vcl_vector<vil_image_view<double> > weights;
  for (unsigned int i = 0; i < frames.size(); i++)
  {
    q[i].set_size(warps[i].dst_ni(), warps[i].dst_nj(), 3);
    q[i].fill(0.0);
    weights.push_back(warps[i].weight_map());
  }

  double ssd;
  vil_image_view<double> u_bar, last;
  u_bar.deep_copy(u);
  last.deep_copy(u);
  unsigned int i = 0;
#ifdef DEBUG
  vil_image_view<vxl_byte> output;
#endif

  do
  {
    vcl_cout << "Iteration: " << i;
    dual_step_p(u, p, srp);
    dual_step_q(frames, warps, weights, u, q, srp);
    primal_step_u(q, warps, p, u, u_bar, srp);
    double minv, maxv;
    vil_math_value_range(u, minv, maxv);
#ifdef DEBUG
    if (!(i % 20))
    {
      vil_image_view<double> outd;
      vil_convert_stretch_range_limited(u, outd, vcl_max(0.0,minv), vcl_min(1.0,maxv), 0.0, 255.0);
      vil_convert_cast(outd, output);
      vil_save(output, config::inst()->get_value<vcl_string>("output_image").c_str());
    }
#endif

    ssd = vil_math_ssd(last, u, double());
    vcl_cout << " SSD: " << ssd << " " << minv << " " << maxv << "\n";
    vcl_swap(last, u);
  } while (++i < iterations);
}

//*****************************************************************************

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

#include "tv_refine_search.h"
#include "tv_utils.h"
#include "super_config.h"

#include <sstream>
#include <iomanip>
#include <limits>

#include <vil/vil_convert.h>
#include <vil/vil_save.h>
#include <vil/algo/vil_median.h>
#include <vil/vil_math.h>

#include <vcl_vector.h>

#include <vnl/vnl_double_3.h>
#include <vnl/vnl_double_2.h>
#include <super3d/image/dual_rof_denoise.h>


namespace super3d
{

void
refine_depth(vil_image_view<double> &cost_volume,
             const vil_image_view<double> &g,
             vil_image_view<double> &d,
             unsigned int iterations,
             double theta0,
             double theta_end,
             double lambda,
             double epsilon,
             depth_refinement_monitor *drm)
{
  vil_image_view<double> sqrt_cost_range(cost_volume.ni(), cost_volume.nj(), 1);
  double a_step = 1.0 / cost_volume.nplanes();

  for (unsigned int j = 0; j < cost_volume.nj(); j++)
  {
    for (unsigned int i = 0; i < cost_volume.ni(); i++)
    {
      double min, max;
      min = max = cost_volume(i,j,0);
      unsigned int min_k = 0;
      for (unsigned int k = 1; k < cost_volume.nplanes(); k++)
      {
        const double &cost = cost_volume(i,j,k);
        if (cost < min) {
          min = cost;
          min_k = k;
        }
        if (cost > max)
          max = cost;
      }
      sqrt_cost_range(i,j) = vcl_sqrt(max - min);
      d(i,j) = (min_k + 0.5) * a_step;
    }
  }

  vil_image_view<double> q(cost_volume.ni(), cost_volume.nj(), 2);
  q.fill(0.0);

  vil_image_view<double> a(cost_volume.ni(), cost_volume.nj(), 1);

  double theta = theta0;
  double denom = log(10.0);
  double orders = log(theta0)/denom - log(theta_end)/denom;
  double beta = orders / static_cast<double>(iterations);

  for (unsigned int iter = 1; iter <= iterations; iter++)
  {
    vcl_cout << "theta: " << theta << "\n";
    min_search_bound(a, d, cost_volume, sqrt_cost_range, theta, lambda);
    huber(q, d, a, g, theta, 0.25/theta, epsilon);
    theta = pow(10.0, log(theta)/denom - beta);

    if (drm)
    {
      if (*drm->interrupted_)
      {
        return;
      }

      if (drm->callback_ && !(iter % drm->interval_))
      {
        depth_refinement_monitor::update_data data;
        data.current_result.deep_copy(d);
        data.num_iterations = iter;
        drm->callback_(data);
      }
    }
  }
}


//semi-implicit gradient ascent on q and descent on d
void huber(vil_image_view<double> &q,
           vil_image_view<double> &d,
           const vil_image_view<double> &a,
           const vil_image_view<double> &g,
           double theta,
           double step,
           double epsilon)
{
  unsigned int ni = d.ni() - 1, nj = d.nj() - 1;
  double stepsilon1 = 1.0 + step*epsilon;
  for (unsigned int j = 0; j < nj; j++)
  {
    for (unsigned int i = 0; i < ni; i++)
    {
      double &x = q(i,j,0), &y = q(i,j,1);
      double dij = d(i,j);
      x = (x + step * g(i,j) * (d(i+1,j) - dij))/stepsilon1;
      y = (y + step * g(i,j) * (d(i,j+1) - dij))/stepsilon1;

      //truncate vectors
      double mag = x*x + y*y;
      if (mag > 1.0f)
      {
        mag = sqrt(mag);
        x /= mag;
        y /= mag;
      }
    }
  }

  q(ni,nj,0) = q(ni-1,nj,0);
  q(ni,nj,1) = q(ni,nj-1,1);

  double theta_inv = 1.0 / theta, denom = (1.0 + (step / theta));
  for (unsigned int j = 0; j < d.nj(); j++)
  {
    for (unsigned int i = 0; i < d.ni(); i++)
    {
      //add scaled divergence
      double divx = q(i,j,0), divy = q(i,j,1);
      if (i > 0)  divx -=  q(i-1,j,0);
      if (j > 0)  divy -=  q(i,j-1,1);

      double &dij = d(i,j);
      dij = (dij + step * (g(i,j) * (divx + divy) + theta_inv * a(i,j)))/denom;
    }
  }
}


//semi-implicit gradient ascent on q and descent on d
void huber_central(vil_image_view<double> &q,
           vil_image_view<double> &d,
           const vil_image_view<double> &a,
           const vil_image_view<double> &g,
           double theta,
           double step,
           double epsilon)
{
  unsigned int ni = d.ni() - 1, nj = d.nj() - 1;
  double stepsilon1 = 1.0 + step*epsilon;
  for (unsigned int j = 1; j < nj; j++)
  {
    for (unsigned int i = 1; i < ni; i++)
    {
      double &x = q(i,j,0), &y = q(i,j,1);
      x = (x + step * g(i,j) * (d(i+1,j) - d(i-1,j)))/stepsilon1;
      y = (y + step * g(i,j) * (d(i,j+1) - d(i,j-1)))/stepsilon1;

      //truncate vectors
      double mag = x*x + y*y;
      if (mag > 1.0f)
      {
        mag = sqrt(mag);
        x /= mag;
        y /= mag;
      }
    }
  }

  double theta_inv = 1.0 / theta, denom = (1.0 + (step / theta));
  for (unsigned int j = 0; j < d.nj(); j++)
  {
    for (unsigned int i = 0; i < d.ni(); i++)
    {
      //add scaled divergence
      double divx = q(i,j,0), divy = q(i,j,1);
      if (i > 0 && i < ni)  divx =  q(i+1,j,0) - q(i-1,j,0);
      if (j > 0 && j < nj)  divy =  q(i,j+1,1) - q(i,j-1,1);

      double &dij = d(i,j);
      dij = (dij + step * (g(i,j) * (divx + divy) + theta_inv * a(i,j)))/denom;
    }
  }
}


//semi-implicit gradient ascent on q and descent on d
void hessian_frob(vil_image_view<double> &q,
                  vil_image_view<double> &d,
                  const vil_image_view<double> &a,
                  const vil_image_view<double> &g,
                  double theta,
                  double step,
                  double epsilon)
{
  unsigned int ni = d.ni()-1, nj = d.nj()-1;
  double stepsilon1 = 1.0 + step*epsilon;

  for (unsigned int j = 1; j < nj; j++)
  {
    for (unsigned int i = 1; i < ni; i++)
    {
      double &xx = q(i,j,0), &xy = q(i,j,1), &yy = q(i,j,2);
      double dij = d(i,j);
      double scale = step * g(i,j);
      int ip1 = vcl_min(ni, i+1), jp1 = vcl_min(nj, j+1);
      int im1 = vcl_max(i-1, (unsigned int)0), jm1 = vcl_max(j-1, (unsigned int)0);
      xx = (xx + scale * (d(ip1,j) - 2.0 *dij + d(im1,j)))/stepsilon1;
      xy = (xy + scale * 0.25 * (d(ip1,jp1) - d(ip1,jm1) - d(im1,jp1) + d(im1,jm1)))/stepsilon1;
      yy = (yy + scale * (d(i,jm1) - 2.0*dij + d(i,jp1)))/stepsilon1;

      //truncate vectors
      double mag = xx*xx + 2*xy*xy + yy*yy;
      if (mag > 1.0f)
      {
        mag = sqrt(mag);
        xx /= mag;
        xy /= mag;
        yy /= mag;
      }
    }
  }

  double theta_inv = 1.0 / theta, denom = (1.0 + (step / theta));
  for (unsigned int j = 0; j < d.nj(); j++)
  {
    for (unsigned int i = 0; i < d.ni(); i++)
    {
      double qxx = q(i,j,0), qxy = q(i,j,1), qyy = q(i,j,2);
      double &dij = d(i,j);
      int ip1 = vcl_min(ni, i+1), jp1 = vcl_min(nj, j+1);
      int im1 = vcl_max(i-1, (unsigned int)0), jm1 = vcl_max(j-1, (unsigned int)0);
      qxx = - q(im1,j,0) + 2.0*qxx - q(ip1,j,0);
      qyy = -q(i,jm1,2) + 2.0*qyy - q(i,jp1,2);
      qxy = -0.25 * (q(ip1,jp1,1) - q(ip1,jm1,1) - q(im1,jp1,1) + q(im1,jm1,1));
      dij = (dij + step * (g(i,j) * (qxx + 2.0*qxy + qyy) + theta_inv * a(i,j)))/denom;
    }
  }
}


double eval_hessian_frob(const vil_image_view<double> &d,
                         const vil_image_view<double> &costvol,
                         double lambda)
{
  double cost = 0;
  unsigned int ni = d.ni() - 1, nj = d.nj() - 1;

  for (unsigned int i = 1; i < ni; i++)
  {
    for (unsigned int j = 1; j < nj; j++)
    {
      double ixx = d(i-1,j) - 2.0*d(i,j) + d(i+1,j);
      double iyy = d(i,j-1) - 2.0*d(i,j) + d(i,j+1);
      double ixy = 0.25 * (d(i+1,j+1) - d(i-1,j+1) - d(i+1,j-1) + d(i-1,j-1));
      double frob = sqrt(ixx*ixx + iyy*iyy + 2.0*ixy*ixy);
      int label = (int)(d(i,j) * costvol.nplanes());
      cost += frob + lambda * costvol(i,j,label);
    }
  }

  return cost;
}

} // end namespace super3d

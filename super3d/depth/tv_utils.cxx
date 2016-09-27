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

#include "tv_utils.h"

#include <algorithm>
#include <limits>


namespace super3d
{

void
min_search_bound(vil_image_view<double> &a,
           const vil_image_view<double> &d,
           const vil_image_view<double> &cost_volume,
           const vil_image_view<double> &sqrt_cost_range,
           double theta,
           double lambda)
{
  double coeff = (1.0/(2.0*theta));
  double range_coeff = std::sqrt(2.0 * theta * lambda);
  unsigned int S = cost_volume.nplanes();
  double a_step = 1.0 / S;

  for (unsigned int j = 0; j < d.nj(); j++)
  {
    for (unsigned int i = 0; i < d.ni(); i++)
    {
      double r = range_coeff * sqrt_cost_range(i,j);
      double dij = d(i,j);

      // compute the search range and clip between 0 and 1
      // note that when dij is outside the volume range [0,1] this
      // range needs to be shifted to the closest edge of the volume
      // for example if dij < 0 then search in [0, r] and
      // if dij > 0 search in [1-r, 1]
      double min_idepth = std::min(std::max(dij - r, 0.0), 1.0 - r);
      double max_idepth = std::min(std::max(dij + r, r), 1.0);

      unsigned int k = static_cast<unsigned int>(min_idepth * S);
      if (k >= S)
      {
        k = S-1;
      }
      unsigned int max_k = static_cast<unsigned int>(max_idepth * S);
      if (max_k >= S)
      {
        max_k = S-1;
      }
      unsigned int bestk = 0;
      double best_idepth = 0, best_e = std::numeric_limits<double>::infinity();
      for (; k <= max_k; ++k)
      {
        double aij = (k + 0.5) * a_step;
        double cost = cost_volume(i,j,k);
        if (cost < 0.0) continue;
        double diff = dij - aij;
        double e = coeff*diff*diff + lambda * cost;
        if (e < best_e)
        {
          best_e = e;
          best_idepth = aij;
          bestk = k;
         }
      }
      a(i,j) = subsample(dij, cost_volume, coeff, lambda, i, j, bestk, a_step, best_idepth);
    }
  }
}


void
min_search(vil_image_view<double> &a,
           const vil_image_view<double> &d,
           const vil_image_view<double> &cost_volume,
           double theta,
           double lambda)
{
  double coeff = (1.0/(2.0*theta));
  unsigned int S = cost_volume.nplanes();
  double a_step = 1.0 / S;
  for (unsigned int j = 0; j < d.nj(); j++)
  {
    for (unsigned int i = 0; i < d.ni(); i++)
    {
      double dij = d(i,j);
      unsigned int bestk = 0;
      double best_idepth = 0, best_e = std::numeric_limits<double>::infinity();
      for (unsigned int k = 0; k < S; ++k)
      {
        double aij = (k + 0.5) * a_step;
        double cost = cost_volume(i,j,k);
        if (cost < 0.0) continue;
        double diff = dij - aij;
        double e = coeff*diff*diff + lambda * cost;
        if (e < best_e)
        {
          best_e = e;
          best_idepth = aij;
          bestk = k;
         }
      }
      a(i,j) = subsample(dij, cost_volume, coeff, lambda, i, j, bestk, a_step, best_idepth);
    }
  }
}


double subsample(double dij,
                 const vil_image_view<double> &cost_volume,
                 double coeff,
                 double lambda,
                 unsigned int i,
                 unsigned int j,
                 unsigned int k,
                 double a_step,
                 double aij)
{
  double eval[3];
  const int S = static_cast<int>(cost_volume.nplanes());

  for (int m = -1; m <= 1; m++)
  {
    double aij_m = a_step * m + aij;

    int loc = k + m;
    if (loc < 0 || loc >= S)
      return aij;

    double diff = dij - aij_m;
    eval[m+1] = coeff*diff*diff + lambda * cost_volume(i,j,loc);
  }

  double dleft = (eval[1] - eval[0])/a_step;
  double dright = (eval[2] - eval[1])/a_step;

  double diff1 = (eval[2] - eval[0])/(2.0*a_step);
  double diff2 = (dright - dleft)/a_step;
  double delta = -diff1/diff2;
  return aij + delta;
}

} // end namespace super3d

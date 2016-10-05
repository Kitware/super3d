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

namespace {

/// Interpolate the offset to the subsampled minimum by fitting a parabola
/// This function fits a parabola to 3 points: (-1, ym1), (0, y0), (1, yp1)
/// and estimates the X location of the minimum.  It is assumed that y0 < ym1
/// and y0 < yp1.
inline
double
interp_offset(double ym1, double y0, double yp1)
{
  const double d1 = yp1 - ym1;
  const double d2 = 2 * y0 - ym1 - yp1;
  return d2 == 0.0 ? 0.0 : d1 / (2 * d2);
}

}


void
min_search_bound(vil_image_view<double> &a,
           const vil_image_view<double> &d,
           const vil_image_view<double> &cost_volume,
           const vil_image_view<double> &sqrt_cost_range,
           double theta,
           double lambda)
{
  const int S = static_cast<int>(cost_volume.nplanes());
  const double a_step = 1.0 / S;
  const int last_plane = S-1;
  const double coeff = (1.0 / (2.0 * theta * lambda * S * S));
  const double range_coeff = std::sqrt(2.0 * theta * lambda);

  const std::ptrdiff_t istep_c = cost_volume.istep();
  const std::ptrdiff_t jstep_c = cost_volume.jstep();
  const std::ptrdiff_t pstep_c = cost_volume.planestep();

  const double* row_c = cost_volume.top_left_ptr();
  for (unsigned int j = 0; j < d.nj(); j++, row_c+=jstep_c)
  {
    const double* col_c = row_c;
    for (unsigned int i = 0; i < d.ni(); i++, col_c+=istep_c)
    {
      const int r = std::min(last_plane, std::max(0, static_cast<int>(S * range_coeff * sqrt_cost_range(i,j))));
      const double dij = d(i,j) * S - 0.5;
      const int init_k = std::min(last_plane, std::max(0, static_cast<int>(dij)));

      // compute the search range and clip between 0 and S-1
      // note that when dij is outside the volume range [0,1] this
      // range needs to be shifted to the closest edge of the volume
      // for example if dij < 0 then search in [0, r] and
      // if dij > 1 search in [S-1-r, 1]
      const int min_k = std::max(0, init_k - r);
      const int max_k = std::min(last_plane, init_k + r);

      int bestk = init_k;
      const double diff = dij - bestk;
      double best_e = coeff*diff*diff + *(col_c + bestk);
      const double* cost = col_c + min_k;
      for (int k = min_k; k <= max_k; ++k, cost+=pstep_c)
      {
        if (k == init_k || *cost < 0.0 || *cost > best_e)
        {
          continue;
        }
        const double diff = dij - k;
        const double e = coeff*diff*diff + (*cost);
        if (e < best_e)
        {
          best_e = e;
          bestk = k;
        }
      }
      // fit a parabola to estimate the subsample offset for the best k
      if (bestk > 0 && bestk < S-1)
      {
        cost = col_c + bestk;
        const double diff2 = 2 * coeff * (dij - static_cast<double>(bestk));
        const double ym1 = *(cost - pstep_c) + diff2 + coeff;
        const double yp1 = *(cost + pstep_c) - diff2 + coeff;
        a(i,j) = (static_cast<double>(bestk) + interp_offset(ym1, *cost, yp1) + 0.5) * a_step;
      }
      else
      {
        a(i,j) = (static_cast<double>(bestk) + 0.5) * a_step;
      }
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
  double coeff = (1.0 / (2.0 * theta * lambda));
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
        double e = coeff*diff*diff + cost;
        if (e < best_e)
        {
          best_e = e;
          best_idepth = aij;
          bestk = k;
         }
      }
      a(i,j) = subsample(dij, cost_volume, coeff, i, j, bestk, a_step, best_idepth);
    }
  }
}


double subsample(double dij,
                 const vil_image_view<double> &cost_volume,
                 double coeff,
                 unsigned int i,
                 unsigned int j,
                 unsigned int k,
                 double a_step,
                 double aij)
{
  const int S = static_cast<int>(cost_volume.nplanes());
  if (k < 1 || k >= static_cast<unsigned int>(S-1))
  {
    return aij;
  }

  double eval[3];

  for (int m = -1; m <= 1; m++)
  {
    const double aij_m = a_step * m + aij;
    const double diff = dij - aij_m;
    eval[m+1] = coeff*diff*diff + cost_volume(i,j,k+m);
  }

  const double dleft = (eval[1] - eval[0]);
  const double dright = (eval[2] - eval[1]);

  const double diff1 = (eval[2] - eval[0]) * 0.5;
  const double diff2 = (dright - dleft);
  const double delta = -diff1/diff2;
  return aij + delta * a_step;
}

} // end namespace super3d

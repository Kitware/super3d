/*ckwg +29
 * Copyright 2012 by Kitware, Inc.
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

//Takes an int4 for dimensions because int3 crashes AMD's compiler
__kernel void init_depth(__global float * cost_vol,
                         __global float *d,
                         __global float * sqrt_cost_range,
                         int4 dims)
{
  int i = 1;
  int2 coord = (int2)(get_global_id(0), get_global_id(1));
  int ij = mad24(coord.y, dims.x, coord.x);
  int index = mul24(ij, dims.z);

  float min_v = cost_vol[index];
  float max_v = min_v;
  float min_k = 0;
  for ( ; i < dims.z; i++)
  {
    float val = cost_vol[index + i];
    if (val < min_v)
    {
      min_v = val;
      min_k =  i;
    }
    if (val > max_v)
    {
      max_v = val;
    }
  }

  sqrt_cost_range[ij] = sqrt(max_v - min_v);
  d[ij] = (min_k + 0.5f)/(float)dims.z;
}

//Takes an int4 for dimensions because int3 crashes AMD's compiler
__kernel void min_search(__global float *a,
                         __global float *d,
                         __global float *cost_volume,
                         __global float *sqrt_cost_range,
                         float theta,
                         float lambda,
                         int4 dims)
{
  float coeff = 1.0f / (2.0f * theta);
  int2 coord = (int2)(get_global_id(0), get_global_id(1));
  int ij = mad24(coord.y, dims.x, coord.x);
  int cv_ij = mul24(ij, dims.z);
  float r = sqrt(2.0f * theta * lambda) * sqrt_cost_range[ij];
  float depth = d[ij];
  float S = (float)dims.z;
  float a_step = 1.0f / S;

  float min_d = fmin(fmax(depth - r, 0.0f), 1.0f - r);
  float max_d = fmin(fmax(depth + r, r), 1.0f);

  int k = (int)(min_d * S);
  int maxk = (int)(max_d * S);

  int best_k = 0;
  float best_ad = 0.0f;
  float best_e = INFINITY;

  k = k >= S ? S-1 : k;
  maxk = maxk >= S ? S-1 : maxk;

  //Search range for least energy depth
  for ( ; k <= maxk; ++k)
  {
    float a_depth = (k + 0.5f) * a_step;
    float cost = cost_volume[cv_ij + k];
    if (cost >= 0.0f)
    {
      float diff = depth - a_depth;
      float e = coeff * diff * diff + lambda * cost;
      if (e < best_e)
      {
        best_e = e;
        best_ad = a_depth;
        best_k = k;
      }
    }
  }

  //Check if we can subsample, if not we are done
  if (best_k - 1 < 0 || best_k + 1 >= dims.z)
  {
    a[ij] = best_ad;
    return;
  }

  float3 eval;
  float3 a_depths = (float3)(-a_step + best_ad, best_ad, a_step + best_ad);
  float3 diffs = (float3)depth - a_depths;
  float3 costs = (float3)(cost_volume[cv_ij + best_k - 1],
                          cost_volume[cv_ij + best_k],
                          cost_volume[cv_ij + best_k + 1]);
  eval = (coeff * diffs * diffs) + (lambda * costs);

  //Fit a parabola
  float dleft = (eval.y - eval.x)/a_step;
  float dright = (eval.z - eval.y)/a_step;

  float diff1 = (eval.z - eval.x)/(2.0f*a_step);
  float diff2 = (dright - dleft)/a_step;
  float delta = -diff1/diff2;

  a[ij] = best_ad + delta;
}

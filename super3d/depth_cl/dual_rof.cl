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

__constant sampler_t image_sampler = CLK_NORMALIZED_COORDS_FALSE |
                                     CLK_ADDRESS_CLAMP_TO_EDGE |
                                     CLK_FILTER_NEAREST;

__constant float2 zero2 = (float2)(0.0f, 0.0f);
__kernel void init_dual(__global float2 *dual, int dimx)
{
  int2 coord = (int2)(get_global_id(0), get_global_id(1));
  dual[mad24(coord.y, dimx, coord.x)] = zero2;
}

__kernel void gradient(__global float2 *dual,
                       __global float *dest,
                       __read_only image2d_t g,
                         float step,
                         float epsilon,
                         int2 dims)
{
  int2 coord = (int2)(get_global_id(0), get_global_id(1));
  int index = mad24(coord.y, dims.x, coord.x);

  float intensity = dest[index];

  float2 vec = dual[index];

  float gij = read_imagef(g, image_sampler, coord).x;
  float stepsilon = 1.0f + step * epsilon;
  vec.x = (vec.x + step * gij * (dest[index + 1] - intensity))/stepsilon;
  vec.y = (vec.y + step * gij * (dest[mad24(coord.y+1, dims.x, coord.x)] - intensity))/stepsilon;

  //truncate vectors
  float mag = vec.x*vec.x + vec.y*vec.y;
  vec = mag > 1.0f ? vec / mag : vec;

  dual[index] = vec;
}

__kernel void divergence(__global float *src,
                         __global float2 *dual,
                         __global float *dest,
                         __read_only image2d_t g,
                           float lambda,
                           float step,
                           int2 dims)
{
  int2 coord = (int2)(get_global_id(0), get_global_id(1));

  int index = mad24(coord.y, dims.x, coord.x);
  float a = src[index];
  float2 vec = dual[index];
  float gij = read_imagef(g, image_sampler, coord).x;

  //add scaled divergence
  float2 div;
  div.x = coord.x > 0 ? vec.x - dual[index - 1].x : vec.x;
  div.y = coord.y > 0 ? vec.y - dual[mad24(coord.y-1, dims.x, coord.x)].y : vec.y;

  dest[index] = (dest[index] + step * ( gij * (div.x + div.y) + lambda * a))/(1.0f + step * lambda);
}

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

#ifndef ROF_CL_
#define ROF_CL_

#define STRINGIFY(A) #A

const char *rof_src = STRINGIFY(
__constant sampler_t image_sampler = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_CLAMP | CLK_FILTER_NEAREST;


/* Copy input 2D image to output 2D image */
__kernel void gradient(__global float *dual, __global float *dest, float scale, int2 dims)
{
  int2 coord = (int2)(get_global_id(0), get_global_id(1));

  float intensity = dest[dims.x * coord.y + coord.x];
  int size = dims.x * dims.y;

  float2 vec = (float2)(dual[dims.x * coord.y + coord.x], dual[dims.x * coord.y + coord.x + size]);

  if (coord.x < dims.x - 1)
    vec.x += scale * (dest[dims.x * coord.y + coord.x + 1] - intensity);
  if (coord.y < dims.y - 1)
    vec.y += scale * (dest[dims.x * (coord.y + 1) + coord.x] - intensity);

  //truncate vectors
  float mag = vec.x*vec.x + vec.y*vec.y;
  if (mag > 1.0f)
  {
    mag = sqrt(mag);
    vec.x /= mag;
    vec.y /= mag;
  }

  //barrier
  dual[dims.x * coord.y + coord.x] = vec.x;
  dual[dims.x * coord.y + coord.x + size] = vec.y;
}

  /* Copy input 2D image to output 2D image */
__kernel void divergence(__global image2d_t src, __global float *dual, __global float *dest, float theta, int2 dims)
{
  int2 coord = (int2)(get_global_id(0), get_global_id(1));
  int size = dims.x * dims.y;

  float original = read_imagef(src, image_sampler, coord).x;
  float2 vec = (float2)(dual[dims.x * coord.y + coord.x], dual[dims.x * coord.y + coord.x + size]);

  //add scaled divergence
  float2 div;
  if (coord.x > 0)
    div.x = vec.x - dual[dims.x * coord.y + coord.x - 1];
  else
    div.x = vec.x;

  if (coord.y > 0)
    div.y = vec.y - dual[dims.x * (coord.y-1) + coord.x + size];
  else
    div.y = vec.y;

  dest[dims.x * coord.y + coord.x] = original + theta*(div.x + div.y);
}


);

#endif

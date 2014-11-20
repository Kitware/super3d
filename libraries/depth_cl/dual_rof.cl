/*ckwg +5
 * Copyright 2012 by Kitware, Inc. All Rights Reserved. Please refer to
 * KITWARE_LICENSE.TXT for licensing information, or contact General Counsel,
 * Kitware, Inc., 28 Corporate Drive, Clifton Park, NY 12065.
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

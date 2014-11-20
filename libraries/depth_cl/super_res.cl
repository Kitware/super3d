/*ckwg +5
 * Copyright 2013 by Kitware, Inc. All Rights Reserved. Please refer to
 * KITWARE_LICENSE.TXT for licensing information, or contact General Counsel,
 * Kitware, Inc., 28 Corporate Drive, Clifton Park, NY 12065.
 */

typedef struct {
  float lambda, epsilon_data, epsilon_reg, tau, sigma;
  int2 sdim, ldim;
  float scale_factor;
  float sensor_sigma;

  float dual_p_denom;
  float dual_q_denom;
  float sf_2;
} super_res_params;

//*****************************************************************************

__kernel void zero(__global float *buf)
{
  buf[get_global_id(0)] = 0.0f;
}

//*****************************************************************************

__kernel void zero2(__global float2 *buf)
{
  buf[get_global_id(0)] = (float2)(0.0f, 0.0f);
}

//*****************************************************************************

__kernel void dual_step_p(__global float2 *p,
                          __global float *u,
                          __global super_res_params *srp)
{
  int2 coord = (int2)(get_global_id(0), get_global_id(1));
  int index = mad24(coord.y, srp->sdim.x, coord.x);

  float intensity = u[index];
  float2 vec = p[index];

  vec.x = vec.x + srp->sigma * (u[index + 1] - intensity);
  vec.y = vec.y + srp->sigma * (u[index + srp->sdim.x] - intensity);
  vec = vec / srp->dual_p_denom;

  //truncate vectors
  float mag = sqrt(dot(vec, vec)) / srp->lambda;
  vec = mag > 1.0f ? vec / mag : vec;

  p[index] = vec;
}

//*****************************************************************************

__kernel void dual_step_q(__global float *qi,
                          __global float *dbwi_u,
                          __global float *framei,
                          __global super_res_params *srp,
                          int2 dims)
{
  int2 coord = (int2)(get_global_id(0), get_global_id(1));
  int index = mad24(coord.y, dims.x, coord.x);
  float qval = qi[index];
  qval = (qval + srp->sigma * srp->sf_2 * (dbwi_u[index] - framei[index]))/srp->dual_q_denom;
  qval = min(qval, srp->sf_2);
  qi[index] = max(qval, -srp->sf_2);
}

//*****************************************************************************

__kernel void primal_step_u(__global float *sum_super_q,
                            __global float *u,
                            __global float2 *p,
                            __global super_res_params *srp,
                            int2 dims)
{
  int2 coord = (int2)(get_global_id(0), get_global_id(1));
  int index = mad24(coord.y, dims.x, coord.x);
  float2 vec = p[index];

  //add scaled divergence
  float2 div;
  div.x = coord.x > 0 ? vec.x - p[index - 1].x : vec.x;
  div.y = coord.y > 0 ? vec.y - p[index - srp->sdim.x].y : vec.y;

  u[index] = u[index] - srp->tau * (-(div.x + div.y) + srp->sf_2 * sum_super_q[index]);
}

//*****************************************************************************

__kernel void blur1D_horiz(__global float *input,
                         __global const float *filter,
                         int radius,
                         int2 dims,
                         __global float *output)
{
  int2 coord = (int2)(get_global_id(0), get_global_id(1));
  int index = mad24(coord.y, dims.x, coord.x);

  float temp = 0.0f;
  int i = 0;

  coord.x = coord.x - radius;
  int cur = index - radius;
  int end = index + radius;

  for( ; cur <= end; cur++, i++, coord.x++)
  {
    if (coord.x >= 0 || coord.x < dims.x)
      temp += filter[i] * input[cur];
  }

  input[index] = temp;
}

//*****************************************************************************

__kernel void blur1D_vert(__global float *input,
                          __global const float *filter,
                          int radius,
                          int2 dims,
                          __global float *output)
{
  int2 coord = (int2)(get_global_id(0), get_global_id(1));
  int index = mad24(coord.y, dims.x, coord.x);

  float temp = 0.0f;
  int i = 0;

  coord.y = coord.y - radius;
  int cur = index - radius * dims.x;
  int end = index + radius * dims.x;

  for( ; cur <= end; cur=cur+dims.x, i++, coord.y++)
  {
    if (coord.y >= 0 || coord.y < dims.y)
      temp += filter[i] * input[cur];
  }

  input[index] = temp;
}

//*****************************************************************************

//indexed on the output size
__kernel void down_sample(__global float *input,
                         int2 dims_s,
                         __global float *output,
                         int2 dims_l,
                         int scale)
{
  int2 coord_l = (int2)(get_global_id(0), get_global_id(1));
  int2 coord_s = scale * coord_l;

  coord_s.x = coord_s.x < dims_s.x ? coord_s.x : dims_s.x - 1;
  coord_s.y = coord_s.x < dims_s.y ? coord_s.y : dims_s.y - 1;

  int index = mad24(coord_s.y, dims_s.x, coord_s.x);
  int outdex = mad24(coord_l.y, dims_l.x, coord_l.x);

  output[outdex] = input[index];
}

//*****************************************************************************

//indexed on the input size, requires zero to be called on output before
//calling this function
__kernel void up_sample(__global float *input,
                       int2 dims_l,
                       __global float *output,
                       int2 dims_s,
                       int scale)
{
  int2 coord_l = (int2)(get_global_id(0), get_global_id(1));
  int2 coord_s = scale * coord_l;

  coord_s.x = coord_s.x < dims_s.x ? coord_s.x : dims_s.x - 1;
  coord_s.y = coord_s.x < dims_s.y ? coord_s.y : dims_s.y - 1;

  int index = mad24(coord_l.y, dims_l.x, coord_l.x);
  int outdex = mad24(coord_s.y, dims_s.x, coord_s.x);

  output[outdex] = input[index];
}

//*****************************************************************************

inline float linear_interp(__global float *src, int2 dims, float2 v)
{
  int i0 = (int)floor(v.x - 0.5f);
  int j0 = (int)floor(v.y - 0.5f);
  int i1 = i0 + 1;
  int j1 = j0 + 1;

  int index = mad24(j0, dims.x, i0);

  float tmp;
  float a = fract(v.x - 0.5f, &tmp);
  float b = fract(v.y - 0.5f, &tmp);
  float t00 = 0.0f, t10 = 0.0f, t01 = 0.0f, t11 = 0.0f;
  if (i0 >= 0 && i0 < dims.x)
  {
    if (j0 >= 0 && j0 < dims.y)
    {
      t00 = src[index];
    }
    if (j1 >= 0 && j1 < dims.y)
    {
      t01 =  src[index + dims.x];
    }
  }
  if (i1 >= 0 && i1 < dims.x)
  {
    if (j0 >= 0 && j0 < dims.y)
    {
      t10 = src[index + 1];
    }
    if (j1 >= 0 && j1 < dims.y)
    {
      t11 = src[index + dims.x + 1];
    }
  }

  return (1.0f - a) * (1.0f - b) * t00
         + a * (1.0f - b) * t10
         + (1.0f - a) * b * t01
         + a * b * t11;
}

//*****************************************************************************

//Given flow from I0 -> I1 takes values from I1 and puts them in I0
__kernel void warp_backward(__global float *I0,
                            __global float *I1,
                            __global float2 *flow,
                            int2 dims_I0,
                            int2 dims_I1)
{
  int2 coord = (int2)(get_global_id(0), get_global_id(1));
  int index = mad24(coord.y, dims_I0.x, coord.x);
  float2 flow_vec = flow[index];

  if (any(isnan(flow_vec)))
  {
    I0[index] = 0.0f;
  }
  else
  {
    float2 v = flow_vec + convert_float2(coord);
    I0[index] = linear_interp(I1, dims_I1, v);
  }
}

//*****************************************************************************

inline float atomic_add_float(volatile __global float *source, const float x)
{
  typedef union
  {
    unsigned int uint_val;
    float float_val;
  } uintflt;

  uintflt val, old;

  do
  {
    old.float_val = *source;
    val.float_val = old.float_val + x;
  } while (atomic_cmpxchg((volatile __global unsigned int *)source, old.uint_val, val.uint_val) != old.uint_val);

  return old.float_val;
}

//*****************************************************************************

//Given flow from I0 -> I1, takes values from I0 and splats them into I1
__kernel void warp_foward(__global float *I0,
                          __global float *I1,
                          __global float2 *flow,
                          int2 dims_I0,
                          int2 dims_I1)
{
  int2 coord = (int2)(get_global_id(0), get_global_id(1));
  int index = mad24(coord.y, dims_I0.x, coord.x);

  float2 flow_vec = flow[index];

  if (any(isnan(flow_vec)))
  {
    return;
  }
  else
  {
    float val = I0[index];
    float2 v = flow_vec + convert_float2(coord);

    int i0 = (int)floor(v.x - 0.5f);
    int j0 = (int)floor(v.y - 0.5f);
    int i1 = i0 + 1;
    int j1 = j0 + 1;

    int outdex = mad24(j0, dims_I1.x, i0);

    float tmp;
    float a = fract(v.x - 0.5f, &tmp);
    float b = fract(v.y - 0.5f, &tmp);

    if (i0 >= 0 && i0 < dims_I1.x)
    {
      if (j0 >= 0 && j0 < dims_I1.y)
      {
        atomic_add_float(&I1[outdex], (1.0f - a) * (1.0f - b) * val);
      }
      if (j1 >= 0 && j1 < dims_I1.y)
      {
        atomic_add_float(&I1[outdex + dims_I1.x], (1.0f - a) * b * val);
      }
    }
    if (i1 >= 0 && i1 < dims_I1.x)
    {
      if (j0 >= 0 && j0 < dims_I1.y)
      {
        atomic_add_float(&I1[outdex + 1], (1.0f - a) * b * val);
      }
      if (j1 >= 0 && j1 < dims_I1.y)
      {
        atomic_add_float(&I1[outdex + dims_I1.x + 1], a * b * val);
      }
    }
  }
}

//*****************************************************************************

__kernel void combine_super_q(__global float *q,
                              __global float *sum)
{
  int index = mad24(get_global_id(1), dims_I0.x, get_global_id(0));
  atomic_add_float(sum[index], q[index]);
}

//*****************************************************************************

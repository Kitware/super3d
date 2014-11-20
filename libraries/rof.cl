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

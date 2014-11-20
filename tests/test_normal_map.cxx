/*ckwg +5
 * Copyright 2012 by Kitware, Inc. All Rights Reserved. Please refer to
 * KITWARE_LICENSE.TXT for licensing information, or contact General Counsel,
 * Kitware, Inc., 28 Corporate Drive, Clifton Park, NY 12065.
 */


#include <test_common.h>

#include <cstdlib>
#include <exception>
#include <string>
#include <limits>


#include "normal_map.h"
#include "depth_map.h"
#include <vil/vil_convert.h>
#include <vil/vil_crop.h>
#include <vil/vil_save.h>



static void run_test(std::string const& test_name);

int
main(int argc, char* argv[])
{
  if (argc != 2)
  {
    TEST_ERROR("Expected one argument");

    return EXIT_FAILURE;
  }

  std::string const test_name = argv[1];

  try
  {
    run_test(test_name);
  }
  catch (std::exception const& e)
  {
    TEST_ERROR("Unexpected exception: " << e.what());

    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

static void test_compare_direct_indirect();

void
run_test(std::string const& test_name)
{
  if (test_name == "compare_direct_indirect")
  {
    test_compare_direct_indirect();
  }
  else
  {
    TEST_ERROR("Unknown test: " << test_name);
  }
}


void render_sphere_depth(double radius,
                         const vgl_point_3d<double>& center,
                         const vpgl_perspective_camera<double>& cam,
                         vil_image_view<double>& depth,
                         vil_image_view<double>& normal)
{
  const unsigned int ni = depth.ni();
  const unsigned int nj = depth.nj();
  normal.set_size(ni, nj, 3);

  vgl_vector_3d<double> offset = vgl_point_3d<double>(cam.camera_center()) - center;
  vgl_vector_3d<double> pa = cam.principal_axis();
  normalize(pa);
  if(dot_product(pa,offset) > 0)
  {
    // sphere behind camera
    return;
  }
  double c = offset.sqr_length() - radius*radius;
  for(unsigned int j=0; j<nj; ++j)
  {
    for(unsigned int i=0; i<ni; ++i)
    {
      vgl_ray_3d<double> ray = cam.backproject_ray(vgl_point_2d<double>(i,j));
      vgl_vector_3d<double> n = ray.direction();
      double b = 2*dot_product(offset,n);
      double disc = b*b - 4*c;
      if (disc < 0.0)
      {
        // ray does not intersect sphere
        continue;
      }
      double d = (-b - sqrt(disc)) / 2.0;
      double dist = d * dot_product(pa,n);
      if( dist < depth(i,j))
      {
        depth(i,j) = dist;
        vgl_vector_3d<double> snorm = normalized(offset + d*n);
        normal(i,j,0) = snorm.x();
        normal(i,j,1) = snorm.y();
        normal(i,j,2) = snorm.z();
      }
    }
  }
}


void
make_test_depth_and_normals(const unsigned int ni, const unsigned int nj,
                            vil_image_view<double>& depth,
                            vil_image_view<double>& normal,
                            vpgl_perspective_camera<double>& camera)
{
  const double foclen = 2*ni;
  depth.set_size(ni, nj);
  normal.set_size(ni, nj, 3);
  depth.fill(25.0);
  normal.fill(0.0);
  vil_plane(normal, 2).fill(-1);
  vpgl_calibration_matrix<double> K(foclen, vgl_point_2d<double>(ni/2.0, nj/2.0));
  vgl_point_3d<double> center(0, 0, -20);
  camera = vpgl_perspective_camera<double>(K, center, vgl_rotation_3d<double>());

  // render a bunch of spheres
  render_sphere_depth(3, vgl_point_3d<double>(0,0,0), camera, depth, normal);
  render_sphere_depth(2, vgl_point_3d<double>(3,3,0), camera, depth, normal);
  render_sphere_depth(1, vgl_point_3d<double>(2,-2,1), camera, depth, normal);
  render_sphere_depth(1, vgl_point_3d<double>(-2,2,-2), camera, depth, normal);
  render_sphere_depth(0.25, vgl_point_3d<double>(-1,-1,-4), camera, depth, normal);
}


void
test_compare_direct_indirect()
{
  const unsigned ni=200, nj=200;
  vil_image_view<double> depth, normal;
  vpgl_perspective_camera<double> cam;
  make_test_depth_and_normals(ni, nj, depth, normal, cam);

  // compute normals in world coordinates
  vil_image_view<double> location;
  vil_image_view<double> normal_wld;
  depth_map_to_location_map(cam, depth, location);
  location_map_to_normal_map(location, normal_wld);

  // compute normals directly in depth coordinates
  vil_image_view<double> normal_img;
  depth_map_to_normal_map_inv_len(cam, depth, normal_img);

  vil_image_view<double> dp;
  dot_product_map(normal_wld, normal_img, dp);

  // ignore a 1 pixel border where the normals are not computed
  dp = vil_crop(dp, 1, ni-2, 1, nj-2);

  vil_math_scale_and_offset_values(dp, -1.0, 1.0);
  double error = 0;
  vil_math_sum(error, dp, 0);
  if( error > 1e-8 )
  {
    TEST_ERROR("normals computations method weighted by inverse length do not agree");
  }

#if 1
  // write out the test images
  double min_d, max_d;
  finite_value_range(depth, min_d, max_d);
  vil_image_view<vxl_byte> byte_img;
  vil_convert_stretch_range_limited(depth, byte_img, min_d, max_d);
  vil_save(byte_img, "depth_true.png");

  byte_normal_map(normal, byte_img);
  vil_save(byte_img, "normal_true.png");

  byte_normal_map(normal_wld, byte_img);
  vil_save(byte_img, "normal_wld_inv_len.png");

  byte_normal_map(normal_img, byte_img);
  vil_save(byte_img, "normal_img_inv_len.png");
#endif
}

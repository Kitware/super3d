// This is brl/bbas/imesh/algo/imesh_render.h
#ifndef imesh_render_h_
#define imesh_render_h_
//:
// \file
// \brief Functions to render the faces of a mesh with projected coordinates
// \author Matt Leotta (mleotta@lems.brown.edu)
// \date June 2, 2008
//
// \verbatim
//  Modifications
//   <none yet>
// \endverbatim

#include "imesh_algo_config.h"
#include "../imesh_mesh.h"
#include "../imesh_operations.h"
#include "imesh_project.h"

#include <vil/vil_image_view.h>
#include <vgl/vgl_point_3d.h>
#include <vgl/vgl_vector_3d.h>
#include <vgl/vgl_triangle_scan_iterator.h>
#include <vpgl/vpgl_proj_camera.h>
#include <vcl_cassert.h>

//: Render a triangle defined by its vertices.
// For each 3d point, X and Y map to image coordinates and Z maps to depth.
// Assign label to each pixel in image such that the pixel is within
// the triangle and the interpolated depth is less than the value in depth_img.
// The depth is also updated for each pixel that is labeled
template <class T>
void imesh_render_triangle_label(const vgl_point_3d<double>& v1,
                                 const vgl_point_3d<double>& v2,
                                 const vgl_point_3d<double>& v3,
                                 const T& label,
                                 vil_image_view<T>& image,
                                 vil_image_view<double>& depth_img)
{
  assert(depth_img.ni() == image.ni());
  assert(depth_img.nj() == image.nj());

  vgl_triangle_scan_iterator<double> tsi;
  tsi.a.x = v1.x();  tsi.a.y = v1.y();
  tsi.b.x = v2.x();  tsi.b.y = v2.y();
  tsi.c.x = v3.x();  tsi.c.y = v3.y();
  vgl_vector_3d<double> b1(v2.x()-v1.x(), v2.y()-v1.y(), v2.z()-v1.z());
  vgl_vector_3d<double> b2(v3.x()-v1.x(), v3.y()-v1.y(), v3.z()-v1.z());
  vgl_vector_3d<double> n = cross_product(b1,b2);
  double A = -n.x()/n.z();
  double B = -n.y()/n.z();
  double C = (v1.x()*n.x() + v1.y()*n.y() + v1.z()*n.z())/n.z();
  for (tsi.reset(); tsi.next(); ) {
    int y = tsi.scany();
    if (y<0 || y>=int(image.nj())) continue;
    int min_x = tsi.startx();
    int max_x = tsi.endx();
    if (min_x >= (int)image.ni() || max_x < 0)
      continue;
    if (min_x < 0) min_x = 0;
    if (max_x >= (int)image.ni()) max_x = image.ni()-1;
    double new_i = B*y+C;
    for (int x = min_x; x <= max_x; ++x) {
      double depth = new_i + A*x;
      if (depth < depth_img(x,y)) {
        depth_img(x,y) = depth;
        image(x,y) = label;
      }
    }
  }
}


//: project the mesh onto the image plane using the camera
//  Set each pixel of the image to the depth to the mesh
template <class T>
void imesh_render_mesh_labels(const imesh_mesh& mesh,
                              const vpgl_proj_camera<double>& camera,
                              vil_image_view<T>& label_image,
                              vil_image_view<double>& depth_image)
{
  if (!mesh.faces().has_groups())
  {
    std::cerr << "Warning: rendering mesh labels, but mesh has no groups"
              << std::endl;
  }

  // clear the images by filling with the default values
  depth_image.fill(std::numeric_limits<double>::infinity());
  label_image.fill(T(0));

  // project the mesh vertices into image space points and depths
  assert(mesh.vertices().dim() == 3);
  std::vector<vgl_point_2d<double> > verts2d;
  std::vector<double> depths;
  imesh_project_verts(mesh.vertices<3>(), camera, verts2d, depths);

  // combine 2D points and depth values into 3D points (in camera coordinates)
  std::vector<vgl_point_3d<double> > verts3d;
  for(unsigned i=0; i<verts2d.size(); ++i)
  {
    const vgl_point_2d<double>& pt = verts2d[i];
    verts3d.push_back(vgl_point_3d<double>(pt.x(), pt.y(), depths[i]));
  }

  // convert mesh faces to triangles, if needed
  const imesh_face_array_base& faces = mesh.faces();
  std::unique_ptr<imesh_regular_face_array<3> > tri_data;
  const imesh_regular_face_array<3>* tris;
  if (faces.regularity() != 3)
  {
    tri_data = imesh_triangulate(faces);
    tris = tri_data.get();
  }
  else
  {
    tris = static_cast<const imesh_regular_face_array<3>*>(&faces);
  }

  // render each triangle into the label and depth images
  // the label is equal to the group number + 1 (zero means no label)
  int group = -1;
  if (tris->has_groups())
    group = 0;
  for (unsigned t=0; t<tris->size(); ++t)
  {
    const imesh_regular_face<3>& tri = (*tris)[t];
    imesh_render_triangle_label(verts3d[tri[0]], verts3d[tri[1]], verts3d[tri[2]],
                                static_cast<T>(group+1), label_image, depth_image);
    if (group >= 0 && t+1 >= tris->groups()[group].second)
    {
      ++group;
    }
  }
}


//: Render a textured triangle defined by its vertices
// \p v1,v2,v3 are coordinates in the projected image (plus depth)
// \p t1,t2,t3 are corresponding texture coordinates (in the unit square)
SUPER3D_IMESH_ALGO_EXPORT
void imesh_render_triangle_texture(const vgl_point_3d<double>& v1,
                                   const vgl_point_3d<double>& v2,
                                   const vgl_point_3d<double>& v3,
                                   const vgl_point_2d<double>& t1,
                                   const vgl_point_2d<double>& t2,
                                   const vgl_point_2d<double>& t3,
                                   const vil_image_view<vxl_byte>& texture,
                                   vil_image_view<vxl_byte>& image,
                                   vil_image_view<double>& depth_img);


//: Render the mesh using the camera and a texture image
//  A depth map is also computed and used for occlusion.
//  Texture mapping uses interpolates from the texture image with no
//  additional lighting calculations.
SUPER3D_IMESH_ALGO_EXPORT
void imesh_render_textured(const imesh_mesh& mesh,
                           const vpgl_proj_camera<double>& camera,
                           const vil_image_view<vxl_byte>& texture,
                           vil_image_view<vxl_byte>& image,
                           vil_image_view<double>& depth_img);

#endif // imesh_render_h_

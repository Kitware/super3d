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

#include "depth_map.h"
#include "world_space.h"

#include <vnl/vnl_double_3.h>
#include <vnl/vnl_inverse.h>
#include <vil/algo/vil_gauss_filter.h>
#include <vil/vil_decimate.h>
#include <vil/vil_resample_bilin.h>
#include <vil/vil_load.h>
#include <vil/vil_save.h>
#include <vil/vil_math.h>
#include <vil/vil_convert.h>
#include <vil/algo/vil_threshold.h>
#include <vil/vil_bilin_interp.h>

#ifdef HAVE_VTK
#include <vtkXMLPolyDataWriter.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkFloatArray.h>
#include <vtkSmartPointer.h>
#include <vtkCellArray.h>
#include <vtkPointData.h>
#endif


namespace super3d
{

/// Compute the average of the finite values in the image.
/// \param img The image. It may contain some infinite values.
/// \return The average of finite depths
double average_finite(const vil_image_view<double>& img)
{
  const unsigned ni = img.ni();
  const unsigned nj = img.nj();
  assert(img.nplanes() == 1);

  std::ptrdiff_t istep=img.istep(),  jstep=img.jstep();
  double avg = 0.0;
  unsigned count = 0;

  const double* row = img.top_left_ptr();
  for (unsigned j=0; j<nj; ++j, row += jstep)
  {
    const double* pixel = row;
    for (unsigned i=0; i<ni; ++i, pixel+=istep)
    {
      if (vnl_math::isfinite(*pixel))
      {
        avg += *pixel;
        ++count;
      }
    }
  }
  return avg / count;
}


/// Fill in infinite values in the image with constant \a value.
/// \retval img The image containing some infinite values
///             (modified in place)
/// \param value The value used to replace infinite values.
void fill_infinite(vil_image_view<double>& img, double value)
{
  const unsigned ni = img.ni();
  const unsigned nj = img.nj();
  assert(img.nplanes() == 1);

  std::ptrdiff_t istep=img.istep(),  jstep=img.jstep();

  double* row = img.top_left_ptr();
  for (unsigned j=0; j<nj; ++j, row += jstep)
  {
    double* pixel = row;
    for (unsigned i=0; i<ni; ++i, pixel+=istep)
    {
      if (!vnl_math::isfinite(*pixel))
      {
        *pixel = value;
      }
    }
  }
}


/// Copy finite values from \a src into \a dest
/// For infinite pixels in \a src, the matching pixels in \a dst are unchanaged
/// \retval src The source image containing finite and infinite values
/// \param dest The destination image.
void copy_finite(const vil_image_view<double>& src,
                 vil_image_view<double>& dest)
{
  const unsigned ni = src.ni();
  const unsigned nj = src.nj();
  assert(src.nplanes() == 1);
  assert(dest.ni() == ni);
  assert(dest.nj() == nj);
  assert(dest.nplanes() == 1);

  std::ptrdiff_t istepS=src.istep(),  jstepS=src.jstep();
  std::ptrdiff_t istepD=dest.istep(), jstepD=dest.jstep();

  const double* rowS = src.top_left_ptr();
  double*       rowD = dest.top_left_ptr();
  for (unsigned j=0; j<nj; ++j, rowS += jstepS, rowD += jstepD)
  {
    const double* pixelS = rowS;
    double*       pixelD = rowD;
    for (unsigned i=0; i<ni; ++i, pixelS+=istepS, pixelD+=istepD)
    {
      if (vnl_math::isfinite(*pixelS))
      {
        *pixelD = *pixelS;
      }
    }
  }
}


/// Fill in infinite values in \a src by diffusing from adjacent finite values.
/// The output \a dest is equal to \a src at all finite values of \a src
/// and has diffused finite values where ever \a src in infinite.
/// If \a dest is not initialized it will be filled with the average finite
/// depth computed from \a src.
/// \param src The source image containing finite and infinite values.
/// \retval dest The destination image (modified in place).
/// \param iterations The number of diffusion iterations
void diffuse_into_infinite(const vil_image_view<double>& src,
                           vil_image_view<double>& dest,
                           unsigned iterations)
{
  vil_gauss_filter_5tap_params params(1.0);
  vil_image_view<double> work;

  if (dest.ni() != src.ni() || dest.nj() != src.nj())
  {
    dest.set_size(src.ni(),src.nj());
    dest.fill(average_finite(src));
  }

  for (unsigned k=0; k<iterations; ++k)
  {
    // This code could be optimized by only computing the gauss filtering
    // at the infinite pixel of the source image.
    // This optimization would skip unnecessary filtering and remove the
    // need to copy the source pixels at each iteration.
    copy_finite(src,dest);
    vil_gauss_filter_5tap(dest, dest, params, work);
  }
  copy_finite(src,dest);
}


/// Fill in missing (infinite) depths in the depth map.
/// Initially missing depths are filled with the average of finite depths.
/// Then the know depths are defused into the unknown depths using
/// a multiscale application of the heat equation.
/// \retval depth_map The depth map containing some infinite values
///                   (modified in place)
/// \param levels The number of 2x pyramid levels to use
/// \param iterations The number of iterations to use at each le
void fill_missing_depths(vil_image_view<double>& depth_map,
                         unsigned levels,
                         unsigned iterations)
{
  const unsigned ni = depth_map.ni();
  const unsigned nj = depth_map.nj();
  assert(depth_map.nplanes() == 1);

  double avg_depth = average_finite(depth_map);

  unsigned scale = 1 << levels;
  vil_image_view<double> smooth(ni/scale,nj/scale);
  smooth.fill(avg_depth);

  for(; scale > 0; scale /= 2)
  {
    diffuse_into_infinite(vil_decimate(depth_map,scale), smooth, iterations);

    if (scale > 1)
    {
      vil_image_view<double> smooth2;
      vil_resample_bilin(smooth, smooth2, ni/(scale/2), nj/(scale/2));
      smooth = smooth2;
    }
  }

  depth_map = smooth;
}


/// Create dense mesh vertices by back projecting pixels to their given depth.
/// \param camera The camera to use for back projection
/// \param depth_map The depths at which to project each pixel
/// \retval index_image Image of indices into the array of mesh vertices
///                     The index is unsigned(-1) for infinite depths.
/// \return An array of mesh vertices covering all finite depth pixels
std::auto_ptr<imesh_vertex_array<3> >
depth_map_to_vertices(const vpgl_perspective_camera<double>& camera,
                      const vil_image_view<double>& depth_map,
                      vil_image_view<unsigned>& index_image)
{
  unsigned ni = depth_map.ni(), nj = depth_map.nj();
  assert(depth_map.nplanes() == 1);

  index_image.set_size(ni,nj,1);

  // a special constant to represent missing pixels
  const unsigned missing = static_cast<unsigned>(-1);

  index_image.fill(missing);

  std::auto_ptr<imesh_vertex_array<3> > verts(new imesh_vertex_array<3>);

  vnl_matrix_fixed<double,3,3> M_inv = vnl_inverse(camera.get_matrix().extract(3,3));
  vnl_double_3 p4 = camera.get_matrix().get_column(3);

  std::ptrdiff_t istepD=depth_map.istep(),   jstepD=depth_map.jstep();
  std::ptrdiff_t istepI=index_image.istep(), jstepI=index_image.jstep();

  // compute vertices
  const double*   rowD = depth_map.top_left_ptr();
  unsigned*       rowI = index_image.top_left_ptr();
  for (unsigned j=0; j<nj; ++j, rowD += jstepD, rowI += jstepI)
  {
    const double* pixelD = rowD;
    unsigned*     pixelI = rowI;
    for (unsigned i=0; i<ni; ++i, pixelD+=istepD, pixelI+=istepI)
    {
      if (vnl_math::isfinite(*pixelD))
      {
        *pixelI = verts->size();
        vnl_double_3 pt = M_inv * ((*pixelD) * vnl_double_3(i,j,1) - p4);
        verts->push_back(vgl_point_3d<double>(pt[0], pt[1], pt[2]));
      }
    }
  }

  return verts;
}


/// Create dense mesh vertices from a height map image (e.g. GeoTIFF LIDAR)
/// \param height_map The heights at which to place each pixel
/// \retval index_image Image of indices into the array of mesh vertices
///                     The index is unsigned(-1) for non-finite heights.
/// \param z_scale Amount to scale the height values
/// \param x_scale Amount to scale the unit horizontal distance between pixels
/// \param y_scale Amount to scale the unit vertical distance between pixels
/// \return An array of mesh vertices covering all finite height pixels
std::auto_ptr<imesh_vertex_array<3> >
height_map_to_vertices(const vil_image_view<double>& height_map,
                       vil_image_view<unsigned>& index_image,
                       double z_scale,
                       double x_scale,
                       double y_scale)
{
  unsigned ni = height_map.ni(), nj = height_map.nj();
  assert(height_map.nplanes() == 1);

  index_image.set_size(ni,nj,1);

  // a special constant to represent missing pixels
  const unsigned missing = static_cast<unsigned>(-1);

  index_image.fill(missing);

  std::auto_ptr<imesh_vertex_array<3> > verts(new imesh_vertex_array<3>);

  std::ptrdiff_t istepH=height_map.istep(),  jstepH=height_map.jstep();
  std::ptrdiff_t istepI=index_image.istep(), jstepI=index_image.jstep();

  // compute vertices
  const double*   rowH = height_map.top_left_ptr();
  unsigned*       rowI = index_image.top_left_ptr();
  for (unsigned j=0; j<nj; ++j, rowH += jstepH, rowI += jstepI)
  {
    const double* pixelH = rowH;
    unsigned*     pixelI = rowI;
    for (unsigned i=0; i<ni; ++i, pixelH+=istepH, pixelI+=istepI)
    {
      if (vnl_math::isfinite(*pixelH))
      {
        *pixelI = verts->size();
        vgl_point_3d<double> pt(x_scale*i, -y_scale*j, z_scale*(*pixelH));
        verts->push_back(pt);
      }
    }
  }

  return verts;
}


/// Create a triangulation from an image of vertex indices
/// \param index_image Image of indices into an array of mesh vertices
///                    The index is unsigned(-1) for missing vertices
/// \return An array of mesh faces (triangles) using the indices
std::auto_ptr<imesh_regular_face_array<3> >
triangulate_index_image(const vil_image_view<unsigned>& index_image)
{
  unsigned ni = index_image.ni(), nj = index_image.nj();
  assert(index_image.nplanes() == 1);

  // a special constant to represent missing pixels
  const unsigned missing = static_cast<unsigned>(-1);

  std::auto_ptr<imesh_regular_face_array<3> > faces(new imesh_regular_face_array<3>);

  std::ptrdiff_t istepI=index_image.istep(), jstepI=index_image.jstep();

  // compute triangles
  const unsigned* rowI = index_image.top_left_ptr();
  for (unsigned j=0; j<nj-1; ++j, rowI += jstepI)
  {
    const unsigned* pixelI = rowI;
    for (unsigned i=0; i<ni-1; ++i, pixelI+=istepI)
    {
      const unsigned& p1 = *pixelI;
      const unsigned& p2 = *(pixelI + istepI);
      const unsigned& p3 = *(pixelI + jstepI);
      const unsigned& p4 = *(pixelI + istepI + jstepI);
      unsigned char code = ((p1 == missing) ? 0 : 1) +
                           ((p2 == missing) ? 0 : 2) +
                           ((p3 == missing) ? 0 : 4) +
                           ((p4 == missing) ? 0 : 8);
      switch (code)
      {
        case 15: // has all points
          faces->push_back(imesh_tri(p1,p3,p2));
          faces->push_back(imesh_tri(p3,p4,p2));
          break;
        case 14: // has p2, p3, p4
          faces->push_back(imesh_tri(p3,p4,p2));
          break;
        case 13: // has p1, p3, p4
          faces->push_back(imesh_tri(p1,p3,p4));
          break;
        case 11: // has p1, p2, p4
          faces->push_back(imesh_tri(p1,p4,p2));
          break;
        case 7: // has p1, p2, p3
          faces->push_back(imesh_tri(p1,p3,p2));
          break;
        default:
          break;
      }
    }
  }

  return faces;
}


/// Create a dense mesh by back projecting pixels to their given depth.
/// \param camera The camera to use for back projection
/// \param depth_map The depths at which to project each pixel
/// \return A dense mesh with one vertex per finite depth pixel
imesh_mesh
depth_map_to_mesh(const vpgl_perspective_camera<double>& camera,
                  const vil_image_view<double>& depth_map)
{
  vil_image_view<unsigned> index_img;

  std::auto_ptr<imesh_vertex_array<3> > verts =
      depth_map_to_vertices(camera, depth_map, index_img);

  std::auto_ptr<imesh_regular_face_array<3> > faces =
      triangulate_index_image(index_img);


  std::auto_ptr<imesh_face_array_base> nf(faces);
  std::auto_ptr<imesh_vertex_array_base > nv(verts);
  return imesh_mesh(nv,nf);
}


/// Create a dense mesh from a height map image (e.g. GeoTIFF LIDAR)
/// \param height_map The heights at which to place each pixel
/// \param z_scale Amount to scale the height values
/// \param x_scale Amount to scale the unit horizontal distance between pixels
/// \param y_scale Amount to scale the unit vertical distance between pixels
/// \return A dense mesh with one vertex per finite height pixel
imesh_mesh
height_map_to_mesh(const vil_image_view<double>& height_map,
                   double z_scale,
                   double x_scale,
                   double y_scale)
{
  vil_image_view<unsigned> index_img;

  std::auto_ptr<imesh_vertex_array<3> > verts =
      height_map_to_vertices(height_map, index_img, z_scale, x_scale, y_scale);

  std::auto_ptr<imesh_regular_face_array<3> > faces =
      triangulate_index_image(index_img);


  std::auto_ptr<imesh_face_array_base> nf(faces);
  std::auto_ptr<imesh_vertex_array_base > nv(verts);
  return imesh_mesh(nv,nf);
}


/// Compute the minimum and maximum value ignoring infinite and NaN values.
void finite_value_range(const vil_image_view<double>& img,
                        double& min_value, double& max_value)
{
  min_value = 0;
  max_value = 0;

  if (img.size()==0)
  {
    return;
  }

  const unsigned ni = img.ni();
  const unsigned nj = img.nj();
  assert(img.nplanes() == 1);

  std::ptrdiff_t istep=img.istep(),  jstep=img.jstep();
  bool first_finite = true;

  const double* row = img.top_left_ptr();
  for (unsigned j=0; j<nj; ++j, row += jstep)
  {
    const double* pixel = row;
    for (unsigned i=0; i<ni; ++i, pixel+=istep)
    {
      if (vnl_math::isfinite(*pixel))
      {
        if (first_finite)
        {
          min_value = *pixel;
          max_value = *pixel;
          first_finite = false;
        }
        else if (*pixel<min_value)
        {
          min_value=*pixel;
        }
        else if (*pixel>max_value)
        {
          max_value=*pixel;
        }
      }
    }
  }
}


void save_depth(const vil_image_view<double> &depth, const char *filename)
{
  std::string ext(filename);
  ext = ext.substr(ext.find_last_of(".") + 1);
  std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
  if (ext == "tif" || ext == "tiff")
  {
    std::cout << "Saving "<<filename<<" as a TIFF"<<std::endl;
    vil_save(depth, filename);
    return;
  }

  FILE *file = fopen(filename, "wb");
  unsigned int ni = depth.ni(), nj = depth.nj();
  fwrite(&ni, sizeof(unsigned int), 1, file);
  fwrite(&nj, sizeof(unsigned int), 1, file);

  for (unsigned int i = 0; i < ni; i++)
  {
    for (unsigned int j = 0; j < nj; j++)
    {
      fwrite(&depth(i,j), sizeof(double), 1, file);
    }
  }

  fclose(file);
}


void load_depth(vil_image_view<double> &depth, const char *filename)
{
  std::string ext(filename);
  ext = ext.substr(ext.find_last_of(".") + 1);
  std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
  if (ext == "tif" || ext == "tiff")
  {
    std::cout << "Loading "<<filename<<" as a TIFF"<<std::endl;
    depth = vil_load(filename);
    return;
  }

  FILE *file = fopen(filename, "rb");
  unsigned int ni, nj;
  fread(&ni, sizeof(unsigned int), 1, file);
  fread(&nj, sizeof(unsigned int), 1, file);
  depth.set_size(ni, nj, 1);

  for (unsigned int i = 0; i < ni; i++)
  {
    for (unsigned int j = 0; j < nj; j++)
    {
      fread(&depth(i,j), sizeof(double), 1, file);
    }
  }

  fclose(file);
}


/// Compares an image (depth) with its ground truth then reports the average
/// error between them.  It also reports the percentage of pixels
/// that by more than a threshold.
void score_vs_gt(const vil_image_view<double> &depth,
                 const vil_image_view<double> &gt,
                 double threshold)
{
  vil_image_view<double> error_img;
  vil_math_image_abs_difference(depth, gt, error_img);

  vil_image_view<bool> bad_img;
  unsigned int numbad = 0;
  unsigned int numpixels = depth.ni() * depth.nj();
  vil_threshold_above(error_img, bad_img, threshold);
  vil_math_sum(numbad, bad_img, 0);

  double mean_error = 0.0;
  vil_math_mean(mean_error, error_img, 0);

  std::cout << "Percentage of bad pixels: " << (double)numbad / (double)numpixels << std::endl;
  std::cout << "Average Error: " << mean_error << std::endl;

  vil_image_view<vxl_byte> bad_img_byte;
  vil_convert_stretch_range(bad_img, bad_img_byte);
  vil_save(bad_img_byte, "error.png");
}


#ifdef HAVE_VTK

/// Writes a vtp (vtk) from a depth image
/// \param filename .vtp file to write to
/// \param depth the depth image
/// \param ref reference image that the depth was computed by, used for texturing
/// \param cam perspective camera associated with the depth map
/// \param world_space the world space class used in the depth estimation
void save_depth_to_vtp(const char *filename,
                       const vil_image_view<double> &depth,
                       const vil_image_view<double> &ref,
                       const vpgl_perspective_camera<double> &cam,
                       world_space *ws)
{
  vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
  vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
  //vtkSmartPointer<vtkFloatArray> normals = vtkSmartPointer<vtkFloatArray>::New();

  //normals->SetName("normals");
  //normals->SetNumberOfComponents(3);

  colors->SetName("colors");
  colors->SetNumberOfComponents(ref.nplanes());

  vtkIdType tri[3];
  vtkIdType width = depth.ni(), height = depth.nj();

  //vil_image_view<double> normal_map(depth.ni(), depth.nj(), 3), bp(depth.ni(),depth.nj(),3);
  //compute_normals_eig(depth, bp, normal_map, ws, 1, 1.0);

  for (vtkIdType j = 0; j < height; j++)
  {
    for (vtkIdType i = 0; i < width; i++)
    {
      vnl_double_3 pt3d = ws->point_at_depth((unsigned int)i, (unsigned int)j, depth(i,j));
      pts->InsertNextPoint(pt3d.data_block());

      double u, v;
      cam.project(pt3d[0], pt3d[1], pt3d[2], u, v);
      if (ref.nplanes() == 1)
      {
        colors->InsertNextValue((unsigned char)vil_bilin_interp_safe(ref, u, v, 0));
      }
      else
      {
        colors->InsertNextTuple3((unsigned char)vil_bilin_interp_safe(ref, u, v, 0),
                                 (unsigned char)vil_bilin_interp_safe(ref, u, v, 1),
                                 (unsigned char)vil_bilin_interp_safe(ref, u, v, 2));
      }
      //normals->InsertNextTuple3(normal_map(i,j,0), normal_map(i,j,1), normal_map(i,j,2));

      if (i != 0 && j != 0)
      {
        tri[0] = j*width + i;
        tri[1] = (j-1)*width + (i-1);
        tri[2] = (j-1)*width + i;
        cells->InsertNextCell(3, tri);

        tri[1] = j*width + (i-1);
        tri[2] = (j-1)*width + (i-1);
        cells->InsertNextCell(3, tri);
      }
    }
  }

  vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
  polydata->SetPoints(pts);
  polydata->SetPolys(cells);
  //polydata->GetPointData()->SetNormals(normals);
  polydata->GetPointData()->AddArray(colors);

  vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  writer->SetInputData(polydata);
  writer->SetFileName(filename);
  writer->Update();
}

//*****************************************************************************

void write_points_to_vtp(std::vector<vnl_double_3> &points, const char *filename)
{
  vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> verts = vtkSmartPointer<vtkCellArray>::New();

  for (vtkIdType i = 0; i < points.size(); i++)
  {
    pts->InsertNextPoint(points[i].data_block());
    verts->InsertNextCell(1, &i);
  }

  vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
  polydata->SetPoints(pts);
  polydata->SetVerts(verts);

  vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  writer->SetInputData(polydata);
  writer->SetFileName(filename);
  writer->Update();
}

#endif

} // end namespace super3d

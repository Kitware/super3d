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


// VXL includes
#include <vil/vil_image_view.h>
#include <vil/vil_convert.h>
#include <vil/vil_save.h>
#include <vul/vul_arg.h>
#include <vul/vul_file.h>
#include <vul/vul_sprintf.h>
#include <vnl/vnl_inverse.h>
#include <super3d/imesh/imesh_mesh.h>
#include <super3d/imesh/imesh_fileio.h>
#include <super3d/imesh/algo/imesh_project.h>
#include <super3d/imesh/algo/imesh_render.h>
#include <super3d/depth/file_io.h>
#include <vpgl/vpgl_perspective_camera.h>

#include <super3d/depth/multiscale.h>
#include <super3d/depth/depth_map.h>


// stl includes
#include <iostream>
#include <map>



/// Load a camera file with sequence of cameras in ASCII format: i K R t
/// where i is the frame number, K is the calibration matrix,
/// R is the rotation matrix, and t is the translation vector
/// \returns A mapping from frame number to camera
std::map<int, vpgl_perspective_camera<double> >
load_cams(const std::string& filename)
{
  std::map<int, vpgl_perspective_camera<double> >  cameras;
  std::ifstream ifs(filename.c_str());
  int frame;
  while (ifs >> frame)
  {
    vpgl_perspective_camera<double> cam;
    ifs >> cam;
    cameras[frame] = cam;
  }

  ifs.close();
  return cameras;
}


/// Save a camera mapping to the same ASCII file format as above
void
save_cams(const std::string& filename,
          const std::map<std::string, vpgl_perspective_camera<double> >& cameras)
{
  typedef std::map<std::string, vpgl_perspective_camera<double> >::const_iterator camera_iterator;
  std::ofstream ofs(filename.c_str());
  for (camera_iterator citr=cameras.begin(); citr!=cameras.end(); ++citr)
  {
    ofs << citr->first << "\n";
    ofs << citr->second << "\n";
  }
  ofs.close();
}


int main(int argc, char* argv[])
{
  vul_arg<std::string> input_mesh( 0, "input mesh file (OBJ)", "" );
  vul_arg<std::string> file_list( 0, "list of image files to match to camera files", "" );
  vul_arg<std::string> camera_dir( "-d", "directory containing input camera files", "" );
  vul_arg<std::string> image_size( "-s", "Image size (WIDTHxHEIGHT)","1920x1080" );
  vul_arg<double>      resolution_scale( "-r", "resolution scale ", 1.0 );
  vul_arg<bool>        byte_images( "-b", "Make byte images where [max, min] = [1,255] and 0 = infinity ", false);

  vul_arg<std::string> output_dir( "-o", "Output directory", "depthmaps" );
  vul_arg<std::string> height_dir( "-h", "Output height map directory", "" );
  vul_arg<std::string> label_dir( "-l", "Output label map directory", "" );
  vul_arg<std::string> output_camera_file( "--output-camera-file", "Save the scaled cameras to this file", "" );

  vul_arg_parse( argc, argv );

  std::stringstream ss(image_size());
  unsigned width, height;
  char sep;
  ss >> width >> sep >> height;
  if ( !ss.eof() || ss.fail() || (sep != 'x' && sep != 'X'))
  {
    std::cerr << "Error: Unable to parse image size string \""
              << image_size() << "\"\n"
              << "       Format should be WIDTHxHEIGHT (e.g. 1024x768)"
              << std::endl;
    return -1;
  }


  imesh_mesh mesh;
  if (!imesh_read(input_mesh(), mesh))
  {
    std::cout << "unable to load mesh file: "<<input_mesh()<<std::endl;
  }
  std::cout << "read mesh: "<<mesh.num_verts()
            << " verts and "<<mesh.num_faces()<<" faces"<<std::endl;

  std::ifstream ifs(file_list());
  std::map<std::string, vpgl_perspective_camera<double> > cameras;
  while (ifs)
  {
    std::string filename;
    ifs >> filename;
    filename = vul_file::strip_extension(vul_file::basename(filename));
    std::cout << "Looking for "<< filename;
    std::string camera_file = camera_dir() + "/" + filename + ".krtd";
    if(vul_file::exists(camera_file))
    {
      std::cout << " -- Loading " << camera_file;
      vpgl_perspective_camera<double> camera = super3d::load_cam(camera_file);
      cameras[filename] = camera;
    }
    std::cout << std::endl;
  }


  vil_image_view<double> depth_map(width, height);
  vil_image_view<double> height_map(width, height);
  vil_image_view<vxl_uint_16> label_map(width, height);

  if (!vul_file::is_directory(output_dir()))
  {
    vul_file::make_directory_path(output_dir());
  }
  if (height_dir() != "" && !vul_file::is_directory(height_dir()))
  {
    vul_file::make_directory_path(height_dir());
  }
  if (label_dir() != "" && !vul_file::is_directory(label_dir()))
  {
    vul_file::make_directory_path(label_dir());
  }

  vxl_uint_16 max_label = mesh.faces().groups().size();

  typedef std::map<std::string, vpgl_perspective_camera<double> >::const_iterator camera_iterator;
  std::map<std::string, vpgl_perspective_camera<double> > scaled_cameras;
  for (camera_iterator citr=cameras.begin(); citr!=cameras.end(); ++citr)
  {
    vpgl_perspective_camera<double> camera = citr->second;
    if (resolution_scale.set())
    {
      camera = super3d::scale_camera(camera, resolution_scale());
    }
    scaled_cameras[citr->first] = camera;
    std::cout << "Rendering depth map from camera " << citr->first << std::endl;
    if (label_dir() != "")
    {
      imesh_render_mesh_labels(mesh, camera, label_map, depth_map);
    }
    else
    {
      imesh_project_depth(mesh, camera, depth_map);
    }
    std::string name = output_dir() + "/" + citr->first + "-depth.tiff";
    if (byte_images())
    {
      double min_d, max_d;
      super3d::finite_value_range(depth_map, min_d, max_d);
      std::cout << "min depth: "<<min_d<<" max depth: "<<max_d<<std::endl;
      std::cout << "Saving byte depth image to " << name << std::endl;
      vil_image_view<vxl_byte> byte_image;
      vil_convert_stretch_range_limited(depth_map, byte_image, min_d, max_d, 1, 255);
      vil_math_scale_and_offset_values(byte_image, -1, 255);
      vil_save(byte_image, name.c_str());
    }
    else
    {
      std::cout << "Saving double depth image to " << name << std::endl;
      vil_save(depth_map, name.c_str());
    }

    if(height_dir() != "")
    {
      super3d::depth_map_to_height_map(camera, depth_map, height_map);
      std::string height_name = height_dir() + "/" + citr->first + "-height.tiff";
      if (byte_images())
      {
        double min_h, max_h;
        super3d::finite_value_range(height_map, min_h, max_h);
        std::cout << "min height: "<<min_h<<" max height: "<<max_h<<std::endl;
        std::cout << "Saving byte height image to " << height_name << std::endl;
        vil_image_view<vxl_byte> byte_image;
        vil_convert_stretch_range_limited(height_map, byte_image, min_h, max_h, 1, 255);
        vil_save(byte_image, height_name.c_str());
      }
      else
      {
        std::cout << "Saving double height image to " << height_name << std::endl;
        vil_save(height_map, height_name.c_str());
      }
    }

    if (label_dir() != "")
    {
      std::string label_name = label_dir() + "/" + citr->first + "-labels.tiff";
      if (byte_images())
      {
        std::cout << "Saving mapped label image to " << label_name << std::endl;
        vil_image_view<vxl_byte> byte_image;
        vil_convert_stretch_range_limited(label_map, byte_image, vxl_uint_16(0), max_label, 0, 255);
        vil_save(byte_image, label_name.c_str());
      }
      else
      {
        std::cout << "Saving unscaled label image to " << label_name << std::endl;
        vil_save(label_map, label_name.c_str());
      }

      imesh_render_mesh_labels(mesh, camera, label_map, depth_map);
    }
  }
  // write out scaled cameras if requested
  if (output_camera_file.set())
  {
    save_cams(output_camera_file(), scaled_cameras);
  }

  return 0;
}

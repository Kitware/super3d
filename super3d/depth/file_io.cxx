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

#include "file_io.h"
#include "multiscale.h"

#include <fstream>
#include <iostream>

#include <vul/vul_file.h>
#include <vul/vul_sequence_filename_map.h>
#include <vil/vil_load.h>
#include <vil/vil_image_view.h>
#include <vil/vil_convert.h>
#include <vnl/vnl_double_3.h>

#include "multiscale.h"


namespace super3d
{

/// Load a camera file with sequence of cameras in ASCII format: i K R t
/// where i is the frame number, K is the calibration matrix,
/// R is the rotation matrix, and t is the translation vector
/// \returns A vector of perspective cameras
std::vector<vpgl_perspective_camera<double> >
load_cams(const std::string& filename, vul_sequence_filename_map frame_seq)
{
  std::vector<vpgl_perspective_camera<double> > cameras;
  std::fstream ifs(filename.c_str());
  const std::vector< int > & indices = frame_seq.get_real_indices();
  int frame;
  unsigned int index = 0;
  while (ifs >> frame && index < indices.size())
  {
    vpgl_perspective_camera<double> cam;
    ifs >> cam;
    if (indices[index] == frame)
    {
      vpgl_calibration_matrix<double> K = cam.get_calibration();
      K.set_focal_length(K.focal_length() * K.x_scale());
      K.set_y_scale(K.y_scale() / K.x_scale());
      K.set_x_scale(1.0);
      cam.set_calibration(K);
      cameras.push_back(cam);
      index++;
    }
  }

  ifs.close();

  std::cout << "Loaded " << cameras.size() << " cameras.\n";
  return cameras;
}


/// Find all frames matching the format string and extract the frame number.
/// \returns A vector of image views
std::vector<vil_image_view<double> >
load_frames(vul_sequence_filename_map frame_seq, std::vector<std::string> &filenames, bool color)
{
  std::vector<vil_image_view<double> > frames;
  const std::vector< int > & indices = frame_seq.get_real_indices();
  for (unsigned i=0; i<indices.size(); ++i)
  {
    std::string file = frame_seq.image_name(i);
    std::cout << frame_seq.get_real_index(i) << " ";
    std::cout << "Loading: " << file << " ";
    if (vul_file_exists(file))
    {
      vil_image_resource_sptr img_rsc = vil_load_image_resource(file.c_str());
      if (img_rsc != NULL)
      {
        vil_image_view<vxl_byte> img = img_rsc->get_view();
        vil_image_view<double> flt;
        if (img.nplanes() == 3 && !color)
          vil_convert_planes_to_grey(img, flt);
        else
          vil_convert_cast(img, flt);
        frames.push_back(flt);

        filenames.push_back(file);
      }
      else
        std::cout << " NOT FOUND.";
    }
    std::cout << "\n";
  }
  return frames;
}


/// Load an exposure file with parameters for linear exposure compensation
/// \returns A vector of (scale, offest) pair
std::vector<std::pair<double,double> >
load_exposure(const std::string& filename, vul_sequence_filename_map frame_seq)
{
  typedef std::pair<double,double> Dpair;
  std::vector<Dpair> exposure;
  std::fstream ifs(filename.c_str());
  const std::vector< int > & indices = frame_seq.get_real_indices();
  int frame;
  unsigned int index = 0;
  while (ifs >> frame && index < indices.size())
  {
    double scale, offset;
    ifs >> scale >> offset;
    if (indices[index] == frame)
    {
      exposure.push_back(Dpair(scale,offset));
      index++;
    }
  }

  ifs.close();
  return exposure;
}


/// Load an exposure file with parameters for linear exposure compensation
/// \returns A vector of (scale, offest) pair
std::vector<std::pair<double,double> >
load_exposure(const std::string& filename, const std::vector<int> &framelist)
{
  typedef std::pair<double,double> Dpair;
  std::vector<Dpair> exposure;
  std::fstream ifs(filename.c_str());
  int frame;
  unsigned int index = 0;
  while (ifs >> frame && index < framelist.size())
  {
    double scale, offset;
    ifs >> scale >> offset;
    if (framelist[index] == frame)
    {
      exposure.push_back(Dpair(scale,offset));
      index++;
    }
  }

  ifs.close();
  return exposure;
}

//Load cameras from a file per camera
//Assume cameras are in directory named %04d.krtd
void
load_krtd_cams(const std::string& directory,
               std::vector<int> &framelist,
               std::vector<vpgl_perspective_camera<double> > &cameras)
{
  std::cout << "Loading Cameras: ";
  char buf[32];
  for (unsigned int i = 0; i < framelist.size(); i++)
  {
    sprintf(buf, "%04d.krtd", framelist[i]);
    std::cout << (directory+buf) << "\n";
    std::fstream ifs((directory + buf).c_str());

    vpgl_perspective_camera<double> cam;
    vnl_double_3x3 K, R;
    vgl_vector_3d<double> t;

    ifs >> K >> R >> t;
    ifs.close();

    cam.set_calibration(vpgl_calibration_matrix<double>(K));
    cam.set_rotation(vgl_rotation_3d<double>(R));
    cam.set_translation(t);

    vpgl_calibration_matrix<double> cal = cam.get_calibration();
    cal.set_focal_length(cal.focal_length() * cal.x_scale());
    cal.set_y_scale(cal.y_scale() / cal.x_scale());
    cal.set_x_scale(1.0);
    cam.set_calibration(cal);

    cameras.push_back(cam);
  }

  std::cout << "\n";
}

//Load cameras from a single camera file with a framelist
void load_cams(const char *camerafile,
               std::vector<int> &framelist,
               std::vector<vpgl_perspective_camera<double> > &cameras)
{
  std::ifstream camstream(camerafile);
  int frame;

  std::cout << "Looking for " << framelist.size() << " frames.\n";

  unsigned int cur = 0;
  while (camstream.good() && cur < framelist.size())
  {
    camstream >> frame;
    vpgl_perspective_camera<double> cam;
    camstream >> cam;

    //assumes frames are sorted
    if (frame != framelist[cur])
      continue;
    cur++;

    vpgl_calibration_matrix<double> K = cam.get_calibration();
    K.set_focal_length(K.focal_length() * K.x_scale());
    K.set_y_scale(K.y_scale() / K.x_scale());
    K.set_x_scale(1.0);
    cam.set_calibration(K);

    cameras.push_back(cam);
  }

  camstream.close();
  std::cout << "Loaded " << cameras.size() << " cameras.\n";
}


/// Loads images from a file list of frames paths
/// \param framefile file that lists frame number and frame paths
/// \param directory the working directory from which the framefile appends its paths to
/// \param filenames vector of frame files that were read
/// \param framelist vector of indices of the read frames
/// \param frames images that were read and converted to greyscale
void load_from_frame_file(const char *framefile,
                          const std::string &directory,
                          std::vector<std::string> &filenames,
                          std::vector<int> &framelist,
                          std::vector<vil_image_view<double> > &frames,
                          bool color,
                          bool rgb12)
{
  std::ifstream framestream(framefile);
  int frame;
  std::string imagename;

  while (framestream >> frame >> imagename && framestream.good())
  {
    framelist.push_back(frame);
    filenames.push_back(imagename);

    std::cout << "Reading frame: " << (directory + imagename) << "\n";
    vil_image_resource_sptr img_rsc = vil_load_image_resource((directory + imagename).c_str());
    if (img_rsc != NULL)
    {
      vil_image_view<double> flt;
      if (rgb12)
      {
        vil_image_view<unsigned short> img = img_rsc->get_view();
        vil_convert_cast(img, flt);
        vil_math_scale_values(flt, 255.0 / 4095.0);

        if (img.nplanes() == 3 && !color)
        {
          vil_image_view<double> grey;
          vil_convert_planes_to_grey(flt, grey);
          flt = grey;
        }
      }
      else
      {
        vil_image_view<vxl_byte> img = img_rsc->get_view();
        if (img.nplanes() == 3 && !color)
          vil_convert_planes_to_grey(img, flt);
        else
          vil_convert_cast(img, flt);
      }
      frames.push_back(flt);
    }
    else
      std::cout << "\n" << (directory + imagename) << " NOT FOUND.\n";
  }

  std::cout << "\n";
  framestream.close();
}

/// Load camera from a file per camera
vpgl_perspective_camera<double>
load_cam(const std::string& filename)
{
  std::fstream ifs((filename).c_str());

  vpgl_perspective_camera<double> cam;

  ifs >> cam;
  ifs.close();

  vpgl_calibration_matrix<double> cal = cam.get_calibration();
  cal.set_focal_length(cal.focal_length() * cal.x_scale());
  cal.set_y_scale(cal.y_scale() / cal.x_scale());
  cal.set_x_scale(1.0);
  cam.set_calibration(cal);

  return cam;
}

/// read a flow file into 2-band image
bool read_flow_file(vil_image_view<double> &flowimg, const char* filename)
{
  std::cout << "Reading " << filename << "\n";
  FILE *stream = fopen(filename, "rb");
  if (stream == 0)
  {
    std::cerr << "ReadFlowFile: could not open file\n";
    return false;
  }

  int width, height;
  float tag;

  if ((int)fread(&tag,    sizeof(float), 1, stream) != 1 ||
  (int)fread(&width,  sizeof(int),   1, stream) != 1 ||
  (int)fread(&height, sizeof(int),   1, stream) != 1)
  {
    std::cerr << "read_flow_file: problem reading file" << std::endl;
    return false;
  }

  int nBands = 2;
  flowimg.set_size(width, height, nBands);

  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++) {
      for (int p = 0; p < nBands; p++) {
        float val;
        if( fread(&val, sizeof(float), 1, stream) != 1)
        {
          std::cerr << "read_flow_file: error reading file" << std::endl;
          return false;
        }
        flowimg(x,y,p) = val;
      }
    }
  }

  fclose(stream);
  return true;
}

} // end namespace super3d

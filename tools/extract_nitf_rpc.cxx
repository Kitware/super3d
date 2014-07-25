/*
 * Copyright 2014 Kitware, Inc.
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

#include <iostream>
#include <sstream>
#include <cstdlib>

#include <vgl/vgl_box_3d.h>
#include <vil/vil_convert.h>
#include <vil/vil_load.h>
#include <vil/vil_save.h>
#include <vpgl/file_formats/vpgl_nitf_rational_camera.h>
#include <vpgl/algo/vpgl_camera_convert.h>

#include <vil_plugins/vil_plugin_loader.h>


int main(int argc, char *argv[])
{

  double lat = std::atof(argv[1]);
  double lng = std::atof(argv[2]);
  double elv = std::atof(argv[3]);

  vgl_box_3d<double> roi(lat-0.01, lng-0.01, elv-50,
                         lat+0.01, lng+0.01, elv+50);

  //load rational camera from image files
  std::vector<vpgl_nitf_rational_camera> nitf_cams;
  std::vector<std::string> nitf_paths;
  for(int i=4; i<argc; ++i)
  {
    vpgl_nitf_rational_camera *nitf_cam = new vpgl_nitf_rational_camera(argv[i]);

    if (!nitf_cam)
    {
      std::cerr << "Error: "<<argv[0] <<" Failed to load NITF camera" <<std::endl;
      return -1;
    }

    nitf_cams.push_back(*nitf_cam);
    nitf_paths.push_back(argv[i]);
    delete nitf_cam;
  }

  /// load the vidtk-vil plugins to get the GDAL reader
  /// Note: vpgl_nitf_rational_camera (above) requires the
  /// nitf2 reader built into VXL to extract camera parameters,
  /// but this fails to load J2K encoded imagery.  The GDAL
  /// loader plugin in vidtk can handle the imagery, but doesn't
  /// integrate with NITF rational camera extraction code.
  /// Therefore, we read the camera parameters first with nitf2 (above),
  /// then load the GDAL plugin (here), and then read the pixels
  /// with GDAL in a second pass (below).
  vidtk::load_vil_plugins();

  vgl_h_matrix_3d<double> trans;
  vpgl_lvcs lvcs_converter( roi.min_y(), roi.min_x(),
                            (roi.min_z() + roi.max_z())/2.0,
                            vpgl_lvcs::wgs84, vpgl_lvcs::DEG );
  for(unsigned int i=0; i<nitf_cams.size(); ++i)
  {
    const vpgl_nitf_rational_camera& nitf_cam = nitf_cams[i];
    double u,v;

    std::cout << "loading pixels from "<< nitf_paths[i] <<std::endl;
    vil_image_resource_sptr im = vil_load_image_resource(nitf_paths[i].c_str());
    if(!im)
    {
      vcl_cerr << "Unable to load image" << vcl_endl;
      return EXIT_FAILURE;
    }
    std::cout << "loaded as format: "<< im->file_format()
              << " size "<< im->ni() <<" x "<<im->nj() <<std::endl;

    unsigned int i1 = 0, j1 = 0;
    unsigned int i0 = im->ni(), j0 = im->nj();
    for( unsigned int x=0; x<2; ++x)
    {
      for( unsigned int y=0; y<2; ++y)
      {
        for( unsigned int z=0; z<2; ++z)
        {
          nitf_cam.project(x?roi.min_x():roi.max_x(),
                           y?roi.min_y():roi.max_y(),
                          z?roi.min_z():roi.max_z(), u, v);
          std::cout << (x?roi.min_x():roi.max_x()) << ", "
                    << (y?roi.min_y():roi.max_y()) << ", "
                    << (z?roi.min_z():roi.max_z()) << std::endl;
          std::cout << "  proj: "<<u<<", "<<v<<std::endl;
          if (u < i0)
          {
            i0 = (u<0) ? 0 : static_cast<unsigned>(u);
          }
          if (u > i1)
          {
            i1 = (u>im->ni()) ? im->ni() : static_cast<unsigned>(u);
          }
          if (v < j0)
          {
            j0 = (v<0) ? 0 : static_cast<unsigned>(v);
          }
          if (v > j1)
          {
            j1 = (v>im->nj()) ? im->nj() : static_cast<unsigned>(v);
          }
        }
      }
    }

    unsigned int ni = i1-i0;
    unsigned int nj = j1-j0;
    std::cout << "ROI: "<<ni<<"x"<<nj<<"+"<<i0<<"+"<<j0<<std::endl;

    if( ni == 0 || nj == 0 )
    {
      std::cerr << "selected ROI has zero intersection with image" << std::endl;
      continue;
    }

    vil_image_view_base_sptr view = im->get_copy_view(i0, ni, j0, nj);
    if( !view )
    {
      std::cerr << "unable to extract pixels in ROI" <<std::endl;
      continue;
    }
    std::cout << "type "<< view->pixel_format()
              << " size " << view->ni() <<", "<< view->nj() <<std::endl;
    std::stringstream ss;
    ss << "crop_"<<i<<".tiff";
    vil_image_view<vxl_byte> byte_img = vil_convert_stretch_range(vxl_byte(), view);
    std::cout << "saving "<< ss.str() << " size "<< byte_img.ni() <<", "<<byte_img.nj() <<std::endl;
    vil_save(byte_img, ss.str().c_str());


    nitf_cam.project(lat, lng, elv, u, v);
    std::cout << "RPC test ("<<lat<<", "<<lng<<", "<<elv<<") --> ("<<u<<", "<<v<<")"<<std::endl;


    vpgl_perspective_camera<double> pcam;
    vpgl_perspective_camera_convert::convert_local(nitf_cam, roi, pcam, trans);
    //vgl_homg_point_3d<double> lla(lat, lng, elv);
    double lx, ly, lz;
    lvcs_converter.global_to_local(lat, lng, elv, vpgl_lvcs::wgs84, lx, ly, lz);
    vgl_homg_point_3d<double> lla(lx, ly, lz);
    vgl_point_3d<double> loc(trans*lla);
    std::cout << "transformed point "<<loc<<std::endl;
    pcam.project(loc.x(), loc.y(), loc.z(), u, v);
    std::cout << "Proj test ("<<lat<<", "<<lng<<", "<<elv<<") --> ("<<u<<", "<<v<<")"<<std::endl;

  }

  std::ofstream ofs("geo_transform.txt");
  ofs << trans << '\n' << lvcs_converter;

  return 0;
}

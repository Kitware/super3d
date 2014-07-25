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

#include <boost/filesystem.hpp>

#include <vgl/vgl_box_3d.h>
#include <vil/vil_convert.h>
#include <vil/vil_load.h>
#include <vil/vil_save.h>
#include <vpgl/file_formats/vpgl_nitf_rational_camera.h>
#include <vpgl/algo/vpgl_camera_convert.h>

#include <vil_plugins/vil_plugin_loader.h>

#include <super3d/depth/multiscale.h>

int main(int argc, char *argv[])
{

  vgl_box_3d<double> roi(std::atof(argv[1]),
                         std::atof(argv[2]),
                         std::atof(argv[3]),
                         std::atof(argv[4]),
                         std::atof(argv[5]),
                         std::atof(argv[6]));

  double lat = roi.centroid_x();
  double lng = roi.centroid_y();
  double elv = roi.centroid_z();

  //load rational camera from image files
  std::vector<vpgl_nitf_rational_camera> nitf_cams;
  std::vector<boost::filesystem::path> nitf_paths;
  for(int i=7; i<argc; ++i)
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

    // Compute the image ROI that contains the world ROI.
    // Project each of the 8 corner points into the image and
    // find the bounding box in the image that contains them
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

    // Decode just the pixels in our ROI
    vil_image_view_base_sptr view = im->get_copy_view(i0, ni, j0, nj);
    if( !view )
    {
      std::cerr << "unable to extract pixels in ROI" <<std::endl;
      continue;
    }
    std::cout << "type "<< view->pixel_format()
              << " size " << view->ni() <<", "<< view->nj() <<std::endl;

    // stretch the range to get a byte representation of the pixels
    vil_image_view<vxl_byte> byte_img = vil_convert_stretch_range(vxl_byte(), view);

    // write the cropped byte image out to a TIFF file
    std::string basename = nitf_paths[i].stem().string();
    std::stringstream ss;
    ss << basename << "_crop_"<<ni<<"x"<<nj<<"+"<<i0<<"+"<<j0;
    std::cout << "saving "<< ss.str() << ".tiff" <<std::endl;
    vil_save(byte_img, (ss.str() + ".tiff").c_str());

    // Estimate a perspective (KRT) camera to approximate the RPC within the ROI.
    // Note 'trans' is estimated here but only depends on 'roi' and should be
    // the same for all images.
    vpgl_perspective_camera<double> pcam;
    vpgl_perspective_camera_convert::convert_local(nitf_cam, roi, pcam, trans);

    // Do a test projection to compare the KRT model with the RPC model
    double lx, ly, lz;
    lvcs_converter.global_to_local(lat, lng, elv, vpgl_lvcs::wgs84, lx, ly, lz);
    vgl_homg_point_3d<double> lla(lx, ly, lz);
    vgl_point_3d<double> loc(trans*lla);
    pcam.project(loc.x(), loc.y(), loc.z(), u, v);
    std::cout << "KRT test ("<<lat<<", "<<lng<<", "<<elv<<") --> ("<<u<<", "<<v<<")"<<std::endl;

    nitf_cam.project(lat, lng, elv, u, v);
    std::cout << "RPC test ("<<lat<<", "<<lng<<", "<<elv<<") --> ("<<u<<", "<<v<<")"<<std::endl;


    // crop camera to correspond to the image crop
    pcam = super3d::crop_camera(pcam, i0, j0);

    // write out the KRTD file (pcam is KRT, D is just 0)
    std::ofstream ofs((ss.str() + ".krtd").c_str());
    ofs << pcam << "\n0\n";
  }

  std::ofstream ofs("geo_transform.txt");
  ofs << trans << '\n' << lvcs_converter;

  return 0;
}

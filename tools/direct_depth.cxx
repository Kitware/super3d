/*
 * Copyright 2012 Kitware, Inc.
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

#include "tv_refine.h"

#include <super3d/depth/file_io.h>
#include <super3d/depth/depth_map.h>
#include <super3d/depth/multiscale.h>

#include <super3d/image/gaussian_pyramid_builder.h>

// VXL includes
#include <vil/vil_image_view.h>
#include <vil/vil_convert.h>
#include <vil/vil_save.h>
#include <vil/vil_bilin_interp.h>
#include <vil/vil_crop.h>
#include <vul/vul_arg.h>
#include <vnl/vnl_double_3.h>
#include <vnl/vnl_inverse.h>
#include <super3d/imesh/imesh_mesh.h>
#include <super3d/imesh/imesh_fileio.h>
#include <super3d/imesh/algo/imesh_project.h>
#include <vpgl/vpgl_perspective_camera.h>

#include <vil/vil_decimate.h>
#include <vil/vil_resample_bilin.h>

// vcl includes
#include <vcl_iostream.h>
#include <vcl_string.h>


void compute_plane_parallax(const vpgl_proj_camera<double>& src_camera,
                            const vpgl_proj_camera<double>& dst_camera,
                            vnl_matrix_fixed<double,3,3>& H,
                            vnl_vector_fixed<double,3>& e);
void compute_gaussian_pyramid(vil_image_view<double> &frame,
                              vcl_pair<double, double> *exposure,
                              const vidtk::gaussian_pyramid_builder &gpb,
                              vcl_vector<vil_image_view<double> > &pyramid);


int main(int argc, char* argv[])
{
  vul_arg<vcl_string> input_mesh( 0, "input mesh file (OBJ)", "" );
  vul_arg<vcl_string> camera_file( 0, "input file containing camera sequence", "" );
  vul_arg<vcl_string> frame_fmt( 0, "frame file format string (e.g. \"tex_%03d.png\")", "" );
  vul_arg<vcl_string> output_mesh( 0, "output mesh file (OBJ)", "");

  vul_arg<unsigned>    reference_frame( "-f", "Reference frame for refining depth", 0 );
  vul_arg<vcl_string> exposure_file( "-e", "Use exposure compensation in the provided file", "");
  vul_arg<unsigned>    output_decimate( "-d", "Decimation factor in output mesh", 1);
  vul_arg<double> camera_scale("-s", "Scaling of camera (e.g. super resolution images)", 1.0);
  vul_arg<vcl_string> crop_window("-cw", "Crop window i0, ni, j0, nj (e.g., \"0 100 0 200\"", "");

  vul_arg_parse( argc, argv );

  vul_sequence_filename_map frame_seq(frame_fmt());

  //Read Mesh
  imesh_mesh mesh;
  if (!imesh_read(input_mesh(), mesh))
    vcl_cout << "unable to load mesh file: "<<input_mesh()<<vcl_endl;

  vcl_cout << "read mesh: "<<mesh.num_verts()
            <<" verts and "<<mesh.num_faces()<<" faces"<<vcl_endl;

  //Read Cameras
  vcl_vector<vpgl_perspective_camera<double> >  cameras = super3d::load_cams(camera_file(), frame_seq);
  if (camera_scale() != 1.0)
    for (unsigned int i = 0; i < cameras.size(); i++)
      cameras[i] = super3d::scale_camera(cameras[i], camera_scale());

  //Read Images
  vcl_vector<vcl_string> filenames;
  vcl_vector<vil_image_view<double> > frames = super3d::load_frames(frame_seq, filenames);
  if (frames.empty())
  {
    vcl_cerr << "No frames found"<<vcl_endl;
    return -1;
  }
  else if (frames.size() < 2)
  {
    vcl_cerr << "At least 2 frames are required"<<vcl_endl;
    return -1;
  }

  unsigned int ref_frame = reference_frame();

  //Compute the window cropping, scale the cropping by the specified scale so that we do not
  //need to recompute cropping input for super resolution.
  int i0, ni, j0, nj;
  if (!crop_window().empty())
  {
    vcl_istringstream cwstream(crop_window());
    cwstream >> i0 >> ni >> j0 >> nj;
    i0 = (int)(i0*camera_scale());
    j0 = (int)(j0*camera_scale());
    ni = (int)(ni*camera_scale());
    nj = (int)(nj*camera_scale());
    vcl_cout << "Crop window: " << i0 << " " << ni << " " << j0 << " " << nj << "\n";
  }
  else
  {
    i0 = j0 = 0;
    ni = frames[ref_frame].ni();
    nj = frames[ref_frame].nj();
  }

  vcl_vector<vcl_pair<double, double> > exposures;
  if (exposure_file.set())
    exposures = super3d::load_exposure(exposure_file(), frame_seq);

  vcl_cout << "Making image pyramids"<<vcl_endl;
  unsigned levels = 4;
  vidtk::gaussian_pyramid_builder gpb(levels+1, 2, 1.0);

  vcl_vector<vil_image_view<double> > pyr_ref;
  compute_gaussian_pyramid(frames[ref_frame], exposure_file.set() ? &exposures[ref_frame] : NULL, gpb, pyr_ref);

  vcl_vector<vcl_vector<vil_image_view<double> > > pyrs(frames.size()-1);
  vcl_vector<vnl_matrix_fixed<double,3,3> > Hs;
  vcl_vector<vnl_vector_fixed<double,3> > es;
  for (unsigned int i = 0, index = 0; i < frames.size(); i++)
  {
    if (i == ref_frame)
      continue;

    vnl_matrix_fixed<double, 3, 3> H;
    vnl_double_3 e;
    compute_plane_parallax(cameras[ref_frame], cameras[i], H, e);
    Hs.push_back(H);
    es.push_back(e);

    compute_gaussian_pyramid(frames[i], exposure_file.set() ? &exposures[i] : NULL, gpb, pyrs[index++]);
  }

  vcl_cout << "Initializing depth map"<<vcl_endl;
  vcl_vector<vil_image_view<double> > pyr_depth(levels+1);
  pyr_depth[levels].set_size(pyr_ref[levels].ni(), pyr_ref[levels].nj());
  imesh_project_depth(mesh,
                      super3d::scale_camera(cameras[ref_frame], 1.0/(1<<levels)),
                      pyr_depth[levels]);
  super3d::fill_missing_depths(pyr_depth[levels],4,20);

  vil_image_view<vxl_byte> depth;
  vil_convert_stretch_range_limited(pyr_depth[levels],depth,0.94,0.96);
  vil_save(depth, "depth_image.png");

  vcl_cout << "Refining depth"<<vcl_endl;
  double theta = 1e-4;
  double lambda = .01;
  for (int s = levels; s > 0; --s)
  {
    double scale = 1.0/(1<<s);
    vcl_vector<vil_image_view<double> > scaled_I1;
    vcl_vector<vnl_matrix_fixed<double,3,3> > scaled_H;
    vcl_vector<vnl_vector_fixed<double,3> > scaled_e;
    for (unsigned i = 0; i < pyrs.size(); ++i)
    {
      scaled_I1.push_back(pyrs[i][s]);
      scaled_H.push_back(super3d::scale_homography(Hs[i],scale));
      scaled_e.push_back(super3d::scale_point(es[i],scale));
    }
    refine_depths(pyr_ref[s],scaled_I1,
                  scaled_H, scaled_e,
                  pyr_depth[s], 10, 10, theta, lambda/*/(scale*scale)*/);

    vil_resample_bilin(pyr_depth[s], pyr_depth[s-1], pyr_ref[s-1].ni(), pyr_ref[s-1].nj());
  }

  vcl_vector<vil_image_view<double> > I1;
  for (unsigned i=0; i< pyrs.size(); ++i)
  {
    I1.push_back(pyrs[i][0]);
  }

  //theta = 5.0;
  //for (unsigned int i = 0; i < 5; i++)
  //{
    refine_depths(pyr_ref[0], I1,
                  Hs, es,
                  pyr_depth[0], 10, 10, theta, lambda);

  //  theta /= 10.0;
  //}
  vpgl_perspective_camera<double> croppedcam = super3d::crop_camera(cameras[ref_frame], i0, j0);
  vil_image_view<double> cropped_depth = vil_crop(pyr_depth[0], i0, ni, j0, nj);

  vcl_cout << "writing mesh"<<vcl_endl;
  imesh_mesh nm = super3d::depth_map_to_mesh(super3d::scale_camera(croppedcam, 1.0/output_decimate()),
                                             vil_decimate(cropped_depth,output_decimate()));
  imesh_write_obj(output_mesh(),nm);

  return 0;
}

void compute_gaussian_pyramid(vil_image_view<double> &frame,
                              vcl_pair<double, double> *exposure,
                              const vidtk::gaussian_pyramid_builder &gpb,
                              vcl_vector<vil_image_view<double> > &pyramid)
{
  // scale and offset used to map the range [0,255] to [-1,1]
  const static double norm_scale = 1.0/255.0;
  const static double norm_offset = 0.0;

  double scale = norm_scale;
  double offset = norm_offset;
  if (exposure)
  {
    offset += scale*exposure->second;
    scale *= exposure->first;
  }

  vil_math_scale_and_offset_values(frame, scale, offset);

  double min, max;
  vil_math_value_range(frame, min, max);
  vcl_cout << "min: " << min << " max: " << max << "\n";
  gpb.build_pyramid<double>(frame, pyramid);
}

/// Compute the plane plus parallax model for two cameras.
/// Produces homography \a H and epipole \a e such that a point p = w*(u,v,1),
/// in the image of \a src_camera maps to H*p + e in the image of \a dst_camera.
/// Here w is the projective depth relative to \a src_camera and \a H is the
/// homography induced by the plane at inifinity.
void compute_plane_parallax(const vpgl_proj_camera<double>& src_camera,
                            const vpgl_proj_camera<double>& dst_camera,
                            vnl_matrix_fixed<double,3,3>& H,
                            vnl_vector_fixed<double,3>& e)
{
  // if src_camera matrix is P1 = [M1 | t1] and dst_camera is P2 = [M2 | t2]
  // the H = M2*(M1^-1) and e = -H*t1 + t2
  const vnl_matrix_fixed<double,3,4>& P1 = src_camera.get_matrix();
  const vnl_matrix_fixed<double,3,4>& P2 = dst_camera.get_matrix();
  H = P2.extract(3,3);
  H *= vnl_inverse(P1.extract(3,3));

  e = P2.get_column(3);
  e -= H*P1.get_column(3);
}

/*
 * Copyright 2011 Kitware, Inc.
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

#include <vcl_iostream.h>
#include <vcl_vector.h>
#include <vcl_sstream.h>

#include <vil/vil_convert.h>
#include <vil/vil_save.h>
#include <vil/vil_bilin_interp.h>
#include <vil/algo/vil_median.h>
#include <vil/vil_decimate.h>
#include <vil/vil_resample_bilin.h>
#include <vil/vil_crop.h>
#include <vul/vul_arg.h>

#include <imesh/imesh_mesh.h>
#include <imesh/imesh_fileio.h>
#include <imesh/algo/imesh_project.h>
#include <imesh/imesh_operations.h>

#include <vnl/vnl_double_3.h>
#include <vnl/vnl_double_2.h>
#include <vnl/vnl_matrix_fixed.h>


#include <tracking/dense_optical_flow.h>
#include <vpgl/vpgl_proj_camera.h>

#include "file_io.h"
#include "depth_map.h"
#include "multiscale.h"


void save_flows(const vcl_vector<vil_image_view<double> > &flows, const char *filename);
void load_flows(vcl_vector<vil_image_view<double> > &flows, const char *filename);
vnl_matrix_fixed<double, 2, 3> jacobian_at_x(const vpgl_perspective_camera<double> &cam, const vnl_double_3 &x);
void project_ref_image(const vpgl_perspective_camera<double> &camref,
                const vpgl_perspective_camera<double> &camtarget,
                imesh_mesh &mesh,
                const vil_image_view<double> &ref,
                vil_image_view<double> &reproj);
void flow_to_colormap(const vil_image_view<double> &flow, vil_image_view<vil_rgb<vxl_byte> > &colormap_flow, double max);
void DifferenceImage(const vil_image_view<double> &I0, const vil_image_view<double> &I1, const vil_image_view<double> &flow, vil_image_view<double> &diff);

//C:/Data/meshes/triangluated_klt.obj C:/Data/ground_medians/cameras.txt C:/Data/ground_medians/full_window15_reset50/frame###.png,100:100:900 c:/Data/meshes/cropout.obj -f 0  -e C:/Data/exposure/exposure.txt -cw "767 1320 757 954"
//C:/Data/meshes/bundler_m.ply C:/Data/ground_medians/output.txt C:/Data/ground_medians/full_window15_reset50/frame###.png,100:100:900 c:/Data/meshes/bundlercropout.obj -f 0  -e C:/Data/exposure/exposure.txt -cw "767 1320 757 954"
int main(int argc, char *argv[])
{
  vul_arg<vcl_string> input_mesh( 0, "input mesh file (OBJ)", "" );
  vul_arg<vcl_string> camera_file( 0, "input file containing camera sequence", "" );
  vul_arg<vcl_string> frame_fmt( 0, "frame file format string (e.g. \"tex_%03d.png\")", "" );
  vul_arg<vcl_string> output_mesh( 0, "output mesh file (OBJ)", "");

  vul_arg<unsigned> reference_frame( "-f", "Reference frame for refining depth", 0 );
  vul_arg<vcl_string> exposure_file( "-e", "Use exposure compensation in the provided file", "");
  vul_arg<vcl_string> crop_window("-cw", "Crop window i0, ni, j0, nj (e.g., \"0 100 0 200\"", "");

  //bool loadflow = true;
  bool loadflow = false;

  vul_arg_parse( argc, argv );

  int i0, ni, j0, nj;
  vcl_istringstream cwstream(crop_window());
  cwstream >> i0 >> ni >> j0 >> nj;

  imesh_mesh mesh;
  imesh_read(input_mesh(), mesh);
  imesh_triangulate(mesh);

  vul_sequence_filename_map frame_seq(frame_fmt());
  vcl_vector<vpgl_perspective_camera<double> >  cameras = load_cams(camera_file(), frame_seq);

  //Read Images
  vcl_vector<vcl_string> filenames;
  vcl_vector<vil_image_view<double> > frames = load_frames(frame_seq, filenames);
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

  for (unsigned int i = 0; i < frames.size(); i++)
  {
    cameras[i] = crop_camera(cameras[i], i0, j0);
    frames[i] = vil_crop(frames[i], i0, ni, j0, nj);
  }


  unsigned int ref_frame = reference_frame();

  for (unsigned int iter = 0; iter < 3; iter++)
  {
  vcl_vector<vil_image_view<double> > flows;
  if (loadflow)
  {
    vcl_cout << "Loading flows.\n";
    load_flows(flows, "flowcrop.dat");
  }
  else
  {
    for (unsigned int i = 0; i < frames.size(); i++)
    {
      if (i == ref_frame) {
        flows.push_back(vil_image_view<double>());
        continue;
      }

      vcl_cout << "Computing optical flow of ref and image " << i << ".\n";

      vil_image_view<double> reproj, flow;
      project_ref_image(cameras[ref_frame], cameras[i], mesh, frames[ref_frame], reproj);

      //vil_image_view<vxl_byte> output;
      //vil_convert_cast(frames[i], output);
      //vil_save(output, "ref.png");
      //vil_convert_cast(reproj, output);
      //vil_save(output, "reproj.png");
      char buf[30];
      sprintf(buf, "flow%d.png", i);

      // vidtk::dense_optical_flow<float>(reproj, frames[i], flow, 4, 5, 2, 5, 1, 10);

      vidtk::dense_optical_flow<double> dof;
      dof.set_num_pyramid_levels(4);
      dof.set_num_warps(5);
      dof.set_num_outer_iterations(2);
      dof.set_num_inner_iterations(5);
      dof.set_num_texture_iterations(10);

      dof.set_use_bicubic_interp(true);
      //dof.set_verbose(flag);
      //dof.set_gradient_blend(val);
      //dof.set_structure_removal(val);

      //dof.set_texture_theta(val);
      //dof.set_theta(val);
      //dof.set_lambda(val);

      dof.compute_flow(reproj, frames[i], flow);

      vil_image_view<double> diff;
      DifferenceImage(reproj, frames[i], flow, diff);
      vil_image_view<vxl_byte> diff_byte;
      vil_convert_cast<double, vxl_byte>(diff, diff_byte);
      vil_save(diff_byte, buf);

      //vil_image_view<vil_rgb<vxl_byte> > colormap_flow;
      //flow_to_colormap(flow, colormap_flow, 10);
      //vil_save(colormap_flow, buf);

      flows.push_back(flow);
    }

    save_flows(flows, "flowcrop.dat");
  }

  const vil_image_view<double> &ref = frames[ref_frame];
  const vpgl_perspective_camera<double> &refcam = cameras[ref_frame];

  vcl_cout << "\nEstimating Depth.\n";

    vil_structuring_element se;
  se.set_to_disk(3.0);

  vil_image_view<double> depth(ref.ni(), ref.nj(), 1);
  imesh_project_depth(mesh, cameras[ref_frame], depth);



    unsigned int n = flows.size();
    for (unsigned int i = 0; i < ref.ni(); i++)
    {
      for (unsigned int j = 0; j < ref.nj(); j++)
      {
        if (depth(i,j,0) == vcl_numeric_limits<double>::infinity())
          continue;

        vgl_line_3d_2_points<double> ray = refcam.backproject(vgl_point_2d<double>((double)i, (double)j));
        vgl_vector_3d<double> d = ray.direction();
        normalize(d);
        vgl_point_3d<double> pt3d = refcam.get_camera_center() + depth(i,j,0)*d;
        vnl_double_3 xj =  vnl_double_3(pt3d.x(), pt3d.y(), pt3d.z());
        vnl_double_3 r(d.x(), d.y(), d.z());

        double numerator = 0.0, denom = 0.0;
        for (unsigned int f = 0; f < flows.size(); f++)
        {
          if (f == ref_frame)
            continue;
          vnl_matrix_fixed<double, 2, 3> J = jacobian_at_x(cameras[f], xj);
          vnl_double_2 K = J * r;
          double u, v;
          cameras[f].project(xj(0), xj(1), xj(2), u, v);
          double flowu = vil_bilin_interp_safe(flows[f], u, v, 0);
          double flowv = vil_bilin_interp_safe(flows[f], u, v, 1);
          numerator += K(0) *flowu + K(1) * flowv;
          denom += K(0) * K(0) + K(1) * K(1);
        }

        double lambda = numerator / denom;
        depth(i,j,0) += lambda;
      }
    }



  vil_image_view<double> cropped_depth;
  //vpgl_perspective_camera<double> croppedcam = crop_camera(cameras[ref_frame], i0, j0);
  //cropped_depth = vil_crop(depth, i0, ni, j0, nj);


    vil_median(depth, cropped_depth, se);
    depth = cropped_depth;

    //mesh = depth_map_to_mesh(croppedcam, depth);
    mesh = depth_map_to_mesh(cameras[ref_frame], depth);




}


vcl_cout << "writing mesh"<<std::endl;
  imesh_write_obj(output_mesh(),mesh);

  return 0;
}


void save_flows(const vcl_vector<vil_image_view<double> > &flows, const char *filename)
{
  FILE *file = fopen(filename, "wb");

  unsigned int nflows = flows.size();
  fwrite(&nflows, sizeof(unsigned int), 1, file);
  for (unsigned int f = 0; f < nflows; f++)
  {
    unsigned int ni = flows[f].ni(), nj = flows[f].nj();
    fwrite(&ni, sizeof(unsigned int), 1, file);
    fwrite(&nj, sizeof(unsigned int), 1, file);
    for (unsigned int i = 0; i < ni; i++)
    {
      for (unsigned int j = 0; j < nj; j++)
      {
        fwrite(&flows[f](i,j,0), sizeof(double), 1, file);
        fwrite(&flows[f](i,j,1), sizeof(double), 1, file);
      }
    }
  }

  fclose(file);
}

void load_flows(vcl_vector<vil_image_view<double> > &flows, const char *filename)
{
  FILE *file = fopen(filename, "rb");

  unsigned int nflows;;
  fread(&nflows, sizeof(unsigned int), 1, file);
  flows.resize(nflows);
  for (unsigned int f = 0; f < nflows; f++)
  {
    unsigned int ni, nj;
    fread(&ni, sizeof(unsigned int), 1, file);
    fread(&nj, sizeof(unsigned int), 1, file);
    flows[f].set_size(ni, nj, 2);
    for (unsigned int i = 0; i < ni; i++)
    {
      for (unsigned int j = 0; j < nj; j++)
      {
        fread(&flows[f](i,j,0), sizeof(double), 1, file);
        fread(&flows[f](i,j,1), sizeof(double), 1, file);
      }
    }
  }

  fclose(file);
}


vnl_matrix_fixed<double, 2, 3> jacobian_at_x(const vpgl_perspective_camera<double> &cam, const vnl_double_3 &x)
{
  vnl_matrix_fixed<double, 2, 3> Jac;
  const vnl_matrix_fixed<double, 3, 4> &P = cam.get_matrix();
  vnl_vector_fixed<double, 4> xh(x(0), x(1), x(2), 1.0);
  vnl_double_3 proj = P * xh;
  double denom = proj(2) * proj(2);
  Jac(0,0) = (P(0,0)*proj(2) - P(2,0)*proj(0))/denom;
  Jac(0,1) = (P(0,1)*proj(2) - P(2,1)*proj(0))/denom;
  Jac(0,2) = (P(0,2)*proj(2) - P(2,2)*proj(0))/denom;
  Jac(1,0) = (P(1,0)*proj(2) - P(2,0)*proj(1))/denom;
  Jac(1,1) = (P(1,1)*proj(2) - P(2,1)*proj(1))/denom;
  Jac(1,2) = (P(1,2)*proj(2) - P(2,2)*proj(1))/denom;

  return Jac;
}

void project_ref_image(const vpgl_perspective_camera<double> &camref,
                const vpgl_perspective_camera<double> &camtarget,
                imesh_mesh &mesh,
                const vil_image_view<double> &ref,
                vil_image_view<double> &reproj)
{
  vil_image_view<double> depth(ref.ni(), ref.nj(), 1);
  imesh_project_depth(mesh, camtarget, depth);

  reproj.set_size(ref.ni(), ref.nj(), 1);
  for (unsigned int i = 0; i < reproj.ni(); i++)
  {
    for (unsigned int j = 0; j < reproj.nj(); j++)
    {
      if (depth(i,j,0) == vcl_numeric_limits<double>::infinity())
        reproj(i,j,0) = 0.0;
      else
      {
        vgl_line_3d_2_points<double> ray = camtarget.backproject(vgl_point_2d<double>((double)i, (double)j));
        vgl_vector_3d<double> d = ray.direction();
        normalize(d);
        vgl_point_3d<double> pt3d = camtarget.get_camera_center() + depth(i,j,0)*d;

        double u, v;
        camref.project(pt3d.x(), pt3d.y(), pt3d.z(), u, v);
        reproj(i,j,0) = vil_bilin_interp_safe(ref, u, v);
      }
    }
  }

  //vil_image_view<vxl_byte> output;
  //vil_convert_cast<double, vxl_byte>(reproj, output);
  //vil_save(output, "testwarp.png");
}

///Converts optical flow 2 plane image to an RGB image
void flow_to_colormap(const vil_image_view<double> &flow, vil_image_view<vil_rgb<vxl_byte> > &colormap_flow, double max)
{
  //If no max flow distance was specified, then use the max for the flow.
  if (max <= 0.0)
  {
    for (unsigned int i = 0; i < flow.ni(); i++)
    {
      for (unsigned int j = 0; j < flow.nj(); j++)
      {
        double len = vnl_double_2(flow(i, j, 0), flow(i, j, 1)).two_norm();
        if (len > max)
          max = len;
      }
    }
  }

  vcl_cout << "Using max flow: " << max << ".\n";

  colormap_flow.set_size(flow.ni(), flow.nj(), 1);
  //Convert flow to HSV to RGB
  for (unsigned int i = 0; i < flow.ni(); i++)
  {
    for (unsigned int j = 0; j < flow.nj(); j++)
    {
       vnl_double_2 uv(flow(i, j, 0), flow(i, j, 1));
       double sat = uv.two_norm() / max;
       if (sat < 1e-6)
         colormap_flow(i, j) = vil_rgb<vxl_byte>(0, 0, 0);

       double hue = (atan2(uv(1), uv(0)) * 57.2957795);
       if (hue < 0.0)
         hue += 360.0;
       double hueprime = hue / 60.0;
       double C = sat;
       double X = C * (1.0 - fabs(fmod(hueprime,2.0) - 1.0));
       double m = 1.0 - C;  //Use white for center
       vxl_byte Cmb = (vxl_byte)((C+m)*255.0);
       vxl_byte Xmb = (vxl_byte)((X+m)*255.0);
       vxl_byte mb = (vxl_byte)(m*255.0);
       if (hueprime < 1.0)
         colormap_flow(i, j) = vil_rgb<vxl_byte>(Cmb, Xmb, mb);
       else if (hueprime < 2.0)
         colormap_flow(i, j) = vil_rgb<vxl_byte>(Xmb, Cmb, mb);
       else if (hueprime < 3.0)
         colormap_flow(i, j) = vil_rgb<vxl_byte>(mb, Cmb, Xmb);
       else if (hueprime < 4.0)
         colormap_flow(i, j) = vil_rgb<vxl_byte>(mb, Xmb, Cmb);
       else if (hueprime < 5.0)
         colormap_flow(i, j) = vil_rgb<vxl_byte>(Xmb, mb, Cmb);
       else if (hueprime < 6.0)
         colormap_flow(i, j) = vil_rgb<vxl_byte>(Cmb, mb, Xmb);
    }
  }
}

void DifferenceImage(const vil_image_view<double> &I0, const vil_image_view<double> &I1, const vil_image_view<double> &flow, vil_image_view<double> &diff)
{
  diff.set_size(I0.ni(), I0.nj(), 1);

  for (unsigned int i = 0; i < diff.ni(); i++)
  {
    for (unsigned int j = 0; j < diff.nj(); j++)
    {
      double x = (double)i + flow(i,j,0);
      double y = (double)j + flow(i,j,1);
      diff(i, j) = fabs(I0(i,j,0) - vil_bilin_interp_safe(I1, x, y, 0));
    }
  }
}

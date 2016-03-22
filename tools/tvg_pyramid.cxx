/*ckwg +29
 * Copyright 2012 by Kitware, Inc.
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

#include "tv_refine.h"

#include "file_io.h"
#include "depth_map.h"
#include "multiscale.h"
#include "super_config.h"

#include <video/gaussian_pyramid_builder.h>

// VXL includes
#include <vil/vil_image_view.h>
#include <vil/vil_convert.h>
#include <vil/vil_save.h>
#include <vil/vil_bilin_interp.h>
#include <vil/vil_crop.h>
#include <vul/vul_arg.h>
#include <vnl/vnl_double_3.h>
#include <vnl/vnl_inverse.h>
#include <imesh/imesh_mesh.h>
#include <imesh/imesh_fileio.h>
#include <imesh/algo/imesh_project.h>
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
                              const super3d::gaussian_pyramid_builder &gpb,
                              vcl_vector<vil_image_view<double> > &pyramid);


int main(int argc, char* argv[])
{
  try  {

  config::inst()->read_config(argv[1]);

  vcl_vector<vil_image_view<double> > frames;
  vcl_vector<vpgl_perspective_camera<double> >  cameras;
  vcl_vector<vcl_string> filenames;

  if (!load_frames_and_cameras(frames, cameras, filenames))
    return -1;

  unsigned int ref_frame = config::inst()->get_value<unsigned int>("ref_frame");
  vpgl_perspective_camera<double> ref_cam = cameras[ref_frame];

  double camera_scale = 1.0;
  if (config::inst()->is_set("camera_scale"))
  {
    camera_scale = config::inst()->get_value<double>("camera_scale");
    for (unsigned int i = 0; i < cameras.size(); i++)
      cameras[i] = scale_camera(cameras[i], camera_scale);
  }

  int i0, ni, j0, nj;

  //Compute the window cropping, scale the cropping by the specified scale so that we do not
  //need to recompute cropping input for super resolution.
  if (config::inst()->is_set("crop_window"))
  {
    vcl_istringstream cwstream(config::inst()->get_value<vcl_string>("crop_window"));
    cwstream >> i0 >> ni >> j0 >> nj;
    i0 = (int)(i0*camera_scale);
    j0 = (int)(j0*camera_scale);
    ni = (int)(ni*camera_scale);
    nj = (int)(nj*camera_scale);
    vcl_cout << "Crop window: " << i0 << " " << ni << " " << j0 << " " << nj << "\n";
    frames[ref_frame] = vil_crop(frames[ref_frame], i0, ni, j0, nj);
    cameras[ref_frame] = crop_camera(cameras[ref_frame], i0, j0);
  }
  else
  {
    i0 = j0 = 0;
    ni = frames[ref_frame].ni();
    nj = frames[ref_frame].nj();
  }

  vcl_cout << "Making image pyramids"<<vcl_endl;
  unsigned int levels = config::inst()->get_value<unsigned int>("levels");
  super3d::gaussian_pyramid_builder gpb(levels+1, 2, 1.0);

  vcl_vector<vil_image_view<double> > pyr_ref;
  gpb.build_pyramid(frames[ref_frame], pyr_ref);

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

    gpb.build_pyramid(frames[i], pyrs[index++]);
  }

  double depth_min = config::inst()->get_value<double>("depth_min");
  double depth_max = config::inst()->get_value<double>("depth_max");

  vcl_cout << "Initializing depth map"<<vcl_endl;
  vil_image_view<double> depth;


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
      scaled_H.push_back(scale_homography(Hs[i],scale));
      scaled_e.push_back(scale_point(es[i],scale));
    }

    refine_depths(pyr_ref[s],scaled_I1,
                  scaled_H, scaled_e,
                  depth, 10, 10, theta, lambda);

    if (s != 0)
    {
      vil_image_view<double> temp;
      vil_resample_bilin(depth, temp, pyr_ref[s-1].ni(), pyr_ref[s-1].nj());
      depth = temp;
    }
  }


  }  catch (const config::cfg_exception &e)  {
    vcl_cout << "Error in config: " << e.what() << "\n";
  }

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

/// compute the depth refinement using multiple images.
/// \note The size of vectors \a I1, \a H, and \a e should be the same.
/// \param I0 The reference image.
/// \param I1 The vector of destination images.
/// \param H The vector of homographies induced by the planes at inifinity.
/// \param e The vector eipoles in the destination images.
/// \retval depth The current depth map image (relative to \a I0)
/// \param outer_iterations The number of outer (warping) iterations
/// \param inner_iterations The nummber of inner (propagating) iterations
void
refine_depths(const vil_image_view<double>& I0,
              const std::vector<vil_image_view<double> >& I1,
              const std::vector<vnl_matrix_fixed<double,3,3> >& H,
              const std::vector<vnl_vector_fixed<double,3> >& e,
              vil_image_view<double>& depth,
              unsigned outer_iterations,
              unsigned inner_iterations,
              double theta,
              double lambda)
{
  const unsigned num_views = I1.size();
  assert(H.size() == num_views);
  assert(e.size() == num_views);

  vil_structuring_element se;
  se.set_to_disk(1.9);

  // compute gradients
  std::vector<vil_image_view<double> > I1xy(num_views);
  for (unsigned i=0; i<num_views; ++i)
  {
    I1xy[i] = vil_image_view<double> (I1[i].ni(),I1[i].nj(),2);
    vil_sobel_3x3(I1[i], I1xy[i]);
  }

  vil_image_view<double> bcc(I0.ni(), I0.nj(), 2*num_views);
  vil_image_view<double> new_depth;
  vil_image_view<double> last_depth;
  last_depth.deep_copy(depth);

  vil_image_view<double> I0xy;
  vil_sobel_3x3(I0, I0xy);
  vil_image_view<double> weights(I0.ni(), I0.nj(), 1);

  for (unsigned int i = 0; i < weights.ni(); i++)
  {
    for (unsigned int j = 0; j < weights.nj(); j++)
    {
      double mag = fabs(I0xy(i,j,0)) + fabs(I0xy(i,j,1));
      if (mag > 0.75)
        mag = 0.75;
      weights(i,j,0) = exp(-15.0 * mag);
    }
  }

  for (unsigned int w = 0; w < 5; w++)
  {
    std::cout << "iteration "<< w;
    for (unsigned i=0; i<num_views; ++i)
    {
      vil_image_view<double> bcc_i = vil_planes(bcc, 2*i, 1, 2);
      compute_linear_bcc(I0, I1[i], I1xy[i], depth, H[i], e[i], bcc_i);
    }

    for (unsigned k=0; k<outer_iterations; ++k)
    {


     apply_bcc_to_depth_all(bcc, lambda, theta, depth);


      //super3d::dual_rof_denoise(depth,new_depth,inner_iterations,theta);

      //vil_median(new_depth,depth,se);
    }
    std::cout << "  ssd "<<vil_math_ssd(depth,last_depth,double())<<std::endl;
    last_depth.deep_copy(depth);
  }

  vil_image_view<vxl_byte> depth_byte;
  vil_convert_stretch_range_limited(depth,depth_byte,0.94,0.96);
  std::stringstream file;
  file << "iteration_"<<depth_byte.ni()<<".png";
  vil_save(depth_byte, file.str().c_str());
}

/// compute the linear brightness constancy constraint at each pixel
/// \param I0 The reference image.
/// \param I1 The destination image.
/// \param I1xy The destination image gradients (x in plane 0, y in plane 1).
/// \param depth The current depth map image.
/// \param H The homography induced by the plane at inifinity.
/// \param e The eipole in the destination image.
/// \retval bcc The linear brighness constance constraints (2-plane image).
///             BCC is satisfied when bcc(i,j,0) + depth(i,j)*bcc(i,j,1) == 0
void
compute_linear_bcc(const vil_image_view<double>& I0,
                   const vil_image_view<double>& I1,
                   const vil_image_view<double>& I1xy,
                   const vil_image_view<double>& depth,
                   const vnl_matrix_fixed<double,3,3>& H,
                   const vnl_vector_fixed<double,3>& e,
                   vil_image_view<double>& bcc)
{
  unsigned ni = I0.ni(), nj = I0.nj();
  assert(I1.ni()==ni && I1.nj()==nj);
  assert(I1xy.ni()==ni && I1xy.nj()==nj);
  assert(depth.ni()==ni && depth.nj()==nj);
  assert(I0.nplanes() == 1);
  assert(I1.nplanes() == 1);
  assert(I1xy.nplanes() == 2);
  assert(depth.nplanes() == 1);
  bcc.set_size(ni,nj,2);

  vnl_vector_fixed<double,3> hx = H.get_column(0);
  vnl_vector_fixed<double,3> hy = H.get_column(1);
  vnl_vector_fixed<double,3> h1 = H.get_column(2);

  vcl_ptrdiff_t istep0=I0.istep(),    jstep0=I0.jstep();
  vcl_ptrdiff_t istepD=depth.istep(), jstepD=depth.jstep();
  vcl_ptrdiff_t istepB=bcc.istep(),   jstepB=bcc.jstep(),   pstepB=bcc.planestep();


  const double* row0 = I0.top_left_ptr();
  const double* rowD = depth.top_left_ptr();
  double*       rowB = bcc.top_left_ptr();
  for (unsigned j=0; j<nj; ++j, row0+=jstep0, rowD+=jstepD, rowB+=jstepB)
  {
    vnl_double_3 hj(hy[0]*j,hy[1]*j,hy[2]*j);
    const double* pixel0 = row0;
    const double* pixelD = rowD;
    double*       pixelB0 = rowB;
    for (unsigned i=0; i<ni; ++i, pixel0+=istep0, pixelD+=istepD, pixelB0+=istepB)
    {
      double*       pixelB1 = pixelB0 + pstepB;

      vnl_double_3 h(hx[0]*i+hj[0]+h1[0],
                     hx[1]*i+hj[1]+h1[1],
                     hx[2]*i+hj[2]+h1[2]);

      // compute the (x,y) coordinates mapped into I1
      double w = (*pixelD)*h[2] + e[2];
      double x = ((*pixelD)*h[0] + e[0]) / w;
      double y = ((*pixelD)*h[1] + e[1]) / w;

      if (x < 0.0 || x > ni-1 || y < 0.0 || y > nj-1)
      {
        *pixelB0 = 0.0;
        *pixelB1 = 0.0;
        continue;
      }

      double I1_interp = vil_bilin_interp_safe(I1,x,y,0);

      // compute the change in pixel intensity between images
      double dI_dt = I1_interp - (*pixel0);

      // interpolated image gradients
      double I1x = vil_bilin_interp_safe(I1xy,x,y,0);
      double I1y = vil_bilin_interp_safe(I1xy,x,y,1);

      // compute the derivative (dx,dy) of the coordinates w.r.t. depth
      double denom = (*pixelD)*h[2]+e[2];
      denom *= denom;
      double dx = (h[0]*e[2] - h[2]*e[0])/denom;
      double dy = (h[1]*e[2] - h[2]*e[1])/denom;

      // derivative of image intensity (I1) w.r.t. depth
      double dI_dD = I1x*dx + I1y*dy;

      *pixelB0 = dI_dt - (*pixelD)*dI_dD;
      *pixelB1 = dI_dD;
    }
  }
}

void
apply_bcc_to_depth_all(const vil_image_view<double>& bcc,
                       double lambda, double theta,
                       vil_image_view<double>& depth)
{
  const unsigned ni = bcc.ni(), nj = bcc.nj();
  const unsigned num_views = bcc.nplanes()/2;
  assert(depth.ni()==ni && depth.nj()==nj);
  assert(bcc.nplanes() == 2*num_views);
  assert(depth.nplanes() == 1);

  vcl_ptrdiff_t istepB=bcc.istep(),   jstepB=bcc.jstep(),   pstepB=bcc.planestep();
  vcl_ptrdiff_t istepD=depth.istep(), jstepD=depth.jstep();

  const double*   rowB = bcc.top_left_ptr();
  double*         rowD = depth.top_left_ptr();

  for (unsigned j=0; j<nj; ++j, rowB+=jstepB, rowD+=jstepD)
  {
    const double* pixelB0 = rowB;
    double*       pixelD = rowD;
    for (unsigned i=0; i<ni; ++i, pixelB0+=istepB, pixelD+=istepD)
    {
      double u = *pixelD;
      const double* B0 = pixelB0;
      const double* B1 = pixelB0 + pstepB;

      double l = lambda;

      double min = -(*B0)/(*B1);
      double beste = eval_vij(bcc, l, theta, u, min, i, j);
      unsigned m = 0;
      double bestl = l;

      B0 = B1 + pstepB;
      B1 = B0 + pstepB;
      for (unsigned k = 1; k < num_views; ++k)
      {
        double val = -(*B0)/(*B1);

        double e = eval_vij(bcc, l, theta, u, val, i, j);
        if (e < beste)
        {
          beste = e;
          min = val;
          m = k;
          bestl = l;
        }
        B0 = B1 + pstepB;
        B1 = B0 + pstepB;
      }

      B0 = pixelB0;
      B1 = pixelB0 + pstepB;
      double offset = 0.0;
      for (unsigned k = 0; k < num_views; ++k)
      {
        if (k != m)
        {
          double pi_min = (*B0) + min * (*B1);
          if (pi_min < 0)
            offset += -bcc(i,j,2*k+1);
          else
            offset += bcc(i,j,2*k+1);
        }
        B0 = B1 + pstepB;
        B1 = B0 + pstepB;
      }

      double Imx = bcc(i,j,2*m+1);
      double pm_u = bcc(i,j,2*m) + u*Imx;
      double coeff = bestl * theta / (double)num_views;

      double newdepth = u;
      if ( pm_u > coeff * Imx * ( Imx + offset))
        newdepth = u - coeff * (Imx + offset);
      else if ( pm_u < -coeff * Imx * (Imx - offset))
        newdepth = u + coeff * (Imx - offset);
      else if (Imx != 0.0)
        newdepth = u - pm_u/Imx;

      depth(i,j) = newdepth;
    }
  }
}

double
eval_vij(const vil_image_view<double>& bcc,
                           double lambda, double theta,
                           double depth,
                           double v, unsigned i, unsigned j)
{
  const unsigned num_views = bcc.nplanes()/2;

  double u = depth;

  double e = 0.0;
  for (unsigned k = 0; k < num_views; k++)
  {
    e += (1.0/(2.0*theta))*(u-v)*(u-v) + lambda * fabs(bcc(i,j,2*k) + v*bcc(i,j,2*k+1));
  }
  return e;
}

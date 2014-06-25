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

#include <vcl_sstream.h>
#include <vcl_fstream.h>
#include <vcl_iomanip.h>

#include <vil/vil_convert.h>
#include <vil/vil_save.h>
#include <vil/vil_bilin_interp.h>
#include <vil/algo/vil_sobel_3x3.h>
#include <vil/algo/vil_median.h>
#include <vnl/vnl_double_3.h>

#include <vil/vil_math.h>

#include <video_transforms/dual_rof_denoise.h>

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
  vil_image_view<double> mask(I0.ni(), I0.nj(), num_views);
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
      vil_image_view<double> mask_i = vil_planes(mask, i, 1, 1);
      compute_linear_bcc(I0, I1[i], I1xy[i], depth, H[i], e[i], mask_i, bcc_i);
    }

//    if (I0.ni() > 2000)
//      output_function(bcc, lambda, theta, depth(1958, 981), 1958, 981);
    for (unsigned k=0; k<outer_iterations; ++k)
    {


     apply_bcc_to_depth_all(bcc, mask, lambda, theta, depth);
     //apply_bcc_to_depth_lst_sqr_two(bcc,lambda*theta,depth);

      //vidtk::dual_rof_denoise(depth,new_depth,inner_iterations,theta);
      vidtk::dual_rof_weighted_denoise(depth, weights, new_depth, inner_iterations, theta);
      vil_median(new_depth,depth,se);
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
                   vil_image_view<double>& mask,
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

      double thresh = 0.96;
      if (I1_interp > thresh || (*pixel0) > thresh)
        mask(i,j) = 1.0;
      else
        mask(i,j) = 0.0;

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


/// Apply the brightness constancy constraint to refine depth
/// \param bcc The linear brighness constance constraints (2-plane image).
///            BCC is satisfied when bcc(i,j,0) + depth(i,j)*bcc(i,j,1) == 0
/// \param step Determines the threshold and ammount to step
/// \retval depth The current depth map image, modified in place.
void
apply_bcc_to_depth(const vil_image_view<double>& bcc,
                   double step,
                   vil_image_view<double>& depth)
{
  unsigned ni = bcc.ni(), nj = bcc.nj();
  assert(depth.ni()==ni && depth.nj()==nj);
  assert(bcc.nplanes() == 2);
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
      const double* pixelB1 = pixelB0 + pstepB;

      double val = (*pixelB0) + (*pixelD)*(*pixelB1);
      double thresh = step*(*pixelB1)*(*pixelB1);
      if (val < -thresh)
      {
        *pixelD += step*(*pixelB1);
      }
      else if (val > thresh)
      {
        *pixelD -= step*(*pixelB1);
      }
      else if (*pixelB1 != 0.0)
      {
        *pixelD -= val / (*pixelB1);
      }
    }
  }
}


/// Apply multiple image BCCs for a least squares depth refinement
/// \param bcc The linear brighness constance constraints (2 planes per image).
///            BCC is satisfied for image k when
///            bcc(i,j,2*k) + depth(i,j)*bcc(i,j,2*k+1) == 0
/// \param lambda Determines the influence of the data (BCC) term.
/// \retval depth The current depth map image, modified in place.
void
apply_bcc_to_depth_lst_sqr_two(const vil_image_view<double>& bcc,
                               const vil_image_view<double>& mask,
                               double lambda,
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
      double sum_ab = 0.0;
      double sum_bb = 0.0;
      const double* planeBa = pixelB0;
      const double* planeBb = pixelB0 + pstepB;
      for (unsigned k=0; k<num_views; ++k)
      {
        sum_ab += (*planeBa) * (*planeBb);
        sum_bb += (*planeBb) * (*planeBb);
        planeBa = planeBb + pstepB;
        planeBb = planeBa + pstepB;
      }

      *pixelD -=  lambda*sum_ab;
      *pixelD /= 1 + lambda*sum_bb;
    }
  }
}


void
apply_bcc_to_depth_all(const vil_image_view<double>& bcc,
                       const vil_image_view<double>& mask,
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

      double l;
      if (mask(i,j,0) == 1.0)
        l = 1e-6;
      else
        l = lambda;

      double min = -(*B0)/(*B1);
      double beste = eval_vij(bcc, l, theta, u, min, i, j);
      unsigned m = 0;
      double bestl = l;

      B0 = B1 + pstepB;
      B1 = B0 + pstepB;
      for (unsigned k = 1; k < num_views; ++k)
      {
        if (mask(i,j,k) == 1.0)
          l = 1e-6;
        else
          l = lambda;

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

void output_function(const vil_image_view<double>& bcc,
                     double lambda, double theta,
                     double depth,
                     unsigned i, unsigned j)
{
  const unsigned num_views = bcc.nplanes()/2;

  vcl_ofstream pts("pts.txt");
  double min = 0; double min_e = 1e10;
  for (unsigned k = 0; k < num_views; ++k)
  {
    double pu = bcc(i,j,2*k) + depth*bcc(i,j,2*k+1);
    double val = depth - pu/bcc(i,j,2*k+1);
    //double val = -bcc(i,j,2*k)/bcc(i,j,2*k+1);

    double e = eval_vij(bcc, lambda, theta, depth, val, i, j);
    pts << val << " " << e << "\n";

    if (e < min_e)
    {
      min_e = e;
      min = val;
    }
  }
  pts.close();

  vcl_ofstream outfile("eval.txt");

  for (double v = min-1.0; v <= min+1.0; v+=0.00001)
  {
    outfile << v << " " << eval_vij(bcc, lambda, theta, depth, v, i, j) << "\n";
  }

  outfile.close();

  min = -(bcc(i,j,0))/(bcc(i,j,1));
  double beste = eval_vij(bcc, lambda, theta, depth, min, i, j);
  unsigned m = 0;
  for (unsigned k = 1; k < num_views; ++k)
  {
    double val = -bcc(i,j,2*k)/bcc(i,j,2*k+1);
    double e = eval_vij(bcc, lambda, theta, depth, val, i, j);
    if (e < beste)
    {
      beste = e;
      min = val;
      m = k;
    }
  }

  double offset = 0.0;
  for (unsigned k = 0; k < num_views; ++k)
  {
    if (k != m)
    {
      double pi_min = bcc(i,j,2*k) + min * bcc(i,j,2*k+1);
      if (pi_min < 0)
        offset += -bcc(i,j,2*k+1);
      else
        offset += bcc(i,j,2*k+1);
    }
  }

  double Imx = bcc(i,j,2*m+1);
  double pm_u = bcc(i,j,2*m) + depth*Imx;
  double coeff = lambda * theta / (double)num_views;

  double final = depth;
  if ( pm_u > coeff * Imx * ( Imx + offset))
    final += -coeff * (Imx + offset);
  else if ( pm_u < -coeff * Imx * (Imx - offset))
    final += coeff * (Imx - offset);
  else
    final += -pm_u/Imx;

  vcl_ofstream minimum("min.txt");
  minimum << final << " " << eval_vij(bcc, lambda, theta, depth, final, i, j);
  minimum.close();

  vcl_cout << "wrote\n";
}

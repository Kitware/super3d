/*ckwg +29
 * Copyright 2011-2015 by Kitware, Inc.
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

#include "dense_optical_flow.h"

#include <sstream>

#include <vil/vil_convert.h>
#include <vil/vil_save.h>
#include <vil/vil_math.h>
#include <vil/vil_bilin_interp.h>
#include <vil/vil_bicub_interp.h>
#include <vil/vil_resample_bilin.h>
#include <vil/vil_resample_bicub.h>
#include <vil/algo/vil_sobel_3x3.h>
#include <vil/algo/vil_median.h>

#include <super3d/image/dual_rof_denoise.h>
#include <super3d/image/gaussian_pyramid_builder.h>


// provide missing instantiation of vil_resample_bicub<double,double>
#include <vil/vil_resample_bicub.txx>
VIL_RESAMPLE_BICUB_INSTANTIATE( double , double );


namespace super3d
{


template <typename T>
dense_optical_flow<T>
::dense_optical_flow()
  : num_pyramid_levels_(4),
    num_warps_(35),
    num_outer_iterations_(5),
    num_inner_iterations_(2),
    num_texture_iterations_(100),
    use_bicubic_interp_(false),
    gradient_blend_(T(0.4)),
    structure_removal_(T(0.95)),
    texture_theta_(T(0.125)),
    theta_(T(0.25)),
    lambda_(T(30.0)),
    verbose_(false)
{
}


/// Normalize the intensity range of a pair of images.
/// Both images are scaled to the range [-1,1] using the same scale factor.
template <typename T>
void
dense_optical_flow<T>
::normalize_intensity_ranges(vil_image_view<T>& I0,
                             vil_image_view<T>& I1) const
{
  T min_0, max_0, min_1, max_1;
  vil_math_value_range(I0, min_0, max_0);
  vil_math_value_range(I1, min_1, max_1);
  min_0 = std::min(min_0, min_1);
  max_0 = std::max(max_0, max_1);

  double scale = 2.0 / (max_0 - min_0);
  T offset = static_cast<T>(-1.0 - scale * min_0);
  vil_math_scale_and_offset_values(I0, scale, offset);
  vil_math_scale_and_offset_values(I1, scale, offset);
}


/// compute the flow optical flow between images \a I0 and \a I1
/// This function uses coarse to fine
/// \param I0 The reference image.
/// \param I1 The target image.
/// \retval flow The current flow image (relative to \a I0)
template <typename T>
void
dense_optical_flow<T>
::compute_flow(const vil_image_view<T>& I0,
               const vil_image_view<T>& I1,
               vil_image_view<T>& flow) const
{
  vil_image_view<T> I0tex, I1tex;
  I0tex.deep_copy(I0);
  I1tex.deep_copy(I1);
  if (verbose_)
  {
    std::cout << "Normalizing images" << std::endl;
  }
  normalize_intensity_ranges(I0tex, I1tex);
  if (structure_removal_ > T(0) && num_texture_iterations_ > 0)
  {
    if (verbose_)
    {
      std::cout << "Computing texture images" << std::endl;
    }
    remove_structure(I0tex, structure_removal_,
                     texture_theta_, num_texture_iterations_);
    remove_structure(I1tex, structure_removal_,
                     texture_theta_, num_texture_iterations_);
    normalize_intensity_ranges(I0tex, I1tex);
  }

  gaussian_pyramid_builder gpb(num_pyramid_levels_, 2, 1.0);

  std::vector<vil_image_view<T> > I0_pyr, I0_grad_pyr, I1_pyr, I1_grad_pyr;
  gpb.build_pyramid<T, T>(I0tex, I0_pyr, I0_grad_pyr);
  gpb.build_pyramid<T, T>(I1tex, I1_pyr, I1_grad_pyr);

  vil_image_view<T> dualx, dualy;

  typedef void (*resample_func)(const vil_image_view<T> &src,
                                vil_image_view<T> &dest,
                                int ni, int nj);
  resample_func resample = &vil_resample_bilin;
  if (use_bicubic_interp_)
  {
    resample = &vil_resample_bicub;
  }

  //coarse to fine
  for (int i = num_pyramid_levels_-1; i >= 0; i--)
  {
    const int ni = I0_pyr[i].ni();
    const int nj = I0_pyr[i].nj();
    vil_image_view<T> new_flow(ni, nj, 2);
    vil_image_view<T> new_dualx(ni, nj, 2);
    vil_image_view<T> new_dualy(ni, nj, 2);

    //Don't scale the initial 0 flow
    if (i < static_cast<int>(num_pyramid_levels_) - 1)
    {
      zero_boundaries(dualx);
      zero_boundaries(dualy);
      (*resample)(flow, new_flow, ni, nj);
      (*resample)(dualx, new_dualx, ni, nj);
      (*resample)(dualy, new_dualy, ni, nj);
      vil_math_scale_values(new_flow, 2.0);
    }
    else
    {
      new_flow.fill(T(0));
      new_dualx.fill(T(0));
      new_dualy.fill(T(0));
    }

    refine_flow(I0_pyr[i], I1_pyr[i], I0_grad_pyr[i], I1_grad_pyr[i],
                new_flow, new_dualx, new_dualy);

    flow = new_flow;
    dualx = new_dualx;
    dualy = new_dualy;
  }
}


/// compute the flow refinement using two images \a I0 and \a I1
/// \param I0 The reference image.
/// \param I1 The target image.
/// \retval flow The current flow image (relative to \a I0)
template <typename T>
void
dense_optical_flow<T>
::refine_flow(const vil_image_view<T>& I0,
             const vil_image_view<T>& I1,
             vil_image_view<T>& I0xy,
             vil_image_view<T>& I1xy,
             vil_image_view<T>& flow,
             vil_image_view<T>& dualx,
             vil_image_view<T>& dualy) const
{

  vil_sobel_3x3(I1, I1xy);
  vil_sobel_3x3(I0, I0xy);
  vil_image_view<T> bcc;
  vil_image_view<T> new_flow(I1.ni(), I1.nj(), 2);

  vil_structuring_element se;
  se.set_to_disk(1.9);

  vil_image_view<T> flow_x = vil_plane(flow, 0);
  vil_image_view<T> flow_y = vil_plane(flow, 1);
  vil_image_view<T> new_flow_x = vil_plane(new_flow, 0);
  vil_image_view<T> new_flow_y = vil_plane(new_flow, 1);


  vil_image_view<T> last_flow;
  if (verbose_)
  {
    last_flow.deep_copy(flow);
  }

  for (unsigned int w = 0; w < num_warps_; w++)
  {
    compute_linear_bcc(I0, I1, I0xy, I1xy, flow, bcc);
    for (unsigned k=0; k<num_outer_iterations_; ++k)
    {
      apply_bcc_to_flow(bcc,lambda_*theta_,flow);

      dual_rof_denoise(flow_x,new_flow_x,dualx,num_inner_iterations_,theta_);
      dual_rof_denoise(flow_y,new_flow_y,dualy,num_inner_iterations_,theta_);
      vil_median(new_flow_x,flow_x,se);
      vil_median(new_flow_y,flow_y,se);
    }

    if (verbose_)
    {
      std::cout << "warp " << w << ",  ssd "
                << vil_math_ssd(flow,last_flow,T()) << std::endl;
      last_flow.deep_copy(flow);
    }
  }
}


/// compute the linear brightness constancy constraint at each pixel
/// \param I0 The reference image.
/// \param I1 The target image.
/// \param I1xy The target image gradients (x in plane 0, y in plane 1).
/// \param flow The current flow image (x in plane 0, y in plane 1).
/// \retval bcc The linear brighness constance constraints (3-plane image).
///             BCC is satisfied when
///             bcc(i,j,0) + flow(i,j,0)*bcc(i,j,1) + flow(i,j,1)*bcc(i,j,2) == 0
template <typename T>
void
dense_optical_flow<T>
::compute_linear_bcc(const vil_image_view<T>& I0,
                     const vil_image_view<T>& I1,
                     const vil_image_view<T>& I0xy,
                     const vil_image_view<T>& I1xy,
                     const vil_image_view<T>& u0,
                     vil_image_view<T>& bcc) const
{
  const unsigned ni = I0.ni(), nj = I0.nj();
  assert(I1.ni()==ni && I1.nj()==nj);
  assert(I1xy.ni()==ni && I1xy.nj()==nj);
  assert(u0.ni()==ni && u0.nj()==nj);
  assert(I0.nplanes() == 1);
  assert(I1.nplanes() == 1);
  assert(I1xy.nplanes() == 2);
  assert(u0.nplanes() == 2);
  bcc.set_size(ni,nj,3);

  const std::ptrdiff_t istep0=I0.istep(),   jstep0=I0.jstep();
  const std::ptrdiff_t istepF=u0.istep(),   jstepF=u0.jstep(), pstepF=u0.planestep();
  const std::ptrdiff_t istepB=bcc.istep(),  jstepB=bcc.jstep(),  pstepB=bcc.planestep();

  typedef double (*interpolator_func)(const vil_image_view<T> &view,
                                      double x, double y, unsigned p);
  unsigned min_i = 0, max_i = ni-1;
  unsigned min_j = 0, max_j = nj-1;
  interpolator_func interp = &vil_bilin_interp;
  if (use_bicubic_interp_)
  {
    interp = &vil_bicub_interp;
    min_i = 1;
    max_i = ni-2;
    min_j = 1;
    max_j = nj-2;
  }



  const T* row0 = I0.top_left_ptr();
  const T* rowF = u0.top_left_ptr();
  T*       rowB = bcc.top_left_ptr();
  for (unsigned j=0; j<nj; ++j, row0+=jstep0, rowF+=jstepF, rowB+=jstepB)
  {
    const T* pixel0 = row0;
    const T* pixelFx = rowF;
    T*       pixelB0 = rowB;
    for (unsigned i=0; i<ni; ++i, pixel0+=istep0, pixelFx+=istepF, pixelB0+=istepB)
    {
      const T* pixelFy = pixelFx + pstepF;
      T*       pixelBx = pixelB0 + pstepB;
      T*       pixelBy = pixelBx + pstepB;

      // At this point the image location is (i,j), and the following are true
      // *pixel0  == I0(i,j)
      // *pixelFx == flow(i,j,0)
      // *pixelFy == flow(i,j,1)
      // *pixelB0 == bcc(i,j,0)
      // *pixelBx == bcc(i,j,1)
      // *pixelBy == bcc(i,j,2)

      // compute the (x,y) coordinates mapped into I1
      const T x = *pixelFx + i;
      const T y = *pixelFy + j;

      // if flow maps outside image then set bcc coefficients to zero
      if (x < min_i || x > max_i || y < min_j || y > max_j)
      {
        *pixelB0 = T(0);
        *pixelBx = T(0);
        *pixelBy = T(0);
        continue;
      }

      // get the corresponding (interpolated) intensity at I1(x,y)
      const double I1_interp = (*interp)(I1,x,y,0);

      // compute the change in pixel intensity between images
      const double dI_dt = I1_interp - (*pixel0);

      // get interpolated image gradients and store
      // into bcc(i,j,1) and bcc(i,j,2)
      T& I1x = *pixelBx;
      T& I1y = *pixelBy;
      I1x = static_cast<T>(gradient_blend_ * I0xy(i,j,0) + (1.0-gradient_blend_) * (*interp)(I1xy,x,y,0));
      I1y = static_cast<T>(gradient_blend_ * I0xy(i,j,1) + (1.0-gradient_blend_) * (*interp)(I1xy,x,y,1));

      // compute the constant term, bcc(i,j,0)
      *pixelB0 = static_cast<T>(dI_dt - I1x*(*pixelFx) - I1y*(*pixelFy));
    }
  }
}


/// Apply the brightness constancy constraint to refine flow
/// \param bcc The linear brighness constance constraints (3-plane image).
///            BCC is satisfied when
///            bcc(i,j,0) + flow(i,j,0)*bcc(i,j,1) + flow(i,j,1)*bcc(i,j,2) == 0
/// \param step Determines the threshold and ammount to step
/// \retval flow The current flow image, modified in place.
template <typename T>
void
dense_optical_flow<T>
::apply_bcc_to_flow(const vil_image_view<T>& bcc,
                    T step,
                    vil_image_view<T>& flow) const
{
  const unsigned ni = bcc.ni(), nj = bcc.nj();
  assert(flow.ni()==ni && flow.nj()==nj);
  assert(bcc.nplanes() == 3);
  assert(flow.nplanes() == 2);

  std::ptrdiff_t istepB=bcc.istep(),  jstepB=bcc.jstep(),  pstepB=bcc.planestep();
  std::ptrdiff_t istepF=flow.istep(), jstepF=flow.jstep(), pstepF=flow.planestep();


  const T*   rowB = bcc.top_left_ptr();
  T*         rowF = flow.top_left_ptr();
  for (unsigned j=0; j<nj; ++j, rowB+=jstepB, rowF+=jstepF)
  {

    const T* pixelB0 = rowB;
    T*       pixelFx = rowF;
    for (unsigned i=0; i<ni; ++i, pixelB0+=istepB, pixelFx+=istepF)
    {
      const T* pixelBx = pixelB0 + pstepB;
      const T* pixelBy = pixelBx + pstepB;
      T*       pixelFy = pixelFx + pstepF;

      const T& I1x = *pixelBx;
      const T& I1y = *pixelBy;

      // evaluate the BCC for the current flow
      // val is the residual of the brightness constancy
      const T val = (*pixelB0) + (*pixelFx)*I1x + (*pixelFy)*I1y;

      // compute L1 magnitude of gradient
      const T grad1 = fabs(I1x) + fabs(I1y);
      // apply 1 of 3 flow updates, see paper for details
      const T thresh = step * grad1 * grad1;
      if (val < -thresh)
      {
        *pixelFx += step * I1x;
        *pixelFy += step * I1y;
      }
      else if (val > thresh)
      {
        *pixelFx -= step * I1x;
        *pixelFy -= step * I1y;
      }
      else if (grad1 != T(0))
      {
        const T step2 = -val / (grad1 * grad1);
        *pixelFx += step2 * I1x;
        *pixelFy += step2 * I1y;
      }
    }
  }
}


/// Remove the structure component of the image using TV optimization.
/// The structure is the result of Rudin, Osher and Fatemi noise removal.
/// The resulting image is \a src - \a alpha * structure.
/// \retval img The image to operate on (modified in place)
/// \param alpha The fraction of the structure to remove.
/// \param theta The tuning parameter for structure computation.
/// \param iterations The number of iterations used in solving TV.
template <typename T>
void
dense_optical_flow<T>
::remove_structure(vil_image_view<T>& img,
                   T alpha,
                   T theta,
                   unsigned iterations)
{
  vil_image_view<T> structure;
  dual_rof_denoise(img, structure, iterations, theta);
  vil_math_scale_values(structure, alpha);
  vil_math_image_difference(img,structure,img);
}


/// Set all boundary pixels to zero in all planes.
/// Used to enforce Dirichlet boundary conditions on dual variables.
/// \param img The image to modify in place
template <typename T>
void
dense_optical_flow<T>
::zero_boundaries(vil_image_view<T>& img) const
{
  const unsigned ni = img.ni();
  const unsigned nj = img.nj();
  assert(img.nplanes() == 2);
  assert(ni >= 1);
  assert(nj >= 1);

  for (unsigned int i = 0; i < ni; i++)
  {
    img(i,0,0) = T(0);
    img(i,0,1) = T(0);
    img(i,nj-1,0) = T(0);
    img(i,nj-1,1) = T(0);
  }
  for (unsigned int j = 0; j < nj; j++)
  {
    img(0,j,0) = T(0);
    img(0,j,1) = T(0);
    img(ni-1,j,0) = T(0);
    img(ni-1,j,1) = T(0);
  }
}


// ======= Setters =======

/// Set number of pyramid level to use
template <typename T>
dense_optical_flow<T>&
dense_optical_flow<T>
::set_num_pyramid_levels(unsigned num)
{
  num_pyramid_levels_ = num;
  return *this;
}

/// Set number of image warps to do at each pyramid level
template <typename T>
dense_optical_flow<T>&
dense_optical_flow<T>
::set_num_warps(unsigned num)
{
  num_warps_ = num;
  return *this;
}

/// Set number of outer iterations (enforcing brightness constancy)
template <typename T>
dense_optical_flow<T>&
dense_optical_flow<T>
::set_num_outer_iterations(unsigned num)
{
  num_outer_iterations_ = num;
  return *this;
}

/// Set number of inner iterations (denoising flow field)
template <typename T>
dense_optical_flow<T>&
dense_optical_flow<T>
::set_num_inner_iterations(unsigned num)
{
  num_inner_iterations_ = num;
  return *this;
}

/// Set number of denoising iterations to extract texture component of images
template <typename T>
dense_optical_flow<T>&
dense_optical_flow<T>
::set_num_texture_iterations(unsigned num)
{
  num_texture_iterations_ = num;
  return *this;
}

/// Set use of bicubic interpolation instead of bilinear
template <typename T>
dense_optical_flow<T>&
dense_optical_flow<T>
::set_use_bicubic_interp(bool bicub)
{
  use_bicubic_interp_ = bicub;
  return *this;
}

/// Set amount to blend source gradient with warped target gradient
template <typename T>
dense_optical_flow<T>&
dense_optical_flow<T>
::set_gradient_blend(T val)
{
  gradient_blend_ = val;
  return *this;
}

/// Set amount of structure to remove to produce texture images
template <typename T>
dense_optical_flow<T>&
dense_optical_flow<T>
::set_structure_removal(T val)
{
  structure_removal_ = val;
  return *this;
}

/// Set the regularization tuning parameter for structure computation
template <typename T>
dense_optical_flow<T>&
dense_optical_flow<T>
::set_texture_theta(T val)
{
  texture_theta_ = val;
  return *this;
}

/// Set the regularization tuning parameter for flow denosing
template <typename T>
dense_optical_flow<T>&
dense_optical_flow<T>
::set_theta(T val)
{
  theta_ = val;
  return *this;
}

/// Set the weight on the data term in the energy functional
template <typename T>
dense_optical_flow<T>&
dense_optical_flow<T>
::set_lambda(T val)
{
  lambda_ = val;
  return *this;
}

/// Set verbose mode
template <typename T>
dense_optical_flow<T>&
dense_optical_flow<T>
::set_verbose(bool verb)
{
  verbose_ = verb;
  return *this;
}


// Template Instantiations
template class dense_optical_flow<float>;
template class dense_optical_flow<double>;

}

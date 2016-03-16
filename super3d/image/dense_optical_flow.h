/*ckwg +5
 * Copyright 2011-2013 by Kitware, Inc. All Rights Reserved. Please refer to
 * KITWARE_LICENSE.TXT for licensing information, or contact General Counsel,
 * Kitware, Inc., 28 Corporate Drive, Clifton Park, NY 12065.
 */

#ifndef dense_optical_flow_h_
#define dense_optical_flow_h_

#include <vil/vil_image_view.h>

/// \file
/// Dense optical flow using total variation with a L1 norm.
///
/// The implementation below is based on the equations in the following paper:
///
/// A. Wedel, T. Pock, C. Zach, H. Bischof, and D. Cremers,
/// "An improved algorithm for TV-L1 optical flow",
/// Statistical and Geometrical Approaches to Visual Motion Analysis
/// pg 23-25, Springer, 2009.
/// http://cvpr.in.tum.de/old/pub/pub/DagstuhlOpticalFlowChapter.pdf


namespace vidtk
{

/// Class for managing dense optical flow computation.
/// This templated class is only intended to be instantiated with types
/// "float" and "double".  These instantiations are provided in the .cxx file.
template <typename T>
class dense_optical_flow
{
public:
  /// Default constructor
  dense_optical_flow();

  /// compute the flow optical flow between images \a I0 and \a I1
  /// This function uses coarse to fine
  /// \param I0 The reference image.
  /// \param I1 The target image.
  /// \retval flow The current flow image (relative to \a I0)
  void compute_flow(const vil_image_view<T>& I0,
                    const vil_image_view<T>& I1,
                    vil_image_view<T>& flow) const;

  /// Remove the structure component of the image using TV optimization.
  /// The structure is the result of Rudin, Osher and Fatemi noise removal.
  /// The resulting image is \a src - \a alpha * structure.
  /// \retval img The image to operate on (modified in place)
  /// \param alpha The fraction of the structure to remove.
  /// \param theta The tuning parameter for structure computation.
  /// \param iterations The number of iterations used in solving TV.
  static void remove_structure(vil_image_view<T>& img,
                               T alpha = 0.95,
                               T theta = 0.125,
                               unsigned iterations = 100);

  // ======= Getters =======
  /// Number of pyramid level to use
  unsigned num_pyramid_levels() const { return num_pyramid_levels_; }

  /// Number of image warps to do at each pyramid level
  unsigned num_warps() const { return num_warps_; }

  /// Number of outer iterations (enforcing brightness constancy)
  unsigned num_outer_iterations() const { return num_outer_iterations_; }

  /// Number of inner iterations (denoising flow field)
  unsigned num_inner_iterations() const { return num_inner_iterations_; }

  /// Number of denoising iterations to extract texture component of images
  unsigned num_texture_iterations() const { return num_texture_iterations_; }

  /// Use bicubic interpolation instead of bilinear
  bool use_bicubic_interp() const { return use_bicubic_interp_; }

  /// Amount to blend source gradient with warped target gradient
  T gradient_blend() const { return gradient_blend_; }

  /// Amount of structure to remove to produce texture images
  T structure_removal() const { return structure_removal_; }

  /// The regularization tuning parameter for structure computation
  T texture_theta() const { return texture_theta_; }

  /// The regularization tuning parameter for flow denosing
  T theta() const { return theta_; }

  /// The weight on the data term in the energy functional
  T lambda() const { return lambda_; }

  /// If true then display status messages while computing
  bool verbose() const { return verbose_; }


  // ======= Setters =======
  /// Set number of pyramid level to use
  dense_optical_flow<T>& set_num_pyramid_levels(unsigned num);

  /// Set number of image warps to do at each pyramid level
  dense_optical_flow<T>& set_num_warps(unsigned num);

  /// Set number of outer iterations (enforcing brightness constancy)
  dense_optical_flow<T>& set_num_outer_iterations(unsigned num);

  /// Set number of inner iterations (denoising flow field)
  dense_optical_flow<T>& set_num_inner_iterations(unsigned num);

  /// Set number of denoising iterations to extract texture component of images
  dense_optical_flow<T>& set_num_texture_iterations(unsigned num);

  /// Set use of bicubic interpolation instead of bilinear
  dense_optical_flow<T>& set_use_bicubic_interp(bool bicub);

  /// Set amount to blend source gradient with warped target gradient
  dense_optical_flow<T>& set_gradient_blend(T val);

  /// Set amount of structure to remove to produce texture images
  dense_optical_flow<T>& set_structure_removal(T val);

  /// Set the regularization tuning parameter for structure computation
  dense_optical_flow<T>& set_texture_theta(T val);

  /// Set the regularization tuning parameter for flow denosing
  dense_optical_flow<T>& set_theta(T val);

  /// Set the weight on the data term in the energy functional
  dense_optical_flow<T>& set_lambda(T val);

  /// Set verbose mode
  dense_optical_flow<T>& set_verbose(bool verb);

private:
  /// Normalize the intensity range of a pair of images.
  /// Both images are scaled to the range [-1,1] using the same scale factor.
  void normalize_intensity_ranges(vil_image_view<T>& I0,
                                  vil_image_view<T>& I1) const;

  /// compute the flow refinement using two images \a I0 and \a I1
  /// \param I0 The reference image.
  /// \param I1 The target image.
  /// \retval flow The current flow image (relative to \a I0)
  void refine_flow(const vil_image_view<T>& I0,
                   const vil_image_view<T>& I1,
                   vil_image_view<T>& I0xy,
                   vil_image_view<T>& I1xy,
                   vil_image_view<T>& flow,
                   vil_image_view<T>& dualx,
                   vil_image_view<T>& dualy) const;

  /// compute the linear brightness constancy constraint at each pixel
  /// \param I0 The reference image.
  /// \param I1 The target image.
  /// \param I1xy The target image gradients (x in plane 0, y in plane 1).
  /// \param flow The current flow image (x in plane 0, y in plane 1).
  /// \retval bcc The linear brighness constance constraints (3-plane image).
  ///             BCC is satisfied when
  ///             bcc(i,j,0) + flow(i,j,0)*bcc(i,j,1) + flow(i,j,1)*bcc(i,j,2) == 0
  void compute_linear_bcc(const vil_image_view<T>& I0,
                          const vil_image_view<T>& I1,
                          const vil_image_view<T>& I0xy,
                          const vil_image_view<T>& I1xy,
                          const vil_image_view<T>& flow,
                          vil_image_view<T>& bcc) const;

  /// Apply the brightness constancy constraint to refine flow
  /// \param bcc The linear brighness constance constraints (3-plane image).
  ///            BCC is satisfied when
  ///            bcc(i,j,0) + flow(i,j,0)*bcc(i,j,1) + flow(i,j,1)*bcc(i,j,2) == 0
  /// \param step Determines the threshold and ammount to step
  /// \retval flow The current flow image, modified in place.
  void apply_bcc_to_flow(const vil_image_view<T>& bcc,
                         T step,
                         vil_image_view<T>& flow) const;

  /// Set all boundary pixels to zero in all planes.
  /// Used to enforce Dirichlet boundary conditions on dual variables.
  /// \param img The image to modify in place
  void zero_boundaries(vil_image_view<T>& img) const;


  /// Number of pyramid level to use
  unsigned num_pyramid_levels_;
  /// Number of image warps to do at each pyramid level
  unsigned num_warps_;
  /// Number of outer iterations (enforcing brightness constancy)
  unsigned num_outer_iterations_;
  /// Number of inner iterations (denoising flow field)
  unsigned num_inner_iterations_;
  /// Number of denoising iterations to extract texture component of images
  unsigned num_texture_iterations_;
  /// Use bicubic interpolation instead of bilinear
  bool use_bicubic_interp_;
  /// amount to blend source gradient with warped target gradient
  T gradient_blend_;
  /// amount of structure to remove to produce texture images
  T structure_removal_;
  /// The regularization tuning parameter for structure computation
  T texture_theta_;
  /// The regularization tuning parameter for flow denosing
  T theta_;
  /// The weight on the data term in the energy functional
  T lambda_;
  /// If true then display status messages while computing
  bool verbose_;
};


/// Warp an input image using an optical flow field.
/// If \p flow was computed from I0 to I1 then setting \p input to I1
/// would warp I1 back to reference image I0.
/// \param input is the image to warp.
/// \param flow is the dense optical flow field (2-planes for x and y)
/// \param warped is the result image.
template <typename srcT, typename flowT, typename destT>
void warp_with_flow(const vil_image_view<srcT> &input,
                    const vil_image_view<flowT> &flow,
                    vil_image_view<destT> &warped)
{
  const unsigned ni = input.ni();
  const unsigned nj = input.nj();
  const unsigned np = input.nplanes();
  assert(flow.ni() == ni);
  assert(flow.nj() == nj);
  assert(flow.nplanes() == 2);
  warped.set_size(ni, nj, np);

  for (unsigned int j = 0; j < nj; ++j)
  {
    for (unsigned int i = 0; i < ni; ++i)
    {
      double x = static_cast<double>(i) + flow(i,j,0);
      double y = static_cast<double>(j) + flow(i,j,1);
      for (unsigned int p = 0; p < np; ++p)
      {
        warped(i, j, p) = static_cast<destT>(vil_bicub_interp_safe(input, x, y, p));
      }
    }
  }
}

}

#endif // dense_optical_flow_h_

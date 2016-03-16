/*ckwg +5
 * Copyright 2014 by Kitware, Inc. All Rights Reserved. Please refer to
 * KITWARE_LICENSE.TXT for licensing information, or contact General Counsel,
 * Kitware, Inc., 28 Corporate Drive, Clifton Park, NY 12065.
 */

#include <stdio.h>
#include <super3d/image/adjoint_image_derivs.h>

#include <vil/algo/vil_gauss_filter.h>
#include <vil/vil_save.h>
#include <vil/vil_math.h>
#include <vil/vil_resample_bilin.h>
#include <vil/vil_resample_bicub.h>
#include <vil/vil_convert.h>
#include <vil/algo/vil_greyscale_erode.h>
#include <vil/algo/vil_median.h>
#include <vil/vil_resample_bicub.txx>
#include <super3d/image/warp_image.h>

#include "super_res.h"
#include "super_config.h"
#include "multiscale.h"
#include "super_res_robust_function.h"

VIL_RESAMPLE_BICUB_INSTANTIATE( double , double );

namespace super3d
{

namespace // anonymous
{

//*****************************************************************************

void dual_step_grad_prior(const vil_image_view<double> &primal,
                          vil_image_view<double> &dual,
                          const double lambda,
                          const double sigma,
                          const double alpha,
                          const double gamma=0.0,
                          int cost_function=super_res_params::HUBER_NORM)
{
  vil_image_view<double> work;
  vidtk::forward_gradient(primal, work);

  double denom=1.0;
  switch( cost_function )
  {
  case super_res_params::HUBER_NORM:
    denom = 1.0 + ( sigma * alpha ) / lambda;
    break;
  case super_res_params::TRUNCATED_QUADRATIC:
    denom = 1.0 + sigma / (2.0 * gamma * lambda );
    break;
  case super_res_params::GENERALIZED_HUBER:
    denom = 1.0 + sigma / (2.0 * gamma * lambda );
    break;
  default:
    break;
  }

  for (unsigned int j = 0; j < dual.nj(); j++)
  {
    for (unsigned int i = 0; i < dual.ni(); i++)
    {
      for (unsigned int k = 0; k < primal.nplanes(); k++)
      {
        const unsigned k2=2*k;
        const unsigned k1=k2+1;
        double &x = dual(i,j,k2);
        double &y = dual(i,j,k1);
        x = ( x + sigma * work(i,j,k2) ) / denom;
        y = ( y + sigma * work(i,j,k1) ) / denom;

        const double mag = sqrt(x*x + y*y) / lambda;
        if (mag > 1.0)
        {
          x /= mag;
          y /= mag;
        }
      }
    }
  }
}

//*****************************************************************************

void dual_step_pr(const vil_image_view<double> &u,
                  vil_image_view<double> &pr,
                  const super3d::super_res_params &srp)
{
  return dual_step_grad_prior(u, pr, srp.lambda_r, srp.sigma_pr, srp.alpha_r, srp.gamma_r, srp.cost_function);
}

//*****************************************************************************
void dual_step_pl(const vcl_vector< vil_image_view<double> >&As,
                  vcl_vector< vil_image_view<double> >&pl,
                  const super3d::super_res_params &srp)
{
  if( srp.illumination_prior )
  {
    // TODO: the following line is for intensity scaling
    // for(unsigned int i=0; i<As.size(); i++)
    for(unsigned int i=0; i<As.size(); i=i+2)
    {
      if( (i/2) != srp.ref_frame )
      {
        dual_step_grad_prior(As[i], pl[i], srp.lambda_l, srp.sigma_pl, srp.alpha_l, srp.gamma_l, srp.cost_function);
      }
    }
  }
}

//*****************************************************************************

void dual_step_qa(
  const vil_image_view<double> &u,
  const vcl_vector<vil_image_view<double> > &frames,
  vcl_vector<vil_image_view<double> > &qa,
  const vcl_vector< vil_image_view<double> > &As,
  const vcl_vector<vidtk::adjoint_image_ops_func<double> > &warps,
  const vcl_vector<vil_image_view<double> > &weights,
  const super3d::super_res_params &srp)
{
  const double sf_2 = 1.0 / (srp.scale_factor * srp.scale_factor);

  double denom=1.0;
  switch( srp.cost_function )
  {
  case super3d::super_res_params::HUBER_NORM:
    denom = 1.0 + ( srp.sigma_qa * srp.alpha_a) / sf_2;
    break;
  case super3d::super_res_params::TRUNCATED_QUADRATIC:
  case super3d::super_res_params::GENERALIZED_HUBER:
    denom = 1.0 + srp.sigma_qa / (2.0 * srp.gamma_a * sf_2 );
    break;
  }

  unsigned number_of_frames=0;
  unsigned start_frame=0;
  if( srp.image_data_N )
  {
    number_of_frames = frames.size();
  }
  else if( srp.image_data_1 )
  {
    number_of_frames = 1;
    start_frame = srp.ref_frame;
  }

  for (unsigned int f = start_frame; f < start_frame + number_of_frames; f++)
  {
    const unsigned int ni = weights[f].ni();
    const unsigned int nj = weights[f].nj();

    vil_image_view<double> l_u;
    warps[f].apply_A(u, l_u);

    vil_image_view<double> low_res_frame;
    if( srp.illumination_prior )
    {
      // TODO: the following lines are for intensity scaling
      // vil_math_image_product( frames[f], As[2*f+1], low_res_frame);
      // vil_math_image_sum( low_res_frame, As[2*f], low_res_frame);
      vil_math_image_sum( frames[f], As[2*f], low_res_frame);
    }
    else
    {
      low_res_frame = frames[f];
    }

    for (unsigned int j = 0; j < nj; j++)
    {
      for (unsigned int i = 0; i < ni; i++)
      {
        double val = srp.sigma_qa * sf_2 * weights[f](i,j);
        for (unsigned int k = 0; k < qa[f].nplanes(); k++)
        {
          double diff = l_u(i, j, k) - low_res_frame(i, j, k);

          switch( srp.cost_function )
          {
          case super3d::super_res_params::HUBER_NORM:
            break;
          case super3d::super_res_params::TRUNCATED_QUADRATIC:
            diff = super3d::rho_truncated_quadratic(diff, srp.alpha_a, srp.gamma_a);
            break;
          case super3d::super_res_params::GENERALIZED_HUBER:
            diff = super3d::rho_generalized_huber(diff, srp.alpha_a, srp.beta_a,  srp.gamma_a);
            break;
          }

          double &qfijk = qa[f](i, j, k);
          qfijk = (qfijk +  diff * val)/denom;
          qfijk = vcl_max(qfijk, -sf_2);
          qfijk = vcl_min(qfijk, sf_2);
        }
      }
    }

  }
}

//*****************************************************************************

void dual_step_qg(
  const vil_image_view<double> &u,
  const vcl_vector<vil_image_view<double> > &gradient_frames,
  vcl_vector<vil_image_view<double> > &qg,
  const vcl_vector<vidtk::adjoint_image_ops_func<double> > &warps,
  const vcl_vector<vil_image_view<double> > &weights,
  const super3d::super_res_params &srp)
{
  const double sf_2 = 1.0 / (srp.scale_factor * srp.scale_factor);

  double denom=1.0;
  switch( srp.cost_function )
  {
  case super3d::super_res_params::HUBER_NORM:
    denom = 1.0 + ( srp.sigma_qg * srp.alpha_g) / ( srp.lambda_g * sf_2);
    break;
  case super3d::super_res_params::TRUNCATED_QUADRATIC:
  case super3d::super_res_params::GENERALIZED_HUBER:
    denom = 1.0 + srp.sigma_qg / (2.0 * srp.gamma_g * srp.lambda_g * sf_2 );
    break;
  }

  for (unsigned int f = 0; f < gradient_frames.size(); f++)
  {
    const unsigned int ni = weights[f].ni();
    const unsigned int nj = weights[f].nj();

    vil_image_view<double> l_u;
    warps[f].apply_A(u, l_u);
    vil_image_view<double> gradient_lu;
    vidtk::forward_gradient(l_u, gradient_lu);

    for (unsigned int j = 0; j < nj; j++)
    {
      for (unsigned int i = 0; i < ni; i++)
      {
        double val = srp.sigma_qg * sf_2 * weights[f](i,j);
        for (unsigned int k = 0; k < qg[f].nplanes(); k++)
        {
          double &qfijk = qg[f](i, j, k);
          double diff = gradient_lu(i, j, k) - gradient_frames[f](i, j, k);
          switch( srp.cost_function )
          {
          case super3d::super_res_params::HUBER_NORM:
            break;
          case super3d::super_res_params::TRUNCATED_QUADRATIC:
            diff = super3d::rho_truncated_quadratic(diff, srp.alpha_g, srp.gamma_g);
            break;
          case super3d::super_res_params::GENERALIZED_HUBER:
            diff = super3d::rho_generalized_huber(diff, srp.alpha_g, srp.beta_g, srp.gamma_g );
            break;
          }

          qfijk = (qfijk + diff * val)/denom;
          qfijk = vcl_max(qfijk, -sf_2);
          qfijk = vcl_min(qfijk, sf_2);
        }
      }
    }
  }
}

//*****************************************************************************
void primal_step_A(
  const vcl_vector<vil_image_view<double> > &frames,
  const vcl_vector<vil_image_view<double> > &qa,
  const vcl_vector< vil_image_view<double> >& pl,
  vcl_vector<vil_image_view<double> > &As,
  const vcl_vector<vidtk::adjoint_image_ops_func<double> > &warps,
  const vcl_vector<vil_image_view<double> > &weights,
  const super_res_params &srp)
{
  if( !srp.illumination_prior )
    return;

  const double sf_2 = 1.0 / (srp.scale_factor * srp.scale_factor);
  const double scale = srp.tau;

  for (unsigned int i = 0; i < qa.size(); i++)
  {
    if( i == srp.ref_frame )
      continue;

    unsigned int i2=i*2;
    vil_image_view<double> super_qa(qa[i].ni(), qa[i].nj(), qa[i].nplanes());
    vil_image_view<double> super_qay(qa[i].ni(), qa[i].nj(), qa[i].nplanes());

    vil_math_image_product( qa[i], weights[i], super_qa );
    // TODO: the following line is for intensity scaling
    // vil_math_image_product( super_qa, frames[i], super_qay );

    vil_image_view<double> work;
    vidtk::backward_divergence(pl[i2], work);
    vil_math_add_image_fraction(work, scale, super_qa, scale * sf_2);
    vil_image_view<double>& A0 = As[i2];
    vil_math_image_sum(A0, work, A0);
    // TODO: the following lines are for intensity scaling
    // vil_math_image_sum(A0, work, work);
    // vil_math_add_image_fraction( A0, -1.0, work, 2.0);
    // vidtk::backward_divergence(pl[i*2+1], work);
    // vil_math_add_image_fraction(work, scale, super_qay, scale * sf_2);
    // vil_image_view<double>& A1 = As[2*i+1];
    // vil_math_image_sum(A1, work, work);
    // vil_math_add_image_fraction( A1, -1.0, work, 2.0);
  }
}

//*****************************************************************************
void primal_step_Y(
  vil_image_view<double> &u,
  const vcl_vector<vil_image_view<double> > &qa,
  const vcl_vector<vil_image_view<double> > &qg,
  const vil_image_view<double> &pr,
  const vcl_vector<vidtk::adjoint_image_ops_func<double> > &warps,
  const vcl_vector<vil_image_view<double> > &weights,
  const super3d::super_res_params &srp)
{
  const double sf_2 = 1.0 / (srp.scale_factor * srp.scale_factor);

  vil_image_view<double> sum_super_qa(srp.s_ni, srp.s_nj, u.nplanes());
  sum_super_qa.fill(0.0);

  vil_image_view<double> super_qa(srp.s_ni, srp.s_nj, 1);

  unsigned number_of_frames=0;
  unsigned start_frame=0;
  if( srp.image_data_N )
  {
    number_of_frames = qa.size();
  }
  else if( srp.image_data_1 )
  {
    number_of_frames = 1;
    start_frame = srp.ref_frame;
  }

  for (unsigned int i = start_frame; i < start_frame+number_of_frames; i++)
  {
    vil_image_view<double> weighted_qa;
    vil_math_image_product( qa[i], weights[i], weighted_qa );
    warps[i].apply_At(weighted_qa, super_qa);
    vil_math_image_sum(sum_super_qa, super_qa, sum_super_qa);
  }

  if( srp.gradient_data )
  {
    vil_image_view<double> sum_super_qg(srp.s_ni, srp.s_nj, u.nplanes());
    sum_super_qg.fill(0.0);
    vil_image_view<double> super_qg( srp.s_ni, srp.s_nj, u.nplanes() );
    for (unsigned int i = 0; i < qg.size(); i++)
    {
      vil_image_view<double> div;
      vidtk::backward_divergence(qg[i], div);

      vil_image_view<double> weighted_qg_div;
      vil_math_image_product( div, weights[i], weighted_qg_div );
      warps[i].apply_At( weighted_qg_div, super_qg);

      vil_math_image_sum(sum_super_qg, super_qg, sum_super_qg);
    }

    vil_math_add_image_fraction(sum_super_qa, 1.0, sum_super_qg, -1.0);
  }

  vil_image_view<double> work;
  vidtk::backward_divergence(pr, work);
  vil_math_add_image_fraction(work, srp.tau, sum_super_qa, -srp.tau * sf_2);

  vil_math_image_sum(u, work, u);
}

} // end anonymous namespace


//*****************************************************************************

void super_resolve_robust(
  const vcl_vector<vil_image_view<double> > &frames,
  const vcl_vector<vidtk::adjoint_image_ops_func<double> > &warps,
  vil_image_view<double> &Y,
  super_res_params srp,
  unsigned int iterations,
  vcl_vector< vil_image_view<double> > &As,
  const vcl_string &output_image,
  super_res_monitor *srm)
{
  if (frames.empty())
    return;

  if( srp.tv_method == super3d::super_res_params::SUPER3D_BASELINE ||
      srp.tv_method == super3d::super_res_params::IMAGEDATA_IMAGEPRIOR_ILLUMINATIONPRIOR )
    return super_resolve( frames, warps, Y, srp, iterations, output_image, srm );

  unsigned int np = frames[0].nplanes();
  for (unsigned int i = 1; i < frames.size(); i++)
  {
    if (frames[i].nplanes() != np)
    {
      vcl_cerr << "All frames must have the same number of planes.\n";
      return;
    }
  }

  vil_image_view<double>& u = Y;
  if (u.ni() != srp.s_ni || u.nj() != srp.s_nj)
  {
    u.set_size(srp.s_ni, srp.s_nj, np);
    u.fill(0.0);
  }

  srp.image_data_1=false;
  srp.image_data_N=true;
  srp.gradient_data=false;
  srp.image_prior=true;
  srp.illumination_prior=false;
  srp.median_residue=false;

  switch( srp.tv_method )
  {
  case super3d::super_res_params::IMAGEDATA_IMAGEPRIOR:
    srp.image_data_N=true;
    srp.image_prior=true;
    break;
  case super3d::super_res_params::MEDIANDATA_IMAGEPRIOR:
    srp.image_data_N=true;
    srp.image_prior=true;
    srp.median_residue=true;
    break;
  case super3d::super_res_params::GRADIENTDATA_IMAGEPRIOR:
    srp.image_data_1=true;
    srp.image_data_N=false;
    srp.gradient_data=true;
    srp.image_prior=true;
    break;
  case super3d::super_res_params::IMAGEDATA_IMAGEPRIOR_ILLUMINATIONPRIOR:
    srp.image_data_N=true;
    srp.image_prior=true;
    srp.illumination_prior=true;
    break;
  case super3d::super_res_params::IMAGEDATA_GRADIENTDATA_IMAGEPRIOR_ILLUMINATIONPRIOR:
    srp.image_data_N=true;
    srp.gradient_data=true;
    srp.image_prior=true;
    srp.illumination_prior=true;
    break;
  default:
    vcl_cerr << "unknown tv method.\n";
    return;
  }

  vil_image_view<double> pr(srp.s_ni, srp.s_nj, 2*np);
  pr.fill(0.0);

  vcl_vector<vil_image_view<double> > qa(frames.size());
  vcl_vector<vil_image_view<double> > hat_qa(frames.size());
  vcl_vector<vil_image_view<double> > qg(frames.size());
  vcl_vector<vil_image_view<double> > pl(frames.size()*2);
  vcl_vector<vil_image_view<double> > weights;

  vil_image_view<double> Yref_median;
  vcl_vector<vil_image_view<double> > frames_median;
  vil_image_view<double> Y_baseline;
  vcl_vector<vil_image_view<double> > frames_residue;

  for (unsigned int i = 0; i < frames.size(); i++)
  {
    unsigned int low_ni = warps[i].dst_ni();
    unsigned int low_nj = warps[i].dst_nj();
    vcl_cout << low_ni << " " << low_nj << "\n";

    qa[i].set_size(low_ni, low_nj, np);
    qa[i].fill(0.0);

    if( srp.gradient_data )
    {
      qg[i].set_size(low_ni, low_nj, 2*np);
      qg[i].fill(0.0);
    }

    vil_image_view<double> erode, wmap;
    vil_structuring_element erosion_element;
    erosion_element.set_to_disk( srp.erosion_radius );
    if( warps[i].weight_map().nplanes() == 1 )
    {
      wmap = warps[i].weight_map();
    }
    else
    {
      vil_math_mean_over_planes(warps[i].weight_map(), wmap);
    }
    vil_greyscale_erode( wmap, erode, erosion_element );
    weights.push_back( erode );

    if( srp.illumination_prior )
    {
      hat_qa[i].set_size(low_ni, low_nj, np);
      hat_qa[i].fill(0.0);

      unsigned int j=2*i;
      pl[j].set_size(low_ni, low_nj, 2*np);
      pl[j].fill(0.0);
      ++j;
      pl[j].set_size(low_ni, low_nj, 2*np);
      pl[j].fill(0.0);

      vil_image_view<double> A0( low_ni, low_nj, np );
      A0.fill(0.0);
      As.push_back( A0 );

      vil_image_view<double> A1( low_ni, low_nj, np );
      A1.fill(1.0);
      As.push_back( A1 );

      vil_structuring_element median_element;
      median_element.set_to_disk( srp.median_radius);

      vil_image_view<double> frame_median( frames[i].ni(), frames[i].nj(), frames[i].nplanes() );
      for(unsigned int k=0; k<frames[i].nplanes(); k++)
      {
        vil_image_view<double> frame_plane = vil_plane( frames[i], k );
        vil_image_view<double> frame_median_plane = vil_plane( frame_median, k );
        vil_median( frame_plane, frame_median_plane, median_element );
      }
      frames_median.push_back( frame_median );

      if( i == srp.ref_frame )
      {
        super3d::upsample( frame_median, Yref_median, srp.scale_factor, vidtk::warp_image_parameters::CUBIC );
      }

      if( srp.debug )
      {
        char buf[50];
        sprintf(buf, "images/frame_median_%03d.png", i);
        vil_image_view<vxl_byte> output;
        vil_convert_stretch_range_limited( frames_median[i], output, 0.0, 1.0);
        vil_save( output, buf );
      }
    }

    if( srp.median_residue )
    {
      vil_structuring_element median_element;
      median_element.set_to_disk( srp.median_radius);

      vil_image_view<double> frame_median( frames[i].ni(), frames[i].nj(), frames[i].nplanes() );
      for(unsigned int k=0; k<frames[i].nplanes(); k++)
      {
        vil_image_view<double> frame_plane = vil_plane( frames[i], k );
        vil_image_view<double> frame_median_plane = vil_plane( frame_median, k );
        vil_median( frame_plane, frame_median_plane, median_element );
      }

      if( i == srp.ref_frame )
      {
        super3d::upsample( frame_median, Y_baseline, srp.scale_factor, vidtk::warp_image_parameters::CUBIC );
      }

      vil_math_image_difference( frames[i], frame_median, frame_median );
      frames_residue.push_back( frame_median );

      if( srp.debug )
      {
        char buf[50];
        sprintf(buf, "images/frame_residue_%03d.png", i);
        vil_image_view<vxl_byte> output;
        vil_convert_stretch_range_limited( frames_residue[i], output, 0.0, 1.0);
        vil_save( output, buf );
      }
    }
  }

  if( srp.debug )
  {
    vil_image_view<vxl_byte> output;
    if( srp.illumination_prior )
    {
      vil_convert_stretch_range_limited(Yref_median, output, 0.0, 1.0);
      vil_save( output, "images/Yref_median.png" );
    }
    if( srp.median_residue )
    {
      vil_convert_stretch_range_limited( Y_baseline, output, 0.0, 1.0);
      vil_save( output, "images/Y_baseline.png");
    }
  }

  double ssd;
  vil_image_view<double> last;
  last.deep_copy(u);
  unsigned int i = 0;

  vcl_vector< vil_image_view<double> > frames_gradient;
  if( srp.gradient_data )
  {
    for( unsigned int i=0; i<frames.size(); i++ )
    {
      vil_image_view<double> work;
      vidtk::forward_gradient(frames[i], work);
      frames_gradient.push_back( work );
    }
  }

  do
  {
    vcl_cout << "Iteration: " << i;

    switch( srp.tv_method )
    {
    case super3d::super_res_params::IMAGEDATA_IMAGEPRIOR:
      dual_step_pr(u, pr, srp);
      dual_step_qa(u, frames, qa, As, warps, weights, srp);
      primal_step_Y(u, qa, qg, pr, warps, weights, srp);
      break;

    case super3d::super_res_params::MEDIANDATA_IMAGEPRIOR:
      dual_step_pr(u, pr, srp);
      dual_step_qa(u, frames_residue, qa, As, warps, weights, srp);
      primal_step_Y(u, qa, qg, pr, warps, weights, srp);
      break;

    case super3d::super_res_params::GRADIENTDATA_IMAGEPRIOR:
      dual_step_pr(u, pr, srp);
      dual_step_qa(u, frames, qa, As, warps, weights, srp);
      dual_step_qg(u, frames_gradient, qg, warps, weights, srp);
      primal_step_Y(u, qa, qg, pr, warps, weights, srp);
      break;

    case super3d::super_res_params::IMAGEDATA_IMAGEPRIOR_ILLUMINATIONPRIOR:
      dual_step_pr(u, pr, srp);
      dual_step_pl(As, pl, srp);
      dual_step_qa(u, frames, qa, As, warps, weights, srp);
      dual_step_qa(Yref_median, frames_median, hat_qa, As, warps, weights, srp);
      primal_step_A(frames, hat_qa, pl, As, warps, weights, srp );
      primal_step_Y(u, qa, qg, pr, warps, weights, srp);
      break;

    case super3d::super_res_params::IMAGEDATA_GRADIENTDATA_IMAGEPRIOR_ILLUMINATIONPRIOR:
    default:
      vcl_cerr << "unknown tv method.\n";
      return;
    }

    if (srm)
    {
      if (*srm->interrupted_)
      {
        return;
      }

      if (srm->callback_ && !(i % srm->interval_))
      {
        super_res_monitor::update_data data;
        
        if( srp.median_residue )
          vil_math_image_sum( u, Y_baseline, data.current_result );
        else
          data.current_result.deep_copy(u);
        data.num_iterations = i;
        srm->callback_(data);
      }
    }
    else
    {
      double minv, maxv;
      vil_math_value_range(u, minv, maxv);
      if (!(i % srp.frame_step))
      {
        vil_image_view<double> outd;
        vil_image_view<vxl_byte> output;
        if( srp.median_residue )
        {
          vil_image_view<double> Yall(u.ni(), u.nj(), u.nplanes());
          vil_math_image_sum( Y_baseline, u, Yall );
          vil_convert_stretch_range_limited(Yall, outd, vcl_max(0.0,minv), vcl_min(1.0,maxv), 0.0, 255.0);
        }
        else
        {
          vil_convert_stretch_range_limited(u, outd, vcl_max(0.0,minv), vcl_min(1.0,maxv), 0.0, 255.0);
        }
        vil_convert_cast(outd, output);
        vil_save(output, output_image.c_str());

        if( srp.debug )
        {
          char buf[50];
          for (unsigned int f = 0; f < warps.size(); f++)
          {
            sprintf(buf,"images/frame_warp_%03d.png",f);
            vil_image_view<double> l_u;
            warps[f].apply_A(u, l_u);
            vil_convert_stretch_range_limited(l_u, output, 0.0, 1.0);
            vil_save( output, buf );

            if( srp.illumination_prior )
            {
              vil_convert_stretch_range_limited(As[2*f], output, -1.0, 1.0);
              sprintf(buf, "images/A0-%03d.png", f);
              vil_save( output, buf );

              // TODO: the following lines are for intensity scaling
              // vil_convert_stretch_range_limited(As[2*f+1], output, 0.0, 1.0);
              // sprintf(buf, "images/A1-%03d.png", f);
              // vil_save( output, buf );

              vil_image_view<double> low_res_frame;
              // TODO: the following lines are for intensity scaling
              // vil_math_image_product( frames[f], As[2*f+1], low_res_frame);
              // vil_math_image_sum( low_res_frame, As[2*f], low_res_frame);
              vil_math_image_sum( frames[f], As[2*f], low_res_frame);
              sprintf(buf, "images/lowres_%03d.png", f );
              vil_convert_stretch_range_limited(low_res_frame, output, 0.0, 1.0);
              vil_save( output, buf );

              sprintf(buf, "images/qa_%03d.png", f );
              vil_convert_stretch_range_limited(qa[f], output, -0.5, 0.5);
              vil_save( output, buf );
            }

          }
        }
      }

      ssd = vil_math_ssd(last, u, double());
      vcl_cout << " SSD: " << ssd << " " << minv << " " << maxv << "\n";
      vcl_swap(last, u);
    }

  } while (++i < iterations);

  if( srp.median_residue )
  {
    vil_math_image_sum( u, Y_baseline, u );
  }
}

void read_super_res_params( const boost::scoped_ptr<config>& cfg,
                            super_res_params &srp)
{
    if (cfg->is_set("debug"))
      srp.debug = cfg->get_value<bool>("debug");
    else
      srp.debug = false;

    //Initilize super resolution parameters
    srp.lambda = cfg->get_value<double>("lambda");
    srp.epsilon_data = cfg->get_value<double>("epsilon_data");
    srp.epsilon_reg = cfg->get_value<double>("epsilon_reg");
    srp.sigma = cfg->get_value<double>("sigma");
    srp.tau = cfg->get_value<double>("tau");

    // additional parameters
    vcl_string tv_method_str = cfg->get_value<vcl_string>("tv_method");
    if( tv_method_str.compare("SUPER3D_BASELINE") == 0 )
    {
      srp.tv_method = super3d::super_res_params::SUPER3D_BASELINE;
    }
    else if ( tv_method_str.compare("IMAGEDATA_IMAGEPRIOR") == 0 )
    {
      srp.tv_method = super3d::super_res_params::IMAGEDATA_IMAGEPRIOR;
    }
    else if ( tv_method_str.compare("GRADIENTDATA_IMAGEPRIOR") == 0 )
    {
      srp.tv_method = super3d::super_res_params::GRADIENTDATA_IMAGEPRIOR;
    }
    else if ( tv_method_str.compare("IMAGEDATA_IMAGEPRIOR_ILLUMINATIONPRIOR") == 0 )
    {
      srp.tv_method = super3d::super_res_params::IMAGEDATA_IMAGEPRIOR_ILLUMINATIONPRIOR;
    }
    else if ( tv_method_str.compare("MEDIANDATA_IMAGEPRIOR") == 0 )
    {
      srp.tv_method = super3d::super_res_params::MEDIANDATA_IMAGEPRIOR;
    }
    else if ( tv_method_str.compare("IMAGEDATA_GRADIENTDATA_IMAGEPRIOR_ILLUMINATIONPRIOR") == 0 )
    {
      srp.tv_method = super3d::super_res_params::IMAGEDATA_GRADIENTDATA_IMAGEPRIOR_ILLUMINATIONPRIOR;
    }
    else
    {
      vcl_cerr << "unknown tv method\n";
      return;
    }
    vcl_cout << "tv method : " << tv_method_str << vcl_endl;

    vcl_string str = cfg->get_value<vcl_string>("cost_function");
    if( str.compare("HUBER_NORM") == 0 )
    {
      srp.cost_function = super3d::super_res_params::HUBER_NORM;
    }
    else if( str.compare("TRUNCATED_QUADRATIC") == 0 )
    {
      srp.cost_function = super3d::super_res_params::TRUNCATED_QUADRATIC;
    }
    else if( str.compare("GENERALIZED_HUBER") == 0 )
    {
      srp.cost_function = super3d::super_res_params::GENERALIZED_HUBER;
    }
    else
    {
      vcl_cerr << "unknown cost function\n";
      return;
    }
    vcl_cout << "cost function : " << str << vcl_endl;

    if (cfg->is_set("alpha_a"))
      srp.alpha_a = cfg->get_value<double>("alpha_a");
    if (cfg->is_set("gamma_a"))
      srp.gamma_a = cfg->get_value<double>("gamma_a");
    if (cfg->is_set("beta_a"))
      srp.beta_a = cfg->get_value<double>("beta_a");

    if (cfg->is_set("lambda_g"))
      srp.lambda_g = cfg->get_value<double>("lambda_g");
    if (cfg->is_set("alpha_g"))
      srp.alpha_g = cfg->get_value<double>("alpha_g");
    if (cfg->is_set("gamma_g"))
      srp.gamma_g = cfg->get_value<double>("gamma_g");
    if (cfg->is_set("beta_g"))
      srp.beta_g = cfg->get_value<double>("beta_g");

    if (cfg->is_set("lambda_r"))
      srp.lambda_r = cfg->get_value<double>("lambda_r");
    if (cfg->is_set("alpha_r"))
      srp.alpha_r = cfg->get_value<double>("alpha_r");
    if (cfg->is_set("gamma_r"))
      srp.gamma_r = cfg->get_value<double>("gamma_r");
    if (cfg->is_set("beta_r"))
      srp.beta_r = cfg->get_value<double>("beta_r");

    if (cfg->is_set("lambda_l"))
      srp.lambda_l = cfg->get_value<double>("lambda_l");
    if (cfg->is_set("alpha_l"))
      srp.alpha_l = cfg->get_value<double>("alpha_l");
    if (cfg->is_set("gamma_l"))
      srp.gamma_l = cfg->get_value<double>("gamma_l");
    if (cfg->is_set("beta_l"))
      srp.beta_l = cfg->get_value<double>("beta_l");

    if (cfg->is_set("sigma_pr"))
      srp.sigma_pr = cfg->get_value<double>("sigma_pr");
    if (cfg->is_set("sigma_pl"))
      srp.sigma_pl = cfg->get_value<double>("sigma_pl");
    if (cfg->is_set("sigma_qa"))
      srp.sigma_qa = cfg->get_value<double>("sigma_qa");
    if (cfg->is_set("sigma_qg"))
      srp.sigma_qg = cfg->get_value<double>("sigma_qg");
    if (cfg->is_set("sigma_A"))
      srp.sigma_A = cfg->get_value<double>("sigma_A");
    if (cfg->is_set("sigma_Y"))
      srp.sigma_Y = cfg->get_value<double>("sigma_Y");

    if (cfg->is_set("erosion_radius"))
    {
      srp.erosion_radius = cfg->get_value<double>("erosion_radius");
    }
    else
    {
      srp.erosion_radius = 3.0;
    }

    if (cfg->is_set("median_radius"))
    {
      srp.median_radius = cfg->get_value<double>("median_radius");
    }
    else
    {
      srp.median_radius = 3.0;
    }

    if (cfg->is_set("frame_step"))
    {
      srp.frame_step = cfg->get_value<int>("frame_step");
    }
    else
    {
      srp.frame_step = 15;
    }

}

//*****************************************************************************
} // end namespace super3d

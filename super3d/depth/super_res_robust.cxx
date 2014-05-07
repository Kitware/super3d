/*ckwg +5
 * Copyright 2014 by Kitware, Inc. All Rights Reserved. Please refer to
 * KITWARE_LICENSE.TXT for licensing information, or contact General Counsel,
 * Kitware, Inc., 28 Corporate Drive, Clifton Park, NY 12065.
 */

#include "super_res.h"
#include "super_config.h"
#include <video_transforms/adjoint_image_derivs.h>

#include <vil/algo/vil_gauss_filter.h>
#include <vil/vil_save.h>
#include <vil/vil_math.h>
#include <vil/vil_resample_bilin.h>
#include <vil/vil_resample_bicub.h>
#include <vil/vil_convert.h>

#include <vil/vil_resample_bicub.txx>
VIL_RESAMPLE_BICUB_INSTANTIATE( double , double );

#define DEBUG

//*****************************************************************************

void dual_step_pr(const vil_image_view<double> &u_bar,
                  vil_image_view<double> &pr,
                  const super_res_params &srp)
{
  vil_image_view<double> work;
  vidtk::forward_gradient(u_bar, work);

  double denom;
  switch( srp.cost_function )
  {
  case super_res_params::HUBER_NORM:
    denom = 1.0 + ( srp.sigma_pr * srp.alpha_r) / srp.lambda_r;
    break;
  case super_res_params::TRUNCATED_QUADRATIC:
    denom = 1.0 + srp.sigma_pr / (2.0 * srp.gamma_r * srp.lambda_r );
//    vil_transform( work, rho_truncated_quadratic_functor(srp.alpha_r, srp.gamma_r) );
    break;
  case super_res_params::GENERALIZED_HUBER:
    denom = 1.0 + srp.sigma_pr / (2.0 * srp.gamma_r * srp.lambda_r );
//    vil_transform( work, rho_generalized_huber_functor(srp.alpha_r, srp.beta_r, srp.gamma_r) );
    break;
  }

//  vil_math_add_image_fraction(pr, 1.0/denom, work, srp.sigma_pr/denom);

  for (unsigned int j = 0; j < srp.s_nj; j++)
  {
    for (unsigned int i = 0; i < srp.s_ni; i++)
    {
      for (unsigned int k = 0; k < u_bar.nplanes(); k++)
      {
        const unsigned k2=2*k;
        const unsigned k1=k2+1;
        double work1 = work(i,j,k1);
        double work2 = work(i,j,k2);

        switch( srp.cost_function )
        {
        case super_res_params::TRUNCATED_QUADRATIC:
          work1 = rho_truncated_quadratic( work1, srp.alpha_r, srp.gamma_r );
          work2 = rho_truncated_quadratic( work2, srp.alpha_r, srp.gamma_r );
          break;
        case super_res_params::GENERALIZED_HUBER:
          work1 = rho_generalized_huber( work1, srp.alpha_r, srp.beta_r, srp.gamma_r );
          work2 = rho_generalized_huber( work2, srp.alpha_r, srp.beta_r, srp.gamma_r );
          break;
        }

        double &x = pr(i,j,k2);
        double &y = pr(i,j,k1);
        x = ( x + srp.sigma_pr * work2) / denom;
        y = ( y + srp.sigma_pr * work1) / denom;

        const double mag = sqrt(x*x + y*y)/srp.lambda_r;
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
void dual_step_pl(const vil_image_view<double> &u_bar,
                  vil_image_view<double> &pl,
                  const super_res_params &srp)
{
  vil_image_view<double> work;
  vidtk::forward_gradient(u_bar, work);  // fix me

  double denom;
  switch( srp.cost_function )
  {
  case super_res_params::HUBER_NORM:
    denom = 1.0 + ( srp.sigma_pl * srp.alpha_l) / srp.lambda_l;
    break;
  case super_res_params::TRUNCATED_QUADRATIC:
    denom = 1.0 + srp.sigma_pl / (2.0 * srp.gamma_l * srp.lambda_l );
//    vil_transform( work, rho_truncated_quadratic_functor(srp.alpha_l, srp.gamma_l) );
    break;
  case super_res_params::GENERALIZED_HUBER:
    denom = 1.0 + srp.sigma_pl / (2.0 * srp.gamma_l * srp.lambda_l );
//    vil_transform( work, rho_generalized_huber_functor(srp.alpha_l, srp.beta_l, srp.gamma_l) );
    break;
  }

//  vil_math_add_image_fraction(pl, 1.0/denom, work, srp.sigma_pl/denom);

  for (unsigned int j = 0; j < srp.s_nj; j++)
  {
    for (unsigned int i = 0; i < srp.s_ni; i++)
    {
      for (unsigned int k = 0; k < u_bar.nplanes(); k++)
      {
        const unsigned k2=2*k;
        const unsigned k1=k2+1;
        double work1 = work(i,j,k1);
        double work2 = work(i,j,k2);

        switch( srp.cost_function )
        {
        case super_res_params::TRUNCATED_QUADRATIC:
          work1 = rho_truncated_quadratic( work1, srp.alpha_l, srp.gamma_l );
          work2 = rho_truncated_quadratic( work2, srp.alpha_l, srp.gamma_l );
          break;
        case super_res_params::GENERALIZED_HUBER:
          work1 = rho_generalized_huber( work1, srp.alpha_l, srp.beta_l, srp.gamma_l );
          work2 = rho_generalized_huber( work2, srp.alpha_l, srp.beta_l, srp.gamma_l );
          break;
        }

        double &x = pl(i,j,k2);
        double &y = pl(i,j,k1);
        x = ( x + srp.sigma_pl * work2) / denom;
        y = ( y + srp.sigma_pl * work1) / denom;

        const double mag = sqrt(x*x + y*y)/srp.lambda_l;
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

void dual_step_qa(const vcl_vector<vil_image_view<double> > &frames,
                  const vcl_vector<vidtk::adjoint_image_ops_func<double> > &warps,
                  const vcl_vector<vil_image_view<double> > &weights,
                  const vil_image_view<double> &u_bar,
                  vcl_vector<vil_image_view<double> > &qa,
                  const super_res_params &srp)
{
  vil_image_view<double> l_u, work, dot_mask;
  const double sf_2 = 1.0 / (srp.scale_factor * srp.scale_factor);

  double denom;
  switch( srp.cost_function )
  {
  case super_res_params::HUBER_NORM:
    denom = 1.0 + ( srp.sigma_qa * srp.alpha_a) / sf_2;
    break;
  case super_res_params::TRUNCATED_QUADRATIC:
  case super_res_params::GENERALIZED_HUBER:
    denom = 1.0 + srp.sigma_qa / (2.0 * srp.gamma_a * sf_2 );
    break;
  }

  unsigned number_of_frames=0;
  if( srp.image_data_N )
  {
    number_of_frames = frames.size();
  }
  else if( srp.image_data_1 )
  {
    number_of_frames = 1;
  }

  for (unsigned int f = 0; f < number_of_frames; f++)
  {
    const unsigned int ni = weights[f].ni();
    const unsigned int nj = weights[f].nj();

    // apply the linear operator to warp, blur, and downsample
    warps[f].apply_A(u_bar, l_u);

//    vil_math_image_difference( l_u, frames[f], work );

//    vcl_vector<double> parameters;
//    switch( srp.cost_function )
//    {
//    case super_res_params::TRUNCATED_QUADRATIC:
//      vil_transform( work, rho_truncated_quadratic_functor(srp.alpha_a, srp.gamma_a) );
//      break;
//    case super_res_params::GENERALIZED_HUBER:
//      vil_transform( work, rho_generalized_huber_functor(srp.alpha_a, srp.beta_a, srp.gamma_a) );
//      break;
//    }

//    vil_math_image_product( work, weights[f], dot_mask );
//    vil_math_add_image_fraction(qa[f], 1.0/denom, dot_mask, srp.sigma_qa * sf_2 / denom);
//    vil_math_truncate_range( qa[f], -sf_2, sf_2 );

    for (unsigned int j = 0; j < nj; j++)
    {
      for (unsigned int i = 0; i < ni; i++)
      {
        for (unsigned int k = 0; k < qa[f].nplanes(); k++)
        {
          double &qfijk = qa[f](i, j, k);

          double diff = l_u(i, j, k) - frames[f](i, j, k);
          switch( srp.cost_function )
          {
          case super_res_params::TRUNCATED_QUADRATIC:
            diff = rho_truncated_quadratic(diff, srp.alpha_a, srp.gamma_a);
            break;
          case super_res_params::GENERALIZED_HUBER:
            diff = rho_generalized_huber(diff, srp.alpha_a, srp.beta_a,  srp.gamma_a);
            break;
          }

          qfijk = (qfijk + srp.sigma_qa * sf_2 * diff * weights[f](i,j))/denom;
          qfijk = vcl_max(qfijk, -sf_2);
          qfijk = vcl_min(qfijk, sf_2);
        }
      }
    }
  }

}

//*****************************************************************************

void dual_step_qg(const vcl_vector<vil_image_view<double> > &gradient_frames,
                  const vcl_vector<vidtk::adjoint_image_ops_func<double> > &warps,
                  const vcl_vector<vil_image_view<double> > &weights,
                  const vil_image_view<double> &u_bar,
                  vcl_vector<vil_image_view<double> > &qg,
                  const super_res_params &srp)
{
  vil_image_view<double> l_u, work, dot_mask;
  const double sf_2 = 1.0 / (srp.scale_factor * srp.scale_factor);

  double denom;
  switch( srp.cost_function )
  {
  case super_res_params::HUBER_NORM:
    denom = 1.0 + ( srp.sigma_qg * srp.alpha_g) / ( srp.lambda_g * sf_2);
    break;
  case super_res_params::TRUNCATED_QUADRATIC:
  case super_res_params::GENERALIZED_HUBER:
    denom = 1.0 + srp.sigma_qg / (2.0 * srp.gamma_g * srp.lambda_g * sf_2 );
    break;
  }

  for (unsigned int f = 0; f < gradient_frames.size(); f++)
  {
    const unsigned int ni = weights[f].ni();
    const unsigned int nj = weights[f].nj();

    // apply the linear operator to warp, blur, and downsample
    warps[f].apply_A(u_bar, l_u);
    vil_image_view<double> gradient_lu;
    vidtk::forward_gradient(l_u, gradient_lu);

//    vil_math_image_difference( gradient_lu, gradient_frames[f], work );
//    vil_math_image_product( work, weights[f], dot_mask );
//    vil_math_add_image_fraction(qg[f], 1.0/denom, dot_mask, srp.sigma_qg * sf_2 / denom);
//    vil_math_truncate_range( qg[f], -sf_2, sf_2 );

    for (unsigned int j = 0; j < nj; j++)
    {
      for (unsigned int i = 0; i < ni; i++)
      {
        for (unsigned int k = 0; k < qg[f].nplanes(); k++)
        {
          double &qfijk = qg[f](i, j, k);
          const double &w = weights[f](i,j);
          double diff = gradient_lu(i, j, k) - gradient_frames[f](i, j, k);
          switch( srp.cost_function )
          {
          case super_res_params::TRUNCATED_QUADRATIC:
            diff = rho_truncated_quadratic(diff, srp.alpha_g, srp.gamma_g);
            break;
          case super_res_params::GENERALIZED_HUBER:
            diff = rho_generalized_huber(diff, srp.alpha_g, srp.beta_g, srp.gamma_g );
            break;
          }

          qfijk = (qfijk + srp.sigma_qg * sf_2 * diff * weights[f](i,j))/denom;
          qfijk = vcl_max(qfijk, -sf_2);
          qfijk = vcl_min(qfijk, sf_2);
        }
      }
    }
  }
}

//*****************************************************************************

void primal_step_Y(const vcl_vector<vil_image_view<double> > &qa,
                   const vcl_vector<vil_image_view<double> > &qg,
                   const vcl_vector<vidtk::adjoint_image_ops_func<double> > &warps,
                   const vil_image_view<double> &pr,
                   vil_image_view<double> &u,
                   vil_image_view<double> &u_bar,
                   const super_res_params &srp)
{
  vil_image_view<double> work, div;
  const double sf_2 = 1.0 / (srp.scale_factor * srp.scale_factor);

  vil_image_view<double> sum_super_qa(srp.s_ni, srp.s_nj, u.nplanes());
  sum_super_qa.fill(0.0);

  vil_image_view<double> super_qa(srp.s_ni, srp.s_nj, u.nplanes());

  unsigned number_of_frames=0;
  if( srp.image_data_N )
  {
    number_of_frames = qa.size();
  }
  else if( srp.image_data_1 )
  {
    number_of_frames = 1;
  }

  for (unsigned int i = 0; i < number_of_frames; i++)
  {
    // apply transpose linear operator to upsample, blur, and warp
    warps[i].apply_At(qa[i], super_qa);
    vil_math_image_sum(sum_super_qa, super_qa, sum_super_qa);
  }

  if( srp.gradient_data )
  {
    vil_image_view<double> sum_super_qg(srp.s_ni, srp.s_nj, u.nplanes());
    sum_super_qg.fill(0.0);
    vil_image_view<double> super_qg( srp.s_ni, srp.s_nj, u.nplanes() );
    for (unsigned int i = 0; i < qg.size(); i++)
    {
      // apply transpose linear operator to upsample, blur, and warp
      vidtk::backward_divergence(qg[i], div);
      warps[i].apply_At(div, super_qg);
      vil_math_image_sum(sum_super_qg, super_qg, sum_super_qg);
    }

    vil_math_add_image_fraction(sum_super_qa, 1.0, sum_super_qg, -1.0);
  }

  vidtk::backward_divergence(pr, work);
  vil_math_add_image_fraction(work, srp.tau, sum_super_qa, -srp.tau * sf_2);
  vil_math_image_sum(u, work, work);
  u_bar.deep_copy(work);
  vil_math_add_image_fraction(u_bar, 2.0, u, -1.0);
  vcl_swap(u, work);
}

//*****************************************************************************

void super_resolve_robust(
  const vcl_vector<vil_image_view<double> > &frames,
  const vcl_vector<vidtk::adjoint_image_ops_func<double> > &warps,
  vil_image_view<double> &Y,
  super_res_params srp,
  unsigned int iterations,
  vcl_vector< vil_image_view<double> > &As)
{

  if( srp.tv_method == super_res_params::SUPER3D_BASELINE )
    return super_resolve( frames, warps, Y, srp, iterations );

  vil_image_view<double>& u = Y;

  if (frames.empty())
    return;

  int np = frames[0].nplanes();
  for (unsigned int i = 1; i < frames.size(); i++)
  {
    if (frames[i].nplanes() != np)
    {
      vcl_cerr << "All frames must have the same number of planes.\n";
      return;
    }
  }

  //If the size was correct assume that it was pre-initialized
  if (u.ni() != srp.s_ni || u.nj() != srp.s_nj)
  {
    u.set_size(srp.s_ni, srp.s_nj, np);
    u.fill(0.0);
  }

  vil_image_view<double> pr(srp.s_ni, srp.s_nj, 2*np);
  pr.fill(0.0);

  vcl_vector<vil_image_view<double> > qa(frames.size());
  vcl_vector<vil_image_view<double> > qg(frames.size());
  vcl_vector<vil_image_view<double> > pl(frames.size());
  vcl_vector<vil_image_view<double> > weights;

  for (unsigned int i = 0; i < frames.size(); i++)
  {
    vcl_cout << warps[i].dst_ni() << " " << warps[i].dst_nj() << "\n";

    qa[i].set_size(warps[i].dst_ni(), warps[i].dst_nj(), np);
    qa[i].fill(0.0);

    qg[i].set_size(warps[i].dst_ni(), warps[i].dst_nj(), 2*np);
    qg[i].fill(0.0);

    pl[i].set_size(srp.s_ni, srp.s_nj, 2*np);
    pl[i].fill(0.0);

    weights.push_back(warps[i].weight_map());
  }

  double ssd;
  vil_image_view<double> u_bar, last;
  u_bar.deep_copy(u);
  last.deep_copy(u);
  unsigned int i = 0;
#ifdef DEBUG
  vil_image_view<vxl_byte> output;
#endif

  switch( srp.tv_method )
  {
  case super_res_params::IMAGEDATA_IMAGEPRIOR:
    srp.image_data_1=false;
    srp.image_data_N=true;
    srp.gradient_data=false;
    srp.image_prior=true;
    srp.illumination_prior=false;
    break;
  case super_res_params::GRADIENTDATA_IMAGEPRIOR:
    srp.image_data_1=true;
    srp.image_data_N=false;
    srp.gradient_data=true;
    srp.image_prior=true;
    srp.illumination_prior=false;
    break;
  case super_res_params::IMAGEDATA_IMAGEPRIOR_ILLUMINATIONPRIOR:
    srp.image_data_1=false;
    srp.image_data_N=true;
    srp.gradient_data=false;
    srp.image_prior=true;
    srp.illumination_prior=true;
    break;
  case super_res_params::IMAGEDATA_GRADIENTDATA_IMAGEPRIOR_ILLUMINATIONPRIOR:
    srp.image_data_1=false;
    srp.image_data_N=true;
    srp.gradient_data=true;
    srp.image_prior=true;
    srp.illumination_prior=true;
    break;
  default:
    vcl_cerr << "unknown tv method.\n";
    return;
  }

  vcl_vector< vil_image_view<double> > frames_gradient;
  vil_image_view<double> work;
  if( srp.gradient_data )
  {
    for( unsigned int i=0; i<frames.size(); i++ )
    {
      vidtk::forward_gradient(frames[i], work);
      frames_gradient.push_back( work );
    }
  }

  do
  {
    vcl_cout << "Iteration: " << i;
    double minv, maxv;
    switch( srp.tv_method )
    {
    case super_res_params::IMAGEDATA_IMAGEPRIOR:
      dual_step_pr(u, pr, srp);
      dual_step_qa(frames, warps, weights, u, qa, srp);
      primal_step_Y(qa, qg, warps, pr, u, u_bar, srp);
      break;

    case super_res_params::GRADIENTDATA_IMAGEPRIOR:
      dual_step_pr(u, pr, srp);
//      vil_math_value_range(pr, minv, maxv);
//      vcl_cout << "  pr minv = " << minv << " maxv = " << maxv << "\n";

      dual_step_qa(frames, warps, weights, u, qa, srp);
//      vil_math_value_range(qa[0], minv, maxv);
//      vcl_cout << "  qa[0] minv = " << minv << " maxv = " << maxv << "\n";

      dual_step_qg(frames_gradient, warps, weights, u, qg, srp);
//      vil_math_value_range(qg[0], minv, maxv);
//      vcl_cout << "  qg[0] minv = " << minv << " maxv = " << maxv << "\n";

      primal_step_Y(qa, qg, warps, pr, u, u_bar, srp);
//      vil_math_value_range(u, minv, maxv);
//      vcl_cout << "  u minv = " << minv << " maxv = " << maxv << "\n";
      break;

    case super_res_params::IMAGEDATA_IMAGEPRIOR_ILLUMINATIONPRIOR:
    case super_res_params::IMAGEDATA_GRADIENTDATA_IMAGEPRIOR_ILLUMINATIONPRIOR:
    default:
      vcl_cerr << "unknown tv method.\n";
      return;
    }

    vil_math_value_range(u, minv, maxv);
#ifdef DEBUG
    if (!(i % 15))
    {
      vil_image_view<double> outd;
      vil_convert_stretch_range_limited(u, outd, vcl_max(0.0,minv), vcl_min(1.0,maxv), 0.0, 255.0);
      vil_convert_cast(outd, output);
      vil_save(output, config::inst()->get_value<vcl_string>("output_image").c_str());
    }
#endif

    ssd = vil_math_ssd(last, u, double());
    vcl_cout << " SSD: " << ssd << " " << minv << " " << maxv << "\n";
    vcl_swap(last, u);
  } while (++i < iterations);
}

//*****************************************************************************

inline double rho_huber_norm( double x, double alpha )
{
  double absx = vcl_fabs(x);
  if( absx <= alpha )
    return x*x / 2.0 / alpha;
  else
    return absx - alpha/2.0;
}

inline double psi_huber_norm( double x, double alpha )
{
  if( vcl_fabs(x) <= alpha )
    return 1.0;
  else
    return -1.0;
}

inline double rho_truncated_quadratic( double x, double alpha, double gamma )
{
  if( vcl_fabs(x) <= vcl_sqrt( alpha / gamma ) )
    return gamma * x * x;
  else
    return 0.0;
}

inline double psi_truncated_quadratic( double x, double alpha, double gamma )
{
  if( vcl_fabs(x) <= vcl_sqrt( alpha / gamma ) )
    return 2.0 * gamma * x;
  else
    return 0.0;
}

inline double rho_generalized_huber( double x, double alpha, double beta, double gamma )
{
  double t = vcl_sqrt( alpha/gamma );
  if( vcl_fabs(x) <= t )
    return gamma * x * x;
  else
    return beta * x + alpha - t * beta;
}

inline double psi_generalized_huber( double x, double alpha, double beta, double gamma )
{
  if( vcl_fabs(x) <= vcl_sqrt( alpha/gamma ) )
    return 2.0 * gamma * x;
  else
    return beta;
}

//*****************************************************************************

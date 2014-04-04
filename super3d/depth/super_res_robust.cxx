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
    vil_transform( work, rho_truncated_quadratic_functor(srp.alpha_r, srp.gamma_r) );
    break;
  case super_res_params::GENERALIZED_HUBER:
    denom = 1.0;  // fix me
    vil_transform( work, rho_generalized_huber_functor(srp.alpha_r, srp.beta_r, srp.gamma_r) );
    break;
  }

  vil_math_add_image_fraction(pr, 1.0/denom, work, srp.sigma_pr/denom);

  for (unsigned int j = 0; j < srp.s_nj; j++)
  {
    for (unsigned int i = 0; i < srp.s_ni; i++)
    {
      for (unsigned int k = 0; k < u_bar.nplanes(); k++)
      {
        double &x = pr(i,j,2*k), &y = pr(i,j,2*k+1);

        //truncate vectors
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
  vidtk::forward_gradient(u_bar, work);

  double denom;
  switch( srp.cost_function )
  {
  case super_res_params::HUBER_NORM:
    denom = 1.0 + ( srp.sigma_pl * srp.alpha_l) / srp.lambda_l;
    break;
  case super_res_params::TRUNCATED_QUADRATIC:
    denom = 1.0 + srp.sigma_pl / (2.0 * srp.gamma_l * srp.lambda_l );
    vil_transform( work, rho_truncated_quadratic_functor(srp.alpha_l, srp.gamma_l) );
    break;
  case super_res_params::GENERALIZED_HUBER:
    denom = 1.0;  // fix me
    vil_transform( work, rho_generalized_huber_functor(srp.alpha_l, srp.beta_l, srp.gamma_l) );
    break;
  }

  vil_math_add_image_fraction(pl, 1.0/denom, work, srp.sigma_pl/denom);

  for (unsigned int j = 0; j < srp.s_nj; j++)
  {
    for (unsigned int i = 0; i < srp.s_ni; i++)
    {
      for (unsigned int k = 0; k < u_bar.nplanes(); k++)
      {
        double &x = pl(i,j,2*k), &y = pl(i,j,2*k+1);

        //truncate vectors
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
    denom = 1.0 + srp.sigma_qa / (2.0 * srp.gamma_a * sf_2 );
    break;
  case super_res_params::GENERALIZED_HUBER:
    denom = 1.0;  // fix me
    break;
  }

  for (unsigned int f = 0; f < frames.size(); f++)
  {
    // apply the linear operator to warp, blur, and downsample
    warps[f].apply_A(u_bar, l_u);
    vil_math_image_difference( l_u, frames[f], work );

    switch( srp.cost_function )
    {
    case super_res_params::TRUNCATED_QUADRATIC:
      vil_transform( work, rho_truncated_quadratic_functor(srp.alpha_a, srp.gamma_a) );
      break;
    case super_res_params::GENERALIZED_HUBER:
      vil_transform( work, rho_generalized_huber_functor(srp.alpha_a, srp.beta_a, srp.gamma_a) );
      break;
    }

    vil_math_image_product( work, weights[f], dot_mask );

    vil_math_add_image_fraction(qa[f], 1.0/denom, dot_mask, srp.sigma_qa * sf_2 / denom);

    vil_math_truncate_range( qa[f], -sf_2, sf_2 );

//          qfijk = (qfijk + srp.sigma / sf_2 * (l_u(i, j, k) - frames[f](i, j, k) * weights[f](i,j)))/denom;
//          qfijk = vcl_max(qfijk, -sf_2);
//          qfijk = vcl_min(qfijk, sf_2);
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
  vil_image_view<double> l_u, gradient_lu, work, dot_mask;
  const double sf_2 = 1.0 / (srp.scale_factor * srp.scale_factor);

  double denom;
  switch( srp.cost_function )
  {
  case super_res_params::HUBER_NORM:
    denom = 1.0 + ( srp.sigma_qg * srp.alpha_g) / sf_2;
    break;
  case super_res_params::TRUNCATED_QUADRATIC:
    denom = 1.0 + srp.sigma_qg / (2.0 * srp.gamma_g * sf_2 );
    break;
  case super_res_params::GENERALIZED_HUBER:
    denom = 1.0;  // fix me
    break;
  }

  for (unsigned int f = 0; f < gradient_frames.size(); f++)
  {
    // apply the linear operator to warp, blur, and downsample
    warps[f].apply_A(u_bar, l_u);
    vidtk::forward_gradient(l_u, gradient_lu);

    vil_math_image_difference( gradient_lu, gradient_frames[f], work );

    switch( srp.cost_function )
    {
    case super_res_params::TRUNCATED_QUADRATIC:
      vil_transform( work, rho_truncated_quadratic_functor(srp.alpha_g, srp.gamma_g) );
      break;
    case super_res_params::GENERALIZED_HUBER:
      vil_transform( work, rho_generalized_huber_functor(srp.alpha_g, srp.beta_g, srp.gamma_g) );
      break;
    }

    vil_math_image_product( work, weights[f], dot_mask );

    vil_math_add_image_fraction(qg[f], 1.0/denom, dot_mask, srp.sigma_qg * sf_2 / denom);

    vil_math_truncate_range( qg[f], -sf_2, sf_2 );

//          qfijk = (qfijk + srp.sigma / sf_2 * (l_u(i, j, k) - gradient_frames[f](i, j, k) * weights[f](i,j)))/denom;
//          qfijk = vcl_max(qfijk, -sf_2);
//          qfijk = vcl_min(qfijk, sf_2);
  }
}

//*****************************************************************************

void primal_step_u(const vcl_vector<vil_image_view<double> > &q,
                   const vcl_vector<vidtk::adjoint_image_ops_func<double> > &warps,
                   const vil_image_view<double> &p,
                   vil_image_view<double> &u,
                   vil_image_view<double> &u_bar,
                   const super_res_params &srp)
{
  const double sf_2 = 1.0 / (srp.scale_factor * srp.scale_factor);
  vil_image_view<double> sum_super_q(srp.s_ni, srp.s_nj, u.nplanes());
  sum_super_q.fill(0.0);
  vil_image_view<double> super_q(srp.s_ni, srp.s_nj, u.nplanes()), temp;
  for (unsigned int i = 0; i < q.size(); i++)
  {
    // apply transpose linear operator to upsample, blur, and warp
    warps[i].apply_At(q[i], super_q);
    vil_math_image_sum(sum_super_q, super_q, sum_super_q);
  }

  vil_image_view<double> work;
  vidtk::backward_divergence(p, work);
  vil_math_add_image_fraction(work, srp.tau, sum_super_q, -srp.tau * sf_2);
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

  if( srp.tv_method.compare("imagedata_imageprior") == 0 && srp.cost_function == super_res_params::HUBER_NORM )
    return super_resolve( frames, warps, Y, srp, iterations );

  // determine tv_method
  if( srp.tv_method.compare("imagedata_imageprior") == 0 )
  {
    srp.image_data_1=false;
    srp.image_data_N=true;
    srp.gradient_data=false;
    srp.image_prior=true;
    srp.illumination_prior=false;
  }
  else if ( srp.tv_method.compare("gradientdata_imageprior") == 0 )
  {
    srp.image_data_1=true;
    srp.image_data_N=false;
    srp.gradient_data=true;
    srp.image_prior=true;
    srp.illumination_prior=false;
  }
  else if ( srp.tv_method.compare("imagedata_imageprior_illuminationprior") == 0 )
  {
    srp.image_data_1=false;
    srp.image_data_N=true;
    srp.gradient_data=false;
    srp.image_prior=true;
    srp.illumination_prior=true;
  }
  else if ( srp.tv_method.compare("imagedata_radientdata_imageprior_illuminationprior") == 0 )
  {
    srp.image_data_1=false;
    srp.image_data_N=true;
    srp.gradient_data=true;
    srp.image_prior=true;
    srp.illumination_prior=true;
  }
  else
  {
    vcl_cerr << "unknow tv method.\n";
    return;
  }

  // determine cost_function
  switch( srp.cost_function )
  {
  case super_res_params::HUBER_NORM:
    srp.rho = &rho_huber_norm;
    srp.psi = &psi_huber_norm;
    break;
  case super_res_params::TRUNCATED_QUADRATIC:
    srp.rho = &rho_truncated_quadratic;
    srp.psi = &psi_truncated_quadratic;
    break;
  case super_res_params::GENERALIZED_HUBER:
    srp.rho = &rho_generalized_huber;
    srp.psi = &psi_generalized_huber;
    break;
  default:
    vcl_cerr << "unknown cost function\n";
    break;
  }

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

  vil_image_view<double> p(srp.s_ni, srp.s_nj, 2*np);
  p.fill(0.0);

  vcl_vector<vil_image_view<double> > q(frames.size());
  vcl_vector<vil_image_view<double> > weights;
  for (unsigned int i = 0; i < frames.size(); i++)
  {
    vcl_cout << warps[i].dst_ni() << " " << warps[i].dst_nj() << "\n";
    q[i].set_size(warps[i].dst_ni(), warps[i].dst_nj(), np);
    q[i].fill(0.0);
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

  do
  {
    vcl_cout << "Iteration: " << i;
    dual_step_pr(u, p, srp);
    dual_step_qa(frames, warps, weights, u, q, srp);
    primal_step_u(q, warps, p, u, u_bar, srp);

    double minv, maxv;
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

double rho_huber_norm( const vcl_vector<double>& v, double x )
{
  double y;
  double absx = vcl_fabs(x);
  if( absx <= v[0] )
    y = x*x / 2.0 / v[0];
  else
    y = absx - v[0]/2.0;
  return y;
}

double psi_huber_norm( const vcl_vector<double>& v, double x )
{
  if( vcl_fabs(x) <= v[0] )
    return 1.0;
  else
    return -1.0;
}

double rho_truncated_quadratic( const vcl_vector<double>& v, double x )
{
  double y;
  if( vcl_fabs(x) <= vcl_sqrt( v[0]/v[1] ) )
    y = v[1] * x * x;
  else
    y = 0.0;
  return y;
}

double psi_truncated_quadratic( const vcl_vector<double>& v, double x )
{
  double y;
  if( vcl_fabs(x) <= vcl_sqrt( v[0]/v[1] ) )
    y = 2.0 * v[1] * x;
  else
    y = 0.0;
  return y;
}

double rho_generalized_huber( const vcl_vector<double>& v, double x )
{
  double y;
  double t = vcl_sqrt( v[0]/v[1] );

  if( vcl_fabs(x) <= t )
    y = v[1] * x * x;
  else
    y = v[2] * x + v[0] - t * v[2];
  return y;
}

double psi_generalized_huber( const vcl_vector<double>& v, double x )
{
  double y;
  if( vcl_fabs(x) <= vcl_sqrt( v[0]/v[1] ) )
    y = 2.0 * v[1] * x;
  else
    y = v[2];
  return y;
}

//*****************************************************************************

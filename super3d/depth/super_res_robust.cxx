/*ckwg +5
 * Copyright 2014 by Kitware, Inc. All Rights Reserved. Please refer to
 * KITWARE_LICENSE.TXT for licensing information, or contact General Counsel,
 * Kitware, Inc., 28 Corporate Drive, Clifton Park, NY 12065.
 */

#include "super_res.h"
#include "super_config.h"
#include <stdio.h>
#include <video_transforms/adjoint_image_derivs.h>

#include <vil/algo/vil_gauss_filter.h>
#include <vil/vil_save.h>
#include <vil/vil_math.h>
#include <vil/vil_resample_bilin.h>
#include <vil/vil_resample_bicub.h>
#include <vil/vil_convert.h>
#include <vil/algo/vil_greyscale_erode.h>
#include <vil/vil_resample_bicub.txx>
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
    for(unsigned int i=0; i<As.size(); i++)
    {
      dual_step_grad_prior(As[i], pl[i], srp.lambda_l, srp.sigma_pl, srp.alpha_l, srp.gamma_l, srp.cost_function);
    }
  }
}

//*****************************************************************************

void dual_step_qa(const vcl_vector<vil_image_view<double> > &frames,
                  const vcl_vector<vidtk::adjoint_image_ops_func<double> > &warps,
                  const vcl_vector<vil_image_view<double> > &weights,
                  const vil_image_view<double> &u,
                  vcl_vector<vil_image_view<double> > &qa,
                  vcl_vector< vil_image_view<double> > &As,
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

    // apply the linear operator to warp, blur, and downsample
    vil_image_view<double> l_u;
    if( srp.illumination_prior )
    {
      vil_image_view<double> work(u.ni(),u.nj(),u.nplanes());
      vil_image_view<double> A0=As[2*f];
      vil_image_view<double> A1=As[2*f+1];
      double a,b;
      for(unsigned int j=0; j<u.nj(); j++)
      {
        for(unsigned int i=0; i<u.ni(); i++)
        {
          a = A1(i,j);
          b = A0(i,j);
          for(unsigned int k=0; k<u.nplanes(); k++ )
          {
            work(i,j,k) = u(i,j,k) * a + b;
          }
        }
      }
//      vil_math_image_product( u, As[2*f+1], work);
//      vil_math_image_sum( work, As[2*f], work);
      warps[f].apply_A( work, l_u );
    }
    else
    {
      warps[f].apply_A(u, l_u);
    }

    for (unsigned int j = 0; j < nj; j++)
    {
      for (unsigned int i = 0; i < ni; i++)
      {
        double val = srp.sigma_qa * sf_2 * weights[f](i,j);
        for (unsigned int k = 0; k < qa[f].nplanes(); k++)
        {
          double &qfijk = qa[f](i, j, k);

          double diff = l_u(i, j, k) - frames[f](i, j, k);
          switch( srp.cost_function )
          {
          case super3d::super_res_params::HUBER_NORM:
            break;
          case super3d::super_res_params::TRUNCATED_QUADRATIC:
            diff = rho_truncated_quadratic(diff, srp.alpha_a, srp.gamma_a);
            break;
          case super3d::super_res_params::GENERALIZED_HUBER:
            diff = rho_generalized_huber(diff, srp.alpha_a, srp.beta_a,  srp.gamma_a);
            break;
          }

          qfijk = (qfijk +  diff * val)/denom;
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
                  const vil_image_view<double> &u,
                  vcl_vector<vil_image_view<double> > &qg,
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

    // apply the linear operator to warp, blur, and downsample
    vil_image_view<double> l_u;
    warps[f].apply_A(u, l_u);
    vil_image_view<double> gradient_lu;
    vidtk::forward_gradient(l_u, gradient_lu);

    if( srp.debug )
    {
      char buf[50];
      vil_image_view<vxl_byte> output;
      sprintf(buf,"images/frame_warp_%03d.png",f);
      vil_convert_stretch_range_limited(l_u, output, 0.0, 1.0);
      vil_save( output, buf );
      sprintf(buf,"images/weights_%03d.png",f);
      vil_convert_stretch_range_limited(weights[f], output, 0.0, 1.0);
      vil_save( output, buf );
    }

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
            diff = rho_truncated_quadratic(diff, srp.alpha_g, srp.gamma_g);
            break;
          case super3d::super_res_params::GENERALIZED_HUBER:
            diff = rho_generalized_huber(diff, srp.alpha_g, srp.beta_g, srp.gamma_g );
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
void primal_step_A(const vcl_vector<vil_image_view<double> > &qa,
                   const vcl_vector<vil_image_view<double> > &qg,
                   const vcl_vector<vidtk::adjoint_image_ops_func<double> > &warps,
                   const vil_image_view<double> &pr,
                   const vil_image_view<double> &u,
                   vcl_vector<vil_image_view<double> > &As,
                   const super_res_params &srp)
{
  if( !srp.illumination_prior )
    return;

  const double sf_2 = 1.0 / (srp.scale_factor * srp.scale_factor);

  vil_image_view<double> super_qa(srp.s_ni, srp.s_nj, u.nplanes());
  vil_image_view<double> super_qay(srp.s_ni, srp.s_nj, u.nplanes());

  for (unsigned int i = 0; i < qa.size(); i++)
  {
    vil_image_view<double>& A0 = As[2*i];
    vil_image_view<double>& A1 = As[2*i+1];

    warps[i].apply_At(qa[i], super_qa);

    if( srp.gradient_data )
    {
      vil_image_view<double> super_qg( srp.s_ni, srp.s_nj, u.nplanes() );
      vil_image_view<double> div;
      vidtk::backward_divergence(qg[i], div);
      warps[i].apply_At(div, super_qg);
      vil_math_add_image_fraction(super_qa, 1.0, super_qg, -1.0);
    }
    vil_math_image_product( u, super_qa, super_qay);

    vil_image_view<double> work, buf;
    vidtk::backward_divergence(pr, work);
    buf.deep_copy( work );

    vil_math_add_image_fraction(buf, srp.tau, super_qa, -srp.tau * sf_2);
//    vil_math_image_sum(A0, buf, buf);
//    vil_math_add_image_fraction( A0, -1.0, buf, 2.0);

    vil_math_add_image_fraction(work, srp.tau, super_qay, -srp.tau * sf_2);
//    vil_math_image_sum(A1, work, work);
//    vil_math_add_image_fraction( A1, -1.0, work, 2.0);

    for(unsigned int j=0; j<buf.nj(); j++)
    {
      for(unsigned int i=0; i<buf.ni(); i++)
      {
        double& a0=A0(i,j);
        double& a1=A1(i,j);
        for(unsigned int k=0; k<buf.nplanes(); k++ )
        {
          a0 = (a0 + buf(i,j,k)) * 2.0 - a0;
          a1 = (a1 + work(i,j,k))* 2.0 - a1;
        }
      }
    }

  }
}

//*****************************************************************************
void primal_step_Y(const vcl_vector<vil_image_view<double> > &qa,
                   const vcl_vector<vil_image_view<double> > &qg,
                   const vcl_vector<vidtk::adjoint_image_ops_func<double> > &warps,
                   const vcl_vector<vil_image_view<double> > &weights,
                   const vil_image_view<double> &pr,
                   vil_image_view<double> &u,
                   const vcl_vector< vil_image_view<double> > &As,
                   const super3d::super_res_params &srp)
{
  const double sf_2 = 1.0 / (srp.scale_factor * srp.scale_factor);

  vil_image_view<double> sum_super_qa(srp.s_ni, srp.s_nj, u.nplanes());
  sum_super_qa.fill(0.0);

  vil_image_view<double> super_qa(srp.s_ni, srp.s_nj, u.nplanes());

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
    // apply transpose linear operator to upsample, blur, and warp
    vil_image_view<double> weighted_qa;
    vil_math_image_product( qa[i], weights[i], weighted_qa );
    warps[i].apply_At(weighted_qa, super_qa);

    if( srp.illumination_prior )
    {
//      vil_math_image_product(super_qa, As[2*i+1], super_qa);
      vil_image_view<double> A1=As[2*i+1];
      for(unsigned int jj=0; jj<super_qa.nj(); jj++)
      {
        for(unsigned int ii=0; ii<super_qa.ni(); ii++)
        {
          const double& a1=A1(ii,jj);
          for(unsigned int kk=0; kk<super_qa.nplanes(); kk++ )
          {
            super_qa(ii,jj,kk) *= a1;
          }
        }
      }
    }

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
      vil_image_view<double> div;
      vidtk::backward_divergence(qg[i], div);

      vil_image_view<double> weighted_qg_div;
      vil_math_image_product( div, weights[i], weighted_qg_div );
      warps[i].apply_At( weighted_qg_div, super_qg);

      if( srp.illumination_prior )
      {
//        vil_math_image_product(super_qg, As[2*i+1], super_qg);
        vil_image_view<double> A1=As[2*i+1];
        for(unsigned int jj=0; jj<super_qg.nj(); jj++)
        {
          for(unsigned int ii=0; ii<super_qg.ni(); ii++)
          {
            const double& a1=A1(ii,jj);
            for(unsigned int kk=0; kk<super_qg.nplanes(); kk++ )
            {
              super_qg(ii,jj,kk) *= a1;
            }
          }
        }
      }

      vil_math_image_sum(sum_super_qg, super_qg, sum_super_qg);
    }

    vil_math_add_image_fraction(sum_super_qa, 1.0, sum_super_qg, -1.0);
  }

  vil_image_view<double> work;
  vidtk::backward_divergence(pr, work);
  vil_math_add_image_fraction(work, srp.tau, sum_super_qa, -srp.tau * sf_2);

  vil_math_image_sum(u, work, work);
  vil_math_add_image_fraction(u, -1.0, work, 2.0);
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
  const vcl_string &output_image)
{

  if( srp.tv_method == super3d::super_res_params::SUPER3D_BASELINE )
    return super_resolve( frames, warps, Y, srp, iterations, output_image );

  vil_image_view<double>& u = Y;

  if (frames.empty())
    return;

  unsigned int np = frames[0].nplanes();
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
  vcl_vector<vil_image_view<double> > pl(frames.size()*2);
  vcl_vector<vil_image_view<double> > weights;

  for (unsigned int i = 0; i < frames.size(); i++)
  {
    vcl_cout << warps[i].dst_ni() << " " << warps[i].dst_nj() << "\n";

    qa[i].set_size(warps[i].dst_ni(), warps[i].dst_nj(), np);
    qa[i].fill(0.0);

    if( srp.gradient_data )
    {
      qg[i].set_size(warps[i].dst_ni(), warps[i].dst_nj(), 2*np);
      qg[i].fill(0.0);
    }

    unsigned int j=2*i;
    pl[j].set_size(srp.s_ni, srp.s_nj, 2);
    pl[j].fill(0.0);
    ++j;
    pl[j].set_size(srp.s_ni, srp.s_nj, 2);
    pl[j].fill(0.0);

    vil_image_view<double> erode, wmap;
    vil_structuring_element element;
    element.set_to_disk( srp.erosion_radius );
    if( warps[i].weight_map().nplanes() == 1 )
    {
      wmap = warps[i].weight_map();
    }
    else
    {
      vil_math_mean_over_planes(warps[i].weight_map(), wmap);
    }
    vil_greyscale_erode( wmap, erode, element );
    weights.push_back( erode );

    if( srp.illumination_prior )
    {
      vil_image_view<double> A0( srp.s_ni, srp.s_nj, 1 );
      A0.fill(0.0);
      As.push_back( A0 );

      vil_image_view<double> A1( srp.s_ni, srp.s_nj, 1 );
      A1.fill(1.0);
      As.push_back( A1 );
    }
  }

  double ssd;
  vil_image_view<double> last;
  last.deep_copy(u);
  unsigned int i = 0;

  switch( srp.tv_method )
  {
  case super3d::super_res_params::IMAGEDATA_IMAGEPRIOR:
    srp.image_data_1=false;
    srp.image_data_N=true;
    srp.gradient_data=false;
    srp.image_prior=true;
    srp.illumination_prior=false;
    break;
  case super3d::super_res_params::GRADIENTDATA_IMAGEPRIOR:
    srp.image_data_1=true;
    srp.image_data_N=false;
    srp.gradient_data=true;
    srp.image_prior=true;
    srp.illumination_prior=false;
    break;
  case super3d::super_res_params::IMAGEDATA_IMAGEPRIOR_ILLUMINATIONPRIOR:
    srp.image_data_1=false;
    srp.image_data_N=true;
    srp.gradient_data=false;
    srp.image_prior=true;
    srp.illumination_prior=true;
    break;
  case super3d::super_res_params::IMAGEDATA_GRADIENTDATA_IMAGEPRIOR_ILLUMINATIONPRIOR:
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
    double minv, maxv;
    switch( srp.tv_method )
    {
    case super3d::super_res_params::IMAGEDATA_IMAGEPRIOR:
      dual_step_pr(u, pr, srp);
      dual_step_qa(frames, warps, weights, u, qa, As, srp);
      primal_step_Y(qa, qg, warps, weights, pr, u, As, srp);
      break;

    case super3d::super_res_params::GRADIENTDATA_IMAGEPRIOR:
      dual_step_pr(u, pr, srp);
      dual_step_qa(frames, warps, weights, u, qa, As, srp);
      dual_step_qg(frames_gradient, warps, weights, u, qg, srp);
      primal_step_Y(qa, qg, warps, weights, pr, u, As, srp);
      break;

    case super3d::super_res_params::IMAGEDATA_IMAGEPRIOR_ILLUMINATIONPRIOR:
      dual_step_pr(u, pr, srp);
      dual_step_pl(As, pl, srp);
      dual_step_qa(frames, warps, weights, u, qa, As, srp);
      primal_step_A( qa, qg, warps, pr, u, As, srp );
      primal_step_Y(qa, qg, warps, weights, pr, u, As, srp);
      break;

    case super3d::super_res_params::IMAGEDATA_GRADIENTDATA_IMAGEPRIOR_ILLUMINATIONPRIOR:
    default:
      vcl_cerr << "unknown tv method.\n";
      return;
    }

    vil_math_value_range(u, minv, maxv);
    if (!(i % 15))
    {
      vil_image_view<double> outd;
      vil_image_view<vxl_byte> output;
      vil_convert_stretch_range_limited(u, outd, vcl_max(0.0,minv), vcl_min(1.0,maxv), 0.0, 255.0);
      vil_convert_cast(outd, output);
      vil_save(output, output_image.c_str());
    }

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
} // end namespace super3d

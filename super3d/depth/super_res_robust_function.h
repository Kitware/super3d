/*ckwg +5
 * Copyright 2014 by Kitware, Inc. All Rights Reserved. Please refer to
 * KITWARE_LICENSE.TXT for licensing information, or contact General Counsel,
 * Kitware, Inc., 28 Corporate Drive, Clifton Park, NY 12065.
 */

#ifndef SUPER_RES_ROBUST_FUNCTION_H_
#define SUPER_RES_ROBUST_FUNCTION_H_

#include <vcl_vector.h>
#include <vcl_cmath.h>

namespace super3d
{

inline double rho_huber_norm( double x, double alpha );
inline double psi_huber_norm( double x );

inline double rho_truncated_quadratic( double x, double alpha, double gamma );
inline double psi_truncated_quadratic( double x, double alpha, double gamma );

inline double rho_generalized_huber( double x, double alpha, double beta, double gamma );
inline double psi_generalized_huber( double x, double alpha, double beta, double gamma );

class rho_truncated_quadratic_functor
{
public:
  rho_truncated_quadratic_functor(double alpha, double gamma )
      : alpha_(alpha), gamma_(gamma) {}

  double operator() ( double x ) const
    {
      if( vcl_fabs(x) <= vcl_sqrt( alpha_/gamma_ ) )
        return( gamma_ * x * x );
      else
        return alpha_;
    }
private:
  double alpha_, gamma_;
};

class psi_truncated_quadratic_functor
{
public:
  psi_truncated_quadratic_functor(double alpha, double gamma )
      : alpha_(alpha), gamma_(gamma) {}

  double operator() ( double x ) const
    {
      if( vcl_fabs(x) <= vcl_sqrt( alpha_/gamma_ ) )
        return( 2.0 * gamma_ * x );
      else
        return 0.0;
    }
private:
  double alpha_, gamma_;
};

class rho_generalized_huber_functor
{
public:
  rho_generalized_huber_functor(double alpha, double beta, double gamma )
      : alpha_(alpha), beta_(beta), gamma_(gamma) {}

  double operator() ( double x ) const
    {
      const double t=vcl_sqrt( alpha_/gamma_ );
      if( vcl_fabs(x) <= t )
        return( gamma_ * x * x );
      else
        return beta_*x + alpha_ - t*beta_;
    }
private:
  double alpha_, beta_, gamma_;
};

class psi_generalized_huber_functor
{
public:
  psi_generalized_huber_functor(double alpha, double beta, double gamma )
      : alpha_(alpha), beta_(beta), gamma_(gamma) {}

  double operator() ( double x ) const
    {
      if( vcl_fabs(x) <= vcl_sqrt( alpha_/gamma_ ) )
        return( 2.0 * gamma_ * x );
      else
        return beta_;
    }
private:
  double alpha_, beta_, gamma_;
};

} // end namespace super3d

#endif

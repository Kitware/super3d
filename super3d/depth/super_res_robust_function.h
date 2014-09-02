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
/// The psi_ functions are the first order derivative of the robust function rho_

/// Huber norm is a robust function with quadratic form around 0 and L-1 elsewhere
double rho_huber_norm( double x, double alpha );
/// First order derivative of Huber norm
double psi_huber_norm( double x );

/// Truncated quadratic is a robust function with quadratic form around 0 and flat elsewhere
double rho_truncated_quadratic( double x, double alpha, double gamma );
/// First order derivative of Truncated quadratic function
double psi_truncated_quadratic( double x, double alpha, double gamma );

/// Generalized Huber is a robust function with quadratic form around 0 and constant slope elsewhere
double rho_generalized_huber( double x, double alpha, double beta, double gamma );
/// First order derivative of Generalized Huber
double psi_generalized_huber( double x, double alpha, double beta, double gamma );

class rho_truncated_quadratic_functor
{
public:
  rho_truncated_quadratic_functor(double alpha, double gamma )
      : alpha_(alpha), gamma_(gamma) {}

  double operator() ( double x ) const
    {
      return super3d::rho_truncated_quadratic( x, alpha_, gamma_ );
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
      return super3d::psi_truncated_quadratic( x, alpha_, gamma_ );
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
      return super3d::rho_generalized_huber( x, alpha_, beta_, gamma_ );
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
      return super3d::psi_generalized_huber( x, alpha_, beta_, gamma_ );
    }
private:
  double alpha_, beta_, gamma_;
};

} // end namespace super3d

#endif

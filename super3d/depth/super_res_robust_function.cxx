/*ckwg +5
 * Copyright 2014 by Kitware, Inc. All Rights Reserved. Please refer to
 * KITWARE_LICENSE.TXT for licensing information, or contact General Counsel,
 * Kitware, Inc., 28 Corporate Drive, Clifton Park, NY 12065.
 */

#include "super_res_robust_function.h"

namespace super3d
{

double rho_huber_norm( double x, double alpha )
{
  const double absx = vcl_fabs(x);
  if( absx <= alpha )
    return x*x / 2.0 / alpha;
  else
    return absx - alpha/2.0;
}

double psi_huber_norm( double x, double alpha )
{
  if( vcl_fabs(x) <= alpha )
    return 1.0;
  else
    return -1.0;
}

double rho_truncated_quadratic( double x, double alpha, double gamma )
{
  const double val = gamma * x * x;
  if( val < alpha )
    return val;
  else
    return alpha;
}

double psi_truncated_quadratic( double x, double alpha, double gamma )
{
  if( (x*x*gamma) <= alpha )
    return 2.0 * gamma * x;
  else
    return 0.0;
}

double rho_generalized_huber( double x, double alpha, double beta, double gamma )
{
  double t = vcl_sqrt( alpha/gamma );
  if( vcl_fabs(x) <= t )
    return gamma * x * x;
  else
    return beta * x + alpha - t * beta;
}

double psi_generalized_huber( double x, double alpha, double beta, double gamma )
{
  if( (x*x*gamma) <= alpha )
    return 2.0 * gamma * x;
  else
    return beta;
}

} // end namespace super3d

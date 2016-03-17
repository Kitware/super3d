/*ckwg +29
 * Copyright 2014 by Kitware, Inc.
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

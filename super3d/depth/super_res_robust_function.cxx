/*ckwg +29
 * Copyright 2014-2016 by Kitware, Inc.
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

#include "super_res_robust_function.h"

namespace super3d
{

double rho_huber_norm( double x, double alpha )
{
  const double absx = std::fabs(x);
  if( absx <= alpha )
    return x*x / 2.0 / alpha;
  else
    return absx - alpha/2.0;
}

double psi_huber_norm( double x, double alpha )
{
  if( std::fabs(x) <= alpha )
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
  double t = std::sqrt( alpha/gamma );
  if( std::fabs(x) <= t )
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

/*
 * Copyright 2012 Kitware, Inc.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of this project nor the names of its contributors
 *       may be used to endorse or promote products derived from this software
 *       without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef SUPER_RES_H_
#define SUPER_RES_H_

#include <boost/scoped_ptr.hpp>

#include "super_config.h"
#include "depth_config.h"
#include "super_res_robust_function.h"

#include <vil/vil_image_view.h>
#include <vcl_vector.h>

#include <super3d/image/adjoint_image_op.h>

#include <boost/function.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/locks.hpp>

namespace super3d
{

struct super_res_params {
  double lambda, epsilon_data, epsilon_reg, tau, sigma;
  unsigned int s_ni, s_nj, l_ni, l_nj;
  double scale_factor;
  unsigned int ref_frame;

  // illumination robust optimization parameters
  enum TV_METHOD
  {
    SUPER3D_BASELINE=0,
    IMAGEDATA_IMAGEPRIOR,
    GRADIENTDATA_IMAGEPRIOR,
    IMAGEDATA_IMAGEPRIOR_ILLUMINATIONPRIOR,
    MEDIANDATA_IMAGEPRIOR,
    IMAGEDATA_GRADIENTDATA_IMAGEPRIOR_ILLUMINATIONPRIOR
  } tv_method;

  bool debug;
  double erosion_radius;
  double median_radius;
  int frame_step;
  bool image_data_1, image_data_N, gradient_data, image_prior, illumination_prior, median_residue;

  // 0:huber_norm; 1:truncated_quadratic; 2:generalized_huber
  enum COST_FUNCTION
  {
    HUBER_NORM=0,
    TRUNCATED_QUADRATIC,
    GENERALIZED_HUBER
  } cost_function;

  // image data term (lambda, alpha, beta, gamma)
  double alpha_a, gamma_a, beta_a;

  // gradient data term (lambda, alpha, beta, gamma)
  double lambda_g, alpha_g, gamma_g, beta_g;

  // image prior term (lambda, alpha, beta, gamma)
  double lambda_r, alpha_r, gamma_r, beta_r;

  // illumination prior term (lambda, alpha, beta, gamma)
  double lambda_l, alpha_l, gamma_l, beta_l;

  // dual space
  double sigma_pr, sigma_pl, sigma_qa, sigma_qg, sigma_A, sigma_Y;
};

class super_res_monitor;

SUPER3D_DEPTH_EXPORT
void super_resolve(
  const vcl_vector<vil_image_view<double> > &frames,
  const vcl_vector<vidtk::adjoint_image_ops_func<double> > &warps,
  vil_image_view<double> &u,
  const super_res_params &srp,
  unsigned int iterations,
  const vcl_string &output_image = "",
  super_res_monitor *srm = NULL);

SUPER3D_DEPTH_EXPORT
void super_resolve_robust(
  const vcl_vector<vil_image_view<double> > &frames,
  const vcl_vector<vidtk::adjoint_image_ops_func<double> > &warps,
  vil_image_view<double> &Y,
  super_res_params srp,
  unsigned int iterations,
  vcl_vector< vil_image_view<double> > &As,
  const vcl_string &output_image = "",
  super_res_monitor *srm = NULL);


class SUPER3D_DEPTH_EXPORT super_res_monitor
{
public:

  struct update_data
  {
    vil_image_view<double> current_result;
    unsigned int num_iterations;
  };

  super_res_monitor(boost::function<void (update_data)> callback,
                    unsigned int interval,
                    boost::shared_ptr<bool> interrupted)
    : callback_(callback),
    interval_(interval),
    interrupted_(interrupted) {}

private:

  friend void super_resolve_robust(
    const vcl_vector<vil_image_view<double> > &frames,
    const vcl_vector<vidtk::adjoint_image_ops_func<double> > &warps,
    vil_image_view<double> &Y,
    super_res_params srp,
    unsigned int iterations,
    vcl_vector< vil_image_view<double> > &As,
    const vcl_string &output_image,
    super_res_monitor *srm);

  friend void super_resolve(
    const vcl_vector<vil_image_view<double> > &frames,
    const vcl_vector<vidtk::adjoint_image_ops_func<double> > &warps,
    vil_image_view<double> &u,
    const super_res_params &srp,
    unsigned int iterations,
    const vcl_string &output_image,
    super_res_monitor *srm);

  boost::function<void (update_data)> callback_;
  unsigned int interval_;
  boost::shared_ptr<bool const> interrupted_;
};



SUPER3D_DEPTH_EXPORT
void compare_to_original(const vil_image_view<double> &ref_img,
                         const vil_image_view<double> &super,
                         const vil_image_view<double> &original,
                         unsigned int scale_factor);

SUPER3D_DEPTH_EXPORT
void read_super_res_params( const boost::scoped_ptr<config>& cfg,
                            super_res_params &srp);

} // end namespace super3d

#endif

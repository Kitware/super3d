/*ckwg +29
 * Copyright 2012-2016 by Kitware, Inc.
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

#include "weighted_dbw.h"
#include "super_config.h"

#include <functional>

#include <super3d/image/adjoint_flow_warp.h>
#include <super3d/image/adjoint_resample.h>

#include <vil/algo/vil_gauss_filter.h>
#include <vil/vil_math.h>


namespace super3d
{

/// Helper class to contain a set of 2 image operators that run in sequence.
/// Intermediate results are cached to avoid re-allocation of memory
template <typename T>
class image_op_2_func
{
public:
  typedef std::function<void (const vil_image_view<T>& src,
                              vil_image_view<T>& dst)> func_t;

  image_op_2_func(func_t op1, func_t op2)
  : op1_(op1),
    op2_(op2) {}

  image_op_2_func(func_t op1, func_t op2,
                  unsigned ni, unsigned nj, unsigned np)
  : op1_(op1),
    op2_(op2),
    tmp_(ni, nj, np) {}

  /// Apply the sequence of 2 image operations
  void operator()(const vil_image_view<T>& src,
                  vil_image_view<T>& dst) const
  {
    op1_(src, tmp_);
    op2_(tmp_, dst);
  }

  /// Set the size of the intermediate result image between op1 and op2
  void set_intermediate_size(unsigned ni, unsigned nj, unsigned np)
  {
    tmp_.set_size(ni, nj, np);
  }

  /// Set the first function
  void set_op1(func_t op1)
  {
    op1_ = op1;
  }

  /// Set the second function
  void set_op2(func_t op2)
  {
    op2_ = op2;
  }

private:
  func_t op1_;
  func_t op2_;
  mutable vil_image_view<T> tmp_;
};


super3d::adjoint_image_ops_func<double>
create_dbw_from_flow(const vil_image_view<double> &flow,
                     const vil_image_view<double> &weights,
                     const unsigned ni, const unsigned nj, const unsigned np,
                     int scale_factor, double sensor_sigma, bool down_sample_averaging,
                     bool bicubic_warping)
{
  const unsigned int sni = flow.ni();
  const unsigned int snj = flow.nj();
  const unsigned int wni = ni * scale_factor;
  const unsigned int wnj = nj * scale_factor;

  typedef std::function<void (const vil_image_view<double>& src,
                              vil_image_view<double>& dst)> func_t;

  func_t weight_f = std::bind(vil_math_image_product<double, double, double>,
                              std::placeholders::_1, weights, std::placeholders::_2);

  func_t warp_fwd = std::bind(super3d::warp_forward_with_flow_bilin<double, double, double>,
                              std::placeholders::_1, flow, std::placeholders::_2);
  func_t warp_back = std::bind(super3d::warp_back_with_flow_bilin<double, double, double>,
                               std::placeholders::_1, flow, std::placeholders::_2);
  if(bicubic_warping)
  {
    warp_fwd = std::bind(super3d::warp_forward_with_flow_bicub<double, double, double>,
                         std::placeholders::_1, flow, std::placeholders::_2);
    warp_back = std::bind(super3d::warp_back_with_flow_bicub<double, double, double>,
                          std::placeholders::_1, flow, std::placeholders::_2);
  }

  typedef void (*filter_func_t)(const vil_image_view<double>&, vil_image_view<double> &,
                                double, unsigned, vil_convolve_boundary_option boundary);
  func_t blur_sf = std::bind(static_cast<filter_func_t>(vil_gauss_filter_2d<double, double>),
                             std::placeholders::_1, std::placeholders::_2, sensor_sigma,
                             static_cast<unsigned int>(3.0 * sensor_sigma), vil_convolve_zero_extend);

  func_t down_s = std::bind(super3d::down_sample<double>, std::placeholders::_1,
                            std::placeholders::_2, scale_factor, 0, 0);
  func_t up_s = std::bind(super3d::up_sample<double>, std::placeholders::_1,
                          std::placeholders::_2, scale_factor, 0, 0);

  if(down_sample_averaging)
  {
    down_s = std::bind(super3d::down_scale<double>, std::placeholders::_1,
                       std::placeholders::_2, scale_factor);
    up_s = std::bind(super3d::up_scale<double>, std::placeholders::_1,
                     std::placeholders::_2, scale_factor);
  }

  image_op_2_func<double> forward_ww(weight_f, warp_fwd, sni, snj, np);
  image_op_2_func<double> forward_bs(blur_sf, down_s, wni, wnj, np);
  image_op_2_func<double> forward(forward_ww, forward_bs, wni, wnj, np);

  image_op_2_func<double> backward_ww(warp_back, weight_f, sni, snj, np);
  image_op_2_func<double> backward_bs(up_s, blur_sf, wni, wnj, np);
  image_op_2_func<double> backward(backward_bs, backward_ww, wni, wnj, np);

  return super3d::adjoint_image_ops_func<double>(forward, backward,
                                               sni, snj, np,
                                               ni, nj, np);
}

} // end namespace super3d

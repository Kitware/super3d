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

#include "weighted_dbw.h"
#include "super_config.h"

#include <video_transforms/adjoint_flow_warp.h>
#include <video_transforms/adjoint_resample.h>

#include <boost/bind.hpp>

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
  typedef boost::function<void (const vil_image_view<T>& src,
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


vidtk::adjoint_image_ops_func<double>
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

  typedef boost::function<void (const vil_image_view<double>& src,
                                vil_image_view<double>& dst)> func_t;

  func_t weight_f = boost::bind(vil_math_image_product<double, double, double>,
                                _1, weights, _2);

  func_t warp_fwd = boost::bind(vidtk::warp_forward_with_flow_bilin<double, double, double>,
                                _1, flow, _2);
  func_t warp_back = boost::bind(vidtk::warp_back_with_flow_bilin<double, double, double>,
                                 _1, flow, _2);
  if(bicubic_warping)
  {
    warp_fwd = boost::bind(vidtk::warp_forward_with_flow_bicub<double, double, double>,
                           _1, flow, _2);
    warp_back = boost::bind(vidtk::warp_back_with_flow_bicub<double, double, double>,
                            _1, flow, _2);
  }

  func_t blur_sf = boost::bind(vil_gauss_filter_2d<double, double>, _1, _2, sensor_sigma,
                               static_cast<unsigned int>(3.0 * sensor_sigma), vil_convolve_zero_extend);

  func_t down_s = boost::bind(vidtk::down_sample<double>, _1, _2, scale_factor, 0, 0);
  func_t up_s = boost::bind(vidtk::up_sample<double>, _1, _2, scale_factor, 0, 0);

  if(down_sample_averaging)
  {
    down_s = boost::bind(vidtk::down_scale<double>, _1, _2, scale_factor);
    up_s = boost::bind(vidtk::up_scale<double>, _1, _2, scale_factor);
  }

  image_op_2_func<double> forward_ww(weight_f, warp_fwd, sni, snj, np);
  image_op_2_func<double> forward_bs(blur_sf, down_s, wni, wnj, np);
  image_op_2_func<double> forward(forward_ww, forward_bs, wni, wnj, np);

  image_op_2_func<double> backward_ww(warp_back, weight_f, sni, snj, np);
  image_op_2_func<double> backward_bs(up_s, blur_sf, wni, wnj, np);
  image_op_2_func<double> backward(backward_bs, backward_ww, wni, wnj, np);

  return vidtk::adjoint_image_ops_func<double>(forward, backward,
                                               sni, snj, np,
                                               ni, nj, np);
}

} // end namespace super3d

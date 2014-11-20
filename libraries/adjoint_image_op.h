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

#ifndef norm_estimation_h_
#define norm_estimation_h_

#include <vcl_vector.h>

#include <vil/vil_image_view.h>
#include <boost/function.hpp>


/// Abstract base class for adjoint image operators.
/// This class represents a pair of image operators A and A^t and provides
/// utilities such as verifying the adjoint property and computing the 2-norm.
template <typename T>
class adjoint_image_op
{
public:
  /// Apply the primary image operator
  virtual void apply_A(const vil_image_view<T>& src,
                       vil_image_view<T>& dst) const = 0;

  /// Apply the adjoint image operator
  virtual void apply_At(const vil_image_view<T>& src,
                        vil_image_view<T>& dst) const = 0;

  /// Apply the primary followed by adjoint operators
  virtual void apply_AtA(const vil_image_view<T>& src,
                         vil_image_view<T>& dst) const;

  /// The width of the image that is operated on.
  virtual unsigned src_ni() const = 0;

  /// The height of the image that is operated on.
  virtual unsigned src_nj() const = 0;

  /// The number of planes in the image that is operated on.
  virtual unsigned src_nplanes() const = 0;

  /// The width of the destination image.
  virtual unsigned dst_ni() const = 0;

  /// The height of the destination image.
  virtual unsigned dst_nj() const = 0;

  /// The number of planes in the destination image.
  virtual unsigned dst_nplanes() const = 0;

  /// Estimate the 2-norm of image operator \p op.
  /// The algorithm uses a power series to estimate the 2-norm.
  /// This version is initialized with a random image
  /// \param tol The tolerance on the relative change between iterations.
  double norm_estimation(double tol=1e-6) const;

  /// Estimate the 2-norm of image operator \p op.
  /// The algorithm uses a power series to estimate the 2-norm.
  /// \param vec Initial image "vector" used and max eigenvector returned.
  /// \param tol The tolerance on the relative change between iterations
  double norm_estimation(vil_image_view<T>& vec,
                         double tol=1e-6) const;

  /// Verify that the operators are in fact adjoint.
  /// Compares "vec2 A^t A vec1" to "vec1 A^t A vec2" for
  /// random images vec1 and vec2
  bool is_adjoint(double tol=1e-8) const;

  /// Verify that the operators are in fact adjoint.
  /// Compares "vec2 A^t A vec1" to "vec1 A^t A vec2" for
  /// specified images vec1 and vec2
  bool is_adjoint(const vil_image_view<T>& vec1,
                  const vil_image_view<T>& vec2,
                  double tol=1e-8) const;
};


/// A collection of adjoint image operators that acts like a single operators.
/// This is for the special case a single primal image with many dual images
/// and should be used primarily for norm estimation.
template <typename T>
class adjoint_multi_image_op
 : public adjoint_image_op<T>
{
public:
  /// Constructor
  /// \note Does not take ownership of image_op memory
  adjoint_multi_image_op(const vcl_vector<adjoint_image_op<T>* >& image_ops);

  /// Primary image operator is undefined
  virtual void apply_A(const vil_image_view<T>& src,
                       vil_image_view<T>& dst) const {}

  /// Adjoint image operator is undefined
  virtual void apply_At(const vil_image_view<T>& src,
                        vil_image_view<T>& dst) const {}

  /// Apply the primary followed by adjoint operators
  /// computed as the sum of apply_AtA on all sub-operations
  virtual void apply_AtA(const vil_image_view<T>& src,
                         vil_image_view<T>& dst) const;

  /// The width of the image that is operated on.
  virtual unsigned src_ni() const;

  /// The height of the image that is operated on.
  virtual unsigned src_nj() const;

  /// The number of planes in the image that is operated on.
  virtual unsigned src_nplanes() const;

  /// The width of the destination image is undefined.
  virtual unsigned dst_ni() const { return 0; }

  /// The height of the destination image is undefined.
  virtual unsigned dst_nj() const { return 0; }

  /// The number of planes in the destination image is undefined.
  virtual unsigned dst_nplanes() const { return 0; }

private:
  vcl_vector<adjoint_image_op<T>* > image_ops_;
};


template <typename T>
class adjoint_image_ops_func
 : public adjoint_image_op<T>
{
public:
  typedef boost::function<void (const vil_image_view<T>& src,
                                vil_image_view<T>& dst)> func_t;

  adjoint_image_ops_func(func_t forward, func_t backward,
                         unsigned ni, unsigned nj, unsigned np)
  : forward_(forward),
    backward_(backward),
    src_ni_(ni),
    src_nj_(nj),
    src_np_(np),
    dst_ni_(ni),
    dst_nj_(nj),
    dst_np_(np) {}

  adjoint_image_ops_func(func_t forward, func_t backward,
                         unsigned sni, unsigned snj, unsigned snp,
                         unsigned dni, unsigned dnj, unsigned dnp)
  : forward_(forward),
    backward_(backward),
    src_ni_(sni),
    src_nj_(snj),
    src_np_(snp),
    dst_ni_(dni),
    dst_nj_(dnj),
    dst_np_(dnp) {}

  /// Apply the primary image operator
  virtual void apply_A(const vil_image_view<T>& src,
                       vil_image_view<T>& dst) const
  {
    dst.set_size(dst_ni_, dst_nj_, dst_np_);
    forward_(src, dst);
  }

  /// Apply the adjoint image operator
  virtual void apply_At(const vil_image_view<T>& src,
                        vil_image_view<T>& dst) const
  {
    dst.set_size(src_ni_, src_nj_, src_np_);
    backward_(src, dst);
  }

  /// generate a weight image by running apply_A on an image of all ones
  virtual vil_image_view<T> weight_map() const
  {
    vil_image_view<T> weights(dst_ni_, dst_nj_, dst_np_);
    vil_image_view<T> ones(src_ni_, src_nj_, src_np_);
    ones.fill(1.0);
    this->apply_A(ones, weights);
    return weights;
  }

  /// The width of the image that is operated on.
  virtual unsigned src_ni() const { return src_ni_; }

  /// The height of the image that is operated on.
  virtual unsigned src_nj() const { return src_nj_; }

  /// The number of planes in the image that is operated on.
  virtual unsigned src_nplanes() const { return src_np_; }

  /// The width of the image that is operated on.
  virtual unsigned dst_ni() const { return dst_ni_; }

  /// The height of the image that is operated on.
  virtual unsigned dst_nj() const { return dst_nj_; }

  /// The number of planes in the image that is operated on.
  virtual unsigned dst_nplanes() const { return dst_np_; }

  /// Set the size of the source image
  void set_src_size(unsigned ni, unsigned nj, unsigned np)
  {
    src_ni_ = ni;
    src_nj_ = nj;
    src_np_ = np;
  }

  /// Set the size of the destination image
  void set_dst_size(unsigned ni, unsigned nj, unsigned np)
  {
    dst_ni_ = ni;
    dst_nj_ = nj;
    dst_np_ = np;
  }

protected:
  func_t forward_;
  func_t backward_;
  unsigned src_ni_;
  unsigned src_nj_;
  unsigned src_np_;
  unsigned dst_ni_;
  unsigned dst_nj_;
  unsigned dst_np_;
};



#endif //norm_estimation_h_

/*ckwg +29
 * Copyright 2013-2016 by Kitware, Inc.
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

#include "warp_and_average.hxx"

#define INSTANTIATE( PIXTYPE ) \
\
template bool super3d::warp_and_average(const std::vector<vil_image_view<PIXTYPE> > &src_frames, \
                                      vil_image_view<PIXTYPE> &dest, \
                                      const std::vector<vgl_h_matrix_2d<double> > &homogs, \
                                      unsigned int ref_frame, unsigned int i0, unsigned int j0, \
                                      unsigned int ni, unsigned int nj, \
                                      const warp_image_parameters &wip, double scale_factor); \
\
template void super3d::warp_and_average(const std::vector<vil_image_view<PIXTYPE> > &src_frames, \
                                      vil_image_view<PIXTYPE> &dest, \
                                      const std::vector<vgl_h_matrix_2d<double> > &homogs, \
                                      unsigned int ref_frame, const warp_image_parameters &wip, \
                                      double scale_factor); \

//Does not work with non float
INSTANTIATE( float );
INSTANTIATE( double );

#undef INSTANTIATE

/*ckwg +29
 * Copyright 2013 by Kitware, Inc.
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

#ifndef refine_homog_h_
#define refine_homog_h_

#include <vil/vil_image_view.h>
#include <vnl/vnl_double_3x3.h>
#include <vnl/vnl_double_2.h>

#include <vcl_vector.h>

/// Refines the planar homography that warps from moving -> fixed images
/// \param fixed the destination image of the homography
/// \param moving the source image of the homography
/// \param H the homography to be refined
/// \param grad_thresh the threshold on gradients for edgel computation
/// \param num_iters the number of ICP iterations
/// \param search_radius the max distance to search for correspondences
/// \param normal_cutoff the threshold on cos(angle) between the warped normal and fixed normal
/// \param sigma the smoothing scale of the edgels, 0 means no smoothing of the images
void refine_homography(const vil_image_view<double> &fixed, const vil_image_view<double> &moving,
                       vnl_double_3x3 &H, double grad_thresh, unsigned int num_iters,
                       unsigned int search_radius, double normal_cutoff = 0.8, double sigma = 0.0);

/// Crops a homography
/// \param H the homography to crop
/// \param i0 the left coordiante of the crop region
/// \param j0 the top coordinate of the crop region
void crop_homography(vnl_double_3x3 &H, unsigned int i0, unsigned int j0);


struct edgel
{
  edgel(double x, double y, double angle, double m);

  vnl_double_2 pos, n;
  double mag;
};

/// Extracts edgels from the moving image into a vector
/// \param img the image to compute edgels on
/// \param grad_thresh the threshold on gradients for edgel computation
/// \param sigma the smoothing scale of the edgels, 0 means no smoothing of the images
/// \param edgels the vector the contains the edgels
void extract_driving_edgels(const vil_image_view<double> &img,
                            double grad_thresh,
                            double sigma,
                            vcl_vector<edgel> &edgels);

/// Extracts edgels from the fixed image into a vector and index map
/// \param img the image to compute edgels on
/// \param grad_thresh the threshold on gradients for edgel computation
/// \param sigma the smoothing scale of the edgels, 0 means no smoothing of the images
/// \param edgels the vector the contains the edgels
/// \param pixel map of indicies into the edgel vector, 0 means no edgel
void extract_matchable_edgels(const vil_image_view<double> &img,
                              double grad_thresh,
                              double sigma,
                              vcl_vector<edgel> &edgels,
                              vil_image_view<unsigned int> &index);

struct match
{
  match(const edgel *m, const edgel *f) : moving(m), fixed(f) {}

  const edgel *moving;
  const edgel *fixed;
};

/// Warps a normal with H
/// \param H the homography to warp with
/// \param pt the location of the normal
/// \param normal the vector to warp
vnl_double_2 warp_normal(const vnl_double_3x3 &H, const vnl_double_2 &pt, const vnl_double_2 &normal);

/// Matches a moving edgel set with a fixed edgel set
/// \param H the homography that warps moving to fixed
/// \param e_fixed the vector of fixed edgels
/// \param e_moving the vector of moving edgels
/// \param index the pixel map of indicies into e_fixed
/// \param search_radius the max distance to search for correspondences in fixed coord system
/// \param matches vector of matches for the computed matches
/// \param normal_cutoff the threshold on cos(angle) between the warped normal and fixed normal
void match_edgels(const vnl_double_3x3 &H,
                  const vcl_vector<edgel> &e_fixed,
                  const vcl_vector<edgel> &e_moving,
                  const vil_image_view<unsigned int> &index,
                  const double search_rad,
                  vcl_vector<match> &matches,
                  double normal_cutoff = 0.8);

/// Uses Levenberg-Marquardt to estimate a new homography using normal distances
/// \param corresp the edgel matches
/// \param H the homography
void estimate_homog_lm(const vcl_vector<match> &corresp, vnl_double_3x3 &H);

vnl_double_2 mult_and_norm(const vnl_double_3x3 &T, const vnl_double_2 &x);

#endif

/*ckwg +29
 * Copyright 2011-2016 by Kitware, Inc.
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


#include <vil/vil_decimate.h>
#include <vil/algo/vil_sobel_1x3.h>

#include "gaussian_pyramid_builder.h"

/// \file
/// Build a Gaussian pyramid

namespace super3d
{

/// Compute a Gaussian pyramid according to the constructed parameters
/// Index 0 is the original image.
template <class PixType>
void
gaussian_pyramid_builder
::build_pyramid(const vil_image_view<PixType>& image,
                std::vector<vil_image_view<PixType> >& gauss_pyr) const
{
  gauss_pyr.clear();    //In case a full pyramid is passed in
  gauss_pyr.reserve(levels);

  vil_image_view<PixType> smooth_image;
  vil_image_view<PixType> subsamp_img = image;
  subsamp_img.deep_copy(image);
  gauss_pyr.push_back(subsamp_img);

  for (unsigned i = 1; i < levels; ++i)
  {
    vil_gauss_filter_5tap<PixType, PixType>(subsamp_img, smooth_image, gauss_params);
    subsamp_img = vil_image_view<PixType>();
    subsamp_img.deep_copy(vil_decimate(smooth_image, decimation));
    gauss_pyr.push_back(subsamp_img);
  }
}


/// Compute a Gaussian and gradient pyramid.
/// Level 0 is the original image.
template <class PixType, class GradType>
void
gaussian_pyramid_builder
::build_pyramid(const vil_image_view<PixType>& image,
                std::vector<vil_image_view<PixType> >& gauss_pyr,
                std::vector<vil_image_view<GradType> >& grad_pyr) const
{
  gauss_pyr.clear();
  grad_pyr.clear();

  gauss_pyr.reserve(levels);
  grad_pyr.reserve(levels);

  vil_image_view<PixType> smooth_image, subsamp_img;
  vil_image_view<GradType> grad_ij;

  subsamp_img.deep_copy(image);
  gauss_pyr.push_back(subsamp_img);

  for (unsigned int i = 1; i < levels; ++i)
  {
    //smooth before gradient
    vil_gauss_filter_5tap(subsamp_img, smooth_image, gauss_params);
    vil_sobel_1x3(smooth_image, grad_ij);
    grad_pyr.push_back(grad_ij);

    //decimate after gradient
    subsamp_img = vil_image_view<PixType>();
    subsamp_img.deep_copy(vil_decimate(smooth_image, decimation));
    gauss_pyr.push_back(subsamp_img);
  }

  vil_gauss_filter_5tap(subsamp_img, smooth_image, gauss_params);
  vil_sobel_1x3(smooth_image, grad_ij);
  grad_pyr.push_back(grad_ij);
}


/// This function flattens a Gaussian pyramid of images to a single image
/// If it is given a 2 plane image it will expand it to 3 and display in RGB
/// \param pyramid is the Gaussian pyramid / gradient of Gaussian pyramid
/// \param tiled is the flattened image
template<typename PixType>
void
tile_pyramid(const std::vector<vil_image_view<PixType> > &pyramid,
             vil_image_view<PixType> &tiled)
{
  unsigned int nplanes = pyramid[0].nplanes();
  unsigned int width = pyramid[0].ni();
  unsigned int height = pyramid[0].nj();
  if (pyramid.size() > 1)
    height += pyramid[1].nj();

  //Gradient pyramids have 2 planes, expand the image to 3
  if (nplanes == 2)
    tiled.set_size(width, height, 3);
  else
    tiled.set_size(width, height, 1);
  tiled.fill(static_cast<PixType>(0));
  unsigned int ei = 0, ej = 0;

  //Place pyramid images next to each other below the first image in output image
  for (unsigned int l = 0; l < pyramid.size(); l++)
  {
    const vil_image_view<PixType> &level = pyramid[l];

    for (unsigned int p = 0; p < nplanes; p++)
    {
      for (unsigned int i = 0; i < level.ni(); i++)
      {
        for (unsigned int j = 0; j < level.nj(); j++)
        {
          tiled(i+ei,j+ej,p) = level(i,j,p);
        }
      }
    }

    if (l == 0)
      ej = level.nj();
    else
      ei += level.ni();
  }
}

}

#ckwg +28
#Copyright 2013 by Kitware, Inc.
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
#  * Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
#
#  * Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
#  * Neither name of Kitware, Inc. nor the names of any contributors may be used
#    to endorse or promote products derived from this software without specific
#    prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS ``AS IS''
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHORS OR CONTRIBUTORS BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


"""Module for calculating various error metrics of images.
"""

import numpy as np

from vidpy.image import filters


__all__ = ['rmse', 'max_error', 'mean_error', 'median_error', 'psnr',
           'binary_flips', 'ssim']


def rmse(image1, image2):
    """Computes the root mean square error."""
    abs_diff = np.abs(np.cast['double'](image1) - image2)
    return np.sqrt(np.sum(abs_diff ** 2) / abs_diff.size)


def max_error(image1, image2):
    """Computes the maximum absolute difference."""
    return np.max(np.abs(np.cast['double'](image1) - image2))


def mean_error(image1, image2):
    """Computes the mean absolute difference."""
    return np.mean(np.abs(np.cast['double'](image1) - image2))


def median_error(image1, image2):
    """Computes the median absolute difference."""
    return np.median(np.abs(np.cast['double'](image1) - image2))


def psnr(image1, image2, value_range=None):
    """Computes the peak signal-to-noise ratio.

    Optional value_range argument specifies the value range, if not
    specified it is determined from the image1 data type.
    """
    type_ranges = {'bool':1,
                   'uint8':255, 'int8':255,
                   'uint16':2**16-1, 'int16':2**16-1,
                   'uint32':2**32-1, 'int32':2**32-1,
                   'uint64':2**64-1, 'int64':2**64-1}
    if value_range is None:
        value_range = type_ranges[image1.dtype.name]
    return 20 * np.log10(value_range / rmse(image1, image2))


def binary_flips(image1, image2):
    """Compute the number pixels that differ between two binary images.

    Returns a tuple (lost, gained) indicating the numbers of true pixels
    that have been lost and gained from the image1 to image2.
    """
    assert image1.dtype == 'bool' and image2.dtype == 'bool'
    lost = np.sum(np.logical_and(image1, np.logical_not(image2)))
    gained = np.sum(np.logical_and(np.logical_not(image1), image2))
    return (lost, gained) # sum these to get total change


def ssim(image1, image2, dynamic_range=255):
    """Compute the structureal similarity between two images.

    Implementation of structural similarity from this: Z. Wang, A. C.
    Bovik, H. R. Sheikh and E. P. Simoncelli, "Image quality
    assessment: From error visibility to structural similarity," IEEE
    Transactions on Image Processing, vol. 13, no. 4, pp. 600-612,
    Apr. 2004.
    """

    # stabilizers
    k1 = 0.01
    k2 = 0.03

    # Constants
    C1 = (k1 * dynamic_range) ** 2
    C2 = (k2 * dynamic_range) ** 2

    # Preliminary Calculations
    img1 = np.cast['double'](image1)
    img2 = np.cast['double'](image2)

    img1_sq = img1 ** 2
    img2_sq = img2 ** 2
    img1_img2 = img1 * img2
    mu1 = filters.gaussian_filter(img1, 1.5)
    mu2 = filters.gaussian_filter(img2, 1.5)
    mu1_mu2 = mu1 * mu2
    mu1_sq = mu1 ** 2
    mu2_sq = mu2 ** 2

    sigma1_sq = filters.gaussian_filter(img1_sq, 1.5)
    sigma1_sq = sigma1_sq - mu1_sq

    sigma2_sq = filters.gaussian_filter(img2_sq, 1.5)
    sigma2_sq = sigma2_sq - mu2_sq

    sigma12 = filters.gaussian_filter(img1_img2, 1.5)
    sigma12 = sigma12 - mu1_mu2

    # Bones of the algorithm
    ssim_map = ((2 * mu1_mu2 + C1) * (2 * sigma12 + C2)) / ((mu1_sq + mu2_sq + C1) * (sigma1_sq + sigma2_sq + C2))
    mssim = np.mean(ssim_map)

    return (mssim, ssim_map)

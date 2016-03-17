!/usr/bin/python
#ckwg +28
#Copyright 2012 by Kitware, Inc.
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

from optparse import OptionParser
import Image
from libtiff import TIFF
import numpy as np

import matplotlib.pyplot as plt


def load_image(fname):
    try:
        return np.array(Image.open(fname))
    except IOError:
        return TIFF.open(fname, mode='r').read_image()


def main():
    usage = "usage: %prog [options] gt_depthmap depthmap1 depthmap2 ...\n\n"
    usage += "  Score differences between depthmaps.\n"
    parser = OptionParser(usage=usage)

    parser.add_option("-b", "--border", default=None, type="int",
                      action="store", dest="border",
                      help="ignore this many pixel around the border")

    parser.add_option("-r", "--range", default=None, type="string",
                      action="store", dest="range",
                      help="and expression that evaluates to an X range")

    parser.add_option("-l", "--semilog", default=False,
                      action="store_true", dest="semilog",
                      help="plot the X range on a log scale")

    (options, args) = parser.parse_args()

    if len(args) < 2:
        exit("requires at least 2 depth images")

    gt_file = args[0]
    depth_files = args[1:]

    gt_img = load_image(gt_file)
    means = []
    medians = []
    for f in depth_files:
        img = load_image(f)
        diff = abs(img - gt_img)
        if options.border:
            b = options.border
            diff = diff[b:-b, b:-b]
        means.append(np.mean(diff))
        medians.append(np.median(diff))

    X = range(len(means))
    if options.range:
        rng = eval(options.range)
        if len(rng) == len(means):
            X = rng

    plt.plot(X, means, label='mean')
    plt.ylabel("Mean Depth Error")
    if options.semilog:
        plt.semilogx()
    plt.savefig("mean_error.pdf")
    plt.savefig("mean_error.png")
    plt.figure()
    plt.plot(X, medians, label='median')
    plt.ylabel("Median Depth Error")
    if options.semilog:
        plt.semilogx()
    plt.savefig("median_error.pdf")
    plt.savefig("median_error.png")
    plt.show()


if __name__ == "__main__":
    main()

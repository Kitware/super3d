#!/usr/bin/python
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
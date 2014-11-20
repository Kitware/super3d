#!/usr/bin/python
# Copyright 2012, Kitware, Inc.
# See `LICENSE' file for licensing information.

"""This script converts EPFL camera files to the vpgl camera file format.

EPFL camera files are those included with the dataset found here:
http://cvlab.epfl.ch/~strecha/multiview/denseMVS.html
"""

import os
import glob
import re
from optparse import OptionParser

import numpy as np

from vidpy.io.camera import write_camera_fkrt_file


def read_epfl_camera(fname):
    data = np.fromfile(fname,sep=' ')
    # img_size = data[-2:]
    data = data[:-2].reshape((8,3))
    K = np.matrix(data[:3,:])
    R = np.matrix(data[4:7,:]).T
    t = -R * np.matrix(data[7,:]).T
    return (K, R, t)


def main():
    usage = "usage: %prog [options] input_glob output_file\n\n"
    usage += "  Convert all cameras in the input glob to a vpgl camera file.\n"
    parser = OptionParser(usage=usage)

    (options, args) = parser.parse_args()

    input_glob = args[0]
    output_file = args[1]

    in_files = sorted(glob.glob(input_glob))
    frame_nums = range(len(in_files))

    # find a numbers in the input files names
    fn_numbers = [re.findall('\d+', os.path.basename(name)) for name in in_files]

    # if each filename has a unique number, use these instead
    if all(fn_numbers):
        nums = [int(L[0]) for L in fn_numbers]
        # check that all numbers are unique
        if len(set(nums)) == len(nums):
            frame_nums = nums

    cameras = [read_epfl_camera(fname) for fname in in_files]

    write_camera_fkrt_file(dict(zip(frame_nums,cameras)), output_file)


if __name__ == "__main__":
    main()

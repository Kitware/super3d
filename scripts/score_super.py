# -*- coding: utf-8 -*-
"""
Created on Tue Jul 02 09:11:00 2013

@author: eric.smith
"""

import Image
import numpy as np
import vidpy_image_metric


def main():
    mfavg = np.array(Image.open("C:/Dev/super3d/bin/tools/mfavg.png"))
    bicub = np.array(Image.open("C:/Dev/super3d/bin/tools/bicub.png"))
    original =  np.array(Image.open("C:/Dev/super3d/bin/tools/original.png"))
    superq =  np.array(Image.open("C:/Dev/super3d/bin/tools/super.png"))

    res = vidpy_image_metric.psnr(superq, original)
    print res
    res = vidpy_image_metric.psnr(mfavg, original)
    print res
    res = vidpy_image_metric.psnr(bicub, original)
    print res
    res = vidpy_image_metric.psnr(original, original)
    print res


if __name__ == '__main__':
    main()

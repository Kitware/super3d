import numpy as np
from PIL import Image
import matplotlib.pyplot as plt
import caffe
import os
import cv2

labelmap = { 0 : 255, 1 : 85, 2 : 43, 3 : 128, 4 : 170 }


imgdir = '/home/eric/tmpdata/telesculptor/seq2'
outdir = '/home/eric/tmpdata/telesculptor/seq2labels'

# load net
net = caffe.Net('deploy.prototxt', 'snapshot_iter_200000.caffemodel', caffe.TEST)

lut = [] 
for i in range(256):
    if (labelmap.has_key(i)):
        lut.append(labelmap[i])
    else:
        lut.append(len(labelmap))
                
for root, dirs, files in os.walk(os.path.abspath(imgdir)):
    for f in files:
        baselen = len(imgdir)
        basename = os.path.join(root[baselen+1:],f)
        imgname = os.path.join(imgdir,basename)       
        outname = os.path.join(outdir, os.path.splitext(basename)[0] + '-labels.png')
        
        print imgname
   
        # load image, switch to BGR, subtract mean, and make dims C x H x W for Caffe
        im = Image.open(imgname)
        in_ = np.array(im, dtype=np.float32)
        in_ = in_[:,:,::-1]
        in_ -= np.array((104.00698793,116.66876762,122.67891434))
        in_ = in_.transpose((2,0,1))

        net.blobs['data'].reshape(1, *in_.shape)
        net.blobs['data'].data[...] = in_
        net.forward()
        out = net.blobs['score-telesculptor'].data[0].argmax(axis=0)
                
        lut = np.array(lut,dtype='uint8')
        outimg = cv2.LUT(out.astype('uint8'), lut)
        cv2.imwrite(outname, outimg)        

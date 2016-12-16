# -*- coding: utf-8 -*-
"""
Created on Thu Jan 14 19:00:55 2016

@author: eric.smith
"""

import os
import numpy as np
import cv2
import random
import math
import lmdb
import sys
import caffe
import shutil

maxboxsize = 750        #largest scale
minboxsize = 300        #smallest scale
outputsize = 500        #training image size
numperimg = 25          #number of training samples per image
rotation = 180          #amount of rotation range in augmentation -rotation to rotation

#map the ground truth color values to labels
labelmap = { 255 : 0, 85 : 1, 43 : 2, 128 : 3, 213 : 4, 170 : 4 }



def add_to_db(img, label, txn_data, txn_label, outputdir, index):
        
    rt2 = math.sqrt(2.0)
    
    for i in range(numperimg):
       
        cropsize = random.randint(minboxsize, maxboxsize)
        boxsize = int(cropsize * rt2)
        ul = np.array([random.randint(0,img.shape[1]-boxsize-1), random.randint(0,img.shape[0]-boxsize-1)])
        lr = np.array([ul[0]+boxsize, ul[1]+boxsize])

        crop = img[ul[1]:lr[1],ul[0]:lr[0]]   
        labelcrop = label[ul[1]:lr[1],ul[0]:lr[0]]  
        
        rot = random.randrange(-rotation, rotation)
        halfboxsize = boxsize / 2.0
        A = cv2.getRotationMatrix2D( (halfboxsize, halfboxsize), rot, 1.0 )
        rcrop = cv2.warpAffine( crop, A, (boxsize,boxsize))
        rlabelcrop = cv2.warpAffine( labelcrop, A, (boxsize,boxsize))
        
        halfcropsize = cropsize / 2.0
        ul = np.array([int(halfboxsize - halfcropsize), int(halfboxsize - halfcropsize)])
        lr = np.array([int(halfboxsize + halfcropsize), int(halfboxsize + halfcropsize)])
        rcrop = rcrop[ul[1]:lr[1],ul[0]:lr[0]]   
        rlabelcrop = rlabelcrop[ul[1]:lr[1],ul[0]:lr[0]]
        
        crop = cv2.resize(rcrop,tuple([outputsize,outputsize]))
        labelcrop = cv2.resize(rlabelcrop,tuple([outputsize,outputsize]))    
        
        #cv2.imwrite(os.path.join(os.path.join(outputdir,'images'),"%06dcrop.png"%index), crop)
        #cv2.imwrite(os.path.join(os.path.join(outputdir,'labels'),"%06dlabel.png"%index), labelcrop)
        
        crop = crop.transpose((2,0,1))
        im_dat = caffe.io.array_to_datum(crop)        
        txn_data.put('{:0>10d}'.format(index), im_dat.SerializeToString())
        labelcrop = np.expand_dims(labelcrop, axis=0)   
        im_label = caffe.io.array_to_datum(labelcrop)        
        txn_label.put('{:0>10d}'.format(index), im_label.SerializeToString()) 
        
        index += 1
        
    return index
        
        
    
def create_db(imgdir, labeldir, outputdir, dbname):

    imgnames = []
    labelimgnames = []

    for root, dirs, files in os.walk(os.path.abspath(imgdir)):
        for f in files:
            baselen = len(imgdir)
            imgname = os.path.join(root[baselen+1:],f)
            imgnames.append(os.path.join(imgdir,imgname))
            labelimgnames.append(os.path.join(labeldir, os.path.splitext(imgname)[0] + '-labels.png'))

    lmdb_path_data = os.path.join(outputdir, dbname + '_data')
    lmdb_path_label = os.path.join(outputdir, dbname + '_label')
    
    if os.path.isdir(lmdb_path_data):
        shutil.rmtree(lmdb_path_data)      
    envdata = lmdb.open(lmdb_path_data, map_size=1e11) #100gb database max size
    if os.path.isdir(lmdb_path_label):
        shutil.rmtree(lmdb_path_label)      
    envlabel = lmdb.open(lmdb_path_label, map_size=1e11) #100gb database max size

    lut = [] 
    for i in range(256):
        if (labelmap.has_key(i)):
            lut.append(labelmap[i])
        else:
            lut.append(len(labelmap))
    lut = np.array(lut,dtype='uint8')
    
    numimages = len(imgnames)
    
    index = 0
    
    imgoutdir = ""#/home/eric/tmpdata/telesculptor/output'
    if (not os.path.isdir(os.path.join(imgoutdir,'images'))):
        os.makedirs(os.path.join(imgoutdir,'images'))
    if (not os.path.isdir(os.path.join(imgoutdir,'labels'))):
        os.makedirs(os.path.join(imgoutdir,'labels'))
    
    with envdata.begin(write=True) as txn_data:  
        with envlabel.begin(write=True) as txn_label:
            for imgname,labelname in zip(imgnames,labelimgnames):
                print imgname
                img = cv2.imread(imgname)
                label = cv2.imread(labelname, cv2.IMREAD_GRAYSCALE)
                labelindex = cv2.LUT(label, lut)
                index = add_to_db(img, labelindex, txn_data, txn_label, imgoutdir, index)
    
    envdata.close()
    envlabel.close()
    

 
        
if __name__ == "__main__":    
    #imagedir labeldir outputdir databasename
    create_db(r'/home/eric/tmpdata/telesculptor/train/frames', r'/home/eric/tmpdata/telesculptor/train/labels',  r'/home/eric/tmpdata/telesculptor', 'train_lmdb')
    create_db(r'/home/eric/tmpdata/telesculptor/val/frames', r'/home/eric/tmpdata/telesculptor/val/labels',  r'/home/eric/tmpdata/telesculptor', 'val_lmdb')
    

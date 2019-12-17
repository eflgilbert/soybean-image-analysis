# Project: python_soybean_pixel_extraction.py
# Author: Erin Gilbert
# Created: ~Sept 2018
# Updated: Dec 2019

#Usage:  python_soybean_pixel_extraction.py rgbimage naive_bayes_pdfs.txt
#Purpose: Using results of a supervised machine learning PlantCV run to create a mask of an image with black being background (non-plant) and white being plant. These masks are then used to generate QCC images to judge the applicability of the PlantCV machine learning on our soybean plants.


############################################
#####    Define Environment/Packages  ######
############################################

#!/usr/bin/python
import sys, traceback
import cv2
import numpy as np
import argparse
import string
sys.path.append("/Users/gilbe952/Documents/plantcv")
import plantcv
from plantcv import plantcv as pcv
import PIL


#global variables
device = 0
debug = 'None'  


############################################
######    Load All the Data    #############
############################################
# Read in a color image

rgb_img, path, filename =pcv.readimage(sys.argv[1])
name = filename

device, mask = pcv.naive_bayes_classifier(rgb_img, pdf_file=sys.argv[2], device=0, debug="print") #plantcv model output
mask_image_plant, path, filename = pcv.readimage('1_naive_bayes_Plant_mask.jpg')




############################################
######    Perform Calculations    ##########
############################################

# combine leaf and labels
mask_image_label, path, filename = pcv.readimage('1_naive_bayes_labels_mask.jpg')
mask_image = mask_image_plant + mask_image_label
device, mask_image = pcv.rgb2gray_lab(mask_image, 'l', device)


#clean the mask up
device, img_binary = pcv.binary_threshold(mask_image, 50, 255, 'light', device)
pcv.print_image(img_binary, 'img_binary.tif')
device, blur_img = pcv.erode(img_binary, 3, 1, device, debug='print') # Erode to remove soil and Dilate so that you don't lose leaves (just in case)

mask = np.copy(blur_img)
device, fill_image = pcv.fill(blur_img, mask, 100, device)
pcv.print_image(fill_image, 'fill_image.tif')
device, binary_image = pcv.median_blur(fill_image, 1, device)
pcv.print_image(binary_image, 'binary_image.tif')
device, masked_image = device, dilate_image = pcv.dilate(fill_image, 3, 3, device)

############################################
###########    Create Output   #############
############################################

#Print grid of images for QC
# from https://stackoverflow.com/questions/30227466/combine-several-images-horizontally-with-python

#Different steps created here
list_masked = ['mask_image.tif', 'img_binary.tif', 'blur_img.tif',  'fill_image.tif', 'binary_image.tif', outfile]
imgs    = [ PIL.Image.open(i) for i in list_masked ]
# pick the image which is the smallest, and resize the others to match it (can be arbitrary image shape here)
min_shape = sorted( [(np.sum(i.size), i.size ) for i in imgs])[0][1]
imgs_comb = np.hstack( (np.asarray( i.resize(min_shape) ) for i in imgs ) )

# for a vertical stacking it is simple: use vstack
imgs_comb = np.vstack( (np.asarray( i.resize(min_shape) ) for i in imgs ) )
imgs_comb = PIL.Image.fromarray( imgs_comb)
imgs_comb.save( "./"+str(name[:-4]) + '_steps.tif' )

#masks created here
list_masked = ['1_naive_bayes_Plant_mask.jpg', '1_naive_bayes_Soil_mask.jpg', '1_naive_bayes_Background_mask.jpg','1_naive_bayes_Plastic_mask.jpg', '1_naive_bayes_labels_mask.jpg']
imgs    = [ PIL.Image.open(i) for i in list_masked ]

# pick the image which is the smallest, and resize the others to match it (can be arbitrary image shape here)
min_shape = sorted( [(np.sum(i.size), i.size ) for i in imgs])[0][1]
imgs_comb = np.hstack( (np.asarray( i.resize(min_shape) ) for i in imgs ) )

# save that beautiful mask picture
imgs_comb = PIL.Image.fromarray( imgs_comb)
imgs_comb.save( 'masked_temp.jpg' )    

# for a vertical mask stacking it is simple: use vstack
imgs_comb = np.vstack( (np.asarray( i.resize(min_shape) ) for i in imgs ) )
imgs_comb = PIL.Image.fromarray( imgs_comb)
imgs_comb.save( "./"+str(name[:-4]) + 'masks.tif' )

#Combined Mask/RGB created here
list_masked = [sys.argv[1], outfile1]
imgs    = [ PIL.Image.open(i) for i in list_masked ]

# pick the image which is the smallest, and resize the others to match it (can be arbitrary image shape here)
min_shape = sorted( [(np.sum(i.size), i.size ) for i in imgs])[0][1]
imgs_comb = np.hstack( (np.asarray( i.resize(min_shape) ) for i in imgs ) )

# for a vertical stacking it is simple: use vstack
imgs_comb = np.vstack( (np.asarray( i.resize(min_shape) ) for i in imgs ) )
imgs_comb = PIL.Image.fromarray( imgs_comb)
imgs_comb.save( "./"+str(name[:-4]) + 'results.tif' )

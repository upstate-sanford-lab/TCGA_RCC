import sys
import os
import numpy as np
from PIL import Image
from PIL import ImageOps
import cv2
import openslide
from openslide import open_slide
from openslide.deepzoom import DeepZoomGenerator
import xml.etree.ElementTree as ET
from xml.dom import minidom
import pandas as pd
from skimage import draw
import matplotlib.pyplot as plt
from shapely.geometry import Polygon, Point, MultiPoint, MultiPolygon, box
import argparse
from utils import *

class makeMap:

    def __init__(self):
        self.image_locations = '/media/tom/WD_BLACK/TCGA_RCC_sorted'
        self.has_xml = False
        self.save_location = '/data/train_val_test/heatmaps'
        self.master_pred = '/data/train_val_test/predictions/val_results/02172020-1738/master.csv'
        self.coded_filenames = False # for TCGA work, this is FALSE because data was already anon

        # this is not working right now, will produce ground truth map of annotations
        self.stride_ratio = 1
        self.high_risk_class = 'Dead' #this is whatever class you are wanting to make a map of
        self.writemask = False
        self.simplifyroi = True  # default = False, set to true for complex shapes
        self.nolabel = False  # if all regions in an annotation file belong to the same class, they are labeled as 'tumor'
        #      nolabel=FALSE should only be used if the "Text" attribute in xml corresponds to label

    def parsePredictions(self):
        all_preds = self.read_preds()
        imglist = findImageFiles(image_locations=self.image_locations, has_xml=self.has_xml)

        #find the base layer and define output map from here
        res_list=list(all_preds.keys())
        res_list_f=[float(s[:-1]) for s in all_preds.keys()]
        res_index=res_list_f.index(max(res_list_f))
        lowest_mag = res_list[res_index] #every base mag prediction will fill ONE voxel
        lowest_preds = all_preds[lowest_mag]
        other_levels = dict(filter(lambda elem: elem[0] != lowest_mag, all_preds.items()))

        for pt in list(set(lowest_preds['patient'])):
            pt_preds = lowest_preds[lowest_preds['patient'] == pt]

            for block in list(set(pt_preds['block'])):
                block_preds = pt_preds[pt_preds['block'] == block]
                image_name = block # THIS SHOULD BE CHANGED DEPENDING ON HOW YOU DID THINGS

                #find actual svs image file
                # COMMENT 2-4-2020 YOU NEED TO CHANGE THE UTILS TO FIND FILE NAMES THAT MATCH BLOCK ID
                block_img = imglist[imglist['savename']==image_name]
                # find the image size that was actually extracted on prediction files
                phys_size = max([int(i) for i in list(set(block_preds['size']))[0].split('-')])
                # now we utilize the original image to create a probability map
                # first grab data from digital header
                oslide = openslide.OpenSlide(list(block_img['image'])[0])
                # this is physical microns per pixel
                acq_mag = 10.0 / float(oslide.properties[openslide.PROPERTY_NAME_MPP_X])
                # this is nearest multiple of 20 for base layer
                base_mag = int(20 * round(float(acq_mag) / 20))
                base_dim = oslide.dimensions
                print(base_dim)

                base_img = self.make_level_img(base_dim=base_dim, phys_size=phys_size, block_preds=block_preds, image_name=image_name, level_name = lowest_mag,lvl_size = phys_size,pt=pt)

                for lvl in list(other_levels.keys()):
                    lvl_preds = other_levels[lvl]
                    lvl_pt = lvl_preds[lvl_preds['patient'] == pt]
                    lvl_block = lvl_pt[lvl_pt['block'] == block]
                    lvl_size = max([int(i) for i in list(set(lvl_preds['size']))[0].split('-')])

                    lvl_img = self.make_level_img(base_dim=base_dim, phys_size=phys_size, block_preds=lvl_block, image_name=image_name, level_name = lvl,lvl_size = lvl_size,pt=pt)
                    base_img +=lvl_img

                base_img = base_img/all_preds.__len__()

                slideimg = Image.fromarray(np.uint8(base_img * 255))
                slideimg = slideimg.convert('L')
                slideimg.save(os.path.join(self.save_location,pt, image_name + '_avg.jpeg'))
        return all_preds, imglist

    def read_preds(self):
        ''' create a dictionary of predictions from all magnifications provided in master csv file '''
        df_csv = pd.read_csv(self.master_pred,sep=',',header=0)
        all_preds = {}
        for index, pred_i in df_csv.iterrows():
            mag = pred_i['mag']
            df_i = pd.read_csv(pred_i['file'],sep=',',header=0)
            # YOU WILL NEED TO CHANGE THIS IF YOU CHANGE YOUR NAMING CONVENTION
            #df_i['img_name'], df_i['mag'], df_i['loc'], df_i['size'], df_i['ws'], df_i['label'] = df_i['img_name'].str.split('_', expand=True)
            new_df = df_i['img_name'].str.split('_', expand=True)
            patient_df= pd.DataFrame(df_i['img_name'].str[0:12])
            patient_df = patient_df.rename(columns={"img_name": 'patient'})
            new_df = new_df.rename(columns={0: 'block', 1: 'mag', 2: 'loc', 3: 'size', 4: 'ws', 5: 'label'})
            df_i = pd.merge(df_i, patient_df, left_index=True, right_index=True)
            df_i = pd.merge(df_i,new_df,left_index=True,right_index=True)
            all_preds[mag] = df_i
        return all_preds

    def make_level_img(self,base_dim,phys_size,block_preds,image_name,level_name,lvl_size,pt):
        x = np.zeros((round(base_dim[1]*self.stride_ratio/phys_size),round(base_dim[0]*self.stride_ratio/phys_size)), np.float)
        x_count = np.ones((round(base_dim[1]*self.stride_ratio/phys_size),round(base_dim[0]*self.stride_ratio/phys_size)), np.float)
        x_mask = np.zeros((round(base_dim[1] * self.stride_ratio / phys_size),round(base_dim[0] * self.stride_ratio / phys_size)), np.float)

        for index, patch in block_preds.iterrows():
            patch_start = [round(int(i)*self.stride_ratio/phys_size) for i in patch['loc'].split('-')]
            x_count[patch_start[1]:patch_start[1]+self.stride_ratio*round(lvl_size/phys_size),patch_start[0]:patch_start[0]+self.stride_ratio*round(lvl_size/phys_size)] += 1
            x_mask[patch_start[1]:patch_start[1] + self.stride_ratio*round(lvl_size/phys_size),patch_start[0]:patch_start[0] + self.stride_ratio*round(lvl_size/phys_size)] = 1
            x[patch_start[1]:patch_start[1]+self.stride_ratio*round(lvl_size/phys_size),patch_start[0]:patch_start[0]+self.stride_ratio*round(lvl_size/phys_size)] += patch[self.high_risk_class]

        #setup file structure for saving
        if not os.path.exists(os.path.join(self.save_location,pt)): os.mkdir(os.path.join(self.save_location,pt))

        x = x/x_count
        slideimg = Image.fromarray(np.uint8(x*255))
        slideimg = slideimg.convert('L')
        slideimg.save(os.path.join(self.save_location,pt,image_name+'_'+str(level_name)+'.jpeg'))
        slidemask = Image.fromarray(np.uint8(x_mask*255))
        slidemask = slidemask.convert('L')
        slidemask.save(os.path.join(self.save_location,pt,image_name+'_'+str(level_name)+'_mask.jpeg'))
        return x

if __name__ == '__main__':
    # parser = argparse.ArgumentParser()
    # parser.add_argument('--image_locations')
    # parser.add_argument('--save_location')
    # parser.add_argument('--has_xml')
    # parser.add_argument('--master_pred')
    # parser.add_argument('--coded_filenames')
    # args = parser.parse_args()
    c = makeMap()
    all_preds, imglist = c.parsePredictions()
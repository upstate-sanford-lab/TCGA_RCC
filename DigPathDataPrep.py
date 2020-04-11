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
from shapely.geometry import Polygon, Point, MultiPoint, MultiPolygon
from wsi_patchextractor import *
from img_decon import *
import random

class ProcessTCGA(ExtractPatches):

    def __init__(self):
        self.path_to_files='/data/TCGA_RCC_sorted'  # path to directory subdirectories containing .svs files
        self.save='/data'                           # root directory for saving
        self.save_n='T3_T4_deconv'                 # name to save files
        self.clinical_data='/data/TCGA_RCC_clinical_data/clinical_data_T3_T4.csv'
        self.mag_extract = [2.5,5,10,20]  # specify which magnifications you wish to pull images from
        self.save_image_size = 300  # specify image size to be saved (note this is the same for all magnifications)
        self.pixel_overlap = 0  # specify the level of pixel overlap in your saved images

    def make_patches_all_files(self):
        '''save patches for all files '''

        c=ExtractPatches() #instantiate this class
        vs_dict=vital_status_d(cdata=self.clinical_data) #from method below

        #patients without patches
        pts_need_patches=find_pts_no_patches()

        problem_list=[]
        for dir in pts_need_patches:
            print("processing files for patient {}".format(dir))
            c.file_location=os.path.join(self.path_to_files,dir)
            img_list=[file for file in os.listdir(os.path.join(self.path_to_files,dir)) if file.endswith('.svs')]
            print('total of {} images for this patient'.format(len(img_list)))
            xml_list = [file for file in os.listdir(os.path.join(self.path_to_files, dir)) if file.endswith('.xml')]
            xml_list_name=[file.split('.xml')[0] for file in xml_list]
            problem_list=[]
            for img in img_list:
                print("processing image {}".format(img))
                img_name='-'.join(img.split('-')[0:3])
                c.image_file=img; c.save_name=img_name  #assign __init__ methods to save appropriately
                c.mag_extract=self.mag_extract; c.save_image_size=self.save_image_size #re-assign in case
                vs=vs_dict[img_name] #extract the vital status from dictionary

                #set up file structure and filenames
                if not os.path.exists(os.path.join(self.save,self.save_n)): os.mkdir(os.path.join(self.save,self.save_n))
                if not os.path.exists(os.path.join(self.save, self.save_n,img_name+'_'+vs)): os.mkdir(os.path.join(self.save, self.save_n,img_name+'_'+vs))
                if not os.path.exists(os.path.join(self.save,self.save_n,img_name+'_'+vs,img.split('.svs')[0]+'_'+vs)): os.mkdir(os.path.join(self.save,self.save_n,img_name+'_'+vs,img.split('.svs')[0]+'_'+vs))
                if not os.path.exists(os.path.join(self.save,self.save_n,img_name+'_'+vs,img.split('.svs')[0]+'_'+vs,'patches')): os.mkdir(os.path.join(self.save,self.save_n,img_name+'_'+vs,img.split('.svs')[0]+'_'+vs,'patches'))

                c.save_location=os.path.join(self.save,self.save_n,img_name+'_'+vs,img.split('.svs')[0]+'_'+vs,'patches')
                if img.split('.svs')[0] in xml_list_name: c.xml_file=img.split('.svs')[0]+'.xml'

                #this saves the files
                try:
                    c.parseMeta_and_pullTiles()
                except:
                    print('problem processing image {}'.format(img))
                    problem_list+=[img]
        print("list of problem images: {}".format(problem_list))


    def deconv(self,ws_cutoff=0.8):
        '''loads patches, performs deconvolution, and saves devonvoled patches
        :param ws_cutoff - amount of white space willing to tolerate (0.8 = OK to use images with up to 80% ws value)
        '''
        p_path=os.path.join(self.save,self.save_n)
        for patient in sorted(find_pts_no_deconv(p_path)):
            print("performing deconvolution for patient {}".format(patient))
            for img in os.listdir(os.path.join(p_path,patient)):
                for res in os.listdir(os.path.join(p_path,patient,img,'patches')):

                    #set up file structure for saving
                    if not os.path.exists(os.path.join(p_path, patient, img, 'deconv')): os.mkdir(os.path.join(p_path, patient, img, 'deconv'))
                    for type in ['_norm','_bw','_decon']:
                        if not os.path.exists(os.path.join(p_path, patient, img, 'deconv', res + type)):    os.mkdir(os.path.join(p_path, patient, img, 'deconv', res + type))

                    #loop over patches and perform deconvolution
                    for patch in os.listdir(os.path.join(p_path,patient,img,'patches',res)):
                        img_loc=os.path.join(p_path,patient,img,'patches',res)
                        save_loc=os.path.join(p_path,patient,img,'deconv',res)
                        if float(patch.split('_')[4].split('-')[1])<=ws_cutoff:
                            try:
                                colornorm(image_location=img_loc, save_location=save_loc, img_file=patch)
                            except:
                                print("problem processing image {}".format(patch))
                                problem_img+=[patch]
        print("problems with the following image".format(problem_img))


### helper functions ###
def vital_status_d(cdata):
    '''develop dictionary of TCGA_ID:vital status'''
    data=pd.read_csv(cdata)
    vs_dict={}
    for id in data['submitter_id']: vs_dict[id]=data.loc[data['submitter_id']==id,'vital_status'].tolist()[0]
    return vs_dict

def find_pts_no_patches():
    '''return list of all patients without patch files'''

    path_to_patches='/data/T3_T4_patches'
    path_to_files='/data/TCGA_RCC_sorted'
    patches_files=[file.split('_')[0] for file in os.listdir(path_to_patches)]
    all_files=os.listdir(path_to_files)
    files=list(set(all_files)-set(patches_files))
    return files

def find_pts_no_deconv(path='/data/T3_T4_deconv'):
    '''check for patients without conv files'''
    num=0
    pts=[]
    for pt in os.listdir(path):
        for slide in os.listdir(os.path.join(path,pt)):
            if not os.path.exists(os.path.join(path,pt,slide,'deconv')):
                num+=1
                pts+=[pt]
    print("total {} patient no deconv in this folder".format(num))
    return(pts)


if __name__=='__main__':
    c=ProcessTCGA()
    c.deconv()
import os
import numpy as np
import pandas as pd
import random
import shutil

class SplitData:

    random.seed(0)

    def __init__(self):
        self.path_to_files='/data/T3_T4_deconv'  # path to patches
        self.save='/data'          # root directory for saving
        self.save_n='train_val_test'                  #
        self.train_val_test='train_val_test'


    def split_ds_by_pt(self,filetype='decon', val=0.25, test=0.25):
        '''splits the data into train/val/test datasets on the patient level but saves each slide.
        :param val (float) percentage of patients in validation dataset
        :param test (float) percentage of patients in test dataset
        '''
        # setup saving path
        path, dir = os.path.split(self.path_to_files)
        for pname in ['train_val_test', 'train_val_test/train_sl', 'train_val_test/val_sl', 'train_val_test/test_sl']:
            if not os.path.exists(os.path.join(path, pname)): os.mkdir(os.path.join(path, pname))

        # split data by patient an copy data
        pt_list = os.listdir(self.path_to_files)
        test = random.sample(pt_list, int((len(pt_list) * test)))
        train_val = list(set(pt_list) - set(test))
        val = random.sample(train_val, int(len(pt_list) * val))
        train = list(set(pt_list) - set(test + val))

        #save data into a dictionary and check size of each dataset
        dict = {'train': train, 'val': val, 'test': test}
        for sz in dict.keys():
            print("size of {} dataset is {}".format(sz, len(dict[sz])))

        # save data into appropriate dataset on the slice level
        for ds in dict.keys():
            ds_pts = dict[ds]
            for patient in sorted(ds_pts):
                print('copying files for patient {} to folder {}'.format(patient, ds))
                for slide in os.listdir(os.path.join(self.path_to_files, patient)):
                    for file in os.listdir(os.path.join(self.path_to_files, patient,slide,'deconv')):
                        if file.split('_')[1]==filetype:
                            i_path = os.path.join(self.path_to_files, patient,slide,'deconv',file)
                            f_path = os.path.join(path, self.train_val_test, ds + '_sl', slide,file)
                            shutil.copytree(i_path, f_path)


    def make_training_ds(self, filetype='decon', res=['2.5x', '5x', '10x', '20x']):
        '''separate data by training and validation datasets'''

        logger=[]
        # setup file structure to save to
        for r in res:
            if not os.path.exists(os.path.join(self.save, self.train_val_test, r)):
                os.mkdir(os.path.join(self.save, self.train_val_test, r))
            for set in ['train_sl', 'val_sl']:
                    if not os.path.exists(os.path.join(self.save, self.train_val_test,r, set.split('_')[0])):
                        os.mkdir(os.path.join(self.save, self.train_val_test,r, set.split('_')[0]))
                    for v in ['Alive', 'Dead']:
                        if not os.path.exists(os.path.join(self.save, self.train_val_test, r, set.split('_')[0], v)):
                            os.mkdir(os.path.join(self.save, self.train_val_test, r, set.split('_')[0], v))

        for ds in ['train_sl','val_sl']:
            for slide in sorted(os.listdir(os.path.join(self.save, self.train_val_test, ds))):
                for vital_status in ['Alive', 'Dead']:
                    if slide.split('_')[1] == vital_status:
                        print("transferring files for slide {}".format(slide))
                        for r in res:
                            try:
                                for file in os.listdir(os.path.join(self.save, self.train_val_test, ds, slide,r + '_' + filetype)):
                                    if file.endswith('.png'):
                                            orig_path = os.path.join(self.save, self.train_val_test, ds, slide,r + '_' + filetype, file)
                                            final_path = os.path.join(self.save, self.train_val_test,r, ds.split('_')[0],vital_status,file)
                                            shutil.copy2(orig_path, final_path)
                            except:
                                print("problem with resolution {}".format(r))


if __name__=='__main__':
    c=SplitData()
    #c.split_ds_by_pt()
    c.make_training_ds()
    #c.make_harmon_dataset()

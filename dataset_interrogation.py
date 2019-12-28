import os
import numpy as np
import pandas as pd


class DatasetInterrogate:

    def __init__(self):
        self.basePATH='/home/tom/Desktop/revision_analysis2/model_dev_indvPIRADS'


    def slice_number_calculation(self,file='train',subfile='PIRADS_2',ext='.jpg'):
        '''recursively looks in folder for files with specific file extension'''

        num=0
        for root, dirnames, filenames in os.walk(os.path.join(self.basePATH,file,subfile)):
            for filename in filenames:
                if filename.endswith(ext):
                    num+=1
        print("total of {} files in the directory {} ".format(num,file))


if __name__=='__main__':
    c=DatasetInterrogate()
    c.slice_number_calculation()
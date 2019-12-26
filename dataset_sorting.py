import os
import shutil
import numpy as np
import pandas as pd



class SortFiles:
    '''class to sort through all the files downloaded from TCIA and place them into unique patient directories'''

    def __init__(self):
        self.basePATH='/Users/sanforth/Desktop/TCGA' #the directory with all the downloaded TCGA files
        self.outPATH=os.path.join(self.basePATH,'sorted')

    def sort(self):
        '''goes through each of the directories and sorts all the files into a unique directory'''

        T3T4_list=self.T3_T4_by_pt()

        for dir in os.listdir(self.basePATH):
            if dir.startswith('.'):  # stoopid hidden files
                continue
            files=os.listdir(os.path.join(self.basePATH,dir))
            for file in files:
                if file.startswith('.'): #stoopid hidden files
                    continue
                if file.endswith('.svs'):
                    filename=file.split('.')[0]
                    pt_e_l=filename.split('-')[0:3]
                    pt_name='-'.join(pt_e_l)
                    if pt_name in T3T4_list:
                        if not os.path.exists(self.outPATH):
                            os.mkdir(os.path.join(self.outPATH))
                        if not os.path.exists(os.path.join(self.outPATH,pt_name)):
                            os.mkdir(os.path.join(self.outPATH,pt_name))
                        shutil.copy2(os.path.join(self.basePATH,dir,file),os.path.join(self.outPATH,pt_name,filename+'.svs'))

    def T3_T4_by_pt(self):
        '''select patients with T3/T4.  Expects a folder 'clinical_data in same folder'''

        data=pd.read_csv(os.path.join(self.basePATH,'clinical_data','clinical.csv'))
        include_T=['T3','T4','T3a','T3b','T3c']
        sel_by_T=data[data['ajcc_pathologic_t'].isin(include_T)]
        sel_by_T_u=sel_by_T.drop_duplicates('submitter_id')
        sel_by_T_u.to_csv(os.path.join(self.basePATH,'clinical_data','clinical_data_T3_T4.csv'))
        return sorted(sel_by_T_u['submitter_id'].tolist())

    def dataset_characteristics(self):
        '''evaluate the characteristics of the dataset'''

        data = pd.read_csv(os.path.join(self.basePATH, 'clinical_data', 'clinical_data_T3_T4.csv'))
        vital_status=data['vital_status'].value_counts()
        gender_status=data['gender'].value_counts()
        age=data['age_at_index']
        d_to_death=data['days_to_death'].dropna()

        print("total number of patients: {}".format(data.shape[0]))
        print("Patients Died: {} ({})%".format(vital_status[0],round(vital_status[0] / (vital_status[0] + vital_status[1]) * 100), 2))
        print("Patients Alive: {} ({})%".format(vital_status[1],round(vital_status[1] / (vital_status[0] + vital_status[1]) * 100), 2))
        print("Patients Men: {} ({})%".format(vital_status[0],round(gender_status[0] / (gender_status[0] + gender_status[1]) * 100), 2))
        print("Patients Women: {} ({})%".format(vital_status[1],round(gender_status[1] / (gender_status[0] + gender_status[1]) * 100), 2))
        print("Median Age: {}, Range {}-{}".format(age.median(),age.min(),age.max()))
        print("Median Days To Death: {}, Range {}-{}".format(d_to_death.median(),d_to_death.min(),d_to_death.max()))



if __name__=='__main__':
    c=SortFiles()
    c.dataset_characteristics()


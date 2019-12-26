import os
import shutil
import numpy as np
import pandas as pd
import statistics


class SortFiles:
    '''class to sort through all the files downloaded from TCIA and place them into unique patient directories'''

    def __init__(self):
        self.basePATH='/home/tom/Desktop/TCGA_RCC' #the directory with all the downloaded TCGA files
        self.outPATH=os.path.join(self.basePATH,'TCGA_RCC_sorted')

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

        data=pd.read_csv(os.path.join(self.basePATH,'TCGA_RCC_clinical_data','clinical.csv'))
        include_T=['T3','T4','T3a','T3b','T3c']
        sel_by_T=data[data['ajcc_pathologic_t'].isin(include_T)]
        sel_by_T_u=sel_by_T.drop_duplicates('submitter_id')
        sel_by_T_u.to_csv(os.path.join(self.basePATH,'TCGA_RCC_clinical_data','clinical_data_T3_T4.csv'))
        return sorted(sel_by_T_u['submitter_id'].tolist())

    def dataset_characteristics(self):
        '''evaluate the characteristics of the dataset'''

        data_all=pd.read_csv(os.path.join(self.basePATH,'TCGA_RCC_clinical_data', 'clinical.csv'))
        data_all_u= data_all.drop_duplicates('submitter_id')
        data = pd.read_csv(os.path.join(self.basePATH,'TCGA_RCC_clinical_data', 'clinical_data_T3_T4.csv'))
        list_downoaded=os.listdir(self.outPATH)
        slide_list=[len(os.listdir(os.path.join(self.outPATH,file))) for file in os.listdir(self.outPATH)]
        data=data[data['submitter_id'].isin(list_downoaded)]
        vital_status=data['vital_status'].value_counts()
        gender_status=data['gender'].value_counts()
        age=data['age_at_index']
        num_patients_w_days_to_death=len(data['days_to_death'].dropna().tolist())
        d_to_death=data['days_to_death'].dropna()
        num_patients_w_fu_len=len(data['days_to_last_follow_up'].dropna().tolist())
        d_to_last_fu=data['days_to_last_follow_up'].dropna()

        print("TCGA Clear Cell RCC Cohort")
        print("-------------")
        print("Total number of patients in TCGA cohort with ccRCC {}".format(data_all_u.shape[0]))
        print("Total number of patients with ccRCC and T3/T4 tumors: {}, {}% of total cohort".format(data.shape[0],round((data.shape[0]/data_all_u.shape[0])*100,0)))
        print("-------")
        print("Description of T3/T4 cohort")
        print("---------")
        print("Total number of patients with ccRCC and T3/T4 tumors: {}".format(data.shape[0]))
        print('Total number of slides for this patient cohort is: {}'.format(sum(slide_list)))
        print("Median slides per patient: {}, range {}-{}".format (statistics.median(slide_list), min(slide_list),max(slide_list)))
        print("Patients died: {} ({})%".format(vital_status[0],round(vital_status[0] / (vital_status[0] + vital_status[1]) * 100), 2))
        print("Patients alive: {} ({})%".format(vital_status[1],round(vital_status[1] / (vital_status[0] + vital_status[1]) * 100), 2))
        print("Patients men: {} ({})%".format(vital_status[0],round(gender_status[0] / (gender_status[0] + gender_status[1]) * 100), 2))
        print("Patients women: {} ({})%".format(vital_status[1],round(gender_status[1] / (gender_status[0] + gender_status[1]) * 100), 2))
        print("Median age: {}, Range {}-{}".format(age.median(),age.min(),age.max()))
        print("Number of patients with days to death {}, {}%".format(num_patients_w_days_to_death,round((num_patients_w_days_to_death/data.shape[0])*100,0)))
        print("Median time to death: {} days ({} years), Range {}days-{}days".format(d_to_death.median(),round((d_to_death.median()/365),2),d_to_death.min(),d_to_death.max()))
        print("Number of patients with days to last fu {}, {}%".format(num_patients_w_fu_len,round((num_patients_w_fu_len/data.shape[0])*100,2)))
        print("Median follow up: {} days ({} years), Range {} days-{} days".format(d_to_last_fu.median(),round((d_to_last_fu.median()/365),2),d_to_last_fu.min(),d_to_last_fu.max()))



if __name__=='__main__':
    c=SortFiles()
    c.dataset_characteristics()

'''
this inferencer model is now flexible to accept classification
predictions from binary or multi-class predictions
it produces a csv with per-image predictions and softmax outputs
for every class, as well as overall class prediction
'''


from fastai.vision import *
import torch.nn.functional as F
import datetime

class Inferencer:

    def __init__(self):
        self.imagedir = '/data/train_val_test'
        self.outdir = '/data/train_val_test/predictions'
        self.res=['2.5x','5x','10x','20x']
        self.testPath = os.path.join(self.imagedir, 'test_sl')
        self.device=0

    def fastai_apply_model(self):
        '''
        applies model on patient level
        :param index(int): where to index in
        :return:
        '''

        #set device
        torch.cuda.set_device(self.device)

        #set time once
        time = str(datetime.datetime.now().strftime("%m%d%Y-%H%M"))

        # set up the output directory
        if not os.path.isdir(self.outdir):os.mkdir(self.outdir)
        if not os.path.isdir(os.path.join(self.outdir, 'val_results')): os.mkdir(os.path.join(self.outdir, 'val_results'))
        if not os.path.isdir(os.path.join(self.outdir, 'exported_models')): os.mkdir(os.path.join(self.outdir, 'exported_models'))


        for res in self.res:
            #load resoltion-specific model
            initial_filename = os.path.join(self.outdir, 'models', res,os.listdir(os.path.join(self.outdir, 'models', res))[0])
            final_filename = os.path.join(self.outdir, 'exported_models', 'export.pkl')
            shutil.copy2(initial_filename, final_filename)

            # load learner
            model_path = os.path.join(self.outdir, 'exported_models')
            learn = load_learner(model_path)

            df_out = pd.DataFrame()
            for patient in os.listdir(os.path.join(self.testPath)):
                print('applying model to the patient {} at resolution {}'.format(patient,res))
                for image in os.listdir(os.path.join(self.testPath, patient,res+'_decon')):
                    img_name=patient.split('_')[0]+'_'+'_'.join(image.split('_')[1:len(image.split('_'))])
                    img = open_image(os.path.join(self.testPath, patient,res+'_decon', image))
                    pred_class, pred_idx, outputs = learn.predict(img)
                    outputs=F.softmax(outputs,dim=-1)
                    all_classes = learn.data.classes
                    #t_df = pd.DataFrame([image, pred_class, outputs_np]).transpose()
                    df_i = pd.DataFrame([img_name, pred_class]).transpose()
                    for classi in range(0,len(all_classes)):
                        p_df = pd.DataFrame([outputs[classi]]).transpose()
                        df_i=pd.concat([df_i,p_df], axis=1)
                    df_out=pd.concat([df_out,df_i],axis=0)

            #write dataframe out to csv
            if not os.path.exists(os.path.join(self.outdir,'val_results',time)):
                os.mkdir(os.path.join(self.outdir,'val_results',time))
            if not os.path.exists(os.path.join(self.outdir,'val_results',time,res)):
                    os.mkdir(os.path.join(self.outdir,'val_results',time,res))
            df_out.columns=['img_name','class_pred']+all_classes
            df_out.to_csv(os.path.join(self.outdir,'val_results',time,res,res+'.csv'))

        #make master file
        df_master=pd.DataFrame()
        for res in ['2.5x','5x','10x','20x']:
            path=os.path.join(self.outdir, 'val_results', time, res, res + '.csv')
            df_l = pd.DataFrame([res,path]).transpose()
            df_master=pd.concat([df_master,df_l],axis=0)
        df_master.columns=['mag','file']
        df_master.to_csv(os.path.join(self.outdir,'val_results',time,'master.csv'))

if __name__ == '__main__':
    d = Inferencer()
    d.fastai_apply_model()

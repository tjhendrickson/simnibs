#!/usr/bin/env python3

# -*- coding: utf-8 -*-
"""
Created on Wed Dec  2 13:23:55 2020

To run subject specific SimNIBS E-field simulations based on 5x5 SMA grid

@author: timothy
"""

from simnibs import sim_struct, run_simnibs, opt_struct
import pandas as pd
from multiprocessing import Pool
import os
from bids import BIDSLayout
import nibabel as nib
import numpy as np
from glob import glob
import json
import pdb
import sys

class SimNIBS_simulation:
    def __init__(self,BIDS_folder,coil_name,derivatives_folder,low_mem,optimization,stimulation_type,subject,session='',bold_file=''):     
        '''
        Parameters
        ----------
        subject : string
            Inputted BIDS subject number without "sub-"
        BIDS_folder : string/path
            Path to BIDS foler
        derivatives_folder : string/path
            Path that will store TMS simulation outputs
        bold_file : string/path
            Path to functional bold output

        Returns
        -------
        1) SimNIBS headmodel for inputted subject
        2) 5x5 SMA grid CSV file witin subject native space
        3) NIFTI files for each of the 25 electric current simulations for subject specfic TMS targeting
        4) A JSON file with suggested target coordinates, and target correlation value 
        '''

        self.BIDS_folder = BIDS_folder
        self.bold_file = bold_file
        self.coil_name = coil_name
        self.derivatives_folder = derivatives_folder
        self.low_mem = low_mem
        self.optimization = optimization
        self.subject = subject
        self.stimulation_type = stimulation_type
        self.session = session
        
        
        self.run_headmodels()

        if self.optimization:
            if self.stimulation_type == 'TMS':
                pass #TMS optimization
            elif self.stimulation_type == 'tDCS':
                pass # tDCS Lead fields simulations

        #self.generate_SMA_grid()
        #if not os.path.isdir(os.path.join(self.derivatives_folder,'efield_simulations','sub-{subject}_ses-{session}'.format(subject=self.subject,session=self.session))):
        #    p=Pool(6)
        #    print('start simulation for subject {subject} and session {session}'.format(subject=self.subject,session=self.session))
        #    p.map(self.run_simulation,self.grid_positions_df['Target Name'].to_list())
        
        #self.generate_TMS_target()

    def run_headmodels(self):
        ''' Generate siminibs head models based on inputted data '''
        
        layout = BIDSLayout(self.BIDS_folder)
        # find T1w and T2w data
        if self.session == '':
            T1ws=[T1_file for T1_file in layout.get(subject=self.subject,suffix='T1w',return_type='file') if "nii.gz" in T1_file]
            T2ws=[T2_file for T2_file in layout.get(subject=self.subject,suffix='T2w',return_type='file') if "nii.gz" in T2_file]
        else:
            T1ws=[T1_file for T1_file in layout.get(subject=self.subject,session=self.session,suffix='T1w',return_type='file') if "nii.gz" in T1_file]
            T2ws=[T2_file for T2_file in layout.get(subject=self.subject,session=self.session,suffix='T2w',return_type='file') if "nii.gz" in T2_file]
        
        # restrict T1w and T2w use to last in last, if there is no T1w, continue to next iteration
        assert len(T1ws) > 0, 'Cannot run head model. No T1ws exist for subject {subject} and session {session}'.format(subject=self.subject,session=self.session)
        if len(T1ws) > 1:
            T1w=T1ws[-1]
        else:
            T1w=T1ws[0]
        if T2ws:
            if len(T2ws) > 1:
                T2w=T2ws[-1]
            else:
                T2w = T2ws[0]
        os.chdir(self.derivatives_folder)
        if self.session == '':
            self.simnibs_file_base = 'sub-{subject}'.format(subject=self.subject)
        else:
            self.simnibs_file_base = 'sub-{subject}_ses-{session}'.format(subject=self.subject,session=self.session)
        self.m2m_folder='m2m_'+ self.simnibs_file_base
        self.head_mesh = self.simnibs_file_base+'.msh'
        if not os.path.isdir(self.m2m_folder):
            if T2w:
                os.system('headreco all {file_base} {T1w} {T2w}'.format(file_base=self.simnibs_file_base,T1w=T1w,T2w=T2w))    
            else:
                os.system('headreco all {file_base} {T1w}'.format(file_base=self.simnibs_file_base,T1w=T1w))
    def run_tDCS_optimization(self):
        
        # Initialize structure
        tdcs_lf = sim_struct.TDCSLEADFIELD()
        
        # head mesh
        tdcs_lf.fnamehead = os.path.join(self.derivatives_folder,self.head_mesh)
        
        # Output folder
        os.makedirs(os.path.join(self.derivatives_folder,'tdcs_leadfield',self.simnibs_file_base))
        tdcs_lf.pathfem = os.path.join(self.derivatives_folder,'tdcs_leadfield',self.simnibs_file_base)

        tdcs_lf.open_in_gmsh = False # do not open efield simulations with gmsh
        tdcs_lf.map_to_vol = True

        if not self.low_mem:
            tdcs_lf.solver_options = 'pardiso' # This solver is faster than the default. However, it requires much more memory (~12 GB)
        
        run_simnibs(tdcs_lf)
    def run_TMS_optimization(self):
        # Initialize structure
        tms_opt = opt_struct.TMSoptimize()
        # Select the head mesh
        tms_opt.fnamehead = os.path.join(self.derivatives_folder,self.head_mesh)
        # Select output folder
        os.makedirs(os.path.join(self.derivatives_folder,'tms_optimization',self.simnibs_file_base))
        tms_opt.pathfem = os.path.join(self.derivatives_folder,'tms_optimization',self.simnibs_file_base)
        
        # Select the coil model
        tms_opt.fnamecoil = self.coil_name + '.nii.gz'
        # Select a target for the optimization
        #tms_opt.target = [-43.4, -20.7, 83.4]

        if not self.low_mem:
            tms_opt.solver_options = 'pardiso' # This solver is faster than the default. However, it requires much more memory (~12 GB)

        # Run optimization to get optimal coil position
        opt_pos=tms_opt.run()

# TODO: this may or may not be needed later
"""       
def run_simulation(self,grid_position):
    ''' run TMS SimNIBS electric current simulation based on inputted grid position '''
    
    # Initialize a session
    s = sim_struct.SESSION()
    
    # Name of head mesh
    s.fnamehead = os.path.join(self.derivatives_folder,'sub-{subject}_ses-{session}.msh'.format(subject=self.subject,session=self.session))
    
    s.open_in_gmsh = False # do not open efield simulations with gmsh
    s.map_to_vol = True
    # Output folder
    os.makedirs(os.path.join(self.derivatives_folder,'efield_simulations','sub-{subject}_ses-{session}'.format(subject=self.subject,session=self.session),grid_position))
    s.pathfem = os.path.join(self.derivatives_folder,'efield_simulations','sub-{subject}_ses-{session}'.format(subject=self.subject,session=self.session),grid_position)
    
    # add a TMSLIST to the SESSION
    tms = s.add_tmslist()
    #tms.fnamecoil='Magstim_70mm_Fig8.nii.gz'
    
    # add a new position
    pos=tms.add_position()
    pos.centre=[self.grid_positions_df[self.grid_positions_df['Target Name']==grid_position]['Loc. X'].to_numpy()[0],
                self.grid_positions_df[self.grid_positions_df['Target Name']==grid_position]['Loc. Y'].to_numpy()[0],
                self.grid_positions_df[self.grid_positions_df['Target Name']==grid_position]['Loc. Z'].to_numpy()[0]]
    grid_direction='_'.join(grid_position.split('_')[0:3])+'_0'
    pos.pos_ydir=[self.grid_directions_df[self.grid_directions_df['Target Name']==grid_direction]['Loc. X'].to_numpy()[0],
                    self.grid_directions_df[self.grid_directions_df['Target Name']==grid_direction]['Loc. Y'].to_numpy()[0],
                    self.grid_directions_df[self.grid_directions_df['Target Name']==grid_direction]['Loc. Z'].to_numpy()[0]]
    
    # run simulation
    run_simnibs(s)


    
def generate_SMA_grid(self):
    ''' generate subject specific SMA grid and direction coordinations if they do not exist already '''
    if not os.path.isfile('{derivatives_folder}/sub-{subject}_ses-{session}_SMA_Grid.csv'.format(derivatives_folder=self.derivatives_folder,subject=self.subject,session=self.session)):
        os.system('mni2subject_coords -m {derivatives_folder}/{m2m_folder} -s {derivatives_folder}/SMA_Grid-MNI.csv -o {derivatives_folder}/sub-{subject}_ses-{session}_SMA_Grid.csv'.format(subject=self.subject,session=self.session,m2m_folder=self.m2m_folder,derivatives_folder=self.derivatives_folder))
    self.grid_positions_df = pd.read_csv('{derivatives_folder}/sub-{subject}_ses-{session}_SMA_Grid.csv'.format(derivatives_folder=self.derivatives_folder,subject=self.subject,session=self.session))
    if not os.path.isfile('{derivatives_folder}/sub-{subject}_ses-{session}_SMA_Grid-dirs.csv'.format(derivatives_folder=self.derivatives_folder,subject=self.subject,session=self.session)):
        os.system('mni2subject_coords -m {derivatives_folder}/{m2m_folder} -s {derivatives_folder}/SMA_Grid-dirs-MNI.csv -o {derivatives_folder}/sub-{subject}_ses-{session}_SMA_Grid-dirs.csv'.format(subject=self.subject,session=self.session,m2m_folder=self.m2m_folder,derivatives_folder=self.derivatives_folder))
    self.grid_directions_df = pd.read_csv('{derivatives_folder}/sub-{subject}_ses-{session}_SMA_Grid-dirs.csv'.format(derivatives_folder=self.derivatives_folder,subject=self.subject,session=self.session))

def generate_TMS_target(self):
    ''' correlate each of the simulation outputs with the fMRI outputs to find optimal target '''
    
    TMS_target_output_dir=os.path.join(self.derivatives_folder,
                                'sub-'+str(self.subject),
                                'ses-'+str(self.session),'target_outputs')
    if not os.path.isdir(os.path.join(TMS_target_output_dir)):
        os.makedirs(TMS_target_output_dir)
    SMA_bold_file = os.path.join(self.derivatives_folder,
                                        'sub-'+str(self.subject),
                                        'ses-'+str(self.session),
                                        os.path.basename(self.bold_file).split('.nii.gz')[0]+"_ROI-SMA.nii.gz")
    #os.system('fslmaths {')
    print('fslmaths {input_one} -mul {input_two} {out}'.format(input_one=self.bold_file,                                                               
                                                                            input_two=os.path.join(self.derivatives_folder,
                                                                                                    'sub-'+str(self.subject),
                                                                                                    'ses-'+str(self.session),
                                                                                                    'harvardoxford-cortical_ROI-SMA_tpl-T1w.nii.gz'),
                                                                            out=SMA_bold_file))
    
    os.system('fslmaths {input_one} -mul {input_two} {out}'.format(input_one=self.bold_file,                                                               
                                                                            input_two=os.path.join(self.derivatives_folder,
                                                                                                    'sub-'+str(self.subject),
                                                                                                    'ses-'+str(self.session),
                                                                                                    'harvardoxford-cortical_ROI-SMA_tpl-T1w.nii.gz'),
                                                                            out=SMA_bold_file))
    try:
        fmri_img_data=nib.load(SMA_bold_file).get_fdata()
    except:
        pass # TODO: incorporate first level analysis script here
    fmri_img_data[fmri_img_data==0]=np.nan # replace all zeros with nan
    final_coords = ('SMA_0_0',0)
    json_data = {}
    pdb.set_trace()
    for efield_data in sorted(glob('/home/lnpi15-raid6/conelea-data/bids/316_CBIT/BIDS_output/derivatives/TMS_Targeting/efield_simulations/sub-10622_ses-56920/SMA_Grid*/subject_volumes/sub-10622_ses-56920_TMS_1-0001_Magstim_70mm_Fig8_nii_scalar_normE.nii.gz')):
    #for efield_data in sorted(glob(os.path.join(self.derivatives_folder,'efield_simulations',
    #                                            '{subj_id}_{ses_id}',
    #                                            'SMA_Grid_*','subject_volumes','{subj_id}_{ses_id}_TMS_1-0001_Magstim_70mm_Fig8_nii_scalar_normE.nii.gz'.format(subj_id='sub-'+str(self.subject),ses_id='ses-'+str(self.session))))):
        efield_data_dir = os.path.dirname(efield_data)
        SMA_loc = efield_data.split('/')[-3]
        efield_xfm_out_file=os.path.join(efield_data_dir,'{subj_id}_{ses_id}_TMS_1-0001_Magstim_70mm_Fig8_nii_scalar_normE.nii.gz'.format(subj_id='sub-'+str(self.subject),ses_id='ses-'+str(self.session)))
        os.system('fslmaths {input_one} -mul {input_two} {out}'.format(input_one=efield_data,                                                               
                                                                        input_two=os.path.join(self.derivatives_folder,'sub-'+str(self.subject),'ses-'+str(self.session),'harvardoxford-cortical_ROI-SMA_tpl-T1w.nii.gz'),
                                                                        out=efield_xfm_out_file))
        efield_img_data = nib.load(efield_xfm_out_file).get_fdata()
        efield_img_data[efield_img_data==0]=np.nan
        dataframe_dict = {'fmri_data':fmri_img_data.flatten(),'efield_data':efield_img_data.flatten()}
        dataframe = pd.DataFrame(data=dataframe_dict) 
        r_val = dataframe.corr(method='pearson')['fmri_data']['efield_data']
        if r_val > final_coords[1]:
            final_coords = (SMA_loc,r_val)
    json_data['Grid Number'] = final_coords[0]
    json_data['R-value'] = final_coords[1]
    print(os.path.join(self.derivatives_folder,'{subj_id}_{ses_id}_SMA_Grid.csv'.format(subj_id='sub-'+str(self.subject),ses_id='ses-'+str(self.session))))
    subj_SMA_coords_df = pd.read_csv(os.path.join(self.derivatives_folder,'{subj_id}_{ses_id}_SMA_Grid.csv'.format(subj_id='sub-'+str(self.subject),ses_id='ses-'+str(self.session))))
    json_data['Target X'] = subj_SMA_coords_df[subj_SMA_coords_df['Target Name']==final_coords[0]]['Loc. X'].to_numpy()[0]
    json_data['Target Y'] = subj_SMA_coords_df[subj_SMA_coords_df['Target Name']==final_coords[0]]['Loc. Y'].to_numpy()[0]            
    json_data['Target Z'] = subj_SMA_coords_df[subj_SMA_coords_df['Target Name']==final_coords[0]]['Loc. Z'].to_numpy()[0]
    with open(os.path.join(self.derivatives_folder,'sub-'+str(self.subject),'ses-'+str(self.session),'target_outputs','{subj_id}_{ses_id}_TMS_target.json'.format(subj_id='sub-'+str(self.subject),ses_id='ses-'+str(self.session))),'w') as json_file:
        json.dump(json_data,json_file,indent=4,sort_keys=True)
    #else:
    #    raise Exception('TMS Target has already been produced for sub-{subject}_ses-{session}. If you would like to re-run, delete folder {TMS_target_dir}, and start again.'.format(subject=self.subject,session=self.session,TMS_target_dir=TMS_target_output_dir))

if __name__ == '__main__':
    arg_list = sys.argv[1:]
    SimNIBS_preprocessing(subject=arg_list[0],session=arg_list[1],
                 BIDS_folder=arg_list[2],derivatives_folder=arg_list[3],
                 bold_file=arg_list[4])
"""    
    
            
            
            
                





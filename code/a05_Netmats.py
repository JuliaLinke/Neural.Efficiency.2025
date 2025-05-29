#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#the script assumes that there is a directory called lists in the rootdir
#that contains a list with the participant IDs (e.g., s24567)

import nibabel
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import imageio

rootdir = os.path.join('/','MyWorkingDirectory')
partial = True
remove_vols = 0
fd_thr = .5
cohort = 'Cohort1'                #Cohort1 or Cohort2
visit = '1'                      #can be 1 for pre-CBT and 2 for post-CBT

def loadGifti(fname, NonSteadyState=0, icres=7): # ============================
    gii = nibabel.load(fname)
    gii_data = [d.data[:,None] for d in gii.darrays]
    gii_data = np.concatenate(gii_data, axis=1).T
    nV = 4**icres*10 + 2
    gii_data = gii_data[:,0:nV]
    return gii_data[NonSteadyState:,:]

def loadNifti(fname, NonSteadyState=0): # =====================================
    n = nibabel.load(fname)
    naff   = n.affine
    img4d  = n.get_fdata();
    imgsiz = img4d.shape
    if len(imgsiz) == 4:
        img4d  = img4d[:,:,:,NonSteadyState:]
        imgsiz = img4d.shape
        img2d  = np.reshape(img4d, (np.prod(imgsiz[0:3]), imgsiz[-1]), order='F').T
    else:
        img2d  = np.reshape(img4d, (np.prod(imgsiz[0:3]), 1), order='F').T
    return img2d, imgsiz, naff

def unWrap(netmat, side='lower'): # ===========================================
    if side.lower() == 'lower':
        uw  = netmat[np.tril_indices(netmat.shape[0], k=-1, m=netmat.shape[1])]
    elif side.lower() == 'upper':
        uw  = netmat[np.triu_indices(netmat.shape[0], k= 1, m=netmat.shape[1])]
    elif side.lower() == 'both':
        uwl = netmat[np.tril_indices(netmat.shape[0], k=-1, m=netmat.shape[1])]
        uwu = netmat[np.triu_indices(netmat.shape[0], k= 1, m=netmat.shape[1])]
        uw  = np.concatenate((uwl, uwu))
    return uw

# ======== [ Main ] ===========================================================
if __name__ == '__main__':
    
    #We could run the next steps for all participants in th fmriprep directory like this:
    #listSubjects = [d for d in os.listdir(os.path.join(rootdir,'derivatives','fmriprep')) if (os.path.isdir(os.path.join(rootdir,'derivatives','fmriprep',d)) and d[0:4] == 'sub-')]
    #However, we run it for the participants in our list as not everyone has T2 data
    listSubjects = open(os.path.join(rootdir,'lists/{}_T{}.txt'.format(cohort, visit)),'r')
    listSubjects = listSubjects.read().split()
    listSubjects.sort()
    fmriprepdir  = os.path.join(rootdir,'derivatives/{}/fmriprep/'.format(cohort))
    confounddir  = os.path.join(rootdir,'derivatives/{}/RegressNuissance'.format(cohort))
    netmatsdir   = os.path.join(rootdir,'derivatives/{}/'.format(cohort),'Ses{}.Netmats_NoGlob/'.format(visit))
    
    # List of subcortical structures of interest
    # 26/58 = Accumbens
    # 18/54 = Amygdala
    # 11/50 = Caudate
    # 17/53 = Hippocampus
    # 13/52 = Pallidum
    # 12/51 = Putamen
    # 10/49 = Thalamus
    # 28/60 = Ventral DC
    aseg_list = {
            'L': [26, 18, 11, 17, 13, 12, 10, 28],
            'R': [58, 54, 50, 53, 52, 51, 49, 60]}
    
    # Load the parcellation "annot" files in fsaverage space
    annot = {}
    ctab  = {}
    names = {}
    [annot['L'], ctab['L'], names['L']] = nibabel.freesurfer.io.read_annot(os.path.join(rootdir ,'atlas', 'Schaefer2018', 'Parcellations', 'FreeSurfer5.3', 'fsaverage', 'label','lh.Schaefer2018_100Parcels_17Networks_order.annot'))
    [annot['R'], ctab['R'], names['R']] = nibabel.freesurfer.io.read_annot(os.path.join(rootdir ,'atlas', 'Schaefer2018', 'Parcellations', 'FreeSurfer5.3', 'fsaverage', 'label','rh.Schaefer2018_100Parcels_17Networks_order.annot'))
    
    # For each subject
    # Maybe we need to select a subset of participants
    # for sidx, subj in enumerate(listSubjects[40:]):
    # Right now we just take all participants
    for sidx, subj in enumerate(listSubjects):
        print(sidx, subj)
        
        # Load the confounds; unfortunately there are 2 different delimiters
        # in the confounds file > need to accomodate that
        Rest_confounds = pd.read_csv(os.path.join(confounddir, subj,'{}_ses-{}_task-rest_confounds.csv'.format(subj,visit)), header=None, sep=',',engine='python')
        TaskRun1_confounds = pd.read_csv(os.path.join(confounddir, subj,'{}_ses-{}_task-{}_run-1_confounds.csv'.format(subj,visit,cohort)), header=None, sep=',',engine='python')
        TaskRun2_confounds = pd.read_csv(os.path.join(confounddir, subj,'{}_ses-{}_task-{}_run-2_confounds.csv'.format(subj,visit,cohort)), header=None, sep=',',engine='python')                             
              
        Rest_confounds = Rest_confounds.values
        TaskRun1_confounds = TaskRun1_confounds.values
        TaskRun2_confounds = TaskRun2_confounds.values
        
        # For each hemisphere
        Rest_surf_func = {}
        Rest_surf_parc = {}
        TaskRun1_surf_func = {}
        TaskRun1_surf_parc = {}
        TaskRun2_surf_func = {}
        TaskRun2_surf_parc = {}
        
        for hemi in ['L', 'R']:
                    
            # Load functional data in fsaverage space (surface)
            Rest_surf_func[hemi] = loadGifti(os.path.join(fmriprepdir, subj, 'ses-{}'.format(visit), 'func', '{}_ses-{}_task-rest_space-fsaverage_hemi-{}_bold.func.gii'.format(subj, visit, hemi)))
            TaskRun1_surf_func[hemi] = loadGifti(os.path.join(fmriprepdir, subj, 'ses-{}'.format(visit), 'func', '{}_ses-{}_task-TAU_run-1_space-fsaverage_hemi-{}_bold.func.gii'.format(subj, visit, hemi)))
            TaskRun2_surf_func[hemi] = loadGifti(os.path.join(fmriprepdir, subj, 'ses-{}'.format(visit), 'func', '{}_ses-{}_task-TAU_run-2_space-fsaverage_hemi-{}_bold.func.gii'.format(subj, visit, hemi)))
             
            # Regress out confounds
            Restb = np.linalg.lstsq(Rest_confounds, Rest_surf_func[hemi], rcond=None)[0]
            Rest_surf_func[hemi] = Rest_surf_func[hemi] - np.matmul(Rest_confounds, Restb)
            TaskRunb = np.linalg.lstsq(TaskRun1_confounds, TaskRun1_surf_func[hemi], rcond=None)[0]
            TaskRun1_surf_func[hemi] = TaskRun1_surf_func[hemi] - np.matmul(TaskRun1_confounds, TaskRun1b)
            TaskRun2b = np.linalg.lstsq(TaskRun2_confounds, TaskRun2_surf_func[hemi], rcond=None)[0]
            TaskRun2_surf_func[hemi] = TaskRun2_surf_func[hemi] - np.matmul(TaskRun2_confounds, TaskRun2b)
            
            # For each cortical parcel, extract the average timecourse
            U = np.unique(annot[hemi])
            Rest_surf_parc[hemi] = np.zeros((Rest_surf_func[hemi].shape[0], U.shape[0]))
            TaskRun1_surf_parc[hemi] = np.zeros((TaskRun1_surf_func[hemi].shape[0], U.shape[0]))
            TaskRun2_surf_parc[hemi] = np.zeros((TaskRun2_surf_func[hemi].shape[0], U.shape[0]))
            
            for parc in U:
                Rest_surf_parc[hemi][:,parc] = np.mean(Rest_surf_func[hemi][:,annot[hemi] == parc], axis=1)
                TaskRun1_surf_parc[hemi][:,parc] = np.mean(TaskRun1_surf_func[hemi][:,annot[hemi] == parc], axis=1)
                TaskRun2_surf_parc[hemi][:,parc] = np.mean(TaskRun2_surf_func[hemi][:,annot[hemi] == parc], axis=1)
                     
        # Load the subcortical segmentation (FS "aseg" files) in MNI space
        [Rest_aseg, Rest_aseg_siz, Rest_aseg_aff]   = loadNifti(os.path.join(fmriprepdir, subj,'ses-{}'.format(visit),'func', '{}_ses-{}_task-rest_space-T1w_desc-aseg_dseg.nii.gz'.format(subj, visit)))
        [TaskRun1_aseg, TaskRun1_aseg_siz, TaskRun1_aseg_aff]   = loadNifti(os.path.join(fmriprepdir, subj,'ses-{}'.format(visit),'func', '{}_ses-{}_task-TAU_run-1_space-T1w_desc-aseg_dseg.nii.gz'.format(subj, visit)))
        [TaskRun2_aseg, TaskRun2_aseg_siz, TaskRun2_aseg_aff]   = loadNifti(os.path.join(fmriprepdir, subj,'ses-{}'.format(visit),'func', '{}_ses-{}_task-TAU_run-2_space-T1w_desc-aseg_dseg.nii.gz'.format(subj, visit)))
        
        # Load their functional data in MNI space (volume)
        [Rest_vol_func, Rest_vol_siz, Rest_vol_aff] = loadNifti(os.path.join(fmriprepdir, subj,'ses-{}'.format(visit),'func', '{}_ses-{}_task-rest_space-T1w_desc-preproc_bold.nii.gz'.format(subj, visit))) 
        [TaskRun1_vol_func, TaskRun1_vol_siz, TaskRun1_vol_aff] = loadNifti(os.path.join(fmriprepdir, subj,'ses-{}'.format(visit),'func', '{}_ses-{}_task-TAU_run-1_space-T1w_desc-preproc_bold.nii.gz'.format(subj, visit))) 
        [TaskRun2_vol_func, TaskRun2_vol_siz, TaskRun2_vol_aff] = loadNifti(os.path.join(fmriprepdir, subj,'ses-{}'.format(visit),'func', '{}_ses-{}_task-TAU_run-2_space-T1w_desc-preproc_bold.nii.gz'.format(subj, visit))) 
        
        # Regress out confounds
        Restb = np.linalg.lstsq(Rest_confounds, Rest_vol_func, rcond=None)[0]
        Rest_vol_func = Rest_vol_func - np.matmul(Rest_confounds, Restb)
        TaskRun1b = np.linalg.lstsq(TaskRun1_confounds, TaskRun1_vol_func, rcond=None)[0]
        TaskRun1_vol_func = TaskRun1_vol_func - np.matmul(TaskRun1_confounds, TaskRun1b)
        TaskRun2b = np.linalg.lstsq(TaskRun2_confounds, TaskRun2_vol_func, rcond=None)[0]
        TaskRun2_vol_func = TaskRun2_vol_func - np.matmul(TaskRun2_confounds, TaskRun2b)
        
        # For each subcortical parcel, extract the average timecourse
        Rest_vol_parc = {}
        for hemi in ['L', 'R']:
            Rest_vol_parc[hemi] = np.zeros((Rest_vol_func.shape[0], len(aseg_list[hemi])))
            for pidx, parc in enumerate(aseg_list[hemi]):
                Rest_vol_parc[hemi][:,pidx] = np.mean(Rest_vol_func[:,np.squeeze(Rest_aseg == parc)], axis=1)
     
        TaskRun1_vol_parc = {}
        for hemi in ['L', 'R']:
            TaskRun1_vol_parc[hemi] = np.zeros((TaskRun1_vol_func.shape[0], len(aseg_list[hemi])))
            for pidx, parc in enumerate(aseg_list[hemi]):
                TaskRun1_vol_parc[hemi][:,pidx] = np.mean(TaskRun1_vol_func[:,np.squeeze(TaskRun1_aseg == parc)], axis=1)

        TaskRun2_vol_parc = {}
        for hemi in ['L', 'R']:
            TaskRun2_vol_parc[hemi] = np.zeros((TaskRun2_vol_func.shape[0], len(aseg_list[hemi])))
            for pidx, parc in enumerate(aseg_list[hemi]):
                TaskRun2_vol_parc[hemi][:,pidx] = np.mean(TaskRun2_vol_func[:,np.squeeze(TaskRun2_aseg == parc)], axis=1)

        # Merge cortical and subcortical timecourses of both hemispheres
        Rest_all_parc = np.concatenate((Rest_surf_parc['L'][:,1:], Rest_surf_parc['R'][:,1:], Rest_vol_parc['L'], Rest_vol_parc['R']), axis=1)
        TaskRun1_all_parc = np.concatenate((TaskRun1_surf_parc['L'][:,1:], TaskRun1_surf_parc['R'][:,1:], TaskRun1_vol_parc['L'], TaskRun1_vol_parc['R']), axis=1)
        TaskRun2_all_parc = np.concatenate((TaskRun2_surf_parc['L'][:,1:], TaskRun2_surf_parc['R'][:,1:], TaskRun2_vol_parc['L'], TaskRun2_vol_parc['R']), axis=1)
        
        # Compute the "netmat" between the regions (cortical and subcortical),
        # and unwrap it      
        Rest_rmat = np.corrcoef(Rest_all_parc, rowvar=False);
        if partial == True:
            Rest_rinv = np.linalg.pinv(Rest_rmat)
            Rest_diag = np.diagonal(Rest_rinv)[:,None]
            Rest_rmat = np.multiply(-Rest_rinv, np.power(np.multiply(Rest_diag, Rest_diag.T), -.5))
            
        TaskRun1_rmat = np.corrcoef(TaskRun1_all_parc, rowvar=False);
        if partial == True:
            TaskRun1_rinv = np.linalg.pinv(TaskRun1_rmat)
            TaskRun1_diag = np.diagonal(TaskRun1_rinv)[:,None]
            TaskRun1_rmat = np.multiply(-TaskRun1_rinv, np.power(np.multiply(TaskRun1_diag, TaskRun1_diag.T), -.5))
        TaskRun1_zmat = np.arctanh(TaskRun1_rmat) 

        TaskRun2_rmat = np.corrcoef(TaskRun2_all_parc, rowvar=False);
        if partial == True:
            TaskRun2_rinv = np.linalg.pinv(TaskRun2_rmat)
            TaskRun2_diag = np.diagonal(TaskRun2_rinv)[:,None]
            TaskRun2_rmat = np.multiply(-TaskRun2_rinv, np.power(np.multiply(TaskRun2_diag, TaskRun2_diag.T), -.5))
        TaskRun2_zmat = np.arctanh(TaskRun2_rmat)

    	# average z-mats for the two runs of tau 
        Task_zmat = (TaskRun1_zmat + TaskRun2_zmat)/2
        Task_rmat = np.tanh(Task_zmat)		

        # Save the wrapped matrices as csv file 
        if not os.path.exists(os.path.join(netmatsdir, subj)):
            os.makedirs(os.path.join(netmatsdir, subj))
        np.savetxt(os.path.join(netmatsdir, subj,'{}_ses-{}_task-rest_netmat-{}_atlas-Schaefer2018-100P+17N_space-T1w.csv'.format(subj, visit, 'partial' if partial else 'full')), Rest_rmat, delimiter=",")
        np.savetxt(os.path.join(netmatsdir, subj,'{}_ses-{}_task-TAU_netmat-{}_atlas-Schaefer2018-100P+17N_space-T1w.csv'.format(subj, visit, 'partial' if partial else 'full')), Task_rmat,  delimiter=",")
        imageio.imwrite(os.path.join(netmatsdir, subj,'{}_ses-{}_task-rest_netmat-{}_atlas-Schaefer2018-100P+17N_space-T1w.png'.format(subj, visit, 'partial' if partial else 'full')), (255*(Rest_rmat+1)/2).astype('uint8'))
        imageio.imwrite(os.path.join(netmatsdir, subj,'{}_ses-{}_task-TAU_netmat-{}_atlas-Schaefer2018-100P+17N_space-T1w.png'.format(subj, visit, 'partial' if partial else 'full')), (255*(Task_rmat +1)/2).astype('uint8'))
        
        # Show the netmats (for each subject)
        plt.imshow(Task_rmat)
        plt.show()
        plt.imshow(Rest_rmat)
        plt.show()


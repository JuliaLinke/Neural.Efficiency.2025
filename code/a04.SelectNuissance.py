#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

For each subject, session, run: read in image, event file, confounds file; then
generate and convolve regressors; create nuisance matrix based on denoise method;
create 1st-level design and contrast files.
"""
import nibabel
import os
import glob
import numpy as np
import pandas as pd
import imageio
import scipy.stats as stats
import json

rootdir = os.path.join('/','MyWorkingDirectory')
hastask = True
hasrest = True
cohort = 'Cohort1'                #Cohort1 or Cohort2
time = 'T1'			  #T1 or T2
microtime = .001

def convolveHRF(regressors, Nvols, TR, demean=False): #========================
    # Convolves HRF with regressors and creates design matrix
    # HRF parameters (from FSL, reparameterized as function of k and theta)
    delay1 = 6.;  sigma1 = 2.449
    delay2 = 16.; sigma2 = 4.
    ratiogammas = 6;
    gammaK1 = delay1**2/sigma1**2; gammaT1 = sigma1**2/delay1;
    gammaK2 = delay2**2/sigma2**2; gammaT2 = sigma2**2/delay2;

    fmri_time = np.arange(0, Nvols*TR, microtime)
    fmri_gamm = stats.gamma.pdf(fmri_time, gammaK1, scale=gammaT1) - \
        stats.gamma.pdf(fmri_time, gammaK2, scale=gammaT2)/ratiogammas
    fmri_gamm = fmri_gamm/sum(fmri_gamm)
    fmri_pgamm = np.concatenate((np.zeros(np.shape(fmri_gamm)), fmri_gamm))
    fmri_fgamm = np.fft.fft(fmri_pgamm)
    
    design = pd.DataFrame()
    for key in regressors.keys():
        fmri_stim = np.zeros(int(round(Nvols*TR/microtime)), dtype=float)
        start = (regressors[key]['onset']/microtime).astype(int)
        end   = (start + regressors[key]['duration']/microtime).astype(int)
        start = start.reset_index(drop=True)
        end   = end.reset_index(drop=True)
        for i, _ in enumerate(start):
            if start[i] == end[i]:
                fmri_stim[start[i]:end[i]+1] = regressors[key]['height'].iloc[i]
            else:
                fmri_stim[start[i]:end[i]] = regressors[key]['height'].iloc[i]
        fmri_pstim = np.concatenate((np.zeros(np.shape(fmri_stim)),fmri_stim))
        fmri_fstim = np.fft.fft(fmri_pstim)
        fmri_conv = np.real(np.fft.ifft(fmri_fstim * fmri_fgamm))
        fmri_conv = fmri_conv[0:np.shape(fmri_gamm)[0]]
        design[key] = fmri_conv
    if demean: 
        design = design - design.mean() 
    idx = np.round(np.arange(0, Nvols*TR, TR)/microtime).astype(int)
    design = design.iloc[idx]
    design.reset_index(inplace=True, drop=True)
    return design

def orthogonalize(design, orth, wrt): #========================================
    # Orthogonalize a specified regressor in relation to another specified regressor
    X = design[orth].values[:,None]
    Z = design[wrt].values[:,None]
    Rz = np.eye(design.shape[0]) - np.matmul(Z, np.linalg.pinv(Z))
    design[orth] = np.matmul(Rz, X)[:,0]
    return design 

def make_nuisance(confounds, method='AROMA2'): # ==============================
    # Returns Z matrix with nuisance regressors (based on different noising methods); always includes cosines and steady-state
    Z_lab = []
    
    for col in confounds.columns:   # Cosines for trends
        if 'cosine' in col or 'non_steady_state' in col:
            Z_lab.append(col)
    if method == 'Satterth24+outliers':
        for col in confounds.columns:
            if  'trans_' in col or 'rot_' in col or \
                'motion_outlier' in col or \
                'csf' == col or 'white_matter' == col:
                Z_lab.append(col)
    elif method == 'AROMA':
        for col in confounds.columns:
            if 'aroma' in col or \
               'csf' == col or 'white_matter' == col:
                Z_lab.append(col)
    elif method == 'AROMA2':
        for col in confounds.columns:
            if 'aroma' in col or \
                'csf' == col or 'white_matter' == col or 'global_signal' == col:
                Z_lab.append(col)
    elif method == 'CompCor':
        for col in confounds.columns:
            if 'comp_cor' in col:
                Z_lab.append(col)
    elif method == '6RP':
        for col in confounds.columns:
            if  col == 'trans_x' or col == 'trans_y' or col == 'trans_z' \
                or col == 'rot_x' or col == 'rot_y' or col == 'rot_z':
                Z_lab.append(col)
    Z = confounds[Z_lab]
    return Z

def writeDesignMatrix(design, savedir, subj, ses, run):#=======================
    # Saves design matrix, and some additional diagnostic metrics.
    
    # Design
    design.to_csv(os.path.join(savedir, '{}{}_task-{}{}_design.csv'.format(subj, ses if not ses else '_{}'.format(ses), TASK, run if not run else '_{}'.format(run))), header=True, index=False, sep=',')
    cols = design.columns.tolist()
    design = np.array(design)
    imageio.imwrite(os.path.join(savedir, '{}{}_task-{}{}_design.png'.format(subj, ses if not ses else '_{}'.format(ses), TASK, run if not run else '_{}'.format(run))), (255*(design+0.5)/2).astype('uint8'))
    
    # Correlation, eigenvalues
    design_correl = np.corrcoef(design[:,0:-1].T)
    imageio.imwrite(os.path.join(savedir, '{}{}_task-{}{}_design_correl.png'.format(subj, ses if not ses else '_{}'.format(ses), TASK, run if not run else '_{}'.format(run))), (255*(design_correl+1)/2).astype('uint8'))
    _, s, _ = np.linalg.svd(design[:,0:-1] - design[:,0:-1].mean())
    np.savetxt(os.path.join(savedir, '{}{}_task-{}{}_eigenvalues.csv'.format(subj, ses if not ses else '_{}'.format(ses), TASK, run if not run else '_{}'.format(run))), s[None], delimiter=",", fmt='%g')
    
    # Amount of usable data, based on df
    J = {}
    J['Number of rows (timepoints)']    = design.shape[0]
    J['Number of columns (regressors)'] = design.shape[1]
    J['Number of motion outliers'] = 0
    for col in cols:
        if type(col) is str and 'motion_outlier' in col:
            J['Number of motion outliers'] = J['Number of motion outliers']+1
    J['Rank']                 = int(np.linalg.matrix_rank(design))
    J['Degrees of freedom']   = design.shape[0] - J['Rank']
    J['Condition number']     = np.linalg.cond(design)
    with open(os.path.join(savedir, '{}{}_task-{}{}_diagnostics.json'.format(subj, ses if not ses else '_{}'.format(ses), TASK, run if not run else '_{}'.format(run))), 'w') as fp:
        json.dump(J, fp, indent=2)
    return

# ======== [ Main ] ===========================================================

listSubjects = open(os.path.join(rootdir,'lists/{}.{}.txt'.format(cohort,time)),'r')
listSubjects = listSubjects.read().split()
listSubjects.sort()
fmriprepdir  = os.path.join(rootdir,'derivatives/{}/fmriprep/'.format(cohort))
timingdir  = os.path.join(rootdir,'derivatives/{}/Timing_Files'.format(cohort))
BIDSdir = os.path.join(rootdir, 'BIDS', cohort)
confounddir  = os.path.join(rootdir,'derivatives/{}/RegressNuissance'.format(cohort))


#1st do things for task data
if hastask == True:
    TASK = 'TAU'
    # For each subject
    for sidx, subj in enumerate(listSubjects):
        print(sidx, subj)
    
    # For each session (if existing)
        sesList = glob.glob(os.path.join(BIDSdir, subj, 'ses-*'))
        if not sesList:
            sesList = ['']
        else:
            if len(sesList) == 1:
                sesList[0] = os.path.split(sesList[0])[1]
            elif len(sesList) > 1:
                nsess = [1,(len(sesList))]
                for idx in nsess:
                    sesList[idx-1] = os.path.join('ses-{}'.format(idx))
    
            for ses in sesList:
        
                # For each run (if existing)
                runList = glob.glob(os.path.join(timingdir, '{}_{}_task-{}_run-*_events.tsv'.format(subj, ses,TASK)))
                if not runList:
                    runList = ['']
                else:
                    if len(runList) == 1:
                        runList[0] = os.path.split(runList[0])[1]
                    elif len(runList) > 1:
                        nrun = [1,(len(runList))]
                        for idx in nrun:
                            runList[idx-1] = os.path.join('run-{}'.format(idx))
        
                for run in runList:
            
                    # These files will be loaded later. Continue only if both exist:
                    event_file = os.path.join(timingdir, '{}_{}_task-{}_{}_events.tsv'.format( subj, ses, TASK, run))
                    image_file = os.path.join(BIDSdir, subj, ses, 'func', '{}{}_task-{}{}_bold.nii.gz'.format(subj, ses if not ses else '_{}'.format(ses), TASK, run if not run else '_{}'.format(run)))
                    confs_file = os.path.join(fmriprepdir, subj, ses, 'func', '{}{}_task-{}{}_desc-confounds_timeseries.tsv'.format(subj, ses if not ses else '_{}'.format(ses), TASK, run if not run else '_{}'.format(run)))
            
                    if os.path.exists(event_file) and os.path.exists(image_file) and os.path.exists(confs_file):

                        # Load the confounds and create Z matrix
                        confounds = pd.read_csv(confs_file, delimiter='\t')
                        confounds = confounds.fillna(0)
                        Z = make_nuisance(confounds, method='AROMA2')
                        # Delete the non-steady state column; for task we need 4 non-steady state regressors irrespective of fMRIprep
                        cols = [c for c in Z.columns if 'non_steady_state_outlier0' not in c.lower()]
                        Z=Z[cols]
                        Z['non-steady_1'] = 0
                        Z.loc[:0, 'non-steady_1'] = 1
                        Z['non-steady_2'] = 0
                        Z.loc[1:1, 'non-steady_2'] = 1
                        Z['non-steady_3'] = 0
                        Z.loc[2:2, 'non-steady_3'] = 1
                        Z['non-steady_4'] = 0
                        Z.loc[3:3, 'non-steady_4'] = 1
                
                
                        # Save variables that need to be regressed out in netmats script
                        if not os.path.exists(os.path.join(confounddir, subj)):
                            os.makedirs(os.path.join(confounddir, subj))
                        np.savetxt(os.path.join(confounddir, subj,'{}{}_task-{}{}_confounds.csv'.format(subj, ses if not ses else '_{}'.format(ses), TASK, run if not run else '_{}'.format(run))), Z, delimiter=",")
        
                else:
                    print("Missing event or image file:")
                    print("- {}".format(event_file))
                    print("- {}".format(image_file))

if hasrest == True:
    TASK = 'rest'
    # For each subject
    for sidx, subj in enumerate(listSubjects):
        print(sidx, subj)
    
    # For each session (if existing)
        sesList = glob.glob(os.path.join(BIDSdir, subj, 'ses-*'))
        if not sesList:
            sesList = ['']
        else:
            if len(sesList) == 1:
                sesList[0] = os.path.split(sesList[0])[1]
            elif len(sesList) > 1:
                nsess = [1,(len(sesList))]
                for idx in nsess:
                    sesList[idx-1] = os.path.join('ses-{}'.format(idx))
    
            for ses in sesList:
                # These files will be loaded later. Continue only if both exist:
                    confs_file = os.path.join(fmriprepdir, subj, ses, 'func', '{}{}_task-{}_desc-confounds_timeseries.tsv'.format(subj, ses if not ses else '_{}'.format(ses), TASK))

                    if os.path.exists(confs_file):

                        # Load the confounds and create Z matrix
                        confounds = pd.read_csv(confs_file, delimiter='\t')
                        confounds = confounds.fillna(0)
                        Z = make_nuisance(confounds, method='AROMA2')

		    else:
                        print(f"Missing file: {confs_file}")

                    # Save variables that need to be regressed out in the netmats script
                    if not os.path.exists(os.path.join(confounddir, subj)):
                        os.makedirs(os.path.join(confounddir, subj))
                    np.savetxt(os.path.join(confounddir, subj,'{}{}_task-{}_confounds.csv'.format(subj, ses if not ses else '_{}'.format(ses), TASK)), Z, delimiter=",")

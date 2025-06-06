# Neural.Efficiency.2025
Analysis files and meta-data to accompany the manuscript "Reduced Threat-Related Neural Efficiency: A Possible Biomarker for Pediatric Anxiety Disorders", which has been accepted for publication by The American Journal of Psychiatry.

This repository contains analysis scripts (bash, Python, Matlab, R) used for preprocessing and analysis of the data found in our OpenNeuro dataset (link will be included, once dataset has been uploaded completely). Although the code used to deface the scans is not included in this repository, we de-identified these images using using the [DSST Defacing Pipeline](https://github.com/nimh-dsst/dsst-defacing-pipeline).

#### Directory Structure for Data Preprocessing

The scripts assume the following directory structure:

```
MyWorkingDirectory
├── atlas
│   └──  Schaefer
├── BIDS
│   ├── Cohort1
│   │   └──  sub-[participant-id]
│   │       ├── ses-1
│   │       │   ├── anat
│   │       │   └── func
│   │       └── ses-2
│   │           ├── anat
│   │           └── func
│   └── Cohort2
│       └──  sub-[participant-id]
│           └── ses-1
│           │   ├── anat
│           │   └── func
│           └── ses-2
│               ├── anat
│               └── func
├── code
├── derivatives
│   ├── Cohort1
│   │   └──  Efficiency
│   │   └──  Fmriprep
│   │       └──  sub-[participant-id]
│   │   └──  MRIQC
│   │   └──  RegressNuissance
│   │       └──  sub-[participant-id]
│   │   └──  Ses-1_Netmats
│   │       └──  sub-[participant-id]
│   │   └──  Ses-2_Netmats
│   │       └──  sub-[participant-id]
│   ├── Cohort2
│   │   └──  Efficiency
│   │   └──  Fmriprep
│   │       └──  sub-[participant-id]
│   │   └──  MRIQC
│   │   └──  RegressNuissance
│   │       └──  sub-[participant-id]
│   │   └──  Ses-1_Netmats
│   │       └──  sub-[participant-id]
│   │   └──  Ses-2_Netmats
│   │       └──  sub-[participant-id]
│   ├── CovBat
├── lists
├── slurm
├── stats
│   ├── Palm
│   ├── Figures
```

## WORKFLOW
### Preprocessing
#### a00.MakeCondaEnv
This script installs a custom Python environment for the project using [Miniconda](https://www.anaconda.com/docs/getting-started/miniconda/main). It sets up all required packages and tools for data processing and analysis.

#### a01.MRIQC
This script automates the quality assessment of structural and functional MRI data using MRIQC. It is designed to run on a high-performance computing (HPC) cluster using the Swarm job scheduler (SLURM-compatible). The script iterates over a list of participant IDs for a given cohort and timepoint. First, it checks whether [MRIQC](https://mriqc.readthedocs.io/en/stable/) outputs already exist. For each participant, it the (a) creates a temporary workspace on local scratch, (b) runs MRIQC with participant-specific settings, (c) copies the output to a shared derivatives folder,and (d) cleans up temporary files. Finally, the script compiles all commands into a Swarm file and submits the job for parallel execution. This setup allows for efficient parallel processing of many participants.

#### a02.MRIQC.Group
This script helps organize participant-level outputs from [MRIQC](https://mriqc.readthedocs.io/en/stable/) and runs a group-level summary analysis. It is intended to be used after all individual MRIQC runs are complete 

#### a03.FMRIPREP
This script automates participant-level preprocessing of BIDS-formatted MRI data using [fMRIPrep](https://fmriprep.org/en/stable/) on a high-performance computing (HPC) cluster with the Swarm job scheduler.

#### a04.SelectNuissance
This is a Python script. Therefore, first, you need to run two lines to activate the environment you made with a00.MakeCondaEnv:
>> source /MyWorkingDirectory/code/conda/etc/profile.d/conda.sh

>> conda activate /MyWorkingDirectory/code/conda/envs/Efficiency

Next you open Spyder

>> spyder

From here you run this script, which prepares first-level regressors and nuisance confound matrices for resting-state and task-based fMRI data.

#### a05.Netmats
This Python script, which you also run via spyder, extracts surface and subcortical time series from fMRI data, regresses out nuisance signals, and computes full or partial correlation matrices (netmats) across brain regions.

#### a06.PrepCovBat
This MATLAB script extracts and formats subject-level functional connectivity matrices (rest and task) for harmonization using methods like [CovBat](https://github.com/andy1764/CovBat_Harmonization).

#### a07.CovBat
This R script harmonizes functional connectivity matrices derived from resting-state and task-based fMRI across cohorts using the CovBat method. It removes unwanted site/scanner effects while preserving biologically meaningful variability.

#### a08.NeuralEfficiency
This MATLAB script estimates neural efficiency by computing the correlation between functional connectivity matrices derived from resting-state and task-based fMRI (dot-probe task). It also computes contribution scores for each connection.

### Analysis
#### b01.PrepMotion
This R script extracts and summarizes motion-related parameters (e.g., framewise displacement, DVARS, translation and rotation along X/Y/Z axes) from [fMRIPrep](https://fmriprep.org/en/stable/) confounds files across different task conditions (resting-state, dot-probe task runs 1 and 2).

#### b02.AnalyseMotion
This R script investigates the relationship between head motion during fMRI and parent-reported anxiety symptoms in youth.

#### b03.Reliability
This R script assesses the **test–retest reliability** of neural efficiency across the two scanning sessions.

#### b04.PrepPermutationTests
This script prepares input matrices and vectors for permutation-based inference using [PALM](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/PALM), a toolbox for nonparametric statistics.

#### b05.PermutationTests
This MATLAB script runs permutation tests using [PALM](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/PALM) to investigate how (a) single edges contribute to neural efficiency, and how (b) neural efficiency relates anxiety as a diagnostic category and symptom dimension, and drift diffusion model (DDM) parameters (i.e., drifte rate, bias).

#### b06.Figures
This R script visualizes results for the main manuscript and the supplementary materials.

#### b07.CBT
This R script analyzes the relationship between neural efficiency and exposure-based cognitive behavioral therapy (EX-CBT). 

#### Circle plots
The circle plots illustrating edges that significantly contribute to neural efficiency or its relationship with other variables were created using [NeuroMarvl](https://immersive.erc.monash.edu/neuromarvl/). You will need the files in "CirclePlots" to recreate them.

## CITATION AND USE
If you use this code or build upon it in your own work, please cite the corresponding paper or contact me directly. This helps ensure appropriate use and attribution.
📩 Feel free to reach out via GitHub or email if you have questions or plan to use the scripts in your own project.

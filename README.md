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
│   ├── ComBat
├── lists
├── slurm
├── stats
│   ├── Palm
│   ├── Figures
```

#### a01.MRIQC

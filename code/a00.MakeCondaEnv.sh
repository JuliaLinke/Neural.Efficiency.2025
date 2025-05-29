#!/bin/bash

# One-time setup script for installing Miniconda and dependencies for the Rest+TAU project

# ==== CONFIGURATION ====
ROOTDIR=/MyWorkingDirectory
cd ${ROOTDIR}/code

# ==== DOWNLOAD AND INSTALL MINICONDA ====
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod +x Miniconda3-latest-Linux-x86_64.sh

mkdir -p "${ROOTDIR}/temp"
TMPDIR="${ROOTDIR}/temp"

./Miniconda3-latest-Linux-x86_64.sh -p ${ROOTDIR}/code/conda -b
rm Miniconda3-latest-Linux-x86_64.sh

# ==== ACTIVATE CONDA AND INSTALL PACKAGES ====
source ${ROOTDIR}/code/conda/etc/profile.d/conda.sh
conda activate base
# Optionally update conda
# conda update -y conda

# Install Python and common packages
conda install -y python=3.7.7
conda install -y numpy scipy pandas scikit-learn statsmodels h5py matplotlib spyder xlrd natsorted

# Install additional tools via pip
pip install --quiet --no-input mriqc fmriprep nda-tools


# Optional: create a separate environment for MRIQC from YAML file
if [[ -f rest+tau.yml ]]; then
    conda env create --file rest+tau.yml
else
    echo "YAML file MRIQC.yml not found. Skipping environment creation."
fi


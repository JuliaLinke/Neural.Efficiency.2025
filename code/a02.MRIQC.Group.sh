#!/bin/bash

# Load MRIQC module (adjust based on your HPC environment)
module load mriqc

# Set cohort and timepoint
Cohort='Cohort1'          # adjust as needed
Timepoint='T1'            # adjust as needed

# Define root directory
ROOTDIR=/MyWorkingDirectory

# Clean up and reorganize individual MRIQC outputs
for s in $(cat "${ROOTDIR}/lists/${Cohort}.${Timepoint}.txt"); do
  cd "${ROOTDIR}/derivatives/${Cohort}/MRIQC" || exit 1

  # Remove log directory if it exists
  rm -rf "${s}.out/logs"

  # Move contents of participant output folder to main MRIQC dir
  mv "${s}.out/"* .

  # Remove now-empty participant folder
  rm -rf "${s}.out"
done

# Run group-level MRIQC analysis
mriqc "${ROOTDIR}/BIDS/${Cohort}" "${ROOTDIR}/derivatives/${Cohort}/MRIQC" group

# Set permissions for all files
chmod -R 770 "${ROOTDIR}/derivatives/${Cohort}/MRIQC"

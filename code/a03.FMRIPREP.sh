#!/bin/bash

# Load Fmriprep (adjust depending on your system/environment)
module load fmriprep

# ==== CONFIGURATION ====
Cohort='Cohort1'              	# adjust as needed
Timepoint='T1'          	# adjust as needed 
ROOTDIR=/MyWorkingDirectory

# Swarm file to store job commands
SWARMFILE=${ROOTDIR}/slurm/swarm_fmriprep.txt
echo -n "" > ${SWARMFILE}

# ==== LOOP OVER PARTICIPANTS ====
for s in $(cat ${ROOTDIR}/lists/${Cohort}.${Timepoint}.txt) ; do
   if [[ ! -d ${ROOTDIR}/derivatives/${Cohort}/Fmriprep/${s} ]] ; then 
      echo "let \"rnd = ${RANDOM} % 300\" ; sleep \${rnd} ; \
            export TMPDIR=/lscratch/\${SLURM_JOBID} ; \
            mkdir -p \${TMPDIR}/${s}.out \${TMPDIR}/${s}.wrk ; \
            fmriprep ${ROOTDIR}/${Cohort} \${TMPDIR}/${s}.out participant --participant_label ${s} -w \${TMPDIR}/${s}.wrk --use-aroma --output-space MNI152NLin2009cAsym:res-2 T1w fsnative  fsaverage fsaverage5 --nthreads 1 --omp-nthreads 1 --skip_bids_validation --notrack; \
            mkdir -p ${ROOTDIR}/derivatives/${Cohort}/Fmriprep ; \
            rsync -a \${TMPDIR}/${s}.out/ ${ROOTDIR}/derivatives/${Cohort}/Fmriprep ; \
            rm -rf \${TMPDIR}/${s}.out \${TMPDIR}/${s}.wrk " >> ${SWARMFILE}
   fi
done

# ==== SUBMIT SWARM JOB ====
swarm -f ${SWARMFILE} -g 64 -t 12 -p 1 --gres=lscratch:100 --logdir ${ROOTDIR}/slurm --time=120:00:00 --job-name Efficiency --module fmriprep



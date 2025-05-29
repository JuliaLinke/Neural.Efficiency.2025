#!/bin/bash

# Load MRIQC (adjust depending on your system/environment)
module load mriqc		# depends on system settings

# ==== CONFIGURATION ====
Cohort='Cohort1'              	# adjust as needed
Timepoint='T1'          	# adjust as needed 
ROOTDIR=/MyWorkingDirectory

# Swarm file to store job commands
SWARMFILE=${ROOTDIR}/slurm/swarm_mriqc.txt
echo -n "" > ${SWARMFILE}

# ==== LOOP OVER PARTICIPANTS ====
for s in $(cat ${ROOTDIR}/lists/${Cohort}.${Timepoint}.txt); do
   if [[ ! -d ${ROOTDIR}/derivatives/${Cohort}/MRIQC/${s} ]]; then
      echo "let \"rnd = \${RANDOM} % 300\"; sleep \${rnd}; \
            export TMPDIR=/lscratch/\${SLURM_JOBID}; \
            mkdir -p \${TMPDIR}/${s}.out \${TMPDIR}/${s}.wrk; \
            mriqc ${ROOTDIR}/${Cohort} \${TMPDIR}/${s}.out participant --participant_label ${s} -w \${TMPDIR}/${s}.wrk --n_procs 1 --no-sub --fd_thres 0.5; \
            mkdir -p ${ROOTDIR}/derivatives/${Cohort}/MRIQC; \
            rsync -a \${TMPDIR}/${s}.out ${ROOTDIR}/derivatives/${Cohort}/MRIQC; \
            rm -rf \${TMPDIR}/${s}.out \${TMPDIR}/${s}.wrk" >> ${SWARMFILE}
   fi
done

# ==== SUBMIT SWARM JOB ====
swarm -f ${SWARMFILE} -g 40 -t auto --gres=lscratch:40 --logdir ${ROOTDIR}/slurm --time=120:00:00 --job-name Efficiency

#!/bin/bash

module load mriqc		# depends on system settings

Cohort='Cohort1'              	# adjust as needed
Timepoint='T1'          	# adjust as needed 

ROOTDIR=/MyWorkingDirectory
SWARMFILE=${ROOTDIR}/slurm/swarm_mriqc.txt
echo -n "" > ${SWARMFILE}

for s in $(cat ${ROOTDIR}/lists/${Cohort}.${Timepoint}.txt); do
   if [[ ! -d ${ROOTDIR}/derivatives/MRIQC/${s} ]]; then
      echo "let \"rnd = \${RANDOM} % 300\"; sleep \${rnd}; \
            export TMPDIR=/lscratch/\${SLURM_JOBID}; \
            mkdir -p \${TMPDIR}/${s}.out \${TMPDIR}/${s}.wrk; \
            mriqc ${ROOTDIR}/${Cohort} \${TMPDIR}/${s}.out participant --participant_label ${s} -w \${TMPDIR}/${s}.wrk --n_procs 1 --no-sub --fd_thres 0.5; \
            mkdir -p ${ROOTDIR}/derivatives/${Cohort}/MRIQC; \
            rsync -a \${TMPDIR}/${s}.out ${ROOTDIR}/derivatives/${Cohort}/MRIQC; \
            rm -rf \${TMPDIR}/${s}.out \${TMPDIR}/${s}.wrk" >> ${SWARMFILE}
   fi
done

swarm -f ${SWARMFILE} -g 40 -t auto --gres=lscratch:40 --logdir ${ROOTDIR}/slurm --time=120:00:00 --job-name Efficiency

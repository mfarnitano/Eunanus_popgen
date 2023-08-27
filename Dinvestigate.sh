#!/bin/bash
#SBATCH --job-name=Dinvestigate		                        # Job name
#SBATCH --partition=batch		                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=10		                            # Number of cores per task
#SBATCH --mem=20gb			                                # Total memory for job
#SBATCH --time=48:00:00  		                            # Time limit hrs:min:sec
#SBATCH -o /home/mcf96392/logs/Brevipes_FINAL/%x.%j.out		# Standard output log
#SBATCH -e /home/mcf96392/logs/Brevipes_FINAL/%x.%j.err		#Standard error log
#SBATCH --mail-user=mcf96392@uga.edu                    # Where to send mail
#SBATCH --mail-type=ALL                            # Mail events (BEGIN, END, FAIL, ALL)

###SETUP

BATCH_NAME=Brevipes_FINAL

SCRIPTS_DIR=~/Git_repos/BrevipesPopgen/Brevipes_submission
LOG_DIR=/home/mcf96392/logs/${BATCH_NAME}
WORKING_DIR=/scratch/mcf96392/${BATCH_NAME}

PROGRAM_PATH=~/apps/Dsuite/Build/Dsuite

VCF=${WORKING_DIR}/VCFs/${BATCH_NAME}.pass.31calledSNPs.vcf.gz

NTHREADS=10
MEM=20

###LOG
cd $SCRIPTS_DIR
echo "${SLURM_JOB_ID},${SLURM_JOB_NAME},$(git rev-parse --short HEAD),$(date)" >> ${LOG_DIR}/joblog.txt

##MODULES
ml GCC/9.2.0

###folder setup
mkdir -p ${WORKING_DIR}/Dinvestigate
cd ${WORKING_DIR}/Dinvestigate

FILENAME=${VCF##*/}
PREFIX=${FILENAME%.vcf.gz}

${PROGRAM_PATH} Dinvestigate -w 100,50 -n ${PREFIX}_BCNA $VCF ${SCRIPTS_DIR}/inputs/Dinv_BCNA.txt ${SCRIPTS_DIR}/inputs/Trio_BCN.txt
${PROGRAM_PATH} Dinvestigate -w 100,50 -n ${PREFIX}_CSFA $VCF ${SCRIPTS_DIR}/inputs/Dinv_CSFA.txt ${SCRIPTS_DIR}/inputs/Trio_CSF.txt

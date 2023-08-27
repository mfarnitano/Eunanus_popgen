#!/bin/bash
#SBATCH --job-name=Dsuite		                        # Job name
#SBATCH --partition=batch		                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=10		                            # Number of cores per task
#SBATCH --mem=20gb			                                # Total memory for job
#SBATCH --time=48:00:00  		                            # Time limit hrs:min:sec
#SBATCH -o /home/mcf96392/logs/Brevipes_submission/%x.%j.out		# Standard output log
#SBATCH -e /home/mcf96392/logs/Brevipes_submission/%x.%j.err		#Standard error log
#SBATCH --mail-user=mcf96392@uga.edu                    # Where to send mail
#SBATCH --mail-type=ALL                            # Mail events (BEGIN, END, FAIL, ALL)

###SETUP

BATCH_NAME=Brevipes_submission

SCRIPTS_DIR=~/Git_repos/BrevipesPopgen/Brevipes_submission
LOG_DIR=/home/mcf96392/logs/${BATCH_NAME}
WORKING_DIR=/scratch/mcf96392/${BATCH_NAME}

PROGRAM_PATH=~/apps/Dsuite/Build/Dsuite

NTHREADS=10
MEM=20

###LOG
cd $SCRIPTS_DIR
echo "${SLURM_JOB_ID},${SLURM_JOB_NAME},$(git rev-parse --short HEAD),$(date)" >> ${LOG_DIR}/joblog.txt

##MODULES
ml GCC/9.2.0

#Memory and threads
MEM=20
NTHREADS=10

###LOG
cd $SCRIPTS_DIR
echo "${PBS_JOBID},${PBS_JOBNAME},$(git rev-parse --short HEAD),$(date)" >> ${LOG_DIR}/joblog.txt

###folder setup
mkdir -p ${WORKING_DIR}/Dsuite

cd ${WORKING_DIR}/Dsuite

for VCF_TYPE in pass.31calledSNPs; do
	VCF=${WORKING_DIR}/VCFs/${BATCH_NAME}.${VCF_TYPE}.vcf.gz
	VCF_FILENAME=${BATCH_NAME}.${VCF_TYPE}.vcf.gz
	printf "\nLooping through VCF file %s...\n" $VCF | tee >(cat >&2)

	FILENAME=${VCF##*/}
	PREFIX=${FILENAME%.vcf.gz}

	printf "\nRunning Dsuite Dtrios on sample %s with outgroup configuration %s...\n" $PREFIX AURA | tee >(cat >&2)
	${PROGRAM_PATH} Dtrios $VCF ${SCRIPTS_DIR}/inputs/Dsuite_popfile_AURA.txt -o ${PREFIX}_AURA -k 100 -t ${SCRIPTS_DIR}/inputs/poplevel_treefile.tre
	${PROGRAM_PATH} Fbranch ${SCRIPTS_DIR}/inputs/poplevel_treefile.tre ${PREFIX}_AURA_tree.txt > ${PREFIX}_AURA.Fbranch
	python ~/apps/Dsuite/utils/dtools.py ${PREFIX}_AURA.Fbranch ${SCRIPTS_DIR}/inputs/poplevel_treefile.tre --outgroup AURA
	mv Fbranch.png ${PREFIX}_AURA.Fbranch.png
	mv Fbranch.svg ${PREFIX}_AURA.Fbranch.svg
done


printf "\nCompleted all runs...\n" | tee >(cat >&2)

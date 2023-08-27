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
ITER=$1
BATCH_NAME=downsample_ITER_${ITER}

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

###folder setup
mkdir -p ${WORKING_DIR}
mkdir -p ${WORKING_DIR}/Dsuite

cd ${WORKING_DIR}/Dsuite


printf "\nStarting ITER %s\n" ${ITER} | tee >(cat >&2)
VCF_TYPE=${BATCH_NAME}_downsample_ITER_${ITER}.pass.11calledSNPs
VCF=${WORKING_DIR}/VCFs/${VCF_TYPE}.vcf.gz

TREE=${SCRIPTS_DIR}/inputs/poplevel_treefile.tre
COMPLETE_POPFILE=${SCRIPTS_DIR}/inputs/Dsuite_popfile_AURA.txt
FASTQ_TABLE=${WORKING_DIR}/subset_fastq_table_ITER_${ITER}.txt

###create popfile
POPFILE=${WORKING_DIR}/Dsuite_popfile_ITER_${ITER}.txt
cut -f2 $FASTQ_TABLE | while read -r line; do
	grep $line $COMPLETE_POPFILE >> $POPFILE
done

printf "\n...Running Dtrios %s\n" ${ITER} | tee >(cat >&2)
${PROGRAM_PATH} Dtrios $VCF $POPFILE -o ${VCF_TYPE}_I${ITER} -k 100 -t $TREE

printf "\nFinished Dtrios, running Fbranch... %s\n" ${ITER} | tee >(cat >&2)
${PROGRAM_PATH} Fbranch $TREE ${VCF_TYPE}_I${ITER}_tree.txt > ${VCF_TYPE}_I${ITER}.Fbranch

printf "\nFinished Fbranch, running dtools %s\n" ${ITER} | tee >(cat >&2)
python ~/apps/Dsuite/utils/dtools.py ${VCF_TYPE}_I${ITER}.Fbranch ${POPFILE} --outgroup AURA
mv fbranch.png ${VCF_TYPE}_I${ITER}.fbranch.png

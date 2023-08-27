#!/bin/bash
#SBATCH --job-name=samplefastqs	                        # Job name
#SBATCH --partition=batch		                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=10		                            # Number of cores per task
#SBATCH --mem=20gb			                                # Total memory for job
#SBATCH --time=72:00:00  		                            # Time limit hrs:min:sec
#SBATCH -o /home/mcf96392/logs/Brevipes_submission/%x.%j.out		# Standard output log
#SBATCH -e /home/mcf96392/logs/Brevipes_submission/%x.%j.err		#Standard error log
#SBATCH --mail-user=mcf96392@uga.edu                    # Where to send mail
#SBATCH --mail-type=ALL                            # Mail events (BEGIN, END, FAIL, ALL)

BATCH_NAME=downsample_ITER_0

SCRIPTS_DIR=~/Git_repos/BrevipesPopgen/Brevipes_submission
LOG_DIR=/home/mcf96392/logs/${BATCH_NAME}
WORKING_DIR=/scratch/mcf96392/${BATCH_NAME}

NTHREADS=10
MEM=20

GENOME=~/Genomes/Maurantiacus/M_aurantiacus_v1_splitline_ordered.fasta

###LOG
cd $SCRIPTS_DIR
echo "${SLURM_JOB_ID},${SLURM_JOB_NAME},$(git rev-parse --short HEAD),$(date)" >> ${LOG_DIR}/joblog.txt

ml seqtk/1.3-GCC-8.3.0

mkdir -p ${WORKING_DIR}/sampled_fastqs
mkdir -p ${WORKING_DIR}/sampled_fastqs

while read -r line; do
	folder=$(echo -n $line | tr -s ' ' | cut -d' ' -f1)
	prefix=$(echo -n $line | tr -s ' ' | cut -d' ' -f2)
	R1=${folder}/${prefix}_R1.fastq.gz
	R2=${folder}/${prefix}_R2.fastq.gz

	#sample
	seqtk sample -s100 ${R1} 12000000 > ${WORKING_DIR}/sampled_fastqs/${prefix}_R1.fastq
	seqtk sample -s100 ${R2} 12000000 > ${WORKING_DIR}/sampled_fastqs/${prefix}_R2.fastq
	gzip ${WORKING_DIR}/sampled_fastqs/${newprefix}_R1.fastq
	gzip ${WORKING_DIR}/sampled_fastqs/${newprefix}_R2.fastq

	#trim and align
	sbatch ${SCRIPTS_DIR}/fastq_align_each.sh ${SCRIPTS_DIR} ${BATCH_NAME} ${GENOME} ${WORKING_DIR}/sampled_fastqs/ ${prefix}

done < ${SCRIPTS_DIR}/inputs/subset_fastq_table_ITER_0.txt

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

###call this script multiple times with ITER=1,2,3... to get unique subsets
ITER=$1
BATCH_NAME=downsample_ITER_$ITER

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

mkdir -p ${WORKING_DIR}
mkdir -p ${WORKING_DIR}/sampled_fastqs/
cd ${WORKING_DIR}/sampled_fastqs/

MASTER=${SCRIPTS_DIR}/inputs/Brevipes_complete_fastq_table.txt
grep "BREV" $MASTER | shuf -n2 > ${WORKING_DIR}/subset_fastq_table_ITER_${ITER}.txt
grep "FREM" $MASTER | shuf -n2 >> ${WORKING_DIR}/subset_fastq_table_ITER_${ITER}.txt
grep "JOHN" $MASTER | shuf -n2 >> ${WORKING_DIR}/subset_fastq_table_ITER_${ITER}.txt
grep "NANU\|SESP\|CONS" $MASTER >> ${WORKING_DIR}/subset_fastq_table_ITER_${ITER}.txt
grep "AURA\|CLEV" $MASTER | shuf -n2 >> ${WORKING_DIR}/subset_fastq_table_ITER_${ITER}.txt

SEED=$(echo "${ITER} * 123" | bc)
while read -r line; do
	folder=$(echo -n $line | tr -s ' ' | cut -d' ' -f1)
	prefix=$(echo -n $line | tr -s ' ' | cut -d' ' -f2)
	R1=${folder}/${prefix}_R1.fastq.gz
	R2=${folder}/${prefix}_R2.fastq.gz

	#sample
	seqtk sample -s${SEED} ${R1} 12000000 > ${WORKING_DIR}/sampled_fastqs/${newprefix}_R1.fastq
	seqtk sample -s${SEED} ${R2} 12000000 > ${WORKING_DIR}/sampled_fastqs/${newprefix}_R2.fastq
	gzip ${WORKING_DIR}/sampled_fastqs/${newprefix}_R1.fastq
	gzip ${WORKING_DIR}/sampled_fastqs/${newprefix}_R2.fastq
	SEED=$(echo "${SEED} + 13" | bc)

	#trim and align
	sbatch ${SCRIPTS_DIR}/fastq_align_each.sh ${SCRIPTS_IDR} ${BATCH_NAME} ${GENOME} ${WORKING_DIR}/ITER_${ITER}/sampled_fastqs/ ${prefix}

done < ${WORKING_DIR}/subset_fastq_table_ITER_${ITER}.txt

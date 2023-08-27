#!/bin/bash
#SBATCH --job-name=vcf_counts	                        # Job name
#SBATCH --partition=batch		                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=10		                            # Number of cores per task
#SBATCH --mem=20gb			                                # Total memory for job
#SBATCH --time=72:00:00  		                            # Time limit hrs:min:sec
#SBATCH -o /home/mcf96392/logs/Brevipes_submission/%x.%j.out		# Standard output log
#SBATCH -e /home/mcf96392/logs/Brevipes_submission/%x.%j.err		#Standard error log
#SBATCH --mail-user=mcf96392@uga.edu                    # Where to send mail
#SBATCH --mail-type=ALL                            # Mail events (BEGIN, END, FAIL, ALL)


###SETUP

BATCH_NAME=Brevipes_submission

SCRIPTS_DIR=~/Git_repos/BrevipesPopgen/Brevipes_submission
LOG_DIR=/home/mcf96392/logs/${BATCH_NAME}
WORKING_DIR=/scratch/mcf96392/${BATCH_NAME}

GENOME=~/Genomes/Maurantiacus/M_aurantiacus_v1_splitline_ordered.fasta
ANNOTATION=~/Genomes/Maurantiacus/MAUR_annotation_functional_submission.gff

NTHREADS=10
MEM=20

###LOG
cd $SCRIPTS_DIR
echo "${SLURM_JOB_ID},${SLURM_JOB_NAME},$(git rev-parse --short HEAD),$(date)" >> ${LOG_DIR}/joblog.txt

#folders
cd $WORKING_DIR/VCFs

touch vcf_counts.txt
for VCF in ${WORKING_DIR}/VCFs/*vcf.gz; do
	if ! grep -q $VCF vcf_counts.txt ; then
		COUNT=$(zcat $VCF | grep -v '^#' | wc -l)
		printf '%s\t%s\n' $VCF $COUNT >> vcf_counts.txt
	fi
done

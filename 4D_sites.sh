#!/bin/bash
#SBATCH --job-name=4D_sites	                        # Job name
#SBATCH --partition=batch		                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=10		                            # Number of cores per task
#SBATCH --mem=20gb			                                # Total memory for job
#SBATCH --time=72:00:00  		                            # Time limit hrs:min:sec
#SBATCH -o /home/mcf96392/logs/Brevipes_submission/%x.%j.out		# Standard output log
#SBATCH -e /home/mcf96392/logs/Brevipes_submission/%x.%j.err		#Standard error log
#SBATCH --mail-user=mcf96392@uga.edu                    # Where to send mail
#SBATCH --mail-type=ALL                            # Mail events (BEGIN, END, FAIL, ALL)

###perl script taken from https://github.com/tsackton/linked-selection/blob/master/misc_scripts/Identify_4D_Sites.pl

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

###MODULES
ml BCFtools/1.10.2-GCC-8.3.0

#folders
cd $WORKING_DIR
mkdir -p ${WORKING_DIR}/4D_sites
cd ${WORKING_DIR}/4D_sites


cp ${SCRIPTS_DIR}/inputs/codon_table.txt ${WORKING_DIR}/4D_sites/codon_table.txt
perl ${SCRIPTS_DIR}/Identify_4D_Sites.pl $ANNOTATION $GENOME > ${WORKING_DIR}/4D_sites/4Dsites.txt

for VCF in ${WORKING_DIR}/VCFs/${BATCH_NAME}.pass.*vcf.gz; do
	PREFIX=${VCF%.vcf.gz}
	if [ ! -f ${PREFIX}.4D.vcf.gz ]; then

		bcftools view -R ${WORKING_DIR}/4D_sites/4Dsites.txt -Oz \
			-o ${PREFIX}.4D.vcf.gz \
			$VCF
		bcftools index -t ${PREFIX}.4D.vcf.gz
	fi
done

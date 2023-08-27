#!/bin/bash
#SBATCH --job-name=phylo		                        # Job name
#SBATCH --partition=batch		                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=56		                            # Number of cores per task
#SBATCH --mem=60gb			                                # Total memory for job
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

NTHREADS=56
MEM=60

###LOG
cd $SCRIPTS_DIR
echo "${SLURM_JOB_ID},${SLURM_JOB_NAME},$(git rev-parse --short HEAD),$(date)" >> ${LOG_DIR}/joblog.txt

###SETUP DIRS
mkdir -p $WORKING_DIR
cd $WORKING_DIR
mkdir -p ${WORKING_DIR}/phylo

#MODULES
ml BCFtools/1.13-GCC-8.3.0
ml seqtk/1.3-GCC-8.3.0
ml ClustalW2/2.1-GCC-10.2.0
ml RAxML/8.2.12-foss-2019b-pthreads-avx

###full genome trees and fastas
VCF_TYPE=${BATCH_NAME}.pass.31calledSNPs
VCF=${WORKING_DIR}/VCFs/${VCF_TYPE}.vcf.gz

for sample in $(bcftools query -l $VCF); do
	printf "\n>%s\n" $sample >> ${WORKING_DIR}/phylo/${VCF_TYPE}.fa
	bcftools consensus -f $GENOME -I -a '?' -M 'N' -s $sample  $VCF |
		grep -v "^>" | tr -d '?' | tr -d '\n' >> ${WORKING_DIR}/phylo/${VCF_TYPE}.fa
done

###write partition file and asc file
variablesites=7728322
invariablesites=32658171
printf "[asc~p1.txt], DNA, SNPs = 1-%s\n" $variablesites > partition.txt
printf "%s\n" $invariablesites > p1.txt
###choose randomly from het alleles
seqtk randbase ${WORKING_DIR}/phylo/${VCF_TYPE}.fa > ${WORKING_DIR}/phylo/${VCF_TYPE}.randbase.fa

###run raxml with 1000 bootstraps in 10 reps
for iteration in $(seq 10); do
	sbatch ${SCRIPTS_DIR}/raxml_iter.sh $WORKING_DIR $VCF_TYPE $iteration
done

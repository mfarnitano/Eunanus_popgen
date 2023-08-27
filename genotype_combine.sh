#!/bin/bash
#SBATCH --job-name=genotype_combine                     # Job name
#SBATCH --partition=batch		                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=16		                            # Number of cores per task
#SBATCH --mem=120gb			                                # Total memory for job
#SBATCH --time=96:00:00  		                            # Time limit hrs:min:sec
#SBATCH -o /home/mcf96392/logs/Brevipes_submission/%x.%j.out		# Standard output log
#SBATCH -e /home/mcf96392/logs/Brevipes_submission/%x.%j.err		#Standard error log
#SBATCH --mail-user=mcf96392@uga.edu                    # Where to send mail
#SBATCH --mail-type=ALL                            # Mail events (BEGIN, END, FAIL, ALL)

###SETUP

BATCH_NAME=Brevipes_submission

SCRIPTS_DIR=~/Git_repos/BrevipesPopgen/Brevipes_submission
LOG_DIR=/home/mcf96392/logs/${BATCH_NAME}
WORKING_DIR=/scratch/mcf96392/${BATCH_NAME}

NTHREADS=16
MEM=120

GENOME=~/Genomes/Maurantiacus/M_aurantiacus_v1_splitline_ordered.fasta
GVCF_MAP=${WORKING_DIR}/gvcf_map.txt

###LOG
cd $SCRIPTS_DIR
echo "${SLURM_JOB_ID},${SLURM_JOB_NAME},$(git rev-parse --short HEAD),$(date)" >> ${LOG_DIR}/joblog.txt

###SETUP DIRS
mkdir -p $WORKING_DIR
cd $WORKING_DIR
mkdir -p VCFs

###MODULES
ml SAMtools/1.10-iccifort-2019.5.281
ml BWA/0.7.17-GCC-8.3.0
ml GATK/4.1.6.0-GCCcore-8.3.0-Java-1.8
module list

###create chr_list
printf "\n...creating chromosome list\n" | tee >(cat >&2)
if [ ! -f ${WORKING_DIR}/chr_positions.list ]; then
	head -n10 ${GENOME}.fai | awk '{print $1 ":1-"$2}' > ${WORKING_DIR}/chr_positions.list
fi

##merge gvcfs into database
printf "\n...merging gvcfs into database, selecting only major (chromosome) linkage groups\n" | tee >(cat >&2)
cd ${WORKING_DIR}/gvcfs
gatk --java-options "-Xmx${MEM}g" GenomicsDBImport \
	--sample-name-map ${GVCF_MAP} \
	--genomicsdb-workspace-path ${WORKING_DIR}/${BATCH_NAME}_vcf_db \
	-L ${WORKING_DIR}/chr_positions.list
cd $WORKING_DIR

###joint genotyping
printf "\n...Genotyping VCFs...outputing all-sites VCF with min call confidence = 30" | tee >(cat >&2)
gatk --java-options "-Xmx${MEM}g" GenotypeGVCFs \
	-R $GENOME \
	-V gendb://${WORKING_DIR}/${BATCH_NAME}_vcf_db \
	-stand-call-conf 30 \
	-all-sites \
	-O ${WORKING_DIR}/VCFs/${BATCH_NAME}.allsites.vcf.gz

printf "\n...Finished genotyping\n" | tee >(cat >&2)

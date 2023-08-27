#!/bin/bash
#SBATCH --job-name=genotype_each		                        # Job name
#SBATCH --partition=batch		                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=10		                            # Number of cores per task
#SBATCH --mem=40gb			                                # Total memory for job
#SBATCH --time=72:00:00  		                            # Time limit hrs:min:sec
#SBATCH -o /home/mcf96392/logs/Brevipes_submission/%x.%j.out		# Standard output log
#SBATCH -e /home/mcf96392/logs/Brevipes_submission/%x.%j.err		#Standard error log
#SBATCH --mail-user=mcf96392@uga.edu                    # Where to send mail
#SBATCH --mail-type=ALL                            # Mail events (BEGIN, END, FAIL, ALL)

###Creates gvcf for single sample from filtered bam, requires input variables
###See submission script genotype_wrapper.sh

###SETUP
SCRIPTS_DIR=$1
BATCH_NAME=$2
GENOME=$3
INPUT_DIR=$4
INPUT_PREFIX=$5
SEQ_BATCH=$6

LOG_DIR=/home/mcf96392/logs/${BATCH_NAME}
WORKING_DIR=/scratch/mcf96392/${BATCH_NAME}

NTHREADS=10
MEM=40

###LOG
cd $SCRIPTS_DIR
echo "${SLURM_JOB_ID},${SLURM_JOB_NAME},$(git rev-parse --short HEAD),$(date)" >> ${LOG_DIR}/joblog.txt

###SETUP DIRS
cd $WORKING_DIR
mkdir -p ${WORKING_DIR}/qualimap
mkdir -p ${WORKING_DIR}/genotyping_temp
mkdir -p ${WORKING_DIR}/gvcfs

###MODULES
ml SAMtools/1.10-iccifort-2019.5.281
ml GATK/4.1.6.0-GCCcore-8.3.0-Java-1.8
ml Qualimap/2.2.1-foss-2019b-R-3.6.2
module list

printf "\nstarting analysis for sample ${INPUT_PREFIX}\n" | tee >(cat >&2)

#qualimap to get coverage info about bam files
printf "\ngetting bam coverage info with qualimap\n" | tee >(cat >&2)
if [ ! -f ${WORKING_DIR}/qualimap/${INPUT_PREFIX}/genome_results.txt ]; then
	qualimap bamqc -bam ${INPUT_DIR}/${INPUT_PREFIX}.ar.fds.bam -c -outdir ${WORKING_DIR}/qualimap/${INPUT_PREFIX}
fi

###individual genotyping with GATK
printf "\n...adding read groups\n" | tee >(cat >&2)
gatk AddOrReplaceReadGroups \
	-I ${INPUT_DIR}/${INPUT_PREFIX}.ar.fds.bam \
	-O $WORKING_DIR/genotyping_temp/${INPUT_PREFIX}.ar.rg.bam \
	-RGID $INPUT_PREFIX \
	-LB ${SEQ_BATCH} \
	-PL illumina \
	-PU $INPUT_PREFIX \
	-SM $INPUT_PREFIX

printf "\n...indexing bam file\n" | tee >(cat >&2)
samtools index $WORKING_DIR/genotyping_temp/${INPUT_PREFIX}.ar.rg.bam

printf "\n...genotyping with gatk HaplotypeCaller\n" | tee >(cat >&2)
gatk --java-options "-Xmx${MEM}g" HaplotypeCaller  \
	-R $GENOME \
	-I $WORKING_DIR/genotyping_temp/${INPUT_PREFIX}.ar.rg.bam \
	-O $WORKING_DIR/gvcfs/${INPUT_PREFIX}.ar.g.vcf.gz \
	-ERC GVCF

printf "...Finished sample ${INPUT_PREFIX}\n" | tee >(cat >&2)

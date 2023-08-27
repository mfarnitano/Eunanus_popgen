#!/bin/bash
#SBATCH --job-name=pixy		                        # Job name
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

POPFILE=${SCRIPTS_DIR}/inputs/pixy_popfile.txt
INDPOPFILE=${SCRIPTS_DIR}/inputs/pixy_popfile_inds.txt

NTHREADS=10
MEM=20

###LOG
cd $SCRIPTS_DIR
echo "${SLURM_JOB_ID},${SLURM_JOB_NAME},$(git rev-parse --short HEAD),$(date)" >> ${LOG_DIR}/joblog.txt

#folders
cd $WORKING_DIR
mkdir -p ${WORKING_DIR}/pixy

###MODULES
ml BCFtools/1.10.2-GCC-8.3.0
ml tabix/0.2.6-GCCcore-8.3.0

#installation of pixy with conda (do once only)
# source ~/.bashrc
# conda create -n pixy
# conda activate pixy
# conda install -c conda-forge pixy

#conda setup
source ~/.bashrc
conda activate pixy
pixy --version | tee >(cat >&2)

printf "\nLooping through variant+invariant VCF files...\n" | tee >(cat >&2)
#for VCF in ${WORKING_DIR}/VCFs/${BATCH_NAME}.pass.allsites.fs.vcf.gz ${WORKING_DIR}/VCFs/${BATCH_NAME}.pass.*called.vcf.gz; do
for VCF_TYPE in pass.31called.4D pass.31called; do
	VCF=${WORKING_DIR}/VCFs/${BATCH_NAME}.${VCF_TYPE}.vcf.gz
	VCF_FILENAME=${BATCH_NAME}.${VCF_TYPE}.vcf.gz
	printf "\n...Checking file %s\n" $VCF_FILENAME | tee >(cat >&2)

	##index
	#tabix -f $VCF

	if [ ! -f ${WORKING_DIR}/pixy/${VCF_FILENAME%.vcf.gz}.10000_pi.txt ]; then
		printf "\n...Running pixy groupwise to calculate pi, fst, and dxy in 10000bp windows on sample %s\n" $VCF_FILENAME | tee >(cat >&2)
		pixy --stats pi fst dxy \
			--vcf $VCF \
			--window_size 10000 \
			--chromosomes "LG1,LG2,LG3,LG4,LG5,LG6,LG7,LG8,LG9,LG10" \
			--populations $POPFILE \
			--n_cores $NTHREADS \
			--output_folder ${WORKING_DIR}/pixy \
			--output_prefix ${VCF_FILENAME%.vcf.gz}.10000
	fi
	if [ ! -f ${WORKING_DIR}/pixy/${VCF_FILENAME%.vcf.gz}.ind.10000_pi.txt ]; then
		printf "\n...Running pixy individual-pairwise to calculate pi, fst, and dxy in 10000bp windows on sample %s\n" $VCF_FILENAME | tee >(cat >&2)
		pixy --stats pi fst dxy \
			--vcf $VCF \
			--window_size 10000 \
			--chromosomes "LG1,LG2,LG3,LG4,LG5,LG6,LG7,LG8,LG9,LG10" \
			--populations $INDPOPFILE \
			--n_cores $NTHREADS \
			--output_folder ${WORKING_DIR}/pixy \
			--output_prefix ${VCF_FILENAME%.vcf.gz}.ind.10000
	fi
	if [ ! -f ${WORKING_DIR}/pixy/${VCF_FILENAME%.vcf.gz}.1Mb_pi.txt ]; then
		printf "\n...Running pixy groupwise to calculate pi, fst, and dxy in 1Mbp windows on sample %s\n" $VCF_FILENAME | tee >(cat >&2)
		pixy --stats pi fst dxy \
			--vcf $VCF \
			--window_size 1000000 \
			--chromosomes "LG1,LG2,LG3,LG4,LG5,LG6,LG7,LG8,LG9,LG10" \
			--populations $POPFILE \
			--n_cores $NTHREADS \
			--output_folder ${WORKING_DIR}/pixy \
			--output_prefix ${VCF_FILENAME%.vcf.gz}.1Mb
	fi
	if [ ! -f ${WORKING_DIR}/pixy/${VCF_FILENAME%.vcf.gz}.ind.1Mb_pi.txt ]; then
		printf "\n...Running pixy individual-pairwise to calculate pi, fst, and dxy in 1Mbp windows on sample %s\n" $VCF_FILENAME | tee >(cat >&2)
		pixy --stats pi fst dxy \
			--vcf $VCF \
			--window_size 1000000 \
			--chromosomes "LG1,LG2,LG3,LG4,LG5,LG6,LG7,LG8,LG9,LG10" \
			--populations $INDPOPFILE \
			--n_cores $NTHREADS \
			--output_folder ${WORKING_DIR}/pixy \
			--output_prefix ${VCF_FILENAME%.vcf.gz}.ind.1Mb
	fi

done

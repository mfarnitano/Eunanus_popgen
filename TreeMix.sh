#!/bin/bash
#SBATCH --job-name=TreeMix                    # Job name
#SBATCH --partition=batch		                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=16		                            # Number of cores per task
#SBATCH --mem=80gb			                                # Total memory for job
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

#INPUT files
CLUSTERS=${SCRIPTS_DIR}/inputs/treemix_popfile.clust #3 columns, first two are sample names and third is pop

#Memory and threads
MEM=80
NTHREADS=16

###LOG
cd $SCRIPTS_DIR
echo "${PBS_JOBID},${PBS_JOBNAME},$(git rev-parse --short HEAD),$(date)" >> ${LOG_DIR}/joblog.txt

###MODULES
ml PLINK/1.9b_5-x86_64
ml VCFtools/0.1.16-GCC-8.3.0-Perl-5.30.0
ml TreeMix/1.13-foss-2019b
ml BCFtools/1.13-GCC-8.3.0

mkdir -p ${WORKING_DIR}/TreeMix
cd ${WORKING_DIR}/TreeMix
mkdir ${WORKING_DIR}/TreeMix/sampling

#get scripts
cp ${SCRIPTS_DIR}/plink2treemix.py $WORKING_DIR/TreeMix
cp ${SCRIPTS_DIR}/vcf2treemix.sh $WORKING_DIR/TreeMix
chmod +x $WORKING_DIR/TreeMix/plink2treemix.py
wget https://github.com/joanam/scripts/raw/master/ldPruning.sh
chmod +x ldPruning.sh

#run treemix for various migration models
for VCF_TYPE in pass.31calledSNPs; do
	VCF_PREFIX=${BATCH_NAME}.${VCF_TYPE}
	#remove 'LG' from chromosome names for PLINK
	zcat ${WORKING_DIR}/VCFs/${VCF_PREFIX}.vcf.gz | sed 's/LG//' |
		gzip -c > ${WORKING_DIR}/TreeMix/${VCF_PREFIX}.noLG.vcf.gz

	#remove missing data
	vcftools --gzvcf ${WORKING_DIR}/TreeMix/${VCF_PREFIX}.noLG.vcf.gz --max-missing 0 --recode --stdout | gzip > ${WORKING_DIR}/TreeMix/${VCF_PREFIX}.noLG.noN.vcf.gz

	###prune by LD, see https://speciationgenomics.github.io/Treemix/
	./ldPruning.sh ${WORKING_DIR}/TreeMix/${VCF_PREFIX}.noLG.noN.vcf.gz
	gzip ${WORKING_DIR}/TreeMix/${VCF_PREFIX}.noLG.noN.LDpruned.vcf

	#get treemix frequencies from vcf
	bash ${WORKING_DIR}/TreeMix/vcf2treemix.sh ${WORKING_DIR}/TreeMix/${VCF_PREFIX}.noLG.noN.LDpruned.vcf.gz $CLUSTERS

	#reset column format to newlines and tabs
	zcat ${WORKING_DIR}/TreeMix/${VCF_PREFIX}.noLG.noN.LDpruned.treemix.frq.gz | awk 'BEGIN {RS="\n\n"; FS="\n"; ORS="\n"; OFS="\t"} {$1=$1;print}' | gzip -c > ${WORKING_DIR}/TreeMix/${VCF_PREFIX}.ready.treemix.frq.gz

	for i in {1..10}; do

		#sample SNPs from file
		bcftools view -h ${WORKING_DIR}/TreeMix/${VCF_PREFIX}.noLG.noN.vcf.gz > ${WORKING_DIR}/TreeMix/${VCF_PREFIX}.${i}.vcf
		zcat ${WORKING_DIR}/TreeMix/${VCF_PREFIX}.noLG.noN.vcf.gz | grep -v "^#" | shuf -n 150000 >> ${WORKING_DIR}/TreeMix/${VCF_PREFIX}.${i}.vcf

		#get treemix frequencies
		bash ${WORKING_DIR}/TreeMix/vcf2treemix.sh ${WORKING_DIR}/TreeMix/${VCF_PREFIX}.${i}.vcf $CLUSTERS

		#reset column format to newlines and tabs
		zcat ${WORKING_DIR}/TreeMix/${VCF_PREFIX}.${i}.treemix.frq.gz | awk 'BEGIN {RS="\n\n"; FS="\n"; ORS="\n"; OFS="\t"} {$1=$1;print}' | gzip -c > ${WORKING_DIR}/TreeMix/${VCF_PREFIX}.${i}.ready.treemix.frq.gz


		for m in {1..8}; do
			treemix -i ${WORKING_DIR}/TreeMix/${VCF_PREFIX}.${i}.ready.treemix.frq.gz \
				-o ${WORKING_DIR}/TreeMix/sampling/${VCF_PREFIX}.i${i}.m${m} \
				-global -root AURA -m $m -k 1000
		done
	done
done

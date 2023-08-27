#!/bin/bash
#SBATCH --job-name=genetrees		                        # Job name
#SBATCH --partition=batch		                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=10		                            # Number of cores per task
#SBATCH --mem=60gb			                                # Total memory for job
#SBATCH --time=120:00:00  		                            # Time limit hrs:min:sec
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
MEM=60

###LOG
cd $SCRIPTS_DIR
echo "${SLURM_JOB_ID},${SLURM_JOB_NAME},$(git rev-parse --short HEAD),$(date)" >> ${LOG_DIR}/joblog.txt

###SETUP DIRS
mkdir -p $WORKING_DIR
cd $WORKING_DIR
mkdir -p ${WORKING_DIR}/trees_by_genes
mkdir -p ${WORKING_DIR}/trees_by_genes/nj
mkdir -p ${WORKING_DIR}/trees_by_genes/raxml

###MODULES
ml SAMtools/1.10-GCC-8.3.0
ml BCFtools/1.13-GCC-8.3.0
ml RAxML/8.2.12-foss-2019b-pthreads-avx
ml seqtk/1.3-GCC-8.3.0

###split fasta
VCF_TYPE=pass.SNPs.fs
VCF=${WORKING_DIR}/VCFs/${BATCH_NAME}.${VCF_TYPE}.vcf.gz

###get gene coordinates
cat $ANNOTATION | grep "^LG" | awk '$3 == "gene" {print $1 ":" $4 "-" $5}' > ${WORKING_DIR}/trees_by_genes/gene_coordinates.txt

###do it with vcf consensus to randomly choose heterozygous alleles
printf 'coordinates\tminapplied\n' > ${WORKING_DIR}/trees_by_genes/minapplied_summary.txt
while read -r coordinates; do
	name=$(echo $coordinates | sed 's/:/_/')

	###write fasta
	for sample in $(bcftools query -l $VCF); do
		printf "\n>%s\n" $sample >> ${WORKING_DIR}/trees_by_genes/${name}.cons.fa
		samtools faidx $GENOME $coordinates | bcftools consensus -I -a '?' -M 'N' -s $sample  $VCF |
			grep -v "^>" | tr -d '?' | tr -d '\n' >> ${WORKING_DIR}/trees_by_genes/${name}.cons.fa
	done 2> ${WORKING_DIR}/trees_by_genes/${name}.applied.txt

	###get minimum number of variants used
	minapplied=$(cat ${WORKING_DIR}/trees_by_genes/${name}.applied.txt | grep -o '[0-9]*' | awk 'NR == 1 {min = $0} NR > 1 && $0 < min {min = $0} END {print min}')
	printf '%s\t%s\n' $coordinates $minapplied >> ${WORKING_DIR}/trees_by_genes/minapplied_summary.txt

	###randomly select alleles
	seqtk randbase ${WORKING_DIR}/trees_by_genes/${name}.cons.fa > ${WORKING_DIR}/trees_by_genes/${name}.randbase.fa

done < ${WORKING_DIR}/trees_by_genes/gene_coordinates.txt

###separate ML tree for each
while read -r coordinates; do
	name=$(echo $coordinates | sed 's/:/_/')
	###raxml ML tree
	raxmlHPC-PTHREADS-AVX -T $NTHREADS -s ${WORKING_DIR}/trees_by_genes/${name}.randbase.fa \
		-w ${WORKING_DIR}/trees_by_genes/raxml -n ${name}.randbase \
		-m GTRGAMMA -f a -N 10 -x 12345 -p 12345
done < ${WORKING_DIR}/trees_by_genes/gene_coordinates.txt

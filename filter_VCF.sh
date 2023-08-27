#!/bin/bash
#SBATCH --job-name=filter_VCF                     # Job name
#SBATCH --partition=batch		                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=8		                            # Number of cores per task
#SBATCH --mem=80gb			                                # Total memory for job
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

NTHREADS=8
MEM=80

###LOG
cd $SCRIPTS_DIR
echo "${SLURM_JOB_ID},${SLURM_JOB_NAME},$(git rev-parse --short HEAD),$(date)" >> ${LOG_DIR}/joblog.txt

###MODULES
ml SAMtools/1.10-iccifort-2019.5.281
ml GATK/4.1.6.0-GCCcore-8.3.0-Java-1.8
ml Qualimap/2.2.1-foss-2019b-R-3.6.2
module list

###folder setup
cd $WORKING_DIR
mkdir -p ${WORKING_DIR}/qualimap
mkdir -p ${WORKING_DIR}/filter_temp

###preparing genome
GENOME=~/Genomes/Maurantiacus/M_aurantiacus_v1_splitline_ordered.fasta
printf "\n...creating sequence dictionary if needed\n" | tee >(cat >&2)
if [ ! -f ${GENOME%.fasta}.dict ]; then
	gatk CreateSequenceDictionary -R $GENOME
fi

###input VCF and prefixes
INPUT_VCF=${WORKING_DIR}/VCFs/${BATCH_NAME}.allsites.vcf.gz

###loop through bam files
for BAM in ${WORKING_DIR}/bams/*.bam; do

	INPUT_FILENAME=${BAM##*/}
	INPUT_PREFIX=${INPUT_FILENAME%.ar.fds.bam}

	###re-index bam files
	if [ ! -f ${BAM}.bai ]; then
		samtools index -b $BAM
	fi

	###get summary genome coverage info
	if [ ! -f ${WORKING_DIR}/qualimap/${INPUT_PREFIX}/genome_results.txt ]; then
		qualimap bamqc -bam $BAM -c -outdir ${WORKING_DIR}/qualimap/${INPUT_PREFIX}
	fi
	meandepth=$(grep "mean coverageData" ${WORKING_DIR}/qualimap/${INPUT_PREFIX}/genome_results.txt | sed -E 's/.*\s([0-9.]*)X/\1/')
	stddepth=$(grep "std coverageData" ${WORKING_DIR}/qualimap/${INPUT_PREFIX}/genome_results.txt | sed -E 's/.*\s([0-9.]*)X/\1/')
	vardepth=$(echo $stddepth "^" 2 | bc)
	maxdepth=$(echo 2 "*" $stddepth + $meandepth | bc)

	printf "%s\t%s\t%s\t%s\t%s\n" $INPUT_PREFIX $meandepth $stddepth $vardepth $maxdepth >> ${WORKING_DIR}/ar.bam_coverage_summary.txt

done

#get global min and max DP for filtering
nsamples=$(cat ${WORKING_DIR}/ar.bam_coverage_summary.txt | wc -l)
cat ${WORKING_DIR}/ar.bam_coverage_summary.txt | awk -v n=$nsamples \
	'{means+=$2; stds+=$3; vars+=$4; maxes+=$5} END{print n, means, stds, vars, maxes, means/n, stds/n, vars/n, maxes/n, sqrt(vars), means/3, 2*sqrt(vars)+means}' > \
	${WORKING_DIR}/ar.bam_coverage_grand_means.txt

minDP=$(cat ${WORKING_DIR}/ar.bam_coverage_grand_means.txt | tr -s ' ' | cut -d' ' -f11) #1/3 of mean total depth
maxDP=$(cat ${WORKING_DIR}/ar.bam_coverage_grand_means.txt | tr -s ' ' | cut -d' ' -f12) #mean + 2*sqrt(sum(var depth))

### Filter SNPs
printf "\nStarting SNP filtration process\n" | tee >(cat >&2)
printf "\n...selecting variants\n" | tee >(cat >&2)
gatk --java-options "-Xmx${MEM}g" SelectVariants \
  -V ${INPUT_VCF} \
  -select-type SNP \
  --restrict-alleles-to BIALLELIC \
  -O ${WORKING_DIR}/filter_temp/${BATCH_NAME}.SNPs.vcf.gz

printf "\n...filtering\n" | tee >(cat >&2)
gatk --java-options "-Xmx${MEM}g" VariantFiltration \
  -V ${WORKING_DIR}/filter_temp/${BATCH_NAME}.SNPs.vcf.gz \
	-filter "DP < $minDP" --filter-name "minDP_$minDP" \
	-filter "DP > $maxDP" --filter-name "maxDP_$maxDP" \
  -filter "QD < 2.0" --filter-name "QD2" \
  -filter "QUAL < 40.0" --filter-name "QUAL40" \
  -filter "SOR > 3.0" --filter-name "SOR4" \
  -filter "FS > 60.0" --filter-name "FS60" \
  -filter "MQ < 40.0" --filter-name "MQ40" \
  -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
  -filter "ReadPosRankSum < -12.5" --filter-name "ReadPosRankSum-10" \
  -filter "ReadPosRankSum > 12.5" --filter-name "ReadPosRankSum10" \
  -O ${WORKING_DIR}/filter_temp/${BATCH_NAME}.SNPs.f.vcf.gz \
  --verbosity ERROR

printf "\n...sorting\n" | tee >(cat >&2)
gatk --java-options "-Xmx${MEM}g" SortVcf \
  -I ${WORKING_DIR}/filter_temp/${BATCH_NAME}.SNPs.f.vcf.gz \
  -SD ${GENOME%.fasta}.dict \
  -O ${WORKING_DIR}/filter_temp/${BATCH_NAME}.SNPs.fs.vcf.gz

printf "\n...selecting sites that pass filter\n" | tee >(cat >&2)
gatk --java-options "-Xmx${MEM}g" SelectVariants \
	-V ${WORKING_DIR}/filter_temp/${BATCH_NAME}.SNPs.fs.vcf.gz \
	--exclude-filtered \
	-O ${WORKING_DIR}/VCFs/${BATCH_NAME}.pass.SNPs.fs.vcf.gz

### Filter Invariant Sites
printf "\nStarting INVARIANT filtration process\n" | tee >(cat >&2)
printf "\n...selecting variants" | tee >(cat >&2)
gatk --java-options "-Xmx${MEM}g" SelectVariants \
  -V ${INPUT_VCF} \
  -select-type NO_VARIATION \
  -O ${WORKING_DIR}/filter_temp/${BATCH_NAME}.INVT.vcf.gz

printf "\n...filtering\n" | tee >(cat >&2)
gatk --java-options "-Xmx${MEM}g" VariantFiltration \
  -V ${WORKING_DIR}/filter_temp/${BATCH_NAME}.INVT.vcf.gz \
	-filter "DP < $minDP" --filter-name "minDP_$minDP" \
	-filter "DP > $maxDP" --filter-name "maxDP_$maxDP" \
  -filter "QD < 2.0" --filter-name "QD2" \
  -filter "SOR > 3.0" --filter-name "SOR4" \
  -filter "MQ < 40.0" --filter-name "MQ40" \
  -O ${WORKING_DIR}/filter_temp/${BATCH_NAME}.INVT.f.vcf.gz \
  --verbosity ERROR

printf "\n...sorting\n" | tee >(cat >&2)
gatk --java-options "-Xmx${MEM}g" SortVcf \
  -I ${WORKING_DIR}/filter_temp/${BATCH_NAME}.INVT.f.vcf.gz \
  -SD ${GENOME%.fasta}.dict \
  -O ${WORKING_DIR}/filter_temp/${BATCH_NAME}.INVT.fs.vcf.gz

printf "\n...selecting sites that pass filter\n" | tee >(cat >&2)
gatk --java-options "-Xmx${MEM}g" SelectVariants \
	-V ${WORKING_DIR}/filter_temp/${BATCH_NAME}.INVT.fs.vcf.gz \
	--exclude-filtered \
	-O ${WORKING_DIR}/VCFs/${BATCH_NAME}.pass.INVT.fs.vcf.gz

### Merge all Sites
printf "\nMerging all sites together\n" | tee >(cat >&2)
gatk --java-options "-Xmx${MEM}g" MergeVcfs \
	-I ${WORKING_DIR}/VCFs/${BATCH_NAME}.pass.SNPs.fs.vcf.gz \
	-I ${WORKING_DIR}/VCFs/${BATCH_NAME}.pass.INVT.fs.vcf.gz \
	-O ${WORKING_DIR}/filter_temp/${BATCH_NAME}.pass.allsites.f.vcf.gz

printf "\n...Sorting\n" | tee >(cat >&2)
gatk --java-options "-Xmx${MEM}g" SortVcf \
	-I ${WORKING_DIR}/filter_temp/${BATCH_NAME}.pass.allsites.f.vcf.gz \
	-SD ${GENOME%.fasta}.dict \
	-O ${WORKING_DIR}/VCFs/${BATCH_NAME}.pass.allsites.fs.vcf.gz

### Select sites called in all samples
printf "/n...Selecting only sites called in all samples\n" | tee >(cat >&2)
gatk --java-options "-Xmx${MEM}g" SelectVariants \
	-V ${WORKING_DIR}/VCFs/${BATCH_NAME}.pass.allsites.fs.vcf.gz \
	--max-nocall-number 0 \
	--exclude-filtered \
	-O $WORKING_DIR/VCFs/${BATCH_NAME}.pass.allcalled.vcf.gz

printf "/n...Selecting only sites called in 80% (31+) of samples\n" | tee >(cat >&2)
gatk --java-options "-Xmx${MEM}g" SelectVariants \
	-V ${WORKING_DIR}/VCFs/${BATCH_NAME}.pass.allsites.fs.vcf.gz \
	--max-nocall-number 7 \
	--exclude-filtered \
	-O $WORKING_DIR/VCFs/${BATCH_NAME}.pass.31called.vcf.gz

printf "/n...Selecting only SNPs called in all samples\n" | tee >(cat >&2)
gatk --java-options "-Xmx${MEM}g" SelectVariants \
	-V ${WORKING_DIR}/VCFs/${BATCH_NAME}.pass.SNPs.fs.vcf.gz \
	--max-nocall-number 0 \
	--exclude-filtered \
	-O $WORKING_DIR/VCFs/${BATCH_NAME}.pass.allcalledSNPs.vcf.gz

printf "/n...Selecting only SNPs called in 80% (31+) of samples\n" | tee >(cat >&2)
gatk --java-options "-Xmx${MEM}g" SelectVariants \
	-V ${WORKING_DIR}/VCFs/${BATCH_NAME}.pass.SNPs.fs.vcf.gz \
	--max-nocall-number 7 \
	--exclude-filtered \
	-O $WORKING_DIR/VCFs/${BATCH_NAME}.pass.31calledSNPs.vcf.gz

printf "\nCompleted VCF filtration\n" | tee >(cat >&2)

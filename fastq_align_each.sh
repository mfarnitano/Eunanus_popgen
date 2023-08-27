#!/bin/bash
#SBATCH --job-name=fastq_align_each		                        # Job name
#SBATCH --partition=batch		                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=10		                            # Number of cores per task
#SBATCH --mem=20gb			                                # Total memory for job
#SBATCH --time=72:00:00  		                            # Time limit hrs:min:sec
#SBATCH -o /home/mcf96392/logs/Brevipes_submission/%x.%j.out		# Standard output log
#SBATCH -e /home/mcf96392/logs/Brevipes_submission/%x.%j.err		#Standard error log
#SBATCH --mail-user=mcf96392@uga.edu                    # Where to send mail
#SBATCH --mail-type=ALL                            # Mail events (BEGIN, END, FAIL, ALL)

###Creates filtered bam file for single sample from raw fastq, requires input variables

###SETUP
SCRIPTS_DIR=$1
BATCH_NAME=$2
GENOME=$3
INPUT_DIR=$4
INPUT_PREFIX=$5

LOG_DIR=/home/mcf96392/logs/${BATCH_NAME}
WORKING_DIR=/scratch/mcf96392/${BATCH_NAME}

NTHREADS=10
MEM=20

###LOG
cd $SCRIPTS_DIR
echo "${SLURM_JOB_ID},${SLURM_JOB_NAME},$(git rev-parse --short HEAD),$(date)" >> ${LOG_DIR}/joblog.txt

###SETUP DIRS
mkdir -p $WORKING_DIR
cd $WORKING_DIR
mkdir -p ${WORKING_DIR}/fastqc
mkdir -p ${WORKING_DIR}/fastq_align_temp
mkdir -p ${WORKING_DIR}/bams

###MODULES
ml FastQC/0.11.9-Java-11
ml Trimmomatic/0.39-Java-1.8.0_144
ml BWA/0.7.17-GCC-8.3.0
ml SAMtools/1.10-iccifort-2019.5.281
ml picard/2.21.6-Java-11
ml Qualimap/2.2.1-foss-2019b-R-3.6.2
module list

printf "\nstarting analysis for sample ${INPUT_PREFIX}\n" | tee >(cat >&2)

###fastqc
printf "...running fastqc" | tee >(cat >&2)
fastqc ${INPUT_DIR}/${INPUT_PREFIX}_R1.fastq.gz ${INPUT_DIR}/${INPUT_PREFIX}_R2.fastq.gz
mv ${INPUT_DIR}/${INPUT_PREFIX}_*_fastqc* ${WORKING_DIR}/fastqc/

###find overrepresented sequences
printf "\n...finding overrepresented sequences\n" | tee >(cat >&2)
for i in ${WORKING_DIR}/fastqc/${INPUT_PREFIX}_*.zip; do
	fdir=${i%.zip}
  unzip -o -q $i -d ${WORKING_DIR}/fastqc
  grep -A 2 ">>Overrepresented sequences" ${fdir}/fastqc_data.txt |
    tail -n1 |
    sed 's/ /_/' >> ${WORKING_DIR}/fastqc/${INPUT_PREFIX}_tempfile.txt
done

paste ${WORKING_DIR}/fastqc/${INPUT_PREFIX}_tempfile.txt ${WORKING_DIR}/fastqc/${INPUT_PREFIX}_tempfile.txt |
  cut -f4,5 |
  sort -u |
  sed -e 's/^/>/' |
  sed 's/\t/\n/' > ${WORKING_DIR}/fastqc/${INPUT_PREFIX}_overrepresented.fa

###Trim
printf "\n...trimming with Trimmomatic\n" | tee >(cat >&2)
java -jar ${EBROOTTRIMMOMATIC}/trimmomatic-0.39.jar \
	PE -threads $NTHREADS \
	${INPUT_DIR}/${INPUT_PREFIX}_R1.fastq.gz ${INPUT_DIR}/${INPUT_PREFIX}_R2.fastq.gz \
	-baseout ${WORKING_DIR}/fastq_align_temp/${INPUT_PREFIX}.fq.gz \
	ILLUMINACLIP:${WORKING_DIR}/fastqc/${INPUT_PREFIX}_overrepresented.fa:1:30:15 \
	SLIDINGWINDOW:5:20 \
	MINLEN:30

###map to reference and sort
printf "\n...mapping to reference with bwa mem and sorting with samtools\n" | tee >(cat >&2)
bwa mem -t $NTHREADS $GENOME \
	${WORKING_DIR}/fastq_align_temp/${INPUT_PREFIX}_1P.fq.gz ${WORKING_DIR}/fastq_align_temp/${INPUT_PREFIX}_2P.fq.gz |
	samtools sort --threads $NTHREADS -O bam -o ${WORKING_DIR}/fastq_align_temp/${INPUT_PREFIX}.ar.s.bam -

###index bam file
printf "\n...indexing bam file\n" | tee >(cat >&2)
samtools index ${WORKING_DIR}/fastq_align_temp/${INPUT_PREFIX}.ar.s.bam

###mark and remove duplicates
printf "\n...marking and removing PCR duplicates with picard\n" | tee >(cat >&2)
java -jar ${EBROOTPICARD}/picard.jar \
	MarkDuplicates \
	I=${WORKING_DIR}/fastq_align_temp/${INPUT_PREFIX}.ar.s.bam \
	O=${WORKING_DIR}/fastq_align_temp/${INPUT_PREFIX}.ar.ds.bam \
	M=${WORKING_DIR}/fastq_align_temp/${INPUT_PREFIX}.dedupe-metrics.txt \
	REMOVE_DUPLICATES=true

###filter bam file and index
printf "\n...filtering bamfile for MAPQ>=29, both reads mapped and properly paired, passes platform QC\n" | tee >(cat >&2)
samtools view --threads $NTHREADS -q 29 -f 2 -F 524 -b ${WORKING_DIR}/fastq_align_temp/${INPUT_PREFIX}.ar.ds.bam > ${WORKING_DIR}/bams/${INPUT_PREFIX}.ar.fds.bam

printf "\n...indexing filtered bam file\n" | tee >(cat >&2)
samtools index ${WORKING_DIR}/bams/${INPUT_PREFIX}.ar.fds.bam

printf "\n...generating coverage summary stats with qualimap\n" | tee >(cat >&2)
qualimap bamqc -bam ${WORKING_DIR}/bams/${INPUT_PREFIX}.ar.fds.bam -c -outdir ${WORKING_DIR}/qualimap/${INPUT_PREFIX}

printf "\n...completed sample ${INPUT_PREFIX}" | tee >(cat >&2)

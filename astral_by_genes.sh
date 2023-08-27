#!/bin/bash
#SBATCH --job-name=ASTRAL                    # Job name
#SBATCH --partition=batch		                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=16		                            # Number of cores per task
#SBATCH --mem=80gb			                                # Total memory for job
#SBATCH --time=12:00:00  		                            # Time limit hrs:min:sec
#SBATCH -o /home/mcf96392/logs/Brevipes_FINAL/%x.%j.out		# Standard output log
#SBATCH -e /home/mcf96392/logs/Brevipes_FINAL/%x.%j.err		#Standard error log
#SBATCH --mail-user=mcf96392@uga.edu                    # Where to send mail
#SBATCH --mail-type=ALL                            # Mail events (BEGIN, END, FAIL, ALL)

BATCH_NAME=Brevipes_FINAL

SCRIPTS_DIR=~/Git_repos/BrevipesPopgen/Brevipes_submission
LOG_DIR=/home/mcf96392/logs/${BATCH_NAME}
WORKING_DIR=/scratch/mcf96392/${BATCH_NAME}

#Memory and threads
MEM=80
NTHREADS=16

###LOG
cd $SCRIPTS_DIR
echo "${PBS_JOBID},${PBS_JOBNAME},$(git rev-parse --short HEAD),$(date)" >> ${LOG_DIR}/joblog.txt

###MODULES
ml ASTRAL/5.6.1-Java-1.8.0_144

###folders
mkdir -p ${WORKING_DIR}/trees_by_genes/ASTRAL

###gene trees contatenated in twisst_by_genes.sh script, run that first!

###Run ASTRAL
printf "\nStarting ASTRAL analyses...\n" | tee >(cat >&2)
java -jar $EBROOTASTRAL/astral.5.6.1.jar -t 1 -i ${WORKING_DIR}/trees_by_genes/twisst/completeMLbygenes.tre -o ${WORKING_DIR}/trees_by_genes/ASTRAL/Q_astral_by_genes.tre
#
java -jar $EBROOTASTRAL/astral.5.6.1.jar -t 2 -i ${WORKING_DIR}/trees_by_genes/twisst/completeMLbygenes.tre -o ${WORKING_DIR}/trees_by_genes/ASTRAL/Qc_astral_by_genes.tre
#
java -jar $EBROOTASTRAL/astral.5.6.1.jar -t 3 -i ${WORKING_DIR}/trees_by_genes/twisst/completeMLbygenes.tre -o ${WORKING_DIR}/trees_by_genes/ASTRAL/P_astral_by_genes.tre

###Run ASTRAL gene-bootstrap analysis
printf "\nStarting gene bootstrapping analysis...\n" | tee >(cat >&2)
java -jar $EBROOTASTRAL/astral.5.6.1.jar --gene-only -r 1000 -i ${WORKING_DIR}/trees_by_genes/twisst/completeMLbygenes.tre -o ${WORKING_DIR}/trees_by_genes/ASTRAL/Gene_BS_astral_by_genes.tre


printf "\nFinished all...\n" | tee >(cat >&2)

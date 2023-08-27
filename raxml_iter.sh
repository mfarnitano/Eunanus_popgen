#!/bin/bash
#SBATCH --job-name=raxml_iter		                        # Job name
#SBATCH --partition=batch		                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=56		                            # Number of cores per task
#SBATCH --mem=60gb			                                # Total memory for job
#SBATCH --time=168:00:00  		                            # Time limit hrs:min:sec
#SBATCH -o /home/mcf96392/logs/Brevipes_submission/%x.%j.out		# Standard output log
#SBATCH -e /home/mcf96392/logs/Brevipes_submission/%x.%j.err		#Standard error log
#SBATCH --mail-user=mcf96392@uga.edu                    # Where to send mail
#SBATCH --mail-type=ALL                            # Mail events (BEGIN, END, FAIL, ALL)

###SETUP

#run as sbatch ${SCRIPTS_DIR}/raxml_iter.sh $WORKING_DIR $VCF_TYPE $iteration
WORKING_DIR=$1
VCF_TYPE=$2
iteration=$3

NTHREADS=56

ml RAxML/8.2.12-foss-2019b-pthreads-avx

seed=$(echo "$iteration * 1234" | bc)
printf "\nstarting run random seed=%s\n" $seed | tee >(cat >&2)
raxmlHPC-PTHREADS-AVX -T $NTHREADS -s ${WORKING_DIR}/phylo/${VCF_TYPE}.randbase.fa \
	-w ${WORKING_DIR}/phylo/ \
	-n ${VCF_TYPE}.randbase.rep${iteration} \
	-m ASC_GTRGAMMA -f a -N 1000 -x $seed -p $seed --asc-corr=felsenstein -q partition.txt

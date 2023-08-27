#!/bin/bash
#SBATCH --job-name=twisst                    # Job name
#SBATCH --partition=batch		                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=16		                            # Number of cores per task
#SBATCH --mem=80gb			                                # Total memory for job
#SBATCH --time=6:00:00  		                            # Time limit hrs:min:sec
#SBATCH -o /home/mcf96392/logs/Brevipes_submission/%x.%j.out		# Standard output log
#SBATCH -e /home/mcf96392/logs/Brevipes_submission/%x.%j.err		#Standard error log
#SBATCH --mail-user=mcf96392@uga.edu                    # Where to send mail
#SBATCH --mail-type=ALL                            # Mail events (BEGIN, END, FAIL, ALL)

BATCH_NAME=Brevipes_submission

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
source ~/.bashrc
conda activate twisst
###conda install ete3
###conda install numpy
###git clone https://github.com/simonhmartin/twisst.git

###folders
mkdir -p ${WORKING_DIR}/trees_by_genes/twisst

rm -f ${WORKING_DIR}/trees_by_genes/twisst/completeMLbygenes.tre
for i in ${WORKING_DIR}/trees_by_genes/raxml/RAxML_bestTree*.randbase; do
	name=${i#*.}
	if [ -s $i ]; then
		if [ ! -f ${WORKING_DIR}/trees_by_genes/${name}.fa.reduced ]; then
			#cat $i | tr -d "\n" | sed 's/;/;\n/' >> ${WORKING_DIR}/trees_by_genes/twisst/allnjbygenes.tre
			cat $i >> ${WORKING_DIR}/trees_by_genes/twisst/completeMLbygenes.tre
		fi
	fi
done

function run_twisst () { #name P1 P2 P3 OUT groupsfile
	python ~/apps/twisst/twisst.py -t ${WORKING_DIR}/trees_by_genes/twisst/completeMLbygenes.tre	\
		--outputTopos ${WORKING_DIR}/trees_by_genes/twisst/${1}_by_genes_ML_topo.txt \
		-w ${WORKING_DIR}/trees_by_genes/twisst/${1}_by_genes_ML_weights.txt \
		-g $2 -g $3 -g $4 -g $5 \
		--groupsFile ${SCRIPTS_DIR}/inputs/${6} --method complete
}

###using AURA as outgroup, all species-level trios
run_twisst BJFA BREV JOHN FREM AURA twisst_groups_main.txt
run_twisst BJSA BREV JOHN SESP AURA twisst_groups_main.txt
run_twisst BJCA BREV JOHN CONS AURA twisst_groups_main.txt
run_twisst BJNA BREV JOHN NANU AURA twisst_groups_main.txt

run_twisst BFSA BREV FREM SESP AURA twisst_groups_main.txt
run_twisst BFCA BREV FREM CONS AURA twisst_groups_main.txt
run_twisst BFNA BREV FREM NANU AURA twisst_groups_main.txt

run_twisst JFSA JOHN FREM SESP AURA twisst_groups_main.txt
run_twisst JFCA JOHN FREM CONS AURA twisst_groups_main.txt
run_twisst JFNA JOHN FREM NANU AURA twisst_groups_main.txt

run_twisst BSCA BREV SESP CONS AURA twisst_groups_main.txt
run_twisst JSCA JOHN SESP CONS AURA twisst_groups_main.txt
run_twisst FSCA FREM SESP CONS AURA twisst_groups_main.txt

run_twisst BSNA BREV SESP NANU AURA twisst_groups_main.txt
run_twisst BCNA BREV CONS NANU AURA twisst_groups_main.txt

run_twisst JSNA JOHN SESP NANU AURA twisst_groups_main.txt
run_twisst JCNA JOHN CONS NANU AURA twisst_groups_main.txt

run_twisst FSNA FREM SESP NANU AURA twisst_groups_main.txt
run_twisst FCNA FREM CONS NANU AURA twisst_groups_main.txt

run_twisst CSNA CONS SESP NANU AURA twisst_groups_main.txt

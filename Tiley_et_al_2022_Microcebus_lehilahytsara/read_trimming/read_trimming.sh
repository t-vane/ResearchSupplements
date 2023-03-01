#!/bin/sh
#SBATCH -p medium40

set -euo pipefail

source /home/nibtve93/.bashrc

################################################################################
#### SET-UP ####
################################################################################
## Software:
TRIMMOMATIC=/home/nibtve93/software/Trimmomatic-0.39/trimmomatic-0.39.jar

## Command-line args:
MODE=$1
IN_FILE=$2
NT=$3
IN_DIR=$4
OUT_DIR=$5
ADAPTERS=$6
INDV=$(sed -n "$SLURM_ARRAY_TASK_ID"p $IN_FILE)

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### read_trimming.sh: Starting script."
echo -e "#### read_trimming.sh: Paired-end or single-end: $MODE"
echo -e "#### read_trimming.sh: File with individuals: $IN_FILE"
echo -e "#### read_trimming.sh: Number of threads: $NT"
echo -e "#### read_trimming.sh: Directory with raw reads: $IN_DIR"
echo -e "#### read_trimming.sh: Directory for trimmed reads: $OUT_DIR"
echo -e "#### read_trimming.sh: Adapter file: $ADAPTERS"
echo -e "#### reference_mapping.sh: Individual: $INDV \n\n"

################################################################################
#### TRIM RAW READS ####
################################################################################
if [[ $MODE == "PE" ]]
then
	echo -e "#### read_trimming.sh: Trimming paired-end reads for individual $INDV ...\n"
	java -jar $TRIMMOMATIC PE -threads $NT -phred33 $IN_DIR/$INDV.1.fq.gz $IN_DIR/$INDV.2.fq.gz $OUT_DIR/$INDV.trimmed.1.fq.gz $OUT_DIR/$INDV.1.rem.fq.gz $OUT_DIR/$INDV.trimmed.2.fq.gz $OUT_DIR/$INDV.2.rem.fq.gz \
		ILLUMINACLIP:$ADAPTERS:2:30:10 AVGQUAL:20 SLIDINGWINDOW:4:15 LEADING:3 TRAILING:3 MINLEN:60 
elif [[ $MODE == "SE" ]]
then
	echo -e "#### read_trimming.sh: Trimming single-end reads for individual $INDV ...\n"
	java -jar $TRIMMOMATIC SE -threads $NT -phred33 $IN_DIR/$INDV.1.fq.gz $OUT_DIR/$INDV.trimmed.1.fq.gz ILLUMINACLIP:$ADAPTERS:2:30:10 AVGQUAL:20 SLIDINGWINDOW:4:15 LEADING:3 TRAILING:3 MINLEN:60
else
	echo -e "#### read_trimming.sh: Invalid sequencing mode provided - only PE and SE allowed. ...\n"
	exit 1
fi

## Report:
echo -e "\n#### read_trimming.sh: Done with script."
date

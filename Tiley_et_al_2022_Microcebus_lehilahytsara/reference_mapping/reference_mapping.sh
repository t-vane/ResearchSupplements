#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
# BWA needs to be included in $PATH (v7.0.17; https://github.com/lh3/bwa)

## Command-line args:
MODE=$1
NT=$2
INDEX=$3
IN_DIR=$4
BAM_DIR=$5
IN_FILE=$6
INDV=$(sed -n "$SLURM_ARRAY_TASK_ID"p $IN_FILE)
READGROUP="@RG\tID:group1\tSM:$INDV\tPL:illumina\tLB:lib1"

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### reference_mapping.sh: Starting script."
echo -e "#### reference_mapping.sh: Paired-end or single-end: $MODE"
echo -e "#### reference_mapping.sh: Number of threads: $NT"
echo -e "#### reference_mapping.sh: Reference genome index: $INDEX"
echo -e "#### reference_mapping.sh: Directory with trimmed reads: $IN_DIR"
echo -e "#### reference_mapping.sh: Directory for BAM files: $BAM_DIR"
echo -e "#### reference_mapping.sh: File with individuals: $IN_FILE"
echo -e "#### reference_mapping.sh: Individual: $INDV \n\n"

################################################################################
#### MAP TO REFERENCE GENOME ####
################################################################################
if [[ $MODE == "PE" ]]
then
	echo -e "#### reference_mapping.sh: Reference mapping for paired-end individual $INDV ...\n"
	bwa mem -aM -R $READGROUP -t $NT $INDEX $IN_DIR/$INDV.trimmed.1.fq.gz $IN_DIR/$INDV.trimmed.2.fq.gz | samtools view -b -h > $BAM_DIR/$INDV.bam
elif [[ $MODE == "SE" ]]
then
	echo -e "#### reference_mapping.sh: Reference mapping for single-end individual $INDV ...\n"
	bwa mem -aM -R $READGROUP -t $NT $INDEX $IN_DIR/$INDV.trimmed.1.fq.gz | samtools view -b -h > $BAM_DIR/$INDV.bam
else
	echo -e "#### reference_mapping.sh: Invalid sequencing mode provided - only PE and SE allowed. ...\n" && exit 1
fi

## Report:
echo -e "\n#### reference_mapping.sh: Done with script."
date






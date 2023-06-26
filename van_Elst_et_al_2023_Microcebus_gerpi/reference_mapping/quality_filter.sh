#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
# SAMtools needs to be included in $PATH (v1.11; http://www.htslib.org/)

## Command-line args:
MODE=$1
NT=$2
BAM_DIR=$3
IN_FILE=$4
MINMAPQ=$5
INDV=$(sed -n "$SLURM_ARRAY_TASK_ID"p $IN_FILE)

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### quality_filter.sh: Starting script."
echo -e "#### quality_filter.sh: Paired-end or single-end: $MODE"
echo -e "#### quality_filter.sh: Number of threads: $NT"
echo -e "#### quality_filter.sh: Directory for BAM files: $BAM_DIR"
echo -e "#### quality_filter.sh: File with individuals: $IN_FILE"
echo -e "#### quality_filter.sh: Minimum mapping quality for filtering: $MINMAPQ"
echo -e "#### quality_filter.sh: Individual: $INDV \n\n"

################################################################################
#### SORT AND FILTER FOR MINIMUM MAPPING QUALITY AND PROPER PAIRING (IF PE) ####
################################################################################
if [[ $MODE == "PE" ]]
then
	echo -e "#### quality_filter.sh: Minimum mapping quality and proper-pair filtering and sorting for paired-end individual $INDV ...\n"
	samtools view -bhu -q $MINMAPQ -f 0x2 -@ $NT $BAM_DIR/$INDV.bam | samtools sort -@ $NT -m 15G -O bam > $BAM_DIR/$INDV.MQ$MINMAPQ.pp.bam
elif [[ $MODE == "SE" ]]
	echo -e "#### quality_filter.sh: Minimum mapping quality filtering and sorting for single-end individual $INDV ...\n"
	samtools view -bhu -q $MINMAPQ -@ $NT $BAM_DIR/$INDV.bam | samtools sort -@ $NT -m 15G -O bam > $BAM_DIR/$INDV.MQ$MINMAPQ.bam
else
	echo -e "#### quality_filter.sh: Invalid sequencing mode provided - only PE and SE allowed. ...\n" && exit 1
fi

## Report:
echo -e "\n#### quality_filter.sh: Done with script."
date
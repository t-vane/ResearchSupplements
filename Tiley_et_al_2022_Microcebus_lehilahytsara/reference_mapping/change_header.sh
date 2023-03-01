#!/bin/sh
#SBATCH -p medium40

set -euo pipefail

source /home/nibtve93/.bashrc

################################################################################
#### SET-UP ####
################################################################################
## Software:
# SAMtools needs to be included in $PATH

## Command-line args:
MODE=$1
OUT_DIR=$2
IN_FILE=$3
MINMAPQ=$4
SUFFIX=$5
EXCLUDE=$6
INDV=$(sed -n "$SLURM_ARRAY_TASK_ID"p $IN_FILE)

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### change_header.sh: Starting script."
echo -e "#### change_header.sh: Paired-end or single-end: $MODE"
echo -e "#### change_header.sh: Directory with BAM files: $OUT_DIR"
echo -e "#### change_header.sh: Minimum mapping quality: $MINMAPQ"
echo -e "#### change_header.sh: Suffix of BAM files: $SUFFIX"
echo -e "#### change_header.sh: String with chromosomes to exclude from header: $EXCLUDE"
echo -e "#### change_header.sh: Individual: $INDV \n\n"

################################################################################
#### REMOVE CHROMOSOMES THAT ARE NOT REPRESENTED IN HEADER ####
################################################################################
if [[ $MODE == "PE" ]]
then
	echo -e "#### change_header.sh: Removing chromosomes $EXCLUDE from header for paired-end individual $INDV ...\n"
	samtools view -H $OUT_DIR/$INDV.MQ$MINMAPQ.pp.dedup.$SUFFIX.bam | egrep -v $EXCLUDE > $OUT_DIR/$INDV.fixedhead
	cat $OUT_DIR/${INDV}.fixedhead <(samtools view $OUT_DIR/$INDV.MQ$MINMAPQ.pp.dedup.$SUFFIX.bam) | samtools view -bo $OUT_DIR/${INDV}.MQ$MINMAPQ.pp.dedup.$SUFFIX.bam
	rm $OUT_DIR/$INDV.fixedhead

elif [[ $MODE == "SE" ]]
then
	echo -e "#### change_header.sh: Removing chromosomes $EXCLUDE from header for single-end individual $INDV ...\n"
	samtools view -H $OUT_DIR/$INDV.MQ$MINMAPQ.$SUFFIX.bam | egrep -v $EXCLUDE > $OUT_DIR/$INDV.fixedhead
	cat $OUT_DIR/${INDV}.fixedhead <(samtools view $OUT_DIR/$INDV.MQ$MINMAPQ.$SUFFIX.bam) | samtools view -bo $OUT_DIR/${INDV}.MQ$MINMAPQ.$SUFFIX.bam
	rm $OUT_DIR/$INDV.fixedhead
else
	echo -e "#### change_header.sh: Invalid sequencing mode provided - only PE and SE allowed. ...\n"
	exit 1
fi

## Report:
echo -e "\n#### change_header: Done with script."
date







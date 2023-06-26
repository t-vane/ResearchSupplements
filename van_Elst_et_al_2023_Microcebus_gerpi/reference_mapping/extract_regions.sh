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
BAM_DIR=$2
IN_FILE=$3
MINMAPQ=$4
BED=$5
EXCLUDE=$6
SUFFIX=$7
INDV=$(sed -n "$SLURM_ARRAY_TASK_ID"p $IN_FILE)

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### extract_regions.sh: Starting script."
echo -e "#### extract_regions.sh: Paired-end or single-end: $MODE"
echo -e "#### extract_regions.sh: Directory for BAM files: $BAM_DIR"
echo -e "#### extract_regions.sh: File with individuals: $IN_FILE"
echo -e "#### extract_regions.sh: Minimum mapping quality for filtering: $MINMAPQ"
echo -e "#### extract_regions.sh: BED file with information on genomic regions: $BED"
echo -e "#### extract_regions.sh: String with chromosomes to exclude from header: $EXCLUDE"
echo -e "#### extract_regions.sh: Suffix for final BAM files: $SUFFIX"
echo -e "#### extract_regions.sh: Individual: $INDV \n\n"

################################################################################
#### EXTRACT SPECIFIC REGIONS FROM BED FILE ####
################################################################################
if [[ $MODE == "PE" ]]
then
	echo -e "#### extract_regions.sh: Extraction of genomic regions provided in $BED for paired-end individual $INDV ...\n"
	samtools view -b -L $BED $BAM_DIR/$INDV.MQ$MINMAPQ.pp.dedup.bam > $BAM_DIR/$INDV.MQ$MINMAPQ.pp.dedup.$SUFFIX.bam
	
	echo -e "#### extract_regions.sh: Removing chromosomes $EXCLUDE from header for paired-end individual $INDV ...\n"
	samtools view -H $BAM_DIR/$INDV.MQ$MINMAPQ.pp.dedup.$SUFFIX.bam | egrep -v $EXCLUDE > $BAM_DIR/$INDV.fixedhead
	cat $BAM_DIR/$INDV.fixedhead <(samtools view $BAM_DIR/$INDV.MQ$MINMAPQ.pp.dedup.$SUFFIX.bam) | samtools view -bo $BAM_DIR/$INDV.MQ$MINMAPQ.pp.dedup.$SUFFIX.bam
	rm $BAM_DIR/$INDV.fixedhead
elif [[ $MODE == "SE" ]]
	echo -e "#### extract_regions.sh: Extraction of genomic regions provided in $BED for single-end individual $INDV ...\n"
	samtools view -b -L $BED $BAM_DIR/$INDV.MQ$MINMAPQ.bam > $BAM_DIR/$INDV.MQ$MINMAPQ.$SUFFIX.bam
	
	echo -e "#### extract_regions.sh: Removing chromosomes $EXCLUDE from header for single-end individual $INDV ...\n"
	samtools view -H $BAM_DIR/$INDV.MQ$MINMAPQ.$SUFFIX.bam | egrep -v $EXCLUDE > $BAM_DIR/$INDV.fixedhead
	cat $BAM_DIR/$INDV.fixedhead <(samtools view $BAM_DIR/$INDV.MQ$MINMAPQ.$SUFFIX.bam) | samtools view -bo $BAM_DIR/$INDV.MQ$MINMAPQ.$SUFFIX.bam
	rm $BAM_DIR/$INDV.fixedhead
else
	echo -e "#### extract_regions.sh: Invalid sequencing mode provided - only PE and SE allowed. ...\n" && exit 1
fi

## Report:
echo -e "\n#### extract_regions.sh: Done with script."
date

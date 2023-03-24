#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
# SAMtools needs to be included in $PATH (v1.11; http://www.htslib.org/)

## Command-line args:
BAM_DIR=$1
IN_FILE=$2
MINMAPQ=$3
INDV=$(sed -n "$SLURM_ARRAY_TASK_ID"p $IN_FILE)

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### deduplicate.sh: Starting script."
echo -e "#### deduplicate.sh: Directory for BAM files: $BAM_DIR"
echo -e "#### deduplicate.sh: File with individuals: $IN_FILE"
echo -e "#### deduplicate.sh: Minimum mapping quality for filtering: $MINMAPQ"
echo -e "#### deduplicate.sh: Individual: $INDV \n\n"

################################################################################
#### DEDUPLICATE ####
################################################################################
echo -e "#### deduplicate.sh: Deduplication for paired-end individual $INDV ...\n"
samtools collate --output-fmt BAM $BAM_DIR/$INDV.MQ$MINMAPQ.pp.bam -O | samtools fixmate - - -r -m -O BAM| samtools sort -m 15G -O BAM | samtools markdup - $BAM_DIR/$INDV.MQ$MINMAPQ.pp.dedup.bam -r -s -O BAM

## Report:
echo -e "\n#### deduplicate.sh: Done with script."
date

#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
# samtools needs to be included in $PATH (v1.11; http://www.htslib.org/)

## Command-line args:
MODE=$1
IN_FILE=$2
BAM_DIR=$3
OUT_DIR=$4
INDV=$(sed -n "$SLURM_ARRAY_TASK_ID"p $IN_FILE)

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### coverage.sh: Starting script."
echo -e "#### coverage.sh: Sequencing mode: $MODE"
echo -e "#### coverage.sh: File with individuals: $IN_FILE"
echo -e "#### coverage.sh: Directory with BAM files: $BAM_DIR"
echo -e "#### coverage.sh: Output directory: $OUT_DIR"
echo -e "#### coverage.sh: Individual: $INDV \n\n"

################################################################################
#### ESTIMATE NUMBER OF SITES AND COVERAGE ####
################################################################################
if [[ $MODE == "PE" ]]
then
	echo -e "#### coverage.sh: Estimating coverage for each SbfI forward strand site of $INDV ...\n"
	samtools view -f 0x40 $BAM_DIR/$INDV.auto.bam | awk '$10 ~ /^TGCAGG/' | cut -f 3,4 | tr '\t' '_' | sort | uniq -c | sort -k1 -nr | sed -E 's/^ *//; s/ /\t/' > $OUT_DIR/$INDV.sbf1.f1.bamhits
elif [[ $MODE == "SE" ]]
then
	echo -e "#### coverage.sh: Estimating coverage for each SbfI forward strand site of $INDV ...\n"
	samtools view $BAM_DIR/$INDV.auto.bam | awk '$10 ~ /^TGCAGG/' | cut -f 3,4 | tr '\t' '_' | sort | uniq -c | sort -k1 -nr | sed -E 's/^ *//; s/ /\t/' > $OUT_DIR/$INDV.sbf1.f1.bamhits
else
	echo -e "#### coverage.sh: Invalid sequencing mode provided - only PE and SE allowed. ...\n" && exit 1
fi

## Report:
NSITES=$(wc -l $OUT_DIR/$INDV.sbf1.f1.bamhits)
NREADS=$(cat $OUT_DIR/$INDV.sbf1.f1.bamhits | cut -f 1 | paste -sd+ | bc)
MEAN_COV=$(( $NSITES / $NREADS ))

echo -e "#### coverage.sh: Number of SbfI sites recovered for $INDV: $NSITES"
echo -e "#### coverage.sh: Total number of reads across these sites: $NREADS"
echo -e "#### coverage.sh: Mean coverage across these sites: $MEAN_COV"

echo -e "\n#### coverage.sh: Done with script."
date




#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################

## Command-line args:
SCRIPTS_DIR=$1
VCF_IN=$2
SAMPLE_FILE=$3
OUT_FILE=$4

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### hierfstat.sh: Starting script."
echo -e "#### hierfstat.sh: Directory with scripts: $SCRIPTS_DIR"
echo -e "#### hierfstat.sh: Input VCF file: $VCF_IN"
echo -e "#### hierfstat.sh: File with individuals and associated populations: $SAMPLE_FILE"
echo -e "#### hierfstat.sh: Output file: $OUT_FILE \n\n"

################################################################################
#### ESTIMATE PAIRWISE F_ST BETWEEN POPULATIONS ####
################################################################################
echo -e "#### hierfstat.sh: Estimating pairwise F_ST between populations ...\n"
Rscript $SCRIPTS_DIR/hierfstat.R $VCF_IN $SAMPLE_FILE $OUT_FILE

## Report:
echo -e "\n#### hierfstat.sh: Done with script."
date


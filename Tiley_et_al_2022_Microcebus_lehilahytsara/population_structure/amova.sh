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
OUT=$4

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### amova.sh: Starting script."
echo -e "#### amova.sh: Directory with scripts: $SCRIPTS_DIR"
echo -e "#### amova.sh: Input VCF file: $VCF_IN"
echo -e "#### amova.sh: File with individuals and associated populations: $SAMPLE_FILE"
echo -e "#### amova.sh: Output file: $OUT \n\n"

################################################################################
#### CONDUCT AMOVA ####
################################################################################
echo -e "#### amova.sh: Conducting AMOVA ...\n"
Rscript $SCRIPTS_DIR/amova.R $VCF_IN $SAMPLE_FILE $OUT

## Report:
echo -e "\n#### amova.sh: Done with script."
date


#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################

## Command-line args
SCRIPT_DIR=$1
VCF_FILE=$2
OUT_FILE=$3

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### vcfr.sh: Starting script."
echo -e "#### vcfr.sh: Directory with scripts: $SCRIPT_DIR"
echo -e "#### vcfr.sh: Input VCF file: $VCF_FILE"
echo -e "#### vcfr.sh: Output file: $OUT_FILE \n\n"

################################################################################
#### ESTIMATE GENETIC DISTANCES BETWEEN INDIVIDUALS ####
################################################################################
echo -e "#### vcfr.sh: Estimating genetic distances between individuals ...\n"
Rscript $SCRIPTS_DIR/vcfr.R $VCF_FILE $OUT_FILE

## Report:
echo -e "\n#### vcfr.sh: Done with script."
date


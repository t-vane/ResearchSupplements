#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
# VCFtools needs to be included in $PATH (v0.1.17; https://vcftools.github.io/index.html)

## Command-line args:
VCF_IN=$1
VCF_OUT=$2
REM_STRING=$3

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### 07_filter_outgroups.sh: Starting script."
echo -e "#### 07_filter_outgroups.sh: Input VCF: $VCF_IN"
echo -e "#### 07_filter_outgroups.sh: Output VCF: $VCF_OUT"
echo -e "#### 07_filter_outgroups.sh: String to remove outgroups: $REM_STRING \n\n"

################################################################################
#### FILTER OUTGROUPS ####
################################################################################
echo -e "#### 07_filter_outgroups.sh: Filtering outgroups ...\n"
vcftools $REM_STRING --vcf $VCF_FILE --recode --recode-INFO-all --stdout > $VCF_OUT

## Report:
echo -e "\n#### 07_filter_outgroups.sh: Listing output VCF:"
ls -lh $VCF_OUT
[[ $(grep -cv "^#" $VCF_OUT) = 0 ]] && echo -e "\n\n#### 07_filter_outgroups.sh: ERROR: VCF is empty\n" >&2 && exit 1

echo -e "\n#### 07_filter_outgroups.sh: Done with script."
date


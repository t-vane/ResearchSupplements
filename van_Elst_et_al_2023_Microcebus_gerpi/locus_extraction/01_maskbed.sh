#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Script adapted and modified from Poelstra et al. 2021, Systematic Biology (https://doi.org/10.1093/sysbio/syaa053)

## Software:
# BEDtools needs to be included in $PATH (v2.30.0; https://bedtools.readthedocs.io/en/latest/)
# BEDOPS needs to be included in $PATH (v2.4.38; https://bedops.readthedocs.io/en/latest/)

## Command-line args:
VCF_ALTREF=$1
VCF_FILT_MASK=$2
BED_REMOVED_SITES=$3
BED_DIR=$4

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### 01_maskbed.sh: Starting script."
echo -e "#### 01_maskbed.sh: Raw VCF: $VCF_ALTREF"
echo -e "#### 01_maskbed.sh: Fully filtered VCF: $VCF_FILT_MASK"
echo -e "#### 01_maskbed.sh: BED file for removed sites: $BED_REMOVED_SITES"
echo -e "#### 01_maskbed.sh: Directory for BED file: $BED_DIR \n\n"

################################################################################
#### CREATE BED FILE WITH MASKED SITES ####
################################################################################
mkdir -p $BED_DIR/tmpdir

echo -e "#### 01_maskbed.sh: Creating $BED_REMOVED_SITES with masked sites ..."
bedtools intersect -v -a <(vcf2bed --sort-tmpdir=$BED_DIR/tmpdir < $VCF_ALTREF | cut -f 1,2,3) -b <(vcf2bed --sort-tmpdir=$BED_DIR/tmpdir < $VCF_FILT_MASK | cut -f 1,2,3) > $BED_REMOVED_SITES

## Report:
echo -e "\n#### 01_maskbed.sh: Line count of removed sites: $(wc -l < $BED_REMOVED_SITES)."
echo -e "#### 01_maskbed.sh: Done with script."
date


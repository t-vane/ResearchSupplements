#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Script adapted and modified from Poelstra et al. 2021, Systematic Biology (https://doi.org/10.1093/sysbio/syaa053)

## Software:
# VCFtools needs to be included in $PATH (v0.1.17; https://vcftools.github.io/index.html)

## Command-line args:
VCF_IN=$1
MIN_DP=$2
MEAN_DP=$3
VCF_OUT=$4

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### 01_filter_min-dp.sh: Starting script."
echo -e "#### 01_filter_min-dp.sh: Input VCF: $VCF_IN"
echo -e "#### 01_filter_min-dp.sh: Minimum depth: $MIN_DP"
echo -e "#### 01_filter_min-dp.sh: Minimum mean depth: $MEAN_DP"
echo -e "#### 01_filter_min-dp.sh: Output VCF: $VCF_OUT \n\n"

################################################################################
#### FILTER FOR MINIMUM DEPTHS ####
################################################################################
## Filter with VCFtools for minimum depth
echo -e "#### 01_filter_min-dp.sh: Filtering for minimum depth ...\n"
vcftools --vcf $VCF_IN --minDP $MIN_DP --min-meanDP $MEAN_DP --recode --recode-INFO-all --stdout > $VCF_OUT

## Report:
NVAR_IN=$(grep -cv "^#" $VCF_IN)
NVAR_OUT=$(grep -cv "^#" $VCF_OUT)
NVAR_FILT=$(( $NVAR_IN - $NVAR_OUT ))

echo -e "\n\n"
echo -e "#### 01_filter_min-dp.sh: Number of SNPs in input VCF: $NVAR_IN"
echo -e "#### 01_filter_min-dp.sh: Number of SNPs in output VCF: $NVAR_OUT"
echo -e "#### 01_filter_min-dp.sh: Number of SNPs filtered: $NVAR_FILT"
echo
echo -e "#### 01_filter_min-dp.sh: Listing output VCF:"
ls -lh $VCF_OUT
[[ $(grep -cv "^#" $VCFOUT) = 0 ]] && echo -e "\n\n#### 01_filter_min-dp.sh: ERROR: VCF is empty\n" >&2 && exit 1

echo -e "\n#### 01_filter_min-dp.sh: Done with script."
date

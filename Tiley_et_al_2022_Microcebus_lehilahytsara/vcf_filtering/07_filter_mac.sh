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
MAC=$3

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### 07_filter_mac.sh: Starting script."
echo -e "#### 07_filter_mac.sh: Input VCF: $VCF_IN"
echo -e "#### 07_filter_mac.sh: Output VCF: $VCF_OUT"
echo -e "#### 07_filter_mac.sh: Minor allele count: $MAC \n\n"

################################################################################
#### FILTER WITH MINOR ALLELE COUNT ####
################################################################################
echo -e "#### 07_filter_mac.sh: Filtering for minor allele count $MAC ...\n"
vcftools --vcf $VCF_IN --mac $MAC --recode --recode-INFO-all --stdout > $VCF_OUT

## Report:
NVAR_IN=$(grep -cv "^#" $VCF_IN || true)
NVAR_OUT=$(grep -cv "^#" $VCF_OUT || true)
NVAR_FILT=$(( $NVAR_IN - $NVAR_OUT ))

echo -e "\n#### 07_filter_mac.sh: Number of SNPs before MAC filtering: $NVAR_IN"
echo -e "#### 07_filter_mac.sh: Number of SNPs filtered: $NVAR_FILT"
echo -e "#### 07_filter_mac.sh: Number of SNPs after MAC filtering: $NVAR_OUT"

echo -e "\n#### 07_filter_mac.sh: Listing output VCF:"
ls -lh $VCF_OUT
[[ $(grep -cv "^#" $VCF_OUT) = 0 ]] && echo -e "\n\n#### 07_filter_mac.sh: ERROR: VCF is empty\n" >&2 && exit 1

echo -e "\n#### 07_filter_mac.sh: Done with script."
date


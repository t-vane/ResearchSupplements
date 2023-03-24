#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
# vcftools needs to be included in $PATH (v0.1.17; https://vcftools.github.io/index.html)
# bcftools needs to be included in $PATH (v1.11; http://www.htslib.org/)

## Command-line args:
VCF_IN=$1
DROP_INDS=$2
VCF_OUT=$3

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### 08_drop_inds.sh: Starting script."
echo -e "#### 08_drop_inds.sh: Input VCF: $VCF_IN"
echo -e "#### 08_drop_inds.sh: List of individuals to remove: $DROP_INDS"
echo -e "#### 08_drop_inds.sh: Output VCF: $VCF_OUT \n\n"

################################################################################
#### REMOVE UNDESIRED INDIVIDUALS ####
################################################################################
echo -e "#### 08_drop_inds.sh: Removing undesired individuals ...\n"
vcftools --remove $DROP_INDS --vcf $VCF_IN --recode --recode-INFO-all --stdout > $VCF_OUT

## Report:
NIND_IN=$(bcftools query -l | wc -l $VCF_IN)
NIND_OUT=$(bcftools query -l | wc -l $VCF_OUT)
NIND_FILT=$(( $NIND_IN - $NIND_OUT ))

echo -e "\n#### 08_drop_inds.sh: Number of individuals before filtering: $NIND_IN"
echo -e "#### 08_drop_inds.sh: Number of individuals filtered: $NIND_FILT"
echo -e "#### 08_drop_inds.sh: Number of individuals after filtering: $NIND_OUT"

echo -e "\n#### 08_drop_inds.sh: Listing output VCF:"
ls -lh $VCF_OUT
[[ $(grep -cv "^#" $VCF_OUT) = 0 ]] && echo -e "\n\n#### 08_drop_inds.sh: ERROR: VCF is empty\n" >&2 && exit 1

echo -e "\n#### 08_drop_inds.sh: Done with script."
date


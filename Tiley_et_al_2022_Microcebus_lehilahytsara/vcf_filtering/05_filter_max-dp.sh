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
VCF_OUT=$2

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### 05_filter_max-dp.sh: Starting script."
echo -e "#### 05_filter_max-dp.sh: Input VCF: $VCF_IN"
echo -e "#### 05_filter_max-dp.sh: Output VCF: $VCF_OUT \n\n"

################################################################################
#### FILTER FOR MAXIMUM DEPTH ####
################################################################################
## Create a file with original site depths per locus
mkdir -p $(dirname $VCF_OUT)/tmp
DP_FILE=$(dirname $VCF_OUT)/tmp/vcf_depth.$RANDOM
echo -e "#### 05_filter_max-dp.sh: Extracting original site depths per locus ...\n"
vcftools --vcf $VCF_IN --site-depth --stdout | tail -n +2 | cut -f 3 > $DP_FILE

## Calculate the mean maximum depth: (mean depth + 2 * standard deviation) / number of individuals
echo -e "#### 05_filter_max-dp.sh: Estimating mean maximum depth ...\n"
N_IND=$(bcftools query -l $VCF_IN | wc -l || true)
DP_MEAN_ALL=$(awk '{ total += $1; count++ } END { print total/count }' $DP_FILE)
DP_SD_ALL=$(awk '{ sum+=$1; sumsq+=$1*$1 } END { print sqrt(sumsq/NR - (sum/NR)^2) }' $DP_FILE)
DP_HI_ALL=$(python -c "print($DP_MEAN_ALL + (2* $DP_SD_ALL))")
DP_MEAN_IND=$(python -c "print($DP_MEAN_ALL / $N_IND)")
DP_MAX_IND=$(python -c "print($DP_HI_ALL / $N_IND)")

## Filter with VCFtools for maximum depth
echo -e "#### 05_filter_max-dp.sh: Filtering for maximum depth ...\n"
vcftools --vcf $VCF_IN --max-meanDP $DP_MAX_IND --recode --recode-INFO-all --stdout > $VCF_OUT

## Save a separate file with SNPs with too high depth
echo -e "#### 05_filter_max-dp.sh: Creating file with SNPs with too high depth ...\n"
vcftools --vcf $VCF_IN --min-meanDP $DP_MAX_IND --recode --recode-INFO-all --stdout > ${VCF_OUT//.vcf/_too-high-DP.vcf}

## Remove temporary file
rm $DP_FILE

## Report:
NVAR_IN=$(grep -cv "^#" $VCF_IN)
NVAR_OUT=$(grep -cv "^#" $VCF_OUT)
NVAR_FILT=$(( $NVAR_IN - $NVAR_OUT ))

echo -e "\n\n"
echo -e "#### 05_filter_max-dp.sh: Number of SNPs in input VCF: $NVAR_IN"
echo -e "#### 05_filter_max-dp.sh: Number of SNPs in output VCF: $NVAR_OUT"
echo -e "#### 05_filter_max-dp.sh: Number of SNPs filtered: $NVAR_FILT"
echo
echo -e "#### 05_filter_max-dp.sh: Listing output VCF:"
ls -lh $VCF_OUT
[[ $(grep -cv "^#" $VCFOUT) = 0 ]] && echo -e "\n\n#### 05_filter_max-dp.sh: ERROR: VCF is empty\n" >&2 && exit 1

echo -e "\n#### 05_filter_max-dp.sh: Done with script."
date


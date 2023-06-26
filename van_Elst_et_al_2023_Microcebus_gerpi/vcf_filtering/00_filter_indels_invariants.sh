#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
# VCFtools needs to be included in $PATH (v0.1.17; https://vcftools.github.io/index.html)
GATK3=/home/nibtve93/software/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar # (v3.8.1; https://gatk.broadinstitute.org/hc/en-us)

## Command-line args:
REFERENCE=$1
VCF_IN=$2
VCF_OUT=$3

## Activate conda environment
conda activate java

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### 00_filter_indels_invariants.sh: Starting script."
echo -e "#### 00_filter_indels_invariants.sh: Reference genome: $REFERENCE"
echo -e "#### 00_filter_indels_invariants.sh: Input VCF: $VCF_IN"
echo -e "#### 00_filter_indels_invariants.sh: Output VCF: $VCF_OUT \n\n"

################################################################################
#### REMOVE INDELS AND INVARIANT SITES ####
################################################################################
echo -e "#### 00_filter_indels_invariants.sh: Removing indels ...\n"
java -jar $GATK3 -T SelectVariants -R $REFERENCE -V $VCF_IN -o $(dirname $VCF_IN)/$(basename $VCF_IN .vcf)/.tmp.vcf -selectType SNP

echo -e "#### 00_filter_indels_invariants.sh: Removing invariant sites ...\n"
vcftools --vcf $(dirname $VCF_IN)/$(basename $VCF_IN .vcf)/.tmp.vcf --recode --recode-INFO-all --max-non-ref-af 0.99 --min-alleles 2 --stdout > $VCF_OUT

## Report:
NVAR_IN=$(grep -cv "^#" $VCF_IN)
NVAR_OUT=$(grep -cv "^#" $VCF_OUT)
NVAR_FILT=$(( $NVAR_IN - $NVAR_OUT ))

echo -e "\n\n"
echo -e "#### 00_filter_indels_invariants.sh: Number of sites in input VCF: $NVAR_IN"
echo -e "#### 00_filter_indels_invariants.sh: Number of sites in output VCF: $NVAR_OUT"
echo -e "#### 00_filter_indels_invariants.sh: Number of sites filtered: $NVAR_FILT"
echo
echo -e "#### 00_filter_indels_invariants.sh: Listing output VCF:"
ls -lh $VCF_OUT
[[ $(grep -cv "^#" $VCFOUT) = 0 ]] && echo -e "\n\n#### 00_filter_indels_invariants.sh: ERROR: VCF is empty\n" >&2 && exit 1

echo -e "\n#### 00_filter_indels_invariants.sh: Done with script."
date



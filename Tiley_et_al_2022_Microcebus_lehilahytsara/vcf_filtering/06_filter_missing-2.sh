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
MAXMISS_GENO=$3
FILTER_INDS=$4
MAXMISS_IND=$5

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### 06_filter_missing-2.sh: Starting script."
echo -e "#### 06_filter_missing-2.sh: Input VCF: $VCF_IN"
echo -e "#### 06_filter_missing-2.sh: Output VCF: $VCF_OUT"
echo -e "#### 06_filter_missing-2.sh: Maximum missingness value for genotypes: $MAXMISS_GENO"
echo -e "#### 06_filter_missing-2.sh: NOTE: Maximum missingness is inverted for genotypes, with a value of 1 indicating no missing data (this is how the '--max-missing' option in VCFtools works)."
echo -e "#### 06_filter_missing-2.sh: Filter individuals: $FILTER_INDS"
[[ $FILTER_IND ]] && echo -e "#### 06_filter_missing-2.sh: Maximum missingness value for individuals: $MAXMISS_IND \n\n"

################################################################################
#### FILTER FOR MISSINGNESS PER GENOTYPE (AND INDIVIDUAL) ####
################################################################################
## Estimate number of indidviduals in input
echo -e "#### 06_filter_missing-2.sh: Estimating number of individuals in $VCF_IN ...\n"
NIND_IN=$(bcftools query -l $VCF_IN | wc -l || true)

## Create temporary filter files
mkdir -p $(dirname $VCF_OUT)/tmp
FILTERFILE_PREFIX=$(dirname $VCF_OUT)/tmp/filterfile.$RANDOM
VCF_TMP=$(dirname $VCF_OUT)/tmp/vcf_indfilter.$RANDOM.vcf

## Run filtering
if [[ $FILTER_INDS ]]
then
	echo -e "#### 06_filter_missing-2.sh: Filtering individuals by missing data ... "
	# Get amount of missing data per individual
	vcftools --vcf $VCF_IN --missing-indv --stdout > $FILTERFILE_PREFIX.imiss
	# Get list of individuals with too much missing data
	tail -n +2 $FILTERFILE_PREFIX.imiss | awk -v var=$MAXMISS_IND '$5 > var' | cut -f1 > $FILTERFILE_PREFIX.HiMissInds
	# Remove individuals with too much missing data
	vcftools --vcf $VCF_IN --remove $FILTERFILE_PREFIX.HiMissInds --recode --recode-INFO-all --stdout > $VCF_TMP	
else
	echo -e "#### 06_filter_missing-2.sh: Only filtering by missing data at genotype level (no individuals will be removed)"
	VCF_TMP=$VCF_IN
fi

echo -e "#### 06_filter_missing-2.sh: Filtering genotypes by missing data... "
vcftools --vcf $VCF_TMP --max-non-ref-af 0.99 --min-alleles 2 --max-missing $MAXMISS_GENO --recode --recode-INFO-all --stdout > $VCF_OUT

## Remove temporary files
rm $FILTERFILE_PREFIX*
rm $VCF_TMP

## Report:
NVAR_IN=$(grep -cv "^#" $VCF_IN || true)
NVAR_OUT=$(grep -cv "^#" $VCF_OUT || true)
NVAR_FILT=$(( $NVAR_IN - $NVAR_OUT ))

NIND_FILT=$(wc -l < $FILTERFILE_PREFIX.HiMissInds || true)
NIND_OUT=$(bcftools query -l $VCF_OUT | wc -l || true)

if [[ $FILTER_INDS ]]
then
	echo -e "\n#### 06_filter_missing-2.sh: Number of indidviduals before individual filtering: $NIND_IN"
	echo -e "#### 06_filter_missing-2.sh: Number of indidviduals filtered: $NIND_FILT"
	echo -e "#### 06_filter_missing-2.sh: Number of indidviduals after individual filtering: $NIND_OUT"
fi

echo -e "\n#### 06_filter_missing-2.sh: Number of SNPs before genotype filtering: $NVAR_IN"
echo -e "#### 06_filter_missing-2.sh: Number of SNPs filtered: $NVAR_FILT"
echo -e "#### 06_filter_missing-2.sh: Number of SNPs after genotype filtering: $NVAR_OUT"

echo -e "\n#### 06_filter_missing-2.sh: Listing output VCF:"
ls -lh $VCF_OUT
[[ $(grep -cv "^#" $VCF_OUT) = 0 ]] && echo -e "\n\n#### 06_filter_missing-2.sh: ERROR: VCF is empty\n" >&2 && exit 1

echo -e "\n#### 06_filter_missing-2.sh: Done with script."
date


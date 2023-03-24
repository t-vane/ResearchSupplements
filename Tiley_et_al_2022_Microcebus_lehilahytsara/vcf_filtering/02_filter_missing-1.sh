#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Script adapted and modified from Poelstra et al. 2021, Systematic Biology (https://doi.org/10.1093/sysbio/syaa053)

## Software:
# VCFtools needs to be included in $PATH (v0.1.17; https://vcftools.github.io/index.html)
# BCFtools needs to be included in $PATH (v1.11; http://www.htslib.org/)

## Command-line args:
VCF_IN=$1
VCF_OUT=$2
MAXMISS_GENO1=$3
MAXMISS_GENO2=$4
MAXMISS_GENO3=$5
FILTER_INDS=$6
MAXMISS_IND1=$7
MAXMISS_IND2=$8
MAXMISS_IND3=$9

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### 02_filter_missing-1.sh: Starting script."
echo -e "#### 02_filter_missing-1.sh: Input VCF: $VCF_IN"
echo -e "#### 02_filter_missing-1.sh: Output VCF: $VCF_OUT"
echo -e "#### 02_filter_missing-1.sh: Maximum missingness value for genotypes (round 1): $MAXMISS_GENO1"
echo -e "#### 02_filter_missing-1.sh: Maximum missingness value for genotypes (round 2): $MAXMISS_GENO2"
echo -e "#### 02_filter_missing-1.sh: Maximum missingness value for genotypes (round 3): $MAXMISS_GENO3"
echo -e "#### 02_filter_missing-1.sh: NOTE: Maximum missingness is inverted for genotypes, with a value of 1 indicating no missing data (this is how the '--max-missing' option in VCFtools works)."
echo -e "#### 02_filter_missing-1.sh: Filter individuals: $FILTER_INDS"
[[ $FILTER_IND ]] && echo -e "#### 02_filter_missing-1.sh: Maximum missingness value for individuals (round 1): $MAXMISS_IND1"
[[ $FILTER_IND ]] && echo -e "#### 02_filter_missing-1.sh: Maximum missingness value for individuals (round 2): $MAXMISS_IND2"
[[ $FILTER_IND ]] && echo -e "#### 02_filter_missing-1.sh: Maximum missingness value for individuals (round 3): $MAXMISS_IND3 \n\n"

################################################################################
#### FILTER FOR MISSINGNESS PER GENOTYPE (AND INDIVIDUAL) ####
################################################################################
## Estimate number of indidviduals in input
echo -e "#### 02_filter_missing-1.sh: Estimating number of individuals in $VCF_IN ...\n"
NIND_IN=$(bcftools query -l $VCF_IN | wc -l || true)

## Create temporary filter files
mkdir -p $(dirname $VCF_OUT)/tmp
FILTERFILE_PREFIX=$(dirname $VCF_OUT)/tmp/filterfile.$RANDOM
VCF_TMP_PREFIX=$(dirname $VCF_OUT)/tmp/vcf.$RANDOM

## Run filtering
if [[ $FILTER_INDS ]]
then
	echo -e "#### 02_filter_missing-1.sh: Filtering by missing data at genotype and individual level in three rounds ...\n"
	for i in 1 2 3 
	do
		[[ $i == 1 ]] && MAXMISS_GENO=$MAXMISS_GENO1 && MAXMISS_IND=$MAXMISS_IND1
		[[ $i == 2 ]] && MAXMISS_GENO=$MAXMISS_GENO2 && MAXMISS_IND=$MAXMISS_IND2
		[[ $i == 3 ]] && MAXMISS_GENO=$MAXMISS_GENO3 && MAXMISS_IND=$MAXMISS_IND3
		
		echo -e "## 02_filter_missing-1.sh: Filtering genotypes by missing data - round $i ..."
		vcftools --vcf $VCF_IN --max-non-ref-af 0.99 --min-alleles 2 --max-missing $MAXMISS_GENO --recode --recode-INFO-all --stdout > $VCF_TMP_PREFIX.R${i}a.vcf
		
		echo -e "## 02_filter_missing-1.sh: Filtering individuals by missing data - round $i ..."
		# Get amount of missing data per individual
		vcftools --vcf $VCF_TMP_PREFIX.R${i}a.vcf --missing-indv --stdout > $FILTERFILE_PREFIX.round$i.imiss
		# Get list of individuals with too much missing data
		tail -n +2 $FILTERFILE_PREFIX.round$i.imiss | awk -v var=$MAXMISS_IND '$5 > var' | cut -f1 > $FILTERFILE_PREFIX.HiMissInds$i
		# Remove individuals with too much missing data
		vcftools --vcf $VCF_TMP_PREFIX.R${i}a.vcf --remove $FILTERFILE_PREFIX.HiMissInds$i --recode --recode-INFO-all --stdout > $VCF_TMP_PREFIX.R${i}b.vcf
		[[ $i == 3 ]] && mv $VCF_TMP_PREFIX.R${i}b.vcf $VCF_OUT		
	done
	
	##Report:
	NVAR_IN=$(grep -cv "^#" $VCF_IN || true)
	NVAR_OUT=$(grep -cv "^#" $VCF_OUT || true)
	
	NVAR_R1=$(grep -cv "^#" $VCF_TMP_PREFIX.R1a.vcf || true)
	NVAR_R2=$(grep -cv "^#" $VCF_TMP_PREFIX.R2a.vcf || true)
	NVAR_R3=$(grep -cv "^#" $VCF_TMP_PREFIX.R3a.vcf || true)
	
	NVAR_FILT_R1=$(( $NVAR_IN - $NVAR_R1 ))
	NVAR_FILT_R2=$(( $NVAR_IN - $NVAR_R2 ))
	NVAR_FILT_R3=$(( $NVAR_IN - $NVAR_R3 ))
	
	NIND_FILT_R1=$(wc -l < $FILTERFILE_PREFIX.HiMissInds1 || true)
	NIND_FILT_R2=$(wc -l < $FILTERFILE_PREFIX.HiMissInds2 || true)
	NIND_FILT_R3=$(wc -l < $FILTERFILE_PREFIX.HiMissInds3 || true)
	
	NIND_OUT=$(bcftools query -l $VCF_OUT | wc -l || true)
	
	echo -e "\n#### 02_filter_missing-1.sh: Number of indidviduals before individual filtering: $NIND_IN"
	echo -e "#### 02_filter_missing-1.sh: Number of indidviduals filtered in round 1: $NIND_FILT_R1"
	echo -e "#### 02_filter_missing-1.sh: Number of indidviduals filtered in round 2: $NIND_FILT_R2"
	echo -e "#### 02_filter_missing-1.sh: Number of indidviduals filtered in round 3: $NIND_FILT_R3"
	echo -e "#### 02_filter_missing-1.sh: Number of indidviduals left after individual filtering: $NIND_OUT \n"
	
	echo -e "#### 02_filter_missing-1.sh: Number of SNPs before genotype filtering: $NVAR_IN"
	echo -e "#### 02_filter_missing-1.sh: Number of SNPs filtered in round 1: $NVAR_FILT_R1"
	echo -e "#### 02_filter_missing-1.sh: Number of SNPs filtered in round 2: $NVAR_FILT_R2"
	echo -e "#### 02_filter_missing-1.sh: Number of SNPs filtered in round 3: $NVAR_FILT_R3"
	echo -e "#### 02_filter_missing-1.sh: Number of SNPs left after genotype filtering: $NVAR_OUT"	
else
	echo -e "#### 02_filter_missing-1.sh: Only filtering by missing data at genotype level in a single round (no individuals will be removed) ..."
	vcftools --vcf $VCF_IN --max-non-ref-af 0.99 --min-alleles 2 --max-missing $MAXMISS_GENO3 --recode --recode-INFO-all --stdout > $VCF_OUT
	
	## Report:
	NVAR_IN=$(grep -cv "^#" $VCF_IN || true)
	NVAR_OUT=$(grep -cv "^#" $VCF_OUT || true)
	NVAR_FILT=$(( $NVAR_IN - $NVAR_OUT ))
	
	echo -e "\n#### 02_filter_missing-1.sh: Number of SNPs before genotype filtering: $NVAR_IN"
	echo -e "#### 02_filter_missing-1.sh: Number of SNPs after genotype filtering: $NVAR_OUT"
	echo -e "#### 02_filter_missing-1.sh: Number of SNPs filtered: $NVAR_FILT"
fi

## Remove temporary files
rm $FILTERFILE_PREFIX*
rm $VCF_TMP_PREFIX*

## Report:
echo -e "\n#### 02_filter_missing-1.sh: Listing output VCF:"
ls -lh $VCF_OUT
[[ $(grep -cv "^#" $VCF_OUT) = 0 ]] && echo -e "\n\n#### 02_filter_missing-1.sh: ERROR: VCF is empty\n" >&2 && exit 1

echo -e "\n#### 02_filter_missing-1.sh: Done with script."
date

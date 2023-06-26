#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Script adapted and modified from Poelstra et al. 2021, Systematic Biology (https://doi.org/10.1093/sysbio/syaa053)

## Software:
# VCFtools needs to be included in $PATH (v0.1.17; https://vcftools.github.io/index.html)
# gatk needs to be included in $PATH (v4.1.9.0; https://gatk.broadinstitute.org/hc/en-us)

## Command-line args:
VCF_IN=$1
VCF_OUT_SOFT=$2
VCF_OUT_HARD=$3
REFERENCE=$4

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### 04_filter_gatk.sh: Starting script."
echo -e "#### 04_filter_gatk.sh: Input VCF: $VCF_IN"
echo -e "#### 04_filter_gatk.sh: Output VCF for soft-filtering: $VCF_OUT_SOFT"
echo -e "#### 04_filter_gatk.sh: Output VCF for hard-filtering: $VCF_OUT_HARD"
echo -e "#### 04_filter_gatk.sh: Reference genome: $REFERENCE \n\n"

################################################################################
#### SOFT-FILTER VCF FILE ####
################################################################################
echo -e "#### 04_filter_gatk.sh: Soft-filtering $VCF_IN for\n"
echo -e "#### 04_filter_gatk.sh: FisherStrand, RMSMappingQuality, MappingQualityRankSumTest, ReadPosRankSumTest and AlleleBalance ...\n"
gatk VariantFiltration -R $REFERENCE -V $VCF_IN -O $VCF_OUT_SOFT \
	--filter-expression "FS > 60.0" --filter-name "FS_gt60" \
	--filter-expression "MQ < 40.0" --filter-name "MQ_lt40" \
	--filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum_ltm12.5" \
	--filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum_ltm8" \
	--filter-expression "ABHet > 0.01 && ABHet < 0.2 || ABHet > 0.8 && ABHet < 0.99" --filter-name "ABhet_filt"
	
################################################################################
#### HARD-FILTER VCF FILE ####
################################################################################
echo -e "#### 04_filter_gatk.sh: Retaining only bi-allelic sites and hard-filtering $VCF_IN for\n"
echo -e "#### 04_filter_gatk.sh: FisherStrand, RMSMappingQuality, MappingQualityRankSumTest, ReadPosRankSumTest and AlleleBalance ...\n"
vcftools --vcf $VCF_OUT_SOFT --remove-filtered-all --max-non-ref-af 0.99 --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --stdout > $VCF_OUT_HARD

## Report:
NVAR_IN=$(grep -cv "^#" $VCF_IN)
NVAR_OUT=$(grep -cv "^#" $VCF_OUT_HARD)
NVAR_FILT=$(( $NVAR_IN - $NVAR_OUT ))

NFILT_FS=$(grep -c "FS_gt60" $VCF_OUT_SOFT || true)
NFILT_MQ=$(grep -c "MQ_lt40" $VCF_OUT_SOFT || true)
NFILT_MQR=$(grep -c "MQRankSum_ltm12" $VCF_OUT_SOFT || true)
NFILT_READPOS=$(grep -c "readposRankSum_ltm8" $VCF_OUT_SOFT || true)
NFILT_ABHET=$(grep -c "ABhet_filt" $VCF_OUT_SOFT || true)

echo -e "\n\n"
echo -e "#### 04_filter_gatk.sh: Number of SNPs in input VCF: $NVAR_IN"
echo -e "#### 04_filter_gatk.sh: Number of SNPs in output VCF: $NVAR_OUT"
echo -e "#### 04_filter_gatk.sh: Number of SNPs filtered: $NVAR_FILT\n"

echo -e "#### 04_filter_gatk.sh: Number of SNPs filtered by FS_gt60: $NFILT_FS"
echo -e "#### 04_filter_gatk.sh: Number of SNPs filtered by MQ_lt40: $NFILT_MQ"
echo -e "#### 04_filter_gatk.sh: Number of SNPs filtered by MQRankSum_ltm12: $NFILT_MQR"
echo -e "#### 04_filter_gatk.sh: Number of SNPs filtered by readposRankSum_ltm8: $NFILT_READPOS"
echo -e "#### 04_filter_gatk.sh: Number of SNPs filtered by ABhet_filt: $NFILT_ABHET"

echo -e "\n#### 04_filter_gatk.sh: Listing output VCF:"
ls -lh $VCF_OUT
[[ $(grep -cv "^#" $VCF_OUT) = 0 ]] && echo -e "\n\n#### 04_filter_gatk.shh: ERROR: VCF is empty\n" >&2 && exit 1

echo -e "\n#### 04_filter_gatk.sh: Done with script."
date

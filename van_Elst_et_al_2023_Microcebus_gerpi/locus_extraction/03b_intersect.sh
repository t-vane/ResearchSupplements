#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Script adapted and modified from Poelstra et al. 2021, Systematic Biology (https://doi.org/10.1093/sysbio/syaa053)

## Software:
# BEDtools needs to be included in $PATH (v2.30.0; https://bedtools.readthedocs.io/en/latest/)

## Command-line args:
LOCUSBED_INTERMED=$1
LOCUSBED_FINAL=$2
VCF_HIGHDEPTH=$3
VCF_FILT_INTERSECT=$4

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### 03b_intersect.sh: Starting script."
echo -e "#### 03b_intersect.sh: Intermediate locus BED file: $LOCUSBED_INTERMED"
echo -e "#### 03b_intersect.sh: Final locus BED file: $LOCUSBED_FINAL"
echo -e "#### 03b_intersect.sh: VCF file with loci with too high depth: $VCF_HIGHDEPTH"
echo -e "#### 03b_intersect.sh: Fully filtered VCF file: $VCF_FILT_INTERSECT \n\n"

################################################################################
#### INTERSECT BED FILE WITH LOCI WITH TOO HIGH DEPTH ####
################################################################################
echo -e "#### 03b_intersect.sh: Intersecting BED file with loci with too high depth ..."
bedtools intersect -v -a $LOCUSBED_INTERMED -b $VCF_HIGHDEPTH > $LOCUSBED_FINAL

## Report:
echo -e "#### 03b_intersect.sh: Number of loci after removing SNPs with too high depth: $(wc -l < $LOCUSBED_FINAL)"

SNPS_IN_LOCI=$(bedtools intersect -u -a $VCF_FILT_INTERSECT -b $LOCUSBED_FINAL | grep -cv "##")
SNPS_IN_VCF=$(grep -cv "##" $VCF_FILT_INTERSECT)
SNPS_LOST=$(( $SNPS_IN_VCF - $SNPS_IN_LOCI ))

echo -e "#### 03b_intersect.sh: Number of SNPs in VCF: $SNPS_IN_VCF"
echo -e "#### 03b_intersect.sh: Number of SNPs in loci: $SNPS_IN_LOCI"
echo -e "#### 03b_intersect.sh: Number of lost SNPs (in VCF but not in loci): $SNPS_LOST"

echo -e "\n#### 03b_intersect.sh: Done with script."
date


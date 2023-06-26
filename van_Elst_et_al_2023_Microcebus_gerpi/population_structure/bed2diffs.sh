#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
# BCFtools needs to be included in $PATH (v0.1.17; https://vcftools.github.io/index.html)
# bed2diffs of EEMS needs to be included in $PATH (https://github.com/dipetkov/eems)
# PLINK needs to be included in $PATH (v1.90b6.22; https://zzz.bwh.harvard.edu/plink/)

## Command-line args:
NT=$1
EEMS_DIR=$2
VCF_FILE=$3
CHROM_FILE=$4

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### bed2diffs.sh: Starting script."
echo -e "#### bed2diffs.sh: Number of threads: $NT"
echo -e "#### bed2diffs.sh: Directory for estimated effective migration surfaces: $EEMS_DIR"
echo -e "#### bed2diffs.sh: Input VCF file: $VCF_IN"
echo -e "#### bed2diffs.sh: Renaming file for chromosomes: $CHROM_FILE \n\n"

################################################################################
#### ESTIMATE AVERAGE GENETIC DISSIMILARITY MATRIX ####
################################################################################
echo -e "#### bed2diffs.sh: Renaming chromosomes ...\n"
bcftools annotate --rename-chrs $CHROM_FILE $VCF_FILE > $VCF_FILE.tmp

echo -e "#### bed2diffs.sh: Converting VCF to BED file ...\n"
plink --vcf $VCF_FILE.tmp --make-bed --double-id --chr-set 32 --out $EEMS_DIR/$(basename $VCF_FILE .vcf)
rm $VCF_FILE.tmp

echo -e "#### bed2diffs.sh: Estimating average genetic dissimilarity matrix ...\n"
bed2diffs_v1 --bfile $EEMS_DIR/$(basename $VCF_FILE .vcf) --nthreads $NT


echo -e "\n#### bed2diffs.sh: Done with script."
date


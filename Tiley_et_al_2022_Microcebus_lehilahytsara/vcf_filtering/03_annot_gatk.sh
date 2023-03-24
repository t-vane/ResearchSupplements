#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Script adapted and modified from Poelstra et al. 2021, Systematic Biology (https://doi.org/10.1093/sysbio/syaa053)

## Software:
# SAMtools needs to be included in $PATH (v1.11; http://www.htslib.org/)
# VCFtools needs to be included in $PATH (v0.1.17; https://vcftools.github.io/index.html)
GATK3=/home/nibtve93/software/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar # (v3.8.1; https://gatk.broadinstitute.org/hc/en-us)

## Command-line args:
VCF_IN=$1
VCF_OUT=$2
BAM_DIR=$3
REFERENCE=$4
SUFFIX=$5

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### 03_annot_gatk.sh: Starting script."
echo -e "#### 03_annot_gatk.sh: Input VCF: $VCF_IN"
echo -e "#### 03_annot_gatk: Output VCF: $VCF_OUT"
echo -e "#### 03_annot_gatk: Directory with BAM files: $BAM_DIR"
echo -e "#### 03_annot_gatk.sh: Reference genome: $REFERENCE"
echo -e "#### 03_annot_gatk.sh: Suffix for BAM files: $SUFFIX \n\n"

################################################################################
#### ANNOTATE VCF FILE ####
################################################################################
## Sort VCF file
echo -e "#### 03_annot_gatk.sh: Sorting $VCF_IN ...\n"
vcf-sort $VCF_IN > $VCF_IN.tmp
mv $VCF_IN.tmp > $VCF_IN

## Create list of BAM files
echo -e "#### 03_annot_gatk.sh: Creating list of BAM files for individuals in $VCF_IN ...\n"
for i in $(bcftools query -l $VCF_IN)
do
	echo "$BAM_DIR/${i}.$SUFFIX.bam"
done > $(dirname $VCF_IN)/bamFiles.txt

## Annotate VCF file
# Process list of BAM files to create command argument
BAM_ARG=""
while read -r BAM
do
  [[ ! -f $BAM.bai ]] && samtools index $BAM
  BAM_ARG=$(echo "$BAM_ARG -I $BAM")
done < $(dirname $VCF_IN)/bamFiles.txt

# Run GATK for annotation
echo -e "#### 03_annot_gatk.sh: Annotating VCF_IN with INFO fields for\n"
echo -e "#### 03_annot_gatk.sh: FisherStrand, RMSMappingQuality, MappingQualityRankSumTest, ReadPosRankSumTest and AlleleBalance ...\n"
java -jar $GATK3 -T VariantAnnotator -R $REFERENCE -V $VCF_IN -o $VCF_OUT $BAM_ARG -A FisherStrand -A RMSMappingQuality -A MappingQualityRankSumTest -A ReadPosRankSumTest -A AlleleBalance

## Report:
echo -e "\n#### 03_annot_gatk.sh: Listing output VCF:"
ls -lh $VCF_OUT
[[ $(grep -cv "^#" $VCF_OUT) = 0 ]] && echo -e "\n\n#### 03_annot_gatk.sh: ERROR: VCF is empty\n" >&2 && exit 1

echo -e "\n#### 03_annot_gatk.sh: Done with script."
date


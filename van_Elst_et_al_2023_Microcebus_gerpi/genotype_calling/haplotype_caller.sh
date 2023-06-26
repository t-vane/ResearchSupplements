#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
# SAMtools needs to be included in $PATH (v1.11; http://www.htslib.org/)
# gatk needs to be included in $PATH (v4.1.9.0; https://gatk.broadinstitute.org/hc/en-us)

## Command-line args:
NT=$1
MEM=$2
REFERENCE=$3
IND_FILE=$4
BAM_DIR=$5
SUFFIX=$6
GVCF_DIR=$7

INDV=$(sed -n "$SLURM_ARRAY_TASK_ID"p $IND_FILE)

## Activate conda environment
conda activate java

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### haplotype_caller.sh: Starting script."
echo -e "#### haplotype_caller.sh: Number of threads: $NT"
echo -e "#### haplotype_caller.sh: Memory: $MEM"
echo -e "#### haplotype_caller.sh: Reference genome: $REFERENCE"
echo -e "#### haplotype_caller.sh: List with individuals: $IND_FILE"
echo -e "#### haplotype_caller.sh: Directory with BAM files: $BAM_DIR"
echo -e "#### haplotype_caller.sh: Suffix of BAM files: $SUFFIX"
echo -e "#### haplotype_caller.sh: Output directory for GVCF files: $GVCF_DIR"
echo -e "#### haplotype_caller.sh: Individual: $INDV \n\n"

################################################################################
#### CREATE GVCF FILE ####
################################################################################
echo -e "#### haplotype_caller.sh: Indexing BAM file for individual $INDV ...\n"
[[ ! -f $BAM_DIR/$INDV.$SUFFIX.bai ]] && samtools index $BAM_DIR/$INDV.$SUFFIX

echo -e "#### haplotype_caller.sh: Creating GVCF file for individual $INDV ...\n"
gatk --java-options "-Xmx${MEM}g" HaplotypeCaller -R $REFERENCE -I $BAM_DIR/$INDV.$SUFFIX -O $GVCF_DIR/$INDV.rawvariants.g.vcf.gz -ERC GVCF --pairHMM AVX_LOGLESS_CACHING_OMP --native-pair-hmm-threads $NT

## Report:
echo -e "\n#### haplotype_caller.sh: Done with script."
date

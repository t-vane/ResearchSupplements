#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Script adapted and modified from Poelstra et al. 2021, Systematic Biology (https://doi.org/10.1093/sysbio/syaa053)

## Software:
# SAMtools needs to be included in $PATH (v1.11; http://www.htslib.org/)

## Command-line args:
LOCUSLIST=$1
LOCUSFASTA_DIR_INTERMED=$2
FASTA_MERGED=$3
LOCUS=$(sed -n "$SLURM_ARRAY_TASK_ID"p $LOCUSLIST)

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### 03d_locusfasta.sh: Starting script."
echo -e "#### 03d_locusfasta.sh: List of loci: $LOCUSLIST"
echo -e "#### 03d_locusfasta.sh: Directory for intermediate locus FASTA files: $LOCUSFASTA_DIR_INTERMED"
echo -e "#### 03d_locusfasta.sh: Merged FASTA file: $FASTA_MERGED"
echo -e "#### 03d_locusfasta.sh: Locus: $LOCUS \n\n"

################################################################################
#### CREATE BY-LOCUS FASTA FILE  ####
################################################################################
echo -e "\n#### 03d_locusfasta.sh: Creating FASTA file for locus $LOCUS ..."
LOCUS_FAIDX=$(echo $LOCUS | sed 's/:/,/g')
FASTA=$LOCUSFASTA_DIR_INTERMED/$LOCUS.fa
faidx --regex $LOCUS_FAIDX $FASTA_MERGED | sed 's/,/:/g' > $FASTA

## Report:
echo -e "\n#### 03d_locusfasta.sh: Done with script."
date


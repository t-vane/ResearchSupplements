#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
# gstacks of Stacks needs to be included in $PATH (v2.53; http://catchenlab.life.illinois.edu/stacks/)

## Command-line args:
NT=$1
IN_DIR=$2
OUT_DIR=$3
POPMAP=$4
SUFFIX=$5

## Activate conda environment
conda activate stacks

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### gstacks.sh: Starting script."
echo -e "#### gstacks.sh: Number of threads: $NT"
echo -e "#### gstacks.sh: Directory with BAM files: $IN_DIR"
echo -e "#### gstacks.sh: Output directory: $OUT_DIR"
echo -e "#### gstacks.sh: Population map: $POPMAP"
echo -e "#### gstacks.sh: Suffix for BAM files: $SUFFIX \n\n"

################################################################################
#### CREATE STACKS WITH GSTACKS ####
################################################################################
echo -e "#### gstacks.sh: Creating stacks with gstacks for individuals in $POPMAP ...\n"
gstacks -t $NT -I $IN_DIR -O $OUT_DIR -M $POPMAP -S $SUFFIX.bam

## Report:
echo -e "\n#### gstacks.sh: Done with script."
date
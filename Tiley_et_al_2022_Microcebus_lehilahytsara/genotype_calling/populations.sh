#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
# populations of Stacks needs to be included in $PATH (v2.53; http://catchenlab.life.illinois.edu/stacks/)

## Command-line args:
NT=$1
IN_DIR=$2
POP_DIR=$3
POPMAP=$4
FILTERS=$5
OUTPUT=$6

## Activate conda environment:
conda activate stacks

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### populations.sh: Starting script."
echo -e "#### populations.sh: Number of threads: $NT"
echo -e "#### populations.sh: Directory with created stacks: $IN_DIR"
echo -e "#### populations.sh: Output directory: $POP_DIR"
echo -e "#### populations.sh: Population map: $POPMAP"
echo -e "#### populations.sh: Desired filters: $FILTERS"
echo -e "#### populations.sh: Desired outputs: $OUTPUT \n\n"

################################################################################
#### EXTRACT VCF WITH POPULATIONS ####
################################################################################
echo -e "#### populations.sh: Extracting VCF file from created stacks with populations ...\n"
populations -t $NT -P $IN_DIR -O $POP_DIR -M $POPMAP $FILTERS $OUTPUT

## Report:
echo -e "\n#### populations.sh: Done with script."
date








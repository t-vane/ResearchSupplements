#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
# iqtree2 needs to be included in $PATH (v2.2.0; http://www.iqtree.org/)

## Command-line args:
NT=$1
ALIGNMENT=$2
TREES=$3
OUT_FILE=$4

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### au_test.sh: Starting script."
echo -e "#### au_test.sh: Number of threads: $NT"
echo -e "#### au_test.sh: Input alignment: $ALIGNMENT"
echo -e "#### au_test.sh: File with trees to compare: $TREES"
echo -e "#### au_test.sh: Output file: $OUT_FILE \n\n"

################################################################################
#### APPROXIMATELY UNBIASED TEST ####
################################################################################
echo -e "#### au_test.sh: Conducting approximately unbiased test ...\n"

iqtree2 -T $NC -s $ALIGNMENT -m GTR+G --seqtype DNA -n 0 -z $TREES -zb 10000 -au --prefix $OUT_FILE

## Report:
echo -e "\n#### au_test.sh: Done with script."
date





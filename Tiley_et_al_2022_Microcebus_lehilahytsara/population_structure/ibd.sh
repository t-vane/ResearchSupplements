#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################

## Command-line args
SCRIPT_DIR=$1
GEO_DIST=$2
GEN_DIST=$3
OUT=$4

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### ibd.sh: Starting script."
echo -e "#### ibd.sh: Directory with scripts: $SCRIPT_DIR"
echo -e "#### ibd.sh: Geographic distance matrix: $GEO_DIST"
echo -e "#### ibd.sh: Genetic distance matrix: $GEN_DIST"
echo -e "#### ibd.sh: Output prefix: $OUT \n\n"

################################################################################
#### CONDUCT MANTEL TEST AND PLOT ISOLATION BY DISTANCE ####
################################################################################
echo -e "#### ibd.sh: Conducting Mantel test and plotting ...\n"
Rscript $SCRIPTS_DIR/ibd.R $GEO_DIST $GEN_DIST $OUT

## Report:
echo -e "\n#### ibd.sh: Done with script."
date


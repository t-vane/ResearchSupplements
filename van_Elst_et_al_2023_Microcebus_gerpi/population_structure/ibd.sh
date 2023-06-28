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
STRING=$4
OUT=$5
GEN_DIST_SD=$6

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### ibd.sh: Starting script."
echo -e "#### ibd.sh: Directory with scripts: $SCRIPT_DIR"
echo -e "#### ibd.sh: Geographic distance matrix: $GEO_DIST"
echo -e "#### ibd.sh: Genetic distance matrix: $GEN_DIST"
echo -e "#### ibd.sh: String to specify IBD script to be used: $STRING"
echo -e "#### ibd.sh: Output prefix: $OUT"
echo -e "#### ibd.sh: Genetic distance standard deviation matrix (if applicable): $GEN_DIST_SD \n\n"

################################################################################
#### CONDUCT MANTEL TEST AND PLOT ISOLATION BY DISTANCE ####
################################################################################
echo -e "#### ibd.sh: Conducting Mantel test and plotting ...\n"
Rscript $SCRIPTS_DIR/ibd$STRING.R $GEO_DIST $GEN_DIST $OUT $GEN_DIST_SD

## Report:
echo -e "\n#### ibd.sh: Done with script."
date


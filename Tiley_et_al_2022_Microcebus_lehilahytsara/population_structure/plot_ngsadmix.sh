#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################

## Command-line args:
SCRIPTS_DIR=$1
OUT_DIR=$2
LIKE_VALUES=$3
IND_FILE=$4
SET_ID=$5

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### plot_ngsadmix.sh: Starting script."
echo -e "#### plot_ngsadmix.sh: Directory with scripts: $SCRIPTS_DIR"
echo -e "#### plot_ngsadmix.sh: Output directory: $OUT_DIR"
echo -e "#### plot_ngsadmix.sh: File with likelihood values: $LIKE_VALUES"
echo -e "#### plot_ngsadmix.sh: File that maps individuals to populations: $IND_FILE"
echo -e "#### plot_ngsadmix.sh: Set ID: $SET_ID \n\n"

################################################################################
#### PLOT ADMIXTURE RESULTS ####
################################################################################
echo -e "#### plot_ngsadmix.sh: Running script to plot admixture results ...\n"
Rscript $SCRIPTS_DIR/plot_ngsadmix.R $OUT_DIR $LIKE_VALUES $IND_FILE $SET_ID

echo -e "\n#### plot_ngsadmix.sh: Done with script."
date


#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################

## Command-line args:
SCRIPTS_DIR=$1
IN_FILE=$2
BAM_HITS=$3
SET_ID=$4

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### cov_plot.sh: Starting script."
echo -e "#### cov_plot.sh: Directory with scripts: $SCRIPTS_DIR"
echo -e "#### cov_plot.sh: File with individuals: $IN_FILE"
echo -e "#### cov_plot.sh: Directory with BAM hits: $BAM_HITS"
echo -e "#### cov_plot: Set ID: $SET_ID \n\n"

################################################################################
#### ESTIMATE AND PLOT COVERAGE DISTRIBUTIONS ####
################################################################################
echo -e "#### cov_plot.sh: Submitting R script to estimate and plot coverage distributions for each individual ...\n"
Rscript $SCRIPTS_DIR/cov_plot.R $IN_FILE $BAM_HITS $SET_ID

echo -e "\n#### cov_plot.sh: Done with script."
date







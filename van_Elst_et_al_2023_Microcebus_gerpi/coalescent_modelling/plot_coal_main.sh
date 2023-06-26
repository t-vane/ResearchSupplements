#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################

## Command-line args:
SCRIPTS_DIR=$1
MCMC_NM=$2
MCMC_M=$3
SUMMARY=$4
OUT_DIR=$5

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### plot_coal_main.sh: Starting script."
echo -e "#### plot_coal_main.sh: Directory with scripts: $SCRIPTS_DIR"
echo -e "#### plot_coal_main.sh: MCMC (no migration): $MCMC_NM"
echo -e "#### plot_coal_main.sh: MCMC (migration): $MCMC_M"
echo -e "#### plot_coal_main.sh: Summary table with divergence times, migration rates and genealogical divergence indices: $MCMC_M"
echo -e "#### plot_coal_main.sh: Output directory: $OUT_DIR \n\n"

################################################################################
#### PLOT MAIN FIGURES ####
################################################################################
echo -e "#### plot_coal_main.sh: Plotting models ..."
Rscript $SCRIPTS_DIR/plot_coal_models.R $MCMC_NM $MCMC_M $OUT_DIR

echo -e "#### plot_coal_main.sh: Plotting parameter estimates ..."
Rscript $SCRIPTS_DIR/plot_coal_params.R $SUMMARY $OUT_DIR

## Report:
echo -e "\n#### plot_coal_main.sh: Done with script."
date



#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################

## Command-line args:
SCRIPTS_DIR=$1
MCMC_NM1=$2
MCMC_NM2=$3
MCMC_NM3=$4
MCMC_NM4=$5
MCMC_M1=$6
MCMC_M2=$7
MCMC_M3=$8
MCMC_M4=$9
M_SCALE=$10
T_SCALE=$11
OUT_DIR=$12

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### plot_coal_posteriors.sh: Starting script."
echo -e "#### plot_coal_posteriors.sh: Directory with scripts: $SCRIPTS_DIR"
echo -e "#### plot_coal_posteriors.sh: First MCMC (no migration): $MCMC_NM1"
echo -e "#### plot_coal_posteriors.sh: Second MCMC (no migration): $MCMC_NM2"
echo -e "#### plot_coal_posteriors.sh: Third MCMC (no migration): $MCMC_NM3"
echo -e "#### plot_coal_posteriors.sh: Fourth MCMC (no migration): $MCMC_NM4"
echo -e "#### plot_coal_posteriors.sh: First MCMC (migration): $MCMC_M1"
echo -e "#### plot_coal_posteriors.sh: First MCMC (migration): $MCMC_M2"
echo -e "#### plot_coal_posteriors.sh: First MCMC (migration): $MCMC_M3"
echo -e "#### plot_coal_posteriors.sh: First MCMC (migration): $MCMC_M4"
echo -e "#### plot_coal_posteriors.sh: Inverse scaling factor used in the G-PhoCS configuration file for migration parameter: $M_SCALE"
echo -e "#### plot_coal_posteriors.sh: Inverse scaling factor used in the G-PhoCS configuration file for tau and theta: $T_SCALE"
echo -e "#### plot_coal_posteriors.sh: Output directory: $OUT_DIR \n\n"

################################################################################
#### PLOT POSTERIORS ####
################################################################################
echo -e "#### plot_coal_posteriors.sh: Plotting posteriors ..."
Rscript $SCRIPTS_DIR/plot_coal_posteriors.R $MCMC_NM1 $MCMC_NM2 $MCMC_NM3 $MCMC_NM4 $MCMC_M1 $MCMC_M2 $MCMC_M3 $MCMC_M4 $M_SCALE $T_SCALE $OUT_DIR

## Report:
echo -e "\n#### plot_coal_posteriors.sh: Done with script."
date

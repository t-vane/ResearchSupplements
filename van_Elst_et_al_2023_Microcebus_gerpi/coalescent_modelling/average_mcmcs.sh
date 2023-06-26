#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################

## Command-line args:
SCRIPTS_DIR=$1
MCMC1=$2
MCMC2=$3
MCMC3=$4
MCMC4=$5
MCMC_OUT=$6

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### average_mcmcs.sh: Starting script."
echo -e "#### average_mcmcs.sh: Directory with scripts: $SCRIPTS_DIR"
echo -e "#### average_mcmcs.sh: First MCMC: $MCMC1"
echo -e "#### average_mcmcs.sh: Second MCMC: $MCMC2"
echo -e "#### average_mcmcs.sh: Third MCMC: $MCMC3"
echo -e "#### average_mcmcs.sh: Fourth MCMC: $MCMC4"
echo -e "#### average_mcmcs.sh: MCMC output: $MCMC_OUT \n\n"

################################################################################
#### AVERAGE MCMCS ####
################################################################################
echo -e "#### average_mcmcs.sh: Averaging MCMCs ..."
Rscript $SCRIPTS_DIR/average_mcmcs.R $MCMC1 $MCMC2 $MCMC3 $MCMC4 $MCMC_OUT

## Report:
echo -e "\n#### average_mcmcs.sh: Done with script."
date



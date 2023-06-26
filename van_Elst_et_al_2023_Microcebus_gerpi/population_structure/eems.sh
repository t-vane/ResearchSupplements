#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################

## Command-line args
CONFIG=$1
SEED=$2

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### eems.sh: Starting script."
echo -e "#### eems.sh: Configuration file: $CONFIG"
echo -e "#### eems.sh: Seed: $SEED \n\n"

################################################################################
#### ESTIMATE EFFECTIVE MIGRATION SURFACE ####
################################################################################
echo -e "#### eems.sh: Estimating effective migration surface ...\n"
runeems_snps --params $CONFIG --seed=$SEED

## Report:
echo -e "\n#### eems.sh: Done with script."
date


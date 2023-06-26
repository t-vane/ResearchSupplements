#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################

## Command-line args
SCRIPTS_DIR=$1
INPUT=$2
OUT_DIR=$3
POP_COORDS=$4
SHAPE=$5

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### plot_eems.R: Starting script."
echo -e "#### plot_eems.R: Directory with scripts: $SCRIPTS_DIR"
echo -e "#### plot_eems.R: Input prefix: $INPUT"
echo -e "#### plot_eems.R: Output directory: $OUT_DIR"
echo -e "#### plot_eems.R: File with population coordinates: $POP_COORDS"
echo -e "#### plot_eems.R: Shape file prefix: $SHAPE \n\n"

################################################################################
#### PLOT ESTIMATED EFFECTIVE MIGRATION SURFACE ####
################################################################################
echo -e "#### plot_eems.R: Plotting estimated effective migration surface ...\n"
Rscript $SCRIPTS_DIR/plot_eems.R $INPUT $OUT_DIR $POP_COORDS $SHAPE

## Report:
echo -e "\n#### plot_eems.R: Done with script."
date


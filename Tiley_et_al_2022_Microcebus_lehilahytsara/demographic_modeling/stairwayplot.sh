#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
STAIRWAYPLOT=/home/nibtve93/software/stairway_plot_v2.1.1/stairway_plot_es # (v2.0; https://github.com/xiaoming-liu/stairway-plot-v2)

## Command-line args:
BLUEPRINT=$1

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### stairwayplot.sh: Starting script."
echo -e "#### stairwayplot.sh: Additional options: $OPTIONS \n\n"

################################################################################
#### ESTIMATE POPULATION SIZE CHANGES ####
################################################################################
echo -e "#### stairwayplot.sh: Creating shell script ...\n"
java -cp $STAIRWAYPLOT Stairbuilder $BLUEPRINT

echo -e "#### stairwayplot.sh: Running Stairway Plot ...\n"
bash $BLUEPRINT.sh 

## Report:
echo -e "\n#### stairwayplot.sh: Done with script."
date

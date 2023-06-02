#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################

## Command-line args:
LOG=$1
LIKE_FILE=$2

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### print_likes.sh: Starting script."
echo -e "#### print_likes.sh: Respective log file: $LOG"
echo -e "#### print_likes.sh: Likelihoods summary file: $LIKE_FILE \n\n"

################################################################################
#### EXTRACT LIKELIHOOD ####
################################################################################
echo -e "#### print_likes.sh: Extracting likelihood ...\n"
grep "best" $LOG | awk '{print $K}' | cut -d'=' -f2- | sort -g | sed "s/after/$K/g" | sed "s/iterations/$SEED/g" >> $LIKE_FILE

echo -e "#### print_likes.sh: Removing fopt.gz files ...\n"
rm $(dirname $LOG)/$(basename $LOG .log).fopt.gz

echo -e "\n#### print_likes.sh: Done with script."
date


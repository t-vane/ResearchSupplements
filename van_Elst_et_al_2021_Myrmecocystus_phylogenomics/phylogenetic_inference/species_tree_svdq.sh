#!/bin/bash

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
# paup needs to be included in $PATH

## Command-line args:
IN_FILE=$1
LOG_FILE=$2

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### species_tree_svdq.sh: Starting script."
echo -e "#### species_tree_svdq.sh: PAUP input file: $IN_FILE"
echo -e "#### species_tree_svdq.sh: Log file: $LOG_FILE \n\n"

################################################################################
#### SPECIES TREE INFERENCE WITH SVDquartets####
################################################################################
cd $(dirname $IN_FILE)

echo -e "#### species_tree_svdq.sh: Species tree inference with SVDquartets... \n"
paup4a168_ubuntu64 -n $IN_FILE $LOG_FILE

## Report:
echo -e "\n#### species_tree_svdq.sh: Done with script."
date

#!/bin/bash

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
# paup needs to be included in $PATH (v4.0a; https://paup.phylosolutions.com/)

## Command-line args:
IN_FILE=$1
LOG_FILE=$2

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### svdq.sh: Starting script."
echo -e "#### svdq.sh: PAUP input file: $IN_FILE"
echo -e "#### svdq.sh: Log file: $LOG_FILE \n\n"

################################################################################
#### PHYLOGENETIC INFERENCE WITH SVDquartets####
################################################################################
cd $(dirname $IN_FILE)

echo -e "#### svdq.sh: Phylogenetic inference with SVDquartets... \n"
paup4a168_ubuntu64 -n $IN_FILE $LOG_FILE

## Report:
echo -e "\n#### svdq.sh: Done with script."
date

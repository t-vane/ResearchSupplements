#!/bin/bash

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
# baseml of PAML needs to be included in $PATH (v4.8; http://abacus.gene.ucl.ac.uk/software/paml.html)

## Command-line args:
CTL_FILE=$1

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### baseml.sh: Starting script."
echo -e "#### baseml.sh: baseml control file: $CTL_FILE \n\n"

################################################################################
#### BASEML ANALYSIS####
################################################################################
echo -e "#### baseml.sh: Running baseml ... \n"
baseml $CTL_FILE

## Report:
echo -e "\n#### baseml.sh: Done with script."
date

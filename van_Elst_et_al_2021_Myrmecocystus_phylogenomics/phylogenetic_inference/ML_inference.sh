#!/bin/bash

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
# iqtree needs to be included in $PATH

## Command-line args:
ALIGNMENT=$1
PARTITIONS=$2

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### ML_inference.sh: Starting script."
echo -e "#### ML_inference.sh: Alignment: $ALIGNMENT"
echo -e "#### ML_inference.sh: Partitioning scheme: $PARTITIONS \n\n"

################################################################################
#### CREATE PARTITIONS####
################################################################################
cd $(dirname $ALIGNMENT)

echo -e "#### ML_inference.sh: Maximum likelihood inference ... \n"
iqtree -s $ALIGNMENT -nt AUTO -ntmax 10 -bb 1000 -wbt -spp $PARTITIONS -nstop 200

## Report:
echo -e "\n#### ML_inference.sh: Done with script."
date



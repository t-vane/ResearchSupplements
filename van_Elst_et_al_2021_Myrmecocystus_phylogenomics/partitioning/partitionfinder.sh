#!/bin/bash

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
PFINDER=/global/homes/jg/t_vane02/software/partitionfinder-2.1.1/PartitionFinder.py # (v2.1.1; https://www.robertlanfear.com/partitionfinder/)

## Command-line args:
CONF_FILE=$1

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### partitionfinder.sh: Starting script."
echo -e "#### partitionfinder.sh: Configuration file: $CONF_FILE \n\n"

################################################################################
#### CREATE PARTITIONS####
################################################################################
echo -e "#### partitionfinder.sh: Running ParitionFinder for file $CONF_FILE ... \n"
python2.7 $PFINDER -r $CONF_FILE

## Report:
echo -e "\n#### partitionfinder.sh: Done with script."
date

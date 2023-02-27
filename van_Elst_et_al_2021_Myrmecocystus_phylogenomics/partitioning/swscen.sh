#!/bin/bash

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
SWSCEN=/global/homes/jg/t_vane02/software/PFinderUCE-SWSC-EN-master/py_script/SWSCEN.py

## Command-line args:
IN_FILE=$1

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### swscen.sh: Starting script."
echo -e "#### swscen.sh: Input file: $IN_FILE \n\n"

################################################################################
#### CREATE TRIPLET PARTITIONS WITH SWSC-EN####
################################################################################
echo -e "#### swscen.sh: Running SWSC-EN ... \n"
python3 $SWSCEN $IN_FILE

## Report:
echo -e "\n#### swscen.sh: Done with script."
date


#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
# SplitsTree needs to be included in $PATH (v4.19.0; https://uni-tuebingen.de/fakultaeten/mathematisch-naturwissenschaftliche-fakultaet/fachbereiche/informatik/lehrstuehle/algorithms-in-bioinformatics/software/splitstree/)

## Command-line args:
IN_FILE=$1
OUT_FILE=$2

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### splitstree.sh: Starting script."
echo -e "#### splitstree.sh: Input alignment file: $IN_FILE"
echo -e "#### splitstree.sh: Output file: $OUT_FILE \n\n"

################################################################################
#### PREPARE NEXUS FILE FOR GUI VERSION OF SPLITSTREE ####
################################################################################
echo -e "#### splitstree.sh: Preparing NEXUS file for GUI version of SplitsTree ..."
SplitsTree -g -i $IN_FILE -x "UPDATE; SAVE REPLACE=yes FILE=$OUT_FILE.tmp; QUIT"

echo -e "#### splitstree.sh: Removing sequence from NEXUS output file ..."
FIRST_LINE=$(grep -n "BEGIN Characters;" $OUT_FILE.tmp | cut -f1 -d:)
LAST_LINE=$(grep -n "END;.*Characters" $OUT_FILE.tmp | cut -f1 -d:)
sed "$FIRST_LINE,${LAST_LINE}d" $OUT_FILE.tmp > $OUT_FILE
rm -f $OUT_FILE.tmp

## Report:
echo -e "#### splitstree.sh: Done with script."
date

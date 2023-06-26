#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
ASCBIAS=/home/nibtve93/software/raxml_ascbias/ascbias.py # (https://github.com/btmartin721/raxml_ascbias)
AMAS=/home/nibtve93/software/AMAS/amas/AMAS.py # (https://github.com/marekborowiec/AMAS)

## Command-line args:
IN_FILE=$1
OUT_FILE=$2

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### remove_invariants.sh: Starting script."
echo -e "#### remove_invariants.sh: Input alignment: $IN_FILE"
echo -e "#### remove_invariants.sh: Output alignment: $OUT_FILE \n\n"

################################################################################
#### REMOVE INVARIANT SITES ####
################################################################################
echo -e "#### remove_invariants.sh: Removing invariant sites ..."
python $ASCBIAS -p $IN_FILE -o $OUT_FILE

## Report:
echo -e "#### remove_invariants.sh: Done with script."
date




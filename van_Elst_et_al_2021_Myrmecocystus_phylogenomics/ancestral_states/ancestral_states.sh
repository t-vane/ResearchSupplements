#!/bin/bash

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
ANC_STATES=/global/homes/jg/t_vane02/scripts/ancestral_states.R

## Command-line args:
TREE=$1
CHARS=$2
OUT_DIR=$3

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### ancestral_states.sh: Starting script."
echo -e "#### ancestral_states.sh: Tree file: $TREE"
echo -e "#### ancestral_states.sh: Character mappings: $CHARS"
echo -e "#### ancestral_states.sh: Output directory: $OUT_DIR \n\n"

################################################################################
#### Reconstruction of ancestral states####
################################################################################
mkdir -p $OUT_DIR

echo -e "#### ancestral_states.sh: Reconstructing ancestral states ... \n"
Rscript $ANC_STATES $TREE $CHARS $OUT_DIR

## Report:
echo -e "\n#### ancestral_states.sh: Done with script."
date

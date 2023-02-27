#!/bin/bash

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
ASTRAL=/global/homes/jg/t_vane02/software/Astral/astral.5.6.3.jar

## Command-line args:
GENE_TREES=$1
MAPPING=$2
OUT_FILE=$3

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### species_tree_astral.sh: Starting script."
echo -e "#### species_tree_astral.sh: File with gene trees: $GENE_TREES"
echo -e "#### species_tree_astral.sh: Mapping of specimens to species: $MAPPING"
echo -e "#### species_tree_astral.sh: Output file: $OUT_FILE \n\n"

################################################################################
#### SPECIES TREE INFERENCE WITH ASTRAL####
################################################################################
echo -e "#### species_tree_astral.sh: Species tree inference with ASTRAL ... \n"
java -jar $ASTRAL -i $GENE_TREES -a $MAPPING -o $OUT_FILE

## Report:
echo -e "\n#### species_tree_astral.sh: Done with script."
date

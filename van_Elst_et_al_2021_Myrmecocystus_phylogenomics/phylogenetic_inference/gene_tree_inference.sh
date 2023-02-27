#!/bin/bash

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
SEQKIT=/global/homes/jg/t_vane02/software/seqkit
# iqtree needs to be included in $PATH

## Command-line args:
ALIGNMENT_FILE=$1
ALIGNMENT=$(sed -n "$SGE_TASK_ID"p $ALIGNMENT_FILE)

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### gene_tree_inference.sh: Starting script."
echo -e "#### gene_tree_inference.sh: File with input alignments: $ALIGNMENT_FILE \n\n"

################################################################################
#### ALIGNMENT CLEANING AND GENE TREE INFERENCE####
################################################################################
cd $(dirname $ALIGNMENT)

echo -e "#### gene_tree_inference.sh: Removing taxa with gaps or missing data only ... \n"
$SEQKIT fx2tab $ALIGNMENT | awk '{if ($2 !~ "^?+$") print $0}' | $SEQKIT tab2fx > $(dirname $ALIGNMENT)/$(basename $ALIGNMENT .fas)_clean.fas

echo -e "#### gene_tree_inference.sh: Infering gene tree ... \n"
iqtree -s $(dirname $ALIGNMENT)/$(basename $ALIGNMENT .fas)_clean.fas -nt 2 -bb 1000 -wbt -nstop 200 -mset GTR+G+I,GTR+G,GTR -merit AICc

## Report:
echo -e "\n#### gene_tree_inference.sh: Done with script."
date


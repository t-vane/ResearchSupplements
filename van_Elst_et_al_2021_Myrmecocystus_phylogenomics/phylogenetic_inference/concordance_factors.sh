#!/bin/bash

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
IQTREE2=/global/homes/jg/t_vane02/software/iqtree-2.0.4-Linux/bin/iqtree # (v2.0.5; http://www.iqtree.org/)
CONCORD=/global/homes/jg/t_vane02/scripts/concordance.R

## Command-line args:
NT=$1
TREE=$2
GENE_TREES=$3
ALIGNMENT=$4
OUT_DIR=$5

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### concordance_factors.sh: Starting script."
echo -e "#### concordance_factors.sh: Number of threads: $NT"
echo -e "#### concordance_factors.sh: Tree file: $TREE"
echo -e "#### concordance_factors.sh: File with gene trees: $GENE_TREES"
echo -e "#### concordance_factors.sh: Alignment file: $ALIGNMENT"
echo -e "#### concordance_factors.sh: Output directory: $OUT_DIR \n\n"

################################################################################
#### CONCORDANCE FACTOR ANALYSIS####
################################################################################
echo -e "#### concordance_factors.sh: Calculating site and gene concordance ... \n"
$IQTREE2 -nt $NT -t $TREE --gcf $GENE_TREES -s $ALIGNMENT --scf 100 --prefix $OUT_DIR/concord

echo -e "#### concordance_factors.sh: Testing incomplete lineage sorting and calculating internode certainty ... \n"
Rscript $CONCORD $OUT_DIR

## Report:
echo -e "\n#### concordance_factors.sh: Done with script."
date

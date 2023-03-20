#!/bin/bash

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
BINNING_SCRIPTS=/global/homes/jg/t_vane02/software/binning-master/ # Pipeline by S. Mirarab is used (https://github.com/smirarab/binning)

## Command-line args:
GENES_DIR=$1
SUPPORT=$2
TREE=$4

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### statistical_binning.sh: Starting script."
echo -e "#### statistical_binning.sh: Directory containing the gene folders with alignments and trees: $GENES_DIR"
echo -e "#### statistical_binning.sh: Support value: $SUPPORT"
echo -e "#### statistical_binning.sh: Name of tree file: $TREE \n\n"

################################################################################
#### RUN BINNING PIPELINE####
################################################################################
mkdir -p $GENES_DIR/output
cd $GENES_DIR/

echo -e "#### statistical_binning.sh.sh: Preparing ... \n"
$BINNING_SCRIPTS/makecommands.compatibility.sh $GENES_DIR $SUPPORT $GENES_DIR/output $TREE
parallel < commands*

echo -e "#### statistical_binning.sh.sh: Building bin definitions ... \n"
ls | grep -v ge | sed -e "s/.95$//g" > genes   
python $BINNING_SCRIPTS/cluster_genetrees.py genes $SUPPORT

echo -e "#### statistical_binning.sh.sh: Concatenating gene alignments for each bin to create supergenes ... \n"
mkdir -p $GENES_DIR/output/supergenes
$BINNING_SCRIPTS/build.supergene.alignments.sh $GENES_DIR/output $GENES_DIR $GENES_DIR/output/supergenes

## Report:
echo -e "\n#### statistical_binning.sh.sh: Done with script."
date



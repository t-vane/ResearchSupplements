#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
# RAxML-NG needs to be included in $PATH (v1.0.2; https://github.com/amkozlov/raxml-ng)

## Command-line args:
NT=$1
BOOTSTRAP=$2
OUTGROUP=$3
IN_FILE=$4
OUT_FILE=$5
COMMAND=$6

STAMATAKIS1=$(awk '{print $1}' $IN_FILE.stamatakis)
STAMATAKIS2=$(awk '{print $2}' $IN_FILE.stamatakis)
STAMATAKIS3=$(awk '{print $3}' $IN_FILE.stamatakis)
STAMATAKIS4=$(awk '{print $4}' $IN_FILE.stamatakis)

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### raxml_asc.sh: Starting script."
echo -e "#### raxml_asc.sh: Number of threads: $NT"
echo -e "#### raxml_asc.sh: Number of bootstrap replicates: $BOOTSTRAP"
echo -e "#### raxml_asc.sh: Outgroup: $OUTGROUP"
echo -e "#### raxml_asc.sh: Input alignment file: $IN_FILE"
echo -e "#### raxml_asc.sh: Output prefix: $OUT_FILE"
echo -e "#### raxml_asc.sh: Additional commmands: $COMMAND \n\n"

################################################################################
#### INFER ML PHYLOGENY WITH ASCERTAINMENT BIAS CORRECTION IN RAXML-NG ####
################################################################################
echo -e "#### raxml_asc.sh: Maximum likelihood phylogenetic inference ...\n"
raxml-ng --all --threads $NT --bs-trees $BOOTSTRAP --model GTR+G+ASC_STAM{$STAMATAKIS1/$STAMATAKIS2/$STAMATAKIS3/$STAMATAKIS4} --outgroup $OUTGROUP --seed 12345 --msa $IN_FILE \
	--msa-format PHYLIP --data-type DNA --prefix $OUT_FILE $COMMAND

## Report:
echo -e "\n#### raxml_asc.sh: Done with script."
date


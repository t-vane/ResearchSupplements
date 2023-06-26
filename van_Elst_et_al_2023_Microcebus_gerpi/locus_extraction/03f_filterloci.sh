#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################

## Command-line args:#
SCRIPTS_DIR=$1
LOCUSSTATS=$2
IN_DIR=$3
OUT_DIR=$4
MAXMISS=$5
MINDIST=$6

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### 03f_filterloci.sh: Starting script."
echo -e "#### 03f_filterloci.sh: Directory with scripts: $SCRIPTS_DIR"
echo -e "#### 03f_filterloci.sh: File with locus statistics: $LOCUSSTATS"
echo -e "#### 03f_filterloci.sh: Input FASTA directory: $IN_DIR"
echo -e "#### 03f_filterloci.sh: Output FASTA directory: $OUT_DIR"
echo -e "#### 03f_filterloci.sh: Maximum percentage of missing data: $MAXMISS"
echo -e "#### 03f_filterloci.sh: Minimum distance (bp) between loci: $MINDIST \n\n"

################################################################################
#### FILTER LOCI FOR MISSING DATA AND DISTANCE ####
################################################################################
echo -e "#### 03f_filterloci.sh: Running script to filter for missing data and distance ..."
Rscript $SCRIPTS_DIR/03f_filterloci.R $LOCUSSTATS $IN_DIR $OUT_DIR $MAXMISS$ $MINDIST

## Report:
echo -e "\n#### 03f_filterloci.sh: Done with script."
date
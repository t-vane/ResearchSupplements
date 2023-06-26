#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
PCANGSD=/home/nibtve93/software/pcangsd/pcangsd.py # (v1.01; https://github.com/Rosemeis/pcangsd)

## Command-line args:
NT=$1
BEAGLE=$2
OUT_DIR=$3
SCRIPTS_DIR=$4
IND_FILE=$5
SET_ID=$6

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### pca.sh: Starting script."
echo -e "#### pca.sh: Number of threads: $NT"
echo -e "#### pca.sh: Genotype likelihood file in beagle format: $BEAGLE"
echo -e "#### pca.sh: Output directory: $OUT_DIR"
echo -e "#### pca.sh: Directory with scripts: $SCRIPTS_DIR"
echo -e "#### pca.sh: File that maps individuals to populations: $IND_FILE"
echo -e "#### pca.sh: Set ID: $SET_ID \n\n"

################################################################################
#### CONDUCT PRINCIPAL COMPONENT ANALYSIS ####
################################################################################
echo -e "#### pca.sh: Estimating covariance matrix ...\n"
python $PCANGSD -threads $NT -beagle $BEAGLE -out $OUT_DIR/cov_matrix.txt -tree

echo -e "#### pca.sh: Obtaining principal components and plotting ...\n"
Rscript $SCRIPTS_DIR/pca.R $OUT_DIR $OUT_DIR/cov_matrix.txt $IND_FILE $SET_ID

echo -e "\n#### pca.sh: Done with script."
date


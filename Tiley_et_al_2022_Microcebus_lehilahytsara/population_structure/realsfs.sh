#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
# realSFS needs to be included in $PATH (http://www.popgen.dk/angsd/index.php/RealSFS)

## Command-line args:
NT=$1
SAF_IN=$2
OUT_FILE=$3
OPTIONS=$4

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### realsfs.sh: Starting script."
echo -e "#### realsfs.sh: Number of threads: $NT"
echo -e "#### realsfs.sh: Site allele frequency file(s) for population(s): $SAF_IN"
echo -e "#### realsfs.sh: Output file: $OUT_FILE"
echo -e "#### realsfs.sh: Additional options: $OPTIONS \n\n"

################################################################################
#### ESTIMATE MINOR ALLELE FREQUENCY SPECTRUM ####
################################################################################
echo -e "#### realsfs.sh: Estimating minor allele frequency spectrum ...\n"
realSFS $SAF_IN -P $NT -fold 1 $OPTIONS > $OUT_FILE

## Report:
echo -e "\n#### realsfs.sh: Done with script."
date


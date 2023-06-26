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
SFS=$3
OUT=$4
OPTIONS=$5

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### realsfs_fst.sh: Starting script."
echo -e "#### realsfs_fst.sh: Number of threads: $NT"
echo -e "#### realsfs_fst.sh: Site allele frequency file(s) for population(s): $SAF_IN"
echo -e "#### realsfs_fst.sh: Joint minor allele frequency spectrum: $SFS"
echo -e "#### realsfs_fst.sh: Output prefix: $OUT"
echo -e "#### realsfs_fst.sh: Additional options: $OPTIONS \n\n"

################################################################################
#### ESTIMATE PAIRWISE F_ST BETWEEN POPULATIONS ####
################################################################################
echo -e "#### realsfs_fst.sh: Estimating pairwise F_ST between populations ...\n"
realSFS fst index $SAF_IN -P $NT -sfs $SFS -fstout $OUT
realSFS fst stats $OUT.fst.idx > $OUT.fst

## Report:
echo -e "\n#### realsfs_fst.sh: Done with script."
date


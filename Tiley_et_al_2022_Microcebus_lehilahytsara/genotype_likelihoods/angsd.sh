#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
# angsd needs to be included in $PATH (v0.934; http://www.popgen.dk/angsd/index.php/ANGSD)

## Command-line args:
NT=$1
REFERENCE=$2
BAMLIST=$3
TODO=$4
FILTERS=$5
OUT=$6

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### angsd.sh: Starting script."
echo -e "#### angsd.sh: Number of threads: $NT"
echo -e "#### angsd.sh: Reference genome: $REFERENCE"
echo -e "#### angsd.sh: List of BAM files: $BAMLIST"
echo -e "#### angsd.sh: Inferences to do: $TODO"
echo -e "#### angsd.sh: Filters to apply: $FILTERS"
echo -e "#### angsd.sh: Output prefix: $OUT \n\n"

################################################################################
#### RUN ANGSD WITH SPECIFIED FILTERS AND INFERENCES ####
################################################################################
echo -e "#### angsd.sh: Running angsd ...\n"
angsd -nThreads $NC -ref $REFERENCE -bam $BAMLIST $TODO $FILTERS -out $OUT 

echo -e "\n#### angsd.sh: Done with script."
date


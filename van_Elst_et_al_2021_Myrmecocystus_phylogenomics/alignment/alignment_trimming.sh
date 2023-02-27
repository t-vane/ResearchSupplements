#!/bin/bash

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
MAFFT=/global/homes/jg/t_vane02/software/mafft-7.429/bin/mafft
# trimal needs to be included in $PATH

## Command-line args:
NT=$1
LOCUS_FILE=$2
OUT_DIR=$3
LOCUS=$(sed -n "$SGE_TASK_ID"p $LOCUS_FILE)

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### alignment_trimming.sh: Starting script."
echo -e "#### alignment_trimming.sh: Number of threads: $NT"
echo -e "#### alignment_trimming.sh: File with loci: $LOCUS_FILE"
echo -e "#### alignment_trimming.sh: Output directory: $OUT_DIR \n\n"

################################################################################
#### ALIGNMENT####
################################################################################
echo -e "#### alignment_trimming.sh: Aligning sequences in $LCOUS with MAFFT ... \n"
$MAFFT --thread $NT $LOCUS > $OUT_DIR/$(basename $LOCUS)-mafft

################################################################################
#### TRIMMING####
################################################################################
echo -e "#### alignment_trimming.sh: Trimming alignment with trimAl ... \n"
trimal -in $OUT_DIR/$(basename $LOCUS)-mafft -out $OUT_DIR/$(basename $LOCUS)-mafft-trimal -gappyout

## Report:
echo -e "\n#### alignment_trimming.sh: Done with script."
date

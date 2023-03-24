#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
# NGSadmix needs to be included in $PATH (v32; http://www.popgen.dk/software/index.php/NgsAdmix)

## Command-line args:
NT=$1
K=$2
BEAGLE=$3
OUT_DIR=$4
MININD=$5
SET_ID=$6
SEED=$SLURM_ARRAY_TASK_ID

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### ngsadmix.sh: Starting script."
echo -e "#### ngsadmix.sh: Number of threads: $NT"
echo -e "#### ngsadmix.sh: Number of clusters: $K"
echo -e "#### ngsadmix.sh: Genotype likelihood file in beagle format: $BEAGLE"
echo -e "#### ngsadmix.sh: Output directory: $OUT_DIR"
echo -e "#### ngsadmix.sh: Minimum number of represented individuals: $MININD"
echo -e "#### ngsadmix.sh: Set ID: $SET_ID"
echo -e "#### ngsadmix.sh: Seed (repetition number): $SEED \n\n"

################################################################################
#### INFER INDIVIDUAL ANCESTRIES ####
################################################################################
echo -e "#### ngsadmix.sh: Estimating individual ancestries for $K clusters (repetition number $SEED)...\n"
NGSadmix -P $NT -likes $BEAGLE -seed $SEED  -K $K -outfiles $OUT_DIR/$SET_ID.K$K.seed$SEED -minMaf 0.05 -minInd $MININD -tol 0.000001

echo -e "\n#### ngsadmix.sh: Done with script."
date



#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################

## Command-line args:
SCRIPTS_DIR=$1
SET_ID=$2
IND_FILE=$3
BED_DIR=$4
LOCUSBED_INTERMED=$5
MIN_ELEM_OVL=$6
MIN_ELEM_OVL_TRIM=$7
MIN_LOCUS_SIZE=$8
MAX_DIST_WITHIN_IND=$9
MAX_DIST_BETWEEN_IND=$10
MIN_ELEM_SIZE=$11
LAST_ROW=$12

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### 03a_makelocusbed.sh: Starting script."
echo -e "#### 03a_makelocusbed.sh: Directory with scripts: $SCRIPTS_DIR"
echo -e "#### 03a_makelocusbed.sh: Set ID: $SET_ID"
echo -e "#### 03a_makelocusbed.sh: File with individuals: $IND_FILE"
echo -e "#### 03a_makelocusbed.sh: Directory with BED files: $BED_DIR"
echo -e "#### 03a_makelocusbed.sh: Output BED file: $LOCUSBED_INTERMED"
echo -e "#### 03a_makelocusbed.sh: Minimum element overlap for locus creation: $MIN_ELEM_OVL"
echo -e "#### 03a_makelocusbed.sh: Minimum element overlap for locus trimming: $MIN_ELEM_OVL_TRIM"
echo -e "#### 03a_makelocusbed.sh: Minimum locus size: $MIN_LOCUS_SIZE"
echo -e "#### 03a_makelocusbed.sh: Maximum distance within individuals: $MAX_DIST_WITHIN_IND"
echo -e "#### 03a_makelocusbed.sh: Maximum distance between individuals: $MAX_DIST_BETWEEN_IND"
echo -e "#### 03a_makelocusbed.sh: Minimum locus size: $MIN_ELEM_SIZE"
echo -e "#### 03a_makelocusbed.sh: Number of loci to process (all if 0): $LAST_ROW \n\n"

################################################################################
#### MAKE BED FILE WITH LOCUS COORDINATES ####
################################################################################
## Run GATK CallableLoci to produce BED file for sites that are (non-)callable for a single sample
echo -e "#### 03a_makelocusbed.sh: Running script to create BED file with locus coordinates ..."
Rscript $SCRIPTS_DIR/03a_makelocusbed.R $SET_ID $IND_FILE $BED_DIR $LOCUSBED_INTERMED \
	$MIN_ELEM_OVL $MIN_ELEM_OVL_TRIM $MIN_LOCUS_SIZE $MAX_DIST_WITHIN_IND $MAX_DIST_BETWEEN_IND $MIN_ELEM_SIZE $LAST_ROW

## Report:
echo -e "\n#### 03a_makelocusbed.sh: Done with script."
date
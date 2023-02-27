#!/bin/bash

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
# phyluce scripts (https://phyluce.readthedocs.io/en/latest/) need to be included in $PATH

## Command-line args:
CONTIG_DIR=$1
PROBES=$2
OUT_DIR=$3

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### UCE_match_contigs.sh: Starting script."
echo -e "#### UCE_match_contigs.sh: Contig directory: $CONTIG_DIR"
echo -e "#### UCE_match_contigs.sh: Probes: $PROBES"
echo -e "#### UCE_match_contigs.sh: Output directory: $OUT_DIR \n\n"

################################################################################
#### MATCH CONTIGS TO PROBES####
################################################################################
echo -e "#### UCE_match_contigs.sh: Matching contigs to probes ... \n"
phyluce_assembly_match_contigs_to_probes --contigs $CONTIG_DIR --probes $PROBES --output $OUT_DIR

## Report:
echo -e "\n#### UCE_match_contigs.sh: Done with script."
date

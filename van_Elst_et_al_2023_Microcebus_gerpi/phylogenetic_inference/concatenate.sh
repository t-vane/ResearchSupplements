#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
AMAS=/home/nibtve93/software/AMAS/amas/AMAS.py # (https://github.com/marekborowiec/AMAS)

## Command-line args:
NT=$1
ALIGNMENT_DIR=$2
SET_ID=$3

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### concatenate.sh: Starting script."
echo -e "#### concatenate.sh: Number of threads: $NT"
echo -e "#### concatenate.sh: Alignment directory: $ALIGNMENT_DIR"
echo -e "#### concatenate.sh: Set ID: $SET_ID \n\n"

################################################################################
#### CONCATENATE LOCUS ALIGNMENTS ####
################################################################################
echo -e "#### concatenate.sh: Concatenating locus alignments ..."
python $AMAS concat -i $ALIGNMENT_DIR/*.muscle.fa -t $ALIGNMENT_DIR/$SET_ID.concatenated.nex -p $ALIGNMENT_DIR/$SET_ID.partitions.txt -f fasta -d dna -u nexus -c $NT

echo -e "#### concatenate.sh: Replacing 'N' by '?' ..."
head -n 6 $ALIGNMENT_DIR/$SET_ID.concatenated.nex
awk 'NR>=7 {gsub("N","?",$2)}1' $ALIGNMENT_DIR/$SET_ID.concatenated.nex > $ALIGNMENT_DIR/$SET_ID.concatenated.nex.tmp
mv $ALIGNMENT_DIR/$SET_ID.concatenated.nex.tmp $ALIGNMENT_DIR/$SET_ID.concatenated.nex

echo -e "#### concatenate.sh: Converting NEXUS to PHYLIP format ..."
python $AMAS convert -i $ALIGNMENT_DIR/$SET_ID.concatenated.nex -f nexus -u phylip -d dna
mv $ALIGNMENT_DIR/$SET_ID.concatenated.nex-out.phy $ALIGNMENT_DIR/$SET_ID.concatenated.phy

## Report:
echo -e "#### concatenate.sh: Done with script."
date




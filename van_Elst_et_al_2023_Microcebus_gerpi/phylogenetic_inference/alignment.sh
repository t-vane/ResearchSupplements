#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
# muscle needs to be included in $PATH (v3.8.31; https://www.drive5.com/muscle/)

## Command-line args:
NT=$1
LOCUS_DIR=$2
ALIGNMENT_DIR=$3

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### alignment.sh: Starting script."
echo -e "#### alignment.sh: Number of threads: $NT"
echo -e "#### alignment.sh: Directory with locus fasta files: $LOCUS_DIR"
echo -e "#### alignment.sh: Output directory for alignments: $ALIGNMENT_DIR \n\n"

################################################################################
#### ALIGN LOCI WITH MUSCLE AND CHANGE HEADERS ####
################################################################################
echo -e "#### alignment.sh: Aligning loci with muscle and changing headers ..."
for LOCUS in $LOCUS_DIR/*.fa
do 
	echo "muscle -in $LOCUS -out $ALIGNMENT_DIR/$(basename $LOCUS .fa).muscle.fa; sed -i 's/__.*//g' $ALIGNMENT_DIR/$(basename $LOCUS .fa).muscle.fa" >> $ALIGNMENT_DIR/alignment_commands.txt
done
parallel --jobs $NT < $ALIGNMENT_DIR/alignment_commands.txt

## Report:
echo -e "#### alignment.sh: Done with script."
date





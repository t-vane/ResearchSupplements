#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
AMAS=/home/nibtve93/software/AMAS/amas/AMAS.py # (https://github.com/marekborowiec/AMAS)

## Command-line args:
ALIGNMENT_DIR=$1
OUT_FILE=$2

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### locus_statistics.sh: Starting script."
echo -e "#### locus_statistics.sh: Alignment directory: $ALIGNMENT_DIR"
echo -e "#### locus_statistics.sh: Output file: $OUT_FILE \n\n"

################################################################################
#### CALCULATE STATISTICS FOR LOCUS ALIGNMENTS ####
################################################################################
echo -e "#### locus_statistics.sh: Calculating statistics for locus alignments ..."
for FASTA in $ALIGNMENT_DIR/*.muscle.fa
do
    echo -e "#### locus_statistics.sh: Locus $FASTA ..."
    python $AMAS summary -f fasta -d dna -i $FASTA -o $(dirname $OUT_FILE)/$(basename "$FASTA" .muscle.fa).stats.tmp
    grep -v "Alignment_name" $(dirname $OUT_FILE)/$(basename "$FASTA" .muscle.fa).stats.tmp >> $OUT_FILE
done

echo -e "#### locus_statistics.sh: Preparing final statistics file ..."
HEADER=$(head -n 1 `dirname $OUT_FILE`/`basename "$FASTA" .muscle.fa`.stats.tmp)
(echo $HEADER && cat $OUT_FILE) > $(dirname $OUT_FILE)/final.tmp
mv $(dirname $OUT_FILE)/final.tmp $OUT_FILE

echo -e "#### locus_statistics.sh: Removing temporary files ..."
find $(dirname $OUT_FILE) -name '.tmp' -delete

## Report:
echo -e "#### locus_statistics.sh: Done with script."
date






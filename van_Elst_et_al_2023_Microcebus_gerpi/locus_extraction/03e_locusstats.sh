#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Script adapted and modified from Poelstra et al. 2021, Systematic Biology (https://doi.org/10.1093/sysbio/syaa053)

## Software:
AMAS=/home/nibtve93/software/AMAS/amas/AMAS.py # (https://github.com/marekborowiec/AMAS)

## Command-line args:
FASTA_DIR=$1
LOCUSSTATS=$2

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### 03e_locusstats.sh: Starting script."
echo -e "#### 03e_locusstats.sh: Input directory with FASTA files: $FASTA_DIR"
echo -e "#### 03e_locusstats.sh: Output file with statistics: $LOCUSTATS \n\n"

################################################################################
#### ESTIMATE STATISTICS FOR ALL LOCI  ####
################################################################################
echo -e "\n#### 03e_locusstats.sh: Estimate statistics for each locus ..."
for FASTA in $FASTA_DIR/*
do
    echo -e "\n## 03e_locusstats.sh: Processing $FASTA ..."
    FASTA_ID=$(basename $FASTA)
    FILE_STATS=$(dirname $LOCUSSTATS)/tmp.$FASTA_ID.stats.txt

    python3 $AMAS summary -f fasta -d dna -i $FASTA -o $FILE_STATS

    grep -v "Alignment_name" $FILE_STATS >> $LOCUSSTATS
done

## Add header
echo -e "\n#### 03e_locusstats.sh: Adding header line ..."
HEADER=$(head -n 1 $FILE_STATS)
(echo $HEADER && cat $LOCUSSTATS) > $(dirname $LOCUSSTATS)/tmp.txt && mv $(dirname $LOCUSSTATS)/tmp.txt $LOCUSSTATS

## Remove temporary files
echo -e "\n#### 03e_locusstats.sh: Removing temporary files ..."
find $(dirname $LOCUSSTATS) -name 'tmp*txt' -delete

## Report:
echo -e "\n#### 03e_locusstats.sh: Done with script."
date



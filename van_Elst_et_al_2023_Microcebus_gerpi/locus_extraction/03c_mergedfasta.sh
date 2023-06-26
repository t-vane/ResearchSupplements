#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Script adapted and modified from Poelstra et al. 2021, Systematic Biology (https://doi.org/10.1093/sysbio/syaa053)

## Software:
# BEDtools needs to be included in $PATH (v2.30.0; https://bedtools.readthedocs.io/en/latest/)
# SAMtools needs to be included in $PATH (v1.11; http://www.htslib.org/)

## Command-line args:
IND_FILE=$1
INDFASTA_DIR=$2
LOCUSBED_FINAL=$$3
LOCUSLIST=$4
FASTA_MERGED=$5

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### 03c_mergedfasta.sh: Starting script."
echo -e "#### 03c_mergedfasta.sh: File with individuals: $IND_FILE"
echo -e "#### 03c_mergedfasta.sh: Directory for FASTA files: $INDFASTA_DIR"
echo -e "#### 03c_mergedfasta.sh: Final locus BED file: $LOCUSBED_FINAL"
echo -e "#### 03c_mergedfasta.sh: List of loci to be created: $LOCUSLIST"
echo -e "#### 03c_mergedfasta.sh: Merged output FASTA file: $FASTA_MERGED \n\n"

################################################################################
#### CREATE MERGED FASTA FILE WITH ALL INDIVIDUALS AND LOCI ####
################################################################################
## Extract loci in locus BED file from masked FASTA for each individual
echo -e "#### 03c_mergedfasta.sh: Extracting loci from locus BED file from masked FASTA for each individual ..."
while read -r INDV
do
	echo -e "## 03c_mergedfasta.sh: Processing $INDV ..."
    FASTA_IN=$INDFASTA_DIR/${INDV}.altrefmasked.fasta
    FASTA_OUT=$INDFASTA_DIR/${INDV}_allloci.fasta
    bedtools getfasta -fi $FASTA_IN -bed $LOCUSBED_FINAL > $FASTA_OUT
done < $IND_FILE

## Create list of loci
echo -e "#### 03c_mergedfasta.sh: Creating list of loci ..."
FASTA_1=$(find $INDFASTA_DIR/*allloci.fasta | head -1)
grep ">" $FASTA_1 | sed 's/>//' > $LOCUSLIST

## Merge by-individual FASTA files
echo -e "#### 03c_mergedfasta.sh: Merging loci of individuals ..."
while read -r INDV
do
	echo -e "## 03c_mergedfasta.sh: Processing $INDV ..."
    FASTA=$INDFASTA_DIR/${INDV}_allloci.fasta
    # Replace ":" by "," for compatibility with SAMtools faidx:
    sed "s/>/>${ind}__/g" $FASTA | sed 's/:/,/g' >> $FASTA_MERGED
done < $IND_FILE

## Index merged FASTA file with SAMtools faidx
echo -e "#### 03c_mergedfasta.sh: Indexing merged FASTA file with SAMtools faidx ..."
samtools faidx $FASTA_MERGED

## Report:
echo -e "\n#### 03c_mergedfasta.sh: Done with script."
date


#!/bin/bash

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
# illumiprocessor needs to be included in $PATH
# phyluce scripts (https://phyluce.readthedocs.io/en/latest/) need to be included in $PATH

## Command-line args:
NT=$1
IN_DIR=$2
OUT_DIR=$3
CONFIG=$4
R1_PATTERN=$5
R2_PATTERN=$6

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### clean_reads.sh: Starting script."
echo -e "#### clean_reads.sh: Number of threads: $NT"
echo -e "#### clean_reads.sh: Input directory: $IN_DIR"
echo -e "#### clean_reads.sh: Output directory: $OUT_DIR"
echo -e "#### clean_reads.sh: Configuration file: $CONFIG"
echo -e "#### clean_reads.sh: Forward read pattern: $R1_PATTERN"
echo -e "#### clean_reads.sh: Reverse read pattern: $R2_PATTERN \n\n"

################################################################################
#### CLEAN READS AND CHECK FOR INTEGRITY OF OUTPUT####
################################################################################
mkdir -p $OUT_DIR

echo -e "#### clean_reads.sh: Cleaning reads in $IN_DIR ... \n"
illumiprocessor --cores $NT --input $IN_DIR --output $OUT_DIR --config $IN_DIR/$CONFIG --r1-pattern $R1_PATTERN --r2-pattern $R2_PATTERN

echo -e "#### clean_reads.sh: Checking for integrity of output... \n"
gunzip -t -r $OUT_DIR

echo -e "#### clean_reads.sh: Calculate clean read statistics... \n"
for i in $OUT_DIR/*
do 
	phyluce_assembly_get_fastq_lengths --input $i/split-adapter-quality-trimmed/ --csv
done > $OUT_DIR/clean_read_statistics.txt

## Report:
echo -e "\n#### clean_reads.sh: Done with script."
date

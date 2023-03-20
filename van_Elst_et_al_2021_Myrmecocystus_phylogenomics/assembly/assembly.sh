#!/bin/bash

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
METASPADES=/global/homes/jg/t_vane02/software/SPAdes-3.13.1-Linux/bin/metaspades.py # (v3.13.1; https://github.com/ablab/spades)
# seqtk needs to be included in $PATH (https://github.com/lh3/seqtk)

## Command-line args:
NT=$1
MEM=$2
INDS=$3
CLEAN_READS=$4
OUT_DIR=$5
DOWNSAMPLE=$6
READ_NUMBER=$7
IND=$(sed -n "$SGE_TASK_ID"p $INDS)

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### assembly.sh: Starting script."
echo -e "#### assembly.sh: Number of threads: $NT"
echo -e "#### assembly.sh: Maximum memory in GB: $MEM"
echo -e "#### assembly.sh: File with individuals: $INDS"
echo -e "#### assembly.sh: Directory with clean reads: $CLEAN_READS"
echo -e "#### assembly.sh: Output directory: $OUT_DIR"
echo -e "#### assembly.sh: Downsampling: $DOWNSAMPLE"
echo -e "#### assembly.sh: Number of reads to retain if downsampling: $READ_NUMBER \n\n"

################################################################################
#### UNZIP CLEAN READ FILES####
################################################################################
echo -e "#### assembly.sh: Unzipping clean read files for individual $IND \n"
gunzip $CLEAN_READS/$IND/split-adapter-quality-trimmed/$IND-*.fastq.gz

################################################################################
#### DOWNSAMPLE####
################################################################################
if [ $DOWNSAMPLE == TRUE ]
then
	echo -e "#### assembly.sh: Downsampling reads for individual $IND to $READ_NUMBER ...\n"
	RNUM=$RANDOM
	READS=$READ_NUMBER
	for i in 1 2
	do
		cp $CLEAN_READS/$IND/split-adapter-quality-trimmed/$IND-READ$i.fastq $CLEAN_READS/$IND/split-adapter-quality-trimmed/$IND-READ$i.full.fastq
		seqtk sample -s $RNUM $CLEAN_READS/$IND/split-adapter-quality-trimmed/$IND-READ$i.fastq $READS > $CLEAN_READS/$IND/split-adapter-quality-trimmed/$IND-READ$i.tmp.fastq
		mv $CLEAN_READS/$IND/split-adapter-quality-trimmed/$IND-READ$i.tmp.fastq $CLEAN_READS/$IND/split-adapter-quality-trimmed/$IND-READ$i.fastq
	done
fi

################################################################################
#### RUN ASSEMBLY FOR SPECIFIC SAMPLE####
################################################################################
echo -e "#### assembly.sh: Running assembly for sample $IND \n"
python3 $METASPADES -t $NT -m $MEM -1 $CLEAN_READS/$IND/split-adapter-quality-trimmed/$IND-READ1.fastq -2 $CLEAN_READS/$IND/split-adapter-quality-trimmed/$IND-READ2.fastq -s $CLEAN_READS/$IND/split-adapter-quality-trimmed/$IND-singleton.fastq -o $OUT_DIR 

## Report:
echo -e "\n#### assembly.sh: Done with script."
date

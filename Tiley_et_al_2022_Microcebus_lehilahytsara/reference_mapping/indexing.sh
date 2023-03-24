#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
# BWA needs to be included in $PATH (v7.0.17; https://github.com/lh3/bwa)
# SAMtools needs to be included in $PATH (v1.11; http://www.htslib.org/)
# gatk needs to be included in $PATH (v4.1.9.0; https://gatk.broadinstitute.org/hc/en-us)

## Command-line args:
REFERENCE=$1
BWA=$2
BWA_INDEX=$3
SAMTOOLS=$4
GATK=$5

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### indexing.sh: Starting script."
echo -e "#### indexing.sh: Reference genome: $REFERENCE"
echo -e "#### indexing.sh: Index with BWA: $NT"
[[ $BWA ]] && echo -e "#### indexing.sh: Name for BWA index file: $BWA_INDEX"
echo -e "#### indexing.sh: Index with SAMtools: $SAMTOOLS"
echo -e "#### indexing.sh: Index with GATK (Picard): $GATK \n\n"

################################################################################
#### CREATE INDICES ####
################################################################################
cd $(dirname $REFERENCE)

## Index with BWA
[[ $BWA ]] && echo -e "#### indexing.sh: Creating index of $REFERENCE with BWA ...\n" && bwa index -p $BWA_INDEX -a bwtsw $REFERENCE

## Index with SAMtools
[[ $SAMTOOLS ]] && echo -e "#### indexing.sh: Creating index of $REFERENCE with SAMtools ...\n" && samtools faidx $REFERENCE

## Index with GATK (Picard)
if [[ $GATK ]]
then
	echo -e "#### indexing.sh: Creating index of $REFERENCE with GATK (Picard) ...\n"
	gatk CreateSequenceDictionary -R $REFERENCE -O $REFERENCE.dict
	OLDNAME=${REFERENCE}.dict
	NEWNAME="${OLDNAME/fna.dict/dict}"
	cp $OLDNAME $NEWNAME
fi

## Report:
echo -e "\n#### indexing.sh: Done with script."
date

#!/bin/sh
#SBATCH -p medium40

set -euo pipefail

source /home/nibtve93/.bashrc

################################################################################
#### SET-UP ####
################################################################################
## Software:
# BWA needs to be included in $PATH
# SAMtools needs to be included in $PATH

## Command-line args:
MODE=$1
NT=$2
INDEX=$3
IN_DIR=$4
OUT_DIR=$5
IN_FILE=$6
QFILTER=$7
MINMAPQ=$8
DEDUPLICATE=$8
EXTRACT_BED=$9
BED=$10
SUFFIX=$11
INDV=$(sed -n "$SLURM_ARRAY_TASK_ID"p $IN_FILE)
READGROUP="@RG\tID:group1\tSM:$INDV\tPL:illumina\tLB:lib1"

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### reference_mapping.sh: Starting script."
echo -e "#### reference_mapping.sh: Paired-end or single-end: $MODE"
echo -e "#### reference_mapping.sh: Number of threads: $NT"
echo -e "#### reference_mapping.sh: Reference genome index: $INDEX"
echo -e "#### reference_mapping.sh: Directory with trimmed reads: $IN_DIR"
echo -e "#### reference_mapping.sh: Directory for BAM files: $OUT_DIR"
echo -e "#### reference_mapping.sh: File with individuals: $IN_FILE"
echo -e "#### reference_mapping.sh: Filter based on minimum mapping quality: $QFILTER"
if [[ $QFILTER ]]
then
	echo -e "#### reference_mapping.sh: Minimum mapping quality for filtering: $MINMAPQ"
fi
echo -e "#### reference_mapping.sh: Deduplicate: $DEDUPLICATE"
echo -e "#### reference_mapping.sh: Extract specific genomic region: $EXTRACT_BED"
if [[ $EXTRACT_BED ]]
then
	echo -e "#### reference_mapping.sh: BED file with information on genomic region: $BED"
fi
echo -e "#### reference_mapping.sh: Suffix for final BAM files: $SUFFIX"
echo -e "#### reference_mapping.sh: File with individuals: $IN_FILE"
echo -e "#### reference_mapping.sh: Individual: $INDV \n\n"

################################################################################
#### MAP TO REFERENCE GENOME ####
################################################################################
if [[ $MODE == "PE" ]]
then
	echo -e "#### reference_mapping.sh: Reference mapping for paired-end individual $INDV ...\n"
	bwa mem -aM -R $READGROUP -t $NT $INDEX $IN_DIR/$INDV.trimmed.1.fq.gz $IN_DIR/$INDV.trimmed.2.fq.gz | samtools view -b -h > $OUT_DIR/$INDV.bam
elif [[ $MODE == "SE" ]]
then
	echo -e "#### reference_mapping.sh: Reference mapping for single-end individual $INDV ...\n"
	bwa mem -aM -R $READGROUP -t $NT $INDEX $IN_DIR/$INDV.trimmed.1.fq.gz | samtools view -b -h > $OUT_DIR/$INDV.bam
else
	echo -e "#### reference_mapping.sh: Invalid sequencing mode provided - only PE and SE allowed. ...\n"
	exit 1
fi

################################################################################
#### SORT AND FILTER FOR MINIMUM MAPPING QUALITY AND PROPER PAIRING (IF PE) ####
################################################################################
if [[ $QFILTER ]]
then
	if [[ $MODE == "PE" ]]
	then
		echo -e "#### reference_mapping.sh: Minimum mapping quality and proper-pair filtering and sorting for paired-end individual $INDV ...\n"
		samtools view -bhu -q $MINMAPQ -f 0x2 -@ $NT $OUT_DIR/$INDV.bam | samtools sort -@ $NT -m 15G -O bam > $OUT_DIR/$INDV.MQ$MINMAPQ.pp.bam
	else [[ $MODE == "SE" ]]
		echo -e "#### reference_mapping.sh: Minimum mapping quality filtering and sorting for single-end individual $INDV ...\n"
		samtools view -bhu -q $MINMAPQ -@ $NT $OUT_DIR/$INDV.bam | samtools sort -@ $NT -m 15G -O bam > $OUT_DIR/$INDV.MQ$MINMAPQ.bam
	fi
fi

################################################################################
#### DEDUPLICATE ####
################################################################################
if [[ $DEDUPLICATE ]]
then
	if [[ $MODE == "PE" ]]
	then
		echo -e "#### reference_mapping.sh: Deduplication for paired-end individual $INDV ...\n"
		samtools collate --output-fmt BAM $OUT_DIR/$INDV.MQ$MINMAPQ.pp.bam -O | samtools fixmate - - -r -m -O BAM| samtools sort -m 15G -O BAM | samtools markdup - $OUT_DIR/$INDV.MQ$MINMAPQ.pp.dedup.bam -r -s -O BAM
	fi
fi

################################################################################
#### EXTRACT SPECIFIC REGIONS FROM BED FILE ####
################################################################################
if [[ $EXTRACT_BED ]]
then
	if [[ $MODE == "PE" ]]
	then
		echo -e "#### reference_mapping.sh: Extraction of genomic regions provided in $BED for paired-end individual $INDV ...\n"
		samtools view -b -L $BED $OUT_DIR/$INDV.MQ$MINMAPQ.pp.dedup.bam > $OUT_DIR/$INDV.MQ$MINMAPQ.pp.dedup.$SUFFIX.bam
	else
		echo -e "#### reference_mapping.sh: Extraction of genomic regions provided in $BED for single-end individual $INDV ...\n"
		samtools view -b -L $BED $OUT_DIR/$INDV.MQ$MINMAPQ.bam > $OUT_DIR/$INDV.MQ$MINMAPQ.$SUFFIX.bam
	fi
fi

## Report:
echo -e "\n#### reference_mapping.sh: Done with script."
date






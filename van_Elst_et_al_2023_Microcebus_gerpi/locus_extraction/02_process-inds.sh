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
GATK3=/home/nibtve93/software/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar (v3.8.1.0; https://gatk.broadinstitute.org/hc/en-us)

## Command-line args:
IND_FILE=$1
VCF_ALTREF=$2
REFERENCE=$3
BAM_DIR=$4
SUFFIX=$5
MIN_DP=$6
INDFASTA_DIR=$7
BED_DIR=$8
BED_REMOVED_SITES=$9
INDV=$(sed -n "$SLURM_ARRAY_TASK_ID"p $IND_FILE)

BAM=$BAM_DIR/$INDV.*$SUFFIX.bam

CALLABLE_SUMMARY=$BED_DIR/$INDV.callableloci.sumtable.txt
BED_OUT=$BED_DIR/$INDV.callablelocioutput.bed
BED_NONCALLABLE=$BED_DIR/$INDV.noncallable.bed
BED_CALLABLE=$BED_DIR/$INDV.callable.bed

FASTA_ALTREF=$INDFASTA_DIR/$INDV.altref.fasta
FASTA_MASKED=$INDFASTA_DIR/$INDV.altrefmasked.fasta

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### 02_process-inds.sh: Starting script."
echo -e "#### 02_process-inds.sh: File with individuals: $IND_FILE"
echo -e "#### 02_process-inds.sh: Raw VCF: $VCF_ALTREF"
echo -e "#### 02_process-inds.sh: Reference genome: $REFERENCE"
echo -e "#### 02_process-inds.sh: Directory with BAM files: $BAM_DIR"
echo -e "#### 02_process-inds.sh: BAM file suffix: $SUFFIX"
echo -e "#### 02_process-inds.sh: Minimum depth for CallableLoci: $MIN_DP"
echo -e "#### 02_process-inds.sh: Directory for FASTA files: $INDFASTA_DIR"
echo -e "#### 02_process-inds.sh: Directory with BED file: $BED_DIR"
echo -e "#### 02_process-inds.sh: BED file with removed sites: $BED_REMOVED_SITES"
echo -e "#### 02_process-inds.sh: Individual: $INDV"
echo -e "#### 02_process-inds.sh: BAM file: $BAM"
echo -e "#### 02_process-inds.sh: CallableLoci summary: $CALLABLE_SUMMARY"
echo -e "#### 02_process-inds.sh: CallableLoci output BED all: $BED_OUT"
echo -e "#### 02_process-inds.sh: CallableLoci output BED non-callable: $BED_NONCALLABLE"
echo -e "#### 02_process-inds.sh: CallableLoci output BED callable: $BED_CALLABLE"
echo -e "#### 02_process-inds.sh: Raw FASTA genome: $FASTA_ALTREF"
echo -e "#### 02_process-inds.sh: Masked FASTA genome: $FASTA_MASKED \n\n"

################################################################################
#### CREATE MASKED FASTA GENOME FOR INDIVIDUAL ####
################################################################################
## Index BAM file if not yet done
[[ ! -e $BAM.bai ]] && echo -e "#### 02_process-inds.sh: Indexing bam ..." && samtools index $BAM

## Run GATK CallableLoci to produce BED file for sites that are (non-)callable for a single sample
echo -e "#### 02_process-inds.sh: Running CallableLoci ..."
java -jar $GATK3 -T CallableLoci -R $REFERENCE -I $BAM -summary $CALLABLE_SUMMARY --minDepth $MIN_DP -o $BED_OUT
# Separate into two files
grep -v "CALLABLE" $BED_OUT > $BED_NONCALLABLE
grep "CALLABLE" $BED_OUT > $BED_CALLABLE

## Produce whole-genome FASTA file for a single sample
echo -e "#### 02_process-inds.sh: Running FastaAlternateReferenceMaker ..."
java -jar $GATK3 -T FastaAlternateReferenceMaker -IUPAC $INDV -R $REFERENCE -V $VCF_ALTREF -o $FASTA_ALTREF
# Edit FASTA headers
sed -i -e 's/:1//g' -e 's/>[0-9]* />/g' $FASTA_ALTREF
## Count bases
N_ACGT=$(grep -Eo "A|C|G|T" $FASTA_ALTREF | wc -l)
N_AMBIG=$(grep -Eo "M|R|W|S|Y|K" $FASTA_ALTREF | wc -l)
echo -e "#### 02_process-inds.sh: Number of A/C/G/T in FASTA_ALTREF: $N_ACGT"
echo -e "#### 02_process-inds.sh: Number of heterozygous sites (counted by ambiguity code) in FASTA_ALTREF: $N_AMBIG"

## Mask sites identified as non-callable or removed during VCF filtering in whole-genome FASTA
echo -e "#### 02_process-inds.sh: Running BEDtools maskfasta ..."
# Mask non-callable sites
echo -e "## 02_process-inds.sh: Non-callable sites ..."
bedtools maskfasta -fi $FASTA_ALTREF -bed $BED_NONCALLABLE -fo $FASTA_MASKED.intermed.fasta
# Mask removed sites
echo -e "## 02_process-inds.sh: Callable sites ..."
bedtools maskfasta -fi $FASTA_MASKED.intermed.fasta -bed $BED_REMOVED_SITES -fo $FASTA_MASKED
# Counting Ns in FASTA files:
NCOUNT_FASTA_ALTREF=$(grep -Fo "N" $FASTA_ALTREF | wc -l)
NCOUNT_FASTA_MASKED_INTERMED=$(grep -Fo "N" $FASTA_MASKED.intermed.fasta | wc -l)
NCOUNT_FASTA_MASKED=$(grep -Fo "N" $FASTA_MASKED | wc -l)
echo -e "#### 02_process-inds.sh: Number of Ns in FASTA_ALTREF: $NCOUNT_FASTA_ALTREF"
echo -e "#### 02_process-inds.sh: Number of Ns in FASTA_MASKED_INTERMED: $NCOUNT_FASTA_MASKED_INTERMED"
echo -e "#### 02_process-inds.sh: Number of Ns in FASTA_MASKED: $NCOUNT_FASTA_MASKED"

## Report:
echo -e "\n#### 02_process-inds.sh: Done with script."
date

#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
# VCFtools needs to be included in $PATH (v0.1.17; https://vcftools.github.io/index.html)
AMAS=/home/nibtve93/software/AMAS/amas/AMAS.py # (https://github.com/marekborowiec/AMAS)
ASCBIAS=/home/nibtve93/software/raxml_ascbias/ascbias.py # (https://github.com/btmartin721/raxml_ascbias/blob/master/ascbias.py)

## Command-line args:
SCRIPTS_DIR=$1
VCF_IN=$2
OUT_DIR=$3
FORMAT=$4

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### vcf_convert.sh: Starting script."
echo -e "#### vcf_convert.sh: Directory with vcf_tab_to_fasta_alignment script: $SCRIPTS_DIR"
echo -e "#### vcf_convert.sh: Input VCF file: $VCF_IN"
echo -e "#### vcf_convert.sh: Output directory: OUT_DIR"
echo -e "#### vcf_convert.sh: Output format of alignment: $FORMAT \n\n"

################################################################################
#### CONVERT VCF TO PHYLIP FORMAT AND CALCULATE STATISTICS ####
################################################################################
PREFIX=$(basename $IN_FILE .vcf)

if [[ $FORMAT == phylip ]]
then
	SUFFIX=phy
elif [[ $FORMAT == nexus ]]
then
	SUFFIX=nex
else
	echo -e "#### vcf_convert.sh: Invalid format provided - only phylip and nexus accepted. ...\n" && exit 1
fi

echo -e "#### vcf_convert.sh: Creating TAB file from VCF ...\n"
vcf-to-tab < $VCF_IN > $OUT_DIR/$PREFIX.tab

echo -e "#### vcf_convert.sh: Converting TAB file to FASTA format ...\n"
perl $SCRIPTS_DIR/vcf_tab_to_fasta_alignment_TVE.pl -i $OUT_DIR/$PREFIX.tab > $OUT_DIR/$PREFIX.fasta

echo -e "#### vcf_convert.sh: Converting FASTA file to $FORMAT format ...\n"
python $AMAS convert -i $OUT_DIR/$PREFIX.fasta -f fasta -u $FORMAT -d dna
mv $OUT_DIR/$PREFIX.fasta-out.$SUFFIX $OUT_DIR/$PREFIX.$SUFFIX

echo -e "#### vcf_convert.sh: Removing invariant sites ...\n"
python $ASCBIAS -p $OUT_DIR/$PREFIX.$SUFFIX -o $OUT_DIR/$PREFIX.noinv.$SUFFIX

echo -e "#### vcf_convert.sh: Calculating alignment statistics ...\n"
python $AMAS summary -f $FORMAT -d dna -i $OUT_DIR/$PREFIX.noinv.$SUFFIX -o $OUT_DIR/$PREFIX.noinv.$SUFFIX.summary

echo -e "#### vcf_convert.sh: Writing partitions file ...\n"
ALIGN_LENGTH=$(sed -n 2p $OUT_DIR/$PREFIX.noinv.$SUFFIX.summary | cut -f3)
echo "[asc~$OUT_DIR/$PREFIX.noinv.$SUFFIX.stamatakis], ASC_DNA, p1=1-$ALIGN_LENGTH" > $OUT_DIR/$PREFIX.noinv.$SUFFIX.partitions

## Report:
echo -e "\n#### vcf_convert.sh: Done with script."
date



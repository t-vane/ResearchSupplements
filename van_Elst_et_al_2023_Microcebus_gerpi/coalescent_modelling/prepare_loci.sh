#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
# faidx of SAMtools needs to be included in $PATH (v1.11; http://www.htslib.org/)

## Command-line args:
IN_DIR=$1
OUT_DIR=$2

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### prepare_loci.sh: Starting script."
echo -e "#### prepare_loci.sh: Input directory with extracted loci: $IN_DIR"
echo -e "#### prepare_loci.sh: Output directory for reformatted loci: $OUT_DIR \n\n"

################################################################################
#### REFORMAT LOCUS FILES FOR GPHOCS ####
################################################################################
echo -e "#### prepare_loci.sh: Reformatting locus files for G-PhoCS ..."
for FASTA in $IN_DIR/*fa
do
	FASTA_ID=$(basename $FASTA)	
	echo -e "#### prepare_loci.sh: Processing locus $FASTA_ID ..."

	# Remove locus name from header
	sed 's/__.*//' $IN_DIR/$FASTA_ID | sed -e 's/.*_\(.*_A[01]\)_.*/>\1/' > $OUT_DIR/$FASTA_ID.tmp1
	# Reformat to have one sample per line; then remove second and third column
	faidx --transform transposed $OUT_DIR/$FASTA_ID.tmp1 | tee $OUT_DIR/$FASTA_ID.tmp2 | cut -f 1,4 > $OUT_DIR/$FASTA_ID.tmp3

	# Get number of bases
	SEQ_LEN=$(cut -f 3 $OUT_DIR/$FASTA_ID.tmp2 | head -n 1)
	# Get number of individuals
	NO_SEQS=$(cat $OUT_DIR/$FASTA_ID.tmp2 | wc -l)

	if [[ NRSEQS != 0 ]] && [[ SEQLEN != 0 ]]
	then
		echo -e "## prepare_loci.sh: Number of individuals: $NO_SEQS ..."
		echo -e "## prepare_loci.sh: Sequence length: $SEQ_LEN ..."
		
		# Print output file
		echo "$FASTA_ID $NRSEQS $SEQLEN" > $OUT_DIR/$FASTA_ID.locus
		cat $OUT_DIR/$FASTA_ID.tmp3 >> $OUT_DIR/$FASTA_ID.locus
	else
		echo -e "## prepare_loci.sh: Skipping empty locus ..."
	fi
done

echo -e "#### prepare_loci.sh: Removing temporary files ..."
rm -f $OUT_DIR/*tmp*
rm -f $OUT_DIR/*fai

echo -e "#### prepare_loci.sh: Counting loci ..."
NLOCI=$(ls $OUT_DIR/*locus | wc -l)
printf "${NLOCI}\n\n" > $OUT_DIR/loci.count
echo -e "## prepare_loci.sh: Number of loci: $NLOCI ..."

## Report:
echo -e "\n#### prepare_loci.sh: Done with script."
date



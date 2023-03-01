################################################################################
#### READ TRIMMING ####
################################################################################
$SCRIPTS_DIR=/home/nibtve93/scripts/readTrimming

NT=6

IN_DIR=$PWORK/rawReads/lehilahytsara
OUT_DIR=$PWORK/trimmedReads/lehilahytsara

mkdir $OUT_DIR/logFiles

## Submit read trimming script for paired- and single-end samples (PE and SE, respectively)
for i in PE SE
do
	IN_FILE=$PWORK/rawReads/lehilahytsara/trim_$i.txt # List of samples for which reads shall be trimmed
	NO_INDS=$(cat $IN_FILE | wc -l)
	ADAPTERS=$PWORK/trimmomaticAdapters_$i.fa # FASTA file with adapter sequences
	
	sbatch --array=1-$NO_INDS --output=$OUT_DIR/logFiles/read_trimming_$i.%A_%a.oe $SCRIPTS_DIR/read_trimming.sh $i $IN_FILE $NT $IN_DIR $OUT_DIR $ADAPTERS
done


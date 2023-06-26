################################################################################
#### READ TRIMMING ####
################################################################################
$SCRIPTS_DIR=/home/nibtve93/scripts/readTrimming

NT=6

SET_ID=gerpi
IN_DIR=$PWORK/rawReads/$SET_ID
OUT_DIR=$PWORK/trimmedReads/$SET_ID
mkdir $OUT_DIR/logFiles

## Submit read trimming script for paired- and single-end samples (PE and SE, respectively)
for i in PE SE
do
	IN_FILE=$IN_DIR/trim_$i.txt # List of samples for which reads shall be trimmed
	NO_INDS=$(cat $IN_FILE | wc -l)
	ADAPTERS=$PWORK/trimmomaticAdapters_$i.fa # FASTA file with adapter sequences
	
	sbatch --array=1-$NO_INDS --output=$OUT_DIR/logFiles/read_trimming_$i.%A_%a.$SET_ID.oe $SCRIPTS_DIR/read_trimming.sh $i $IN_FILE $NT $IN_DIR $OUT_DIR $ADAPTERS
done


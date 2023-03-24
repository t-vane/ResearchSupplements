################################################################################
#### GENOTYPE CALLING WITH STACKS####
################################################################################
SCRIPTS_DIR=/home/nibtve93/scripts/stacks

NT=16

SET_ID=lehilahytsara
BAM_DIR=$PWORK/bamFiles/$SET_ID
OUT_DIR=$PWORK/$SET_ID/stacks
POPMAP=$OUT_DIR/$SET_ID.popmap # List of samples included in genotyping
SUFFIX=auto	# Suffix for final BAM files (see scripts for reference mapping)

mkdir -p $OUT_DIR/logFiles
mkdir -p $OUT_DIR/populations

## Create stacks with gstacks
sbatch --job-name=stacks --output=$OUT_DIR/logFiles/gstacks.$SET_ID.oe $SCRIPTS_DIR/gstacks.sh $NT $BAM_DIR $OUT_DIR $POPMAP 

## Extract VCF file with populations
FILTERS="-p 1 -r 0.75"
OUTPUT="--vcf"
sbatch --job-name=stacks --dependency=singleton --output=$OUT_DIR/logFiles/populations.$SET_ID.oe $SCRIPTS_DIR/populations.sh $OUT_DIR $OUT_DIR/populations $POPMAP $NT "$FILTERS" "$OUTPUT"

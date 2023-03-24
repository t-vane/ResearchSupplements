################################################################################
#### REFERENCE MAPPING AND FILTERING ####
################################################################################
SCRIPTS_DIR=/home/nibtve93/scripts/referenceMapping

SET_ID=lehilahytsara
REFERENCE_DIR=$PWORK/references/mmur3
REFERENCE=$REFERENCE_DIR/GCF_000165445.2_Mmur_3.0_genomic.fna # Reference genome in fasta format

#################################################################
#### 0 INDEX REFERENCE GENOME IF NOT DONE YET ####
#################################################################
mkdir -p $REFERENCE_DIR/logFiles

BWA=true # Boolean specifying whether to create index with BWA
BWA_INDEX=mmur3 # Output name of bwa index
SAMTOOLS=true # Boolean specifying whether to create index with SAMtools
GATK=true # Boolean specifying whether to create index with GATK (Picard)

## Submit indexing script
sbatch --wait --output=$REFERENCE_DIR/logFiles/indexing.$SET_ID.oe $SCRIPTS_DIR/indexing.sh $REFERENCE $BWA $BWA_INDEX $SAMTOOLS $GATK

#################################################################
#### 1 ALIGN TRIMMED READS TO REFERENCE GENOME AND FILTER ####
#################################################################
NT=6

IN_DIR=$PWORK/trimmedReads/$SET_ID # The following scripts assume that trimmed read files are named *.trimmed.1.fq.gz and *.trimmed.2.fq.gz 
BAM_DIR=$PWORK/bamFiles/$SET_ID

mkdir -p $OUT_DIR/logFiles

## Submit scripts for reference mapping, sorting, filtering, extraction of genomic regions and header cleaning
for i in PE SE
do
	IN_FILE=$BAM_DIR/map_$i.txt # List of samples (without file extensions) for which reference alignment shall be conducted
	NO_INDS=$(cat $IN_FILE | wc -l)
	MINMAPQ=20
	
	# Reference mapping, filtering and extraction of genomic regions
	sbatch --job-name=map_filter_pip --array=1-$NO_INDS --output=$OUT_DIR/logFiles/reference_mapping_$i.%A_%a.$SET_ID.oe $SCRIPTS_DIR/reference_mapping.sh $i $NT $REFERENCE_DIR/$BWA_INDEX $IN_DIR $BAM_DIR $IN_FILE
	
	# Sort and quality filter
	sbatch --job-name=map_filter_pip --dependency=singleton --array=1-$NO_INDS --output=$OUT_DIR/logFiles/quality_filter_$i.%A_%a.$SET_ID.oe $SCRIPTS_DIR/quality_filter.sh $i $NT $BAM_DIR $IN_FILE $MINMAPQ
	
	# Deduplicate 
	[[ $i == PE ]] && sbatch --job-name=map_filter_pip --dependency=singleton --array=1-$NO_INDS --output=$OUT_DIR/logFiles/deduplicate_$i.%A_%a.$SET_ID.oe $SCRIPTS_DIR/deduplicate.sh $BAM_DIR $IN_FILE $MINMAPQ
	
	# Extract genomic regions
	BED=$BAM_DIR/regionFile_autosomes.bed # BED file with genomic regions that shall be extracted
	EXCLUDE="NW_|NC_028718.1|NC_033692.1" # String of chromosomes that are no longer represented (separator: "|")
	SUFFIX=auto # Suffix for naming of final BAM files
	sbatch --job-name=map_filter_pip --dependency=singleton --array=1-$NO_INDS --output=$OUT_DIR/logFiles/extract_regions_$i.%A_%a.$SET_ID.oe $SCRIPTS_DIR/extract_regions.sh $i $BAM_DIR $IN_FILE $MINMAPQ $BED "$EXCLUDE" $SUFFIX
done

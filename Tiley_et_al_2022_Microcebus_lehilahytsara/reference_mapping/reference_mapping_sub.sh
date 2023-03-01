################################################################################
#### REFERENCE MAPPING AND FILTERING ####
################################################################################
SCRIPTS_DIR=/home/nibtve93/scripts/referenceMapping

REFERENCE_DIR=$PWORK/references/mmur3
REFERENCE=$REFERENCE_DIR/GCF_000165445.2_Mmur_3.0_genomic.fna # Reference genome in fasta format

#################################################################
#### 0 INDEX REFERENCE GENOME IF NOT DONE YET####
#################################################################
mkdir -p $REFERENCE_DIR/logFiles

BWA=TRUE # Boolean specifying whether to create index with BWA
BWA_INDEX=mmur3 # Output name of bwa index
SAMTOOLS=TRUE # Boolean specifying whether to create index with SAMtools
GATK=TRUE # Boolean specifying whether to create index with GATK (Picard)

## Submit indexing script
sbatch --output=$REFERENCE_DIR/logFiles/indexing.oe $SCRIPTS_DIR/indexing.sh $REFERENCE $BWA $BWA_INDEX $SAMTOOLS $GATK

#################################################################
#### 1 ALIGN TRIMMED READS TO REFERENCE GENOME AND FILTER ####
#################################################################
NT=6

IN_DIR=$PWORK/trimmedReads/lehilahytsara # Assumes that trimmed read files are named *.trimmed.1.fq.gz and *.trimmed.2.fq.gz 
OUT_DIR=$PWORK/bamFiles/lehilahytsara

QFILTER=TRUE # Boolean specifying whether to sort and filter for proper pairing (only paired-end) and minimum mapping quality 
MINMAPQ=20
DEDUPLICATE=TRUE # Boolean specifying whether to deduplicate (only paired-end)
EXTRACT_BED=TRUE # Boolean specifying whether to extract specific chro
BED=$PWORK/bamFiles/lehilahytsara/regionFile_autosomes.bed # BED file with genomic regions that shall be extracted
SUFFIX=auto # Suffix for naming of final BAM files

mkdir -p $OUT_DIR/logFiles

## Submit scripts for reference mapping, sorting, filtering, extraction of genomic regions and header cleaning
for i in PE SE
do
	IN_FILE=$PWORK/bamFiles/lehilahytsara/map_$i.txt # List of samples (without file extensions) for which reference alignment shall be conducted
	NO_INDS=$(cat $IN_FILE | wc -l)
	
	# Reference mapping, filtering and extraction of genomic regions
	sbatch --job-name=map_filter_pip --array=1-$NO_INDS --output=$OUT_DIR/logFiles/reference_mapping_$i.%A_%a.oe $SCRIPTS_DIR/reference_mapping.sh $i $NT $REFERENCE_DIR/$BWA_INDEX $IN_DIR $OUT_DIR $IN_FILE \
		$QFILTER $MINMAPQ $DEDUPLICATE $EXTRACT_BED $BED $SUFFIX
		
	# Removing chromosome names that are no longer represented in headers
	EXCLUDE="NW_|NC_028718.1|NC_033692.1" # String of chromosomes that are no longer represented (separator: "|")
	sbatch --job-name=map_filter_pip --array=1-$NO_INDS --dependency=singleton --output=$OUT_DIR/logFiles/change_header_$i.%A_%a.oe $SCRIPTS_DIR/change_header.sh $i $OUT_DIR $IN_FILE $MINMAPQ $SUFFIX "$EXCLUDE"
done



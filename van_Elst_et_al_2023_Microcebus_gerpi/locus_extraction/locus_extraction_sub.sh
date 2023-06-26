################################################################################
#### LOCUS EXTRACTION WITH CUSTOM PIPELINE ####
################################################################################
## Pipeline adapted and modified from Poelstra et al. 2021, Systematic Biology (https://doi.org/10.1093/sysbio/syaa053)

SCRIPTS_DIR=/home/nibtve93/scripts/locusExtraction

SET_ID=gerpi
REFERENCE=$PWORK/mmur3/GCF_000165445.2_Mmur_3.0_genomic_reduced.fna
BAM_DIR=$PWORK/bamFiles/$SET_ID
VCF_DIR=$PWORK/$SET_ID/vcf
VCF_ALTREF=$VCF_DIR/allScaffolds.snps.vcf # Raw VCF file without indels and invariants
VCF_FILT_MASK=$VCF_DIR/allScaffolds.snps.05filt.part.vcf
VCF_FILT_INTERSECT=$VCF_DIR/allScaffolds.snps.06filt.vcf # Fully filtered VCF file
VCF_HIGHDEPTH=$VCF_DIR/allScaffolds.snps.05filt_too-high-DP.vcf
OUT_DIR=$PWORK/$SET_ID/locusExtraction

FASTA_DIR=$OUT_DIR/fasta
INDFASTA_DIR=$FASTA_DIR/$SET_ID.byind
BED_DIR=$OUT_DIR/bed
BED_REMOVED_SITES=$BED_DIR/$SET_ID.sitesinvcfremovedbyfilters.bed

LOCUSFASTA_DIR_INTERMED=$FASTA_DIR/$SET_ID}_bylocus_intermed
LOCUSFASTA_DIR_FINAL=$FASTA_DIR/${SET_ID}_bylocus_final

LOCUSBED_INTERMED=$BED_DIR/${SET_ID}_loci_intermed.bed
LOCUSBED_FINAL=$BED_DIR/${SET_ID}_loci_all.bed
LOCUSLIST=$BED_DIR/locuslist.txt
FASTA_MERGED=$INDFASTA_DIR/${SET_ID}_merged.fasta

LOCUSSTATS_INTERMED=$BED_DIR/$SET_ID/${SET_ID}_locusstats_all.txt
LOCUSSTATS_FINAL=$BED_DIR/$SET_ID/${SET_ID}_locusstats_filt.txt

mkdir -p $FASTA_DIR
mkdir -p $INDFASTA_DIR
mkdir -p $BED_DIR

mkdir $LOCUSFASTA_DIR_INTERMED
mkdir $LOCUSFASTA_DIR_FINAL
mkdir $LOCUSFASTA_DIR_FINAL_GPHOCS
mkdir $BED_DIR"/"$SET_ID

mkdir -p $OUT_DIR/logFiles/

## Get individuals present in VCF file
IND_FILE=$OUT_DIR/slurm.indfile.$SET_ID.tmp
bcftools query -l $VCF_FILT_INTERSECT > $IND_FILE
NO_INDS=$(cat $IND_FILE | wc -l)

#################################################################
#### 1 CREATE MASKED REFERENCE GENOME PER INDIVIDUAL ####
#################################################################
## Create BED file with masked sites from VCF
sbatch --job-name=locus_extract_pip --account=nib00015 --output=$OUT_DIR/logFiles/01_maskbed.$SET_ID.oe $SCRIPTS_DIR/01_maskbed.sh $VCF_ALTREF $VCF_FILT_MASK $BED_REMOVED_SITES $BED_DIR

## Produce masked FASTA file per individual
SUFFIX=auto
MIN_DP=3
sbatch --job-name=locus_extract_pip --dependency=singleton --account=nib00015 --array=1-$NO_INDS --output=$OUT_DIR/logFiles/02_process-inds.%A_%a.$SET_ID.oe $SCRIPTS_DIR/02_process-inds.sh \
	$IND_FILE $VCF_ALTREF $REFERENCE $BAM_DIR $SUFFIX $MIN_DP $INDFASTA_DIR $BED_DIR $BED_REMOVED_SITES

#################################################################
#### 2 EXTRACT AND FILTER LOCI ACROSS INDIVIDUALS ####
#################################################################
## Make BED file with desired locus coordinates
MIN_ELEM_OVL=0.9 # Minimum element overlap for locus creation
MIN_ELEM_OVL_TRIM=0.8 # Minimum element overlap for locus trimming
MIN_LOCUS_SIZE=100 # Minimum locus size
MAX_DIST_WITHIN_IND=10 # Maximum distance within individuals
MAX_DIST_BETWEEN_IND=0 # Maximum distance between individuals
MIN_ELEM_SIZE=25 # Minimum locus size
LAST_ROW=0 # Number of loci to process (all if 0)
sbatch --job-name=locus_extract_pip --dependency=singleton --account=nib00015 --output=$OUT_DIR/logFiles/03a_makelocusbed.$SET_ID.oe $SCRIPTS_DIR/03a_makelocusbed.sh $SCRIPTS_DIR $SET_ID $IND_FILE $BED_DIR $LOCUSBED_INTERMED \
	$MIN_ELEM_OVL $MIN_ELEM_OVL_TRIM $MIN_LOCUS_SIZE $MAX_DIST_WITHIN_IND $MAX_DIST_BETWEEN_IND $MIN_ELEM_SIZE $LAST_ROW

## Intersect BED file with loci with too high depth
sbatch --job-name=locus_extract_pip --dependency=singleton --account=nib00015 --output=$OUT_DIR/logFiles/03b_intersect.$SET_ID.oe $SCRIPTS_DIR/03b_intersect.sh $LOCUSBED_INTERMED $LOCUSBED_FINAL $VCF_HIGHDEPTH $VCF_FILT_INTERSECT

## Get merged FASTA file with all individuals and loci
sbatch --job-name=locus_extract_pip --dependency=singleton --account=nib00015 --output=$OUT_DIR/logFiles/03c_mergedfasta.$SET_ID.oe $SCRIPTS_DIR/03c_mergedfasta.sh $IND_FILE $INDFASTA_DIR $LOCUSBED_FINAL $LOCUSLIST $FASTA_MERGED

## Create by-locus FASTA files
NLOCI=$(cat $LOCUSLIST | wc -l)
sbatch --job-name=locus_extract_pip --dependency=singleton --account=nib00015 --array=1-$NLOCI --output=$OUT_DIR/logFiles/03d_locusfasta.%A_%a.$SET_ID.oe $SCRIPTS_DIR/03d_locusfasta.sh $LOCUSLIST $LOCUSFASTA_DIR_INTERMED $FASTA_MERGED

## Estimate statistics for intermediate loci
sbatch --job-name=locus_extract_pip --dependency=singleton --account=nib00015 --output=$OUT_DIR/logFiles/03e_locusstats_intermed.$SET_ID.oe $SCRIPTS_DIR/03e_locusstats.sh $LOCUSFASTA_DIR_INTERMED $LOCUSSTATS_INTERMED

## Filter loci for maximum proportion of missing data and minimum distance between loci (submitted twice with different MAXMISS because a reduced locus set is used for G-PhoCS)
MAXMISS=5 # Maximum percentage of missing data in percent
MINDIST=10000 # Minimum distance (bp) between loci
sbatch --job-name=locus_extract_pip --dependency=singleton --account=nib00015 --output=$OUT_DIR/logFiles/03f_filterloci.$SET_ID.oe $SCRIPTS_DIR/03f_filterloci.sh $LOCUSSTATS_INTERMED $LOCUSFASTA_DIR_INTERMED $LOCUSFASTA_DIR_FINAL $MAXMISS $MINDIST

## Estimate statistics for final loci (submitted twice with because a reduced locus set is used for G-PhoCS)
sbatch --job-name=locus_extract_pip --dependency=singleton --account=nib00015 --output=$OUT_DIR/logFiles/03e_locusstats_final.$SET_ID.oe $SCRIPTS_DIR/03e_locusstats.sh $LOCUSFASTA_DIR_FINAL $LOCUSSTATS_FINAL

## Archive and remove intermediate locus files
for i in 03366 03367 03368 03369
do
tar -vcf $LOCUSFASTA_DIR_INTERMED/${SET_ID}_loci_intermed_$i.tar $LOCUSFASTA_DIR_INTERMED/NC_$i*fa && rm $LOCUSFASTA_DIR_INTERMED/NC_$i*fa
done




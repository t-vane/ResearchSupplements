################################################################################
#### GENOTYPE CALLING WITH GATK####
################################################################################
SCRIPTS_DIR=/home/nibtve93/scripts/gatk

SET_ID=gerpi
BAM_DIR=$PWORK/bamFiles/$SET_ID
IND_FILE=$PWORK/$SET_ID/gatk/individuals.txt # Contains individual IDs in list format (".bam" will be appended in haplotypeCaller.sh)
REFERENCE_DIR=$PWORK/mmur3
REFERENCE=$REFERENCE_DIR/GCF_000165445.2_Mmur_3.0_genomic.fna # Reference genome in fasta format
REGION_FILE=$REFERENCE_DIR/regionFileAutosomes_modified.bed # Regions (i.e., scaffolds) for which genotyping should be conducted jointly; should start at 1 and not at 0

GVCF_DIR=$PWORK/$SET_ID/gatk/gvcfFiles
DB_DIR=$PWORK/$SET_ID/gatk/DBs
VCF_SCAFFOLD_DIR=$PWORK/$SET_ID/gatk/vcfFilesScaffolds
TMP_DIR=$PWORK/$SET_ID/gatk/tmp

mkdir -p $PWORK/$SET_ID/gatk/logFiles
mkdir -p $TMP_DIR
mkdir -p $GVCF_DIR
mkdir -p $DB_DIR
mkdir -p $VCF_SCAFFOLD_DIR

## Variant discovery with haplotype caller
NT=40
MEM=100
NO_INDS=$(cat $IND_FILE | wc -l)
SUFFIX=auto.bam
sbatch --array=1-$NO_INDS --job-name=gatk -c $NT --account=nib00015 --output=$PWORK/$SET_ID/gatk/logFiles/haplotype_caller.%A_%a.oe $SCRIPTS_DIR/haplotype_caller.sh $NT $MEM $REFERENCE $IND_FILE $BAM_DIR $GVCF_DIR $SUFFIX

## Joint genotyping per scaffold
NO_REGIONS=$(cat $REGION_FILE | wc -l)
sbatch --array=1-$NO_REGIONS --job-name=gatk --dependency=singleton -c $NT --account=nib00015 --output=$PWORK/$SET_ID/gatk/logFiles/joint_genotyping.%A_%a.oe $SCRIPTS_DIR/joint_genotyping.sh $NT $REFERENCE $IND_FILE $REGION_FILE $GVCF_DIR $DB_DIR $VCF_SCAFFOLD_DIR $TMP_DIR

## Merge per-scaffold VCFs
ls $VCF_SCAFFOLD_DIR/*vcf.gz > $VCF_SCAFFOLD_DIR/allScaffolds.vcflist
sbatch --job-name=gatk --dependency=singleton --account=nib00015 --output=$PWORK/$SET_ID/gatk/logFiles/merge_vcfs.oe $SCRIPTS_DIR/merge_vcfs.sh $VCF_SCAFFOLD_DIR/allScaffolds.vcflist $PWORK/$SET_ID/gatk/allScaffolds.vcf.gz




################################################################################
#### VCF FILTERING ####
################################################################################
## Pipeline adapted and modified from Poelstra et al. 2021, Systematic Biology (https://doi.org/10.1093/sysbio/syaa053)

SCRIPTS_DIR=/home/nibtve93/scripts/vcfFiltering

SET_ID=gerpi
VCF_RAW=$PWORK/$SET_ID/gatk/allScaffolds.vcf
VCF_DIR=$PWORK/$SET_ID/vcf
BAM_DIR=$PWORK/bamFiles/$SET_ID
REFERENCE=$PWORK/mmur3/GCF_000165445.2_Mmur_3.0_genomic_reduced.fna # Reference genome in fasta format, not containing unlocalized chromosomal scaffolds between chromosomal ones

mkdir -p $VCF_DIR/logFiles

## Softlink raw VCF to VCF directory
ln -s $VCF_RAW $VCF_DIR

#################################################################
#### 1 MAIN FILTERING PIPELINE ####
#################################################################
## Remove indels and invariant sites
sbatch --job-name=vcf_filter_pip --account=nib00015 --output=$VCF_DIR/logFiles/00_filter_indels_invariants.$SET_ID.oe $SCRIPTS_DIR/00_filter_indels_invariants.sh $REFERENCE $VCF_RAW $VCF_DIR/allScaffolds.snps.vcf

## Filter for minimum depth
MIN_DP=5
MEAN_DP=5
JID1=$(sbatch --job-name=vcf_filter_pip --dependency=singleton --account=nib00015 --output=$VCF_DIR/logFiles/01_filter_min-dp.$SET_ID.oe $SCRIPTS_DIR/01_filter_min-dp.sh $VCF_DIR/allScaffolds.snps.vcf $MIN_DP $MEAN_DP $VCF_DIR/$SET_ID.allScaffolds.snps.01filt.vcf)

## Apply three rounds of filtering for missing data across individuals and genotypes
MAXMISS_GENO1=0.5 # One minus maxmimum missingness across genotypes (filtering round 1), i.e., maximum missingness = 1-$MAXMISS_GENO1
MAXMISS_GENO2=0.6 # One minus maxmimum missingness across genotypes (filtering round 2), i.e., maximum missingness = 1-$MAXMISS_GENO2
MAXMISS_GENO3=0.7 # One minus maxmimum missingness across genotypes (filtering round 2), i.e., maximum missingness = 1-$MAXMISS_GENO3
FILTER_INDS=TRUE # Boolean specifying whether to filter for missingness across individuals
MAXMISS_IND1=0.9 # Maxmimum missingness across individuals (filtering round 1)
MAXMISS_IND2=0.7 # Maxmimum missingness across individuals (filtering round 2)
MAXMISS_IND3=0.5 # Maxmimum missingness across individuals (filtering round 2)
sbatch --job-name=vcf_filter_pip --dependency=singleton --account=nib00015 --output=$VCF_DIR/logFiles/02_filter_missing-1.$SET_ID.oe $SCRIPTS_DIR/02_filter_missing-1.sh \
	$VCF_DIR/$SET_ID.allScaffolds.snps.01filt.vcf $VCF_DIR/$SET_ID.allScaffolds.snps.02filt.vcf $MAXMISS_GENO1 $MAXMISS_GENO2 $MAXMISS_GENO3 $FILTER_INDS $MAXMISS_IND1 $MAXMISS_IND2 $MAXMISS_IND3

## Annotate with INFO fields FisherStrand, RMSMappingQuality, MappingQualityRankSumTest, ReadPosRankSumTest and AlleleBalance
SUFFIX=auto
sbatch --job-name=vcf_filter_pip --dependency=singleton --account=nib00015 --output=$VCF_DIR/logFiles/03_annot_gatk.$SET_ID.oe $SCRIPTS_DIR/03_annot_gatk.sh \
	$VCF_DIR/$SET_ID.allScaffolds.snps.02filt.vcf $VCF_DIR/$SET_ID.allScaffolds.snps.03filt.vcf $BAM_DIR $REFERENCE $SUFFIX

## Retain only bi-allelic sites and filter for INFO fields FisherStrand, RMSMappingQuality, MappingQualityRankSumTest, ReadPosRankSumTest and AlleleBalance
sbatch --job-name=vcf_filter_pip --dependency=singleton --account=nib00015 --output=$VCF_DIR/logFiles/04_filter_gatk.$SET_ID.oe $SCRIPTS_DIR/04_filter_gatk.sh \
	$VCF_DIR/$SET_ID.allScaffolds.snps.03filt.vcf $VCF_DIR/$SET_ID.allScaffolds.snps.04filt-soft.vcf $VCF_DIR/$SET_ID.allScaffolds.snps.04filt-hard.vcf $REFERENCE

## Filter for maximum depth as (mean depth + 2 * standard deviation) / number of individuals
sbatch --job-name=vcf_filter_pip --dependency=singleton --account=nib00015 --output=$VCF_DIR/logFiles/05_filter_max-dp.$SET_ID.oe $SCRIPTS_DIR/05_filter_max-dp.sh \
	$VCF_DIR/$SET_ID.allScaffolds.snps.04filt-hard.vcf $VCF_DIR/$SET_ID.allScaffolds.snps.05filt.vcf 

## Apply final round of filtering for missing data across individuals and genotypes
MAXMISS_GENO=0.9 # One minus maxmimum missingness across genotypes, i.e., maximum missingness = 1-$MAXMISS_GENO
FILTER_INDS=TRUE # Boolean specifying whether to filter for missingness across individuals
MAXMISS_IND=0.5 # Maxmimum missingness across individuals
sbatch --job-name=vcf_filter_pip --dependency=singleton --account=nib00015 --output=$VCF_DIR/logFiles/06_filter_missing-2.$SET_ID.oe $SCRIPTS_DIR/06_filter_missing-2.sh \
	$VCF_DIR/$SET_ID.allScaffolds.snps.05filt.vcf $VCF_DIR/$SET_ID.allScaffolds.snps.06filt.vcf $MAXMISS_GENO $FILTER_INDS $MAXMISS_IND
	
## Create VCF file without outgroups
REM_STRING="--remove-indv Mmur_RMR44 --remove-indv Mmur_RMR45 --remove-indv Mmur_RMR49 --remove-indv Mjol_LAKI5.20a --remove-indv Mjol_LAKI5.24 --remove-indv Mmaro_RMR131"
sbatch --job-name=vcf_filter_pip --dependency=singleton --account=nib00015 --output=$VCF_DIR/logFiles/07_filter_outgroups.$SET_ID.oe $SCRIPTS_DIR/07_filter_outgroups.sh \
	$VCF_DIR/$SET_ID.allScaffolds.snps.06filt.vcf $VCF_DIR/$SET_ID.allScaffolds.snps.07filt.vcf "$REM_STRING"

#################################################################
#### 2 CREATING PARTIALLY FILTERED VCF (NECESSARY FOR LOCUS EXTRACTION) ####
#################################################################
## Here, we skip scripts that filter for missing data based on individuals and genotypes, i.e., 02_filter_missing-1.sh and 06_filter_missing-2.sh

## Annotate with INFO fields FisherStrand, RMSMappingQuality, MappingQualityRankSumTest, ReadPosRankSumTest and AlleleBalance
SUFFIX=auto
sbatch --job-name=vcf_filter_part_pip --dependency=afterany:${JID1##* } --account=nib00015 --output=$VCF_DIR/logFiles/03_annot_gatk.part.$SET_ID.oe $SCRIPTS_DIR/03_annot_gatk.sh \
	$VCF_DIR/$SET_ID.allScaffolds.snps.01filt.vcf $VCF_DIR/$SET_ID.allScaffolds.snps.03filt.part.vcf $BAM_DIR $REFERENCE $SUFFIX

## Filter for INFO fields FisherStrand, RMSMappingQuality, MappingQualityRankSumTest, ReadPosRankSumTest and AlleleBalance
sbatch --job-name=vcf_filter_part_pip --dependency=singleton --account=nib00015 --output=$VCF_DIR/logFiles/04_filter_gatk.part.$SET_ID.oe $SCRIPTS_DIR/04_filter_gatk.sh \
	$VCF_DIR/$SET_ID.allScaffolds.snps.03filt.part.vcf $VCF_DIR/$SET_ID.allScaffolds.snps.04filt-soft.part.vcf $VCF_DIR/$SET_ID.allScaffolds.snps.04filt-hard.part.vcf $REFERENCE

## Filter for maximum depth as (mean depth + 2 * standard deviation) / number of individuals
sbatch --job-name=vcf_filter_part_pip --dependency=singleton --account=nib00015 --output=$VCF_DIR/logFiles/05_filter_max-dp.part.$SET_ID.oe $SCRIPTS_DIR/05_filter_max-dp.sh \
	$VCF_DIR/$SET_ID.allScaffolds.snps.04filt-hard.part.vcf $VCF_DIR/$SET_ID.allScaffolds.snps.05filt.part.vcf 


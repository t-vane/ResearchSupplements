################################################################################
#### VCF FILTERING ####
################################################################################
## Pipeline adapted and modified from Poelstra et al. 2021, Systematic Biology (https://doi.org/10.1093/sysbio/syaa053)

## Software:
# bgzip needs to be included in $PATH (v1.11; http://www.htslib.org/)

SCRIPTS_DIR=/home/nibtve93/scripts/vcfFiltering

SET_ID=lehilahytsara
VCF_RAW=$PWORK/$SET_ID/stacks/populations/populations.snps.vcf
VCF_DIR=$PWORK/$SET_ID/vcf
BAM_DIR=$PWORK/bamFiles/$SET_ID
REFERENCE=$PWORK/references/mmur3/GCF_000165445.2_Mmur_3.0_genomic_reduced.fna # Reference genome in fasta format, not containing unlocalized chromosomal scaffolds between chromosomal ones

mkdir -p $VCF_DIR/logFiles

## Softlink raw VCF to VCF directory
ln -s $IN_FILE $VCF_DIR

#################################################################
#### 1 MAIN FILTERING PIPELINE ####
#################################################################
## Filter for minimum depth
MIN_DP=5
MEAN_DP=5
sbatch --job-name=vcf_filter_pip --output=$VCF_DIR/logFiles/01_filter_min-dp.$SET_ID.oe $SCRIPTS_DIR/01_filter_min-dp.sh $VCF_RAW $MIN_DP $MEAN_DP $VCF_DIR/$SET_ID.populations.snps.01filt.vcf

## Apply three rounds of filtering for missing data across individuals and genotypes
MAXMISS_GENO1=0.5 # One minus maxmimum missingness across genotypes (filtering round 1), i.e., maximum missingness = 1-$MAXMISS_GENO1
MAXMISS_GENO2=0.6 # One minus maxmimum missingness across genotypes (filtering round 2), i.e., maximum missingness = 1-$MAXMISS_GENO2
MAXMISS_GENO3=0.7 # One minus maxmimum missingness across genotypes (filtering round 3), i.e., maximum missingness = 1-$MAXMISS_GENO3
FILTER_INDS=TRUE # Boolean specifying whether to filter for missingness across individuals
MAXMISS_IND1=0.875 # Maxmimum missingness across individuals (filtering round 1)
MAXMISS_IND2=0.7 # Maxmimum missingness across individuals (filtering round 2)
MAXMISS_IND3=0.5 # Maxmimum missingness across individuals (filtering round 3)
sbatch --job-name=vcf_filter_pip --dependency=singleton --output=$VCF_DIR/logFiles/02_filter_missing-1.$SET_ID.oe $SCRIPTS_DIR/02_filter_missing-1.sh \
	$VCF_DIR/$SET_ID.populations.snps.01filt.vcf $VCF_DIR/$SET_ID.populations.snps.02filt.vcf $MAXMISS_GENO1 $MAXMISS_GENO2 $MAXMISS_GENO3 $FILTER_INDS $MAXMISS_IND1 $MAXMISS_IND2 $MAXMISS_IND3

## Annotate with INFO fields FisherStrand, RMSMappingQuality, MappingQualityRankSumTest, ReadPosRankSumTest and AlleleBalance
SUFFIX=auto
sbatch --job-name=vcf_filter_pip --dependency=singleton --output=$VCF_DIR/logFiles/03_annot_gatk.$SET_ID.oe $SCRIPTS_DIR/03_annot_gatk.sh \
	$VCF_DIR/$SET_ID.populations.snps.02filt.vcf $VCF_DIR/$SET_ID.populations.snps.03filt.vcf $BAM_DIR $REFERENCE $SUFFIX

## Filter for INFO fields FisherStrand, RMSMappingQuality, MappingQualityRankSumTest, ReadPosRankSumTest and AlleleBalance
sbatch --job-name=vcf_filter_pip --dependency=singleton --output=$VCF_DIR/logFiles/04_filter_gatk.$SET_ID.oe $SCRIPTS_DIR/04_filter_gatk.sh \
	$VCF_DIR/$SET_ID.populations.snps.03filt.vcf $VCF_DIR/$SET_ID.populations.snps.04filt-soft.vcf $VCF_DIR/$SET_ID.populations.snps.04filt-hard.vcf $REFERENCE

## Filter for maximum depth as (mean depth + 2 * standard deviation) / number of individuals
sbatch --job-name=vcf_filter_pip --dependency=singleton --output=$VCF_DIR/logFiles/05_filter_max-dp.$SET_ID.oe $SCRIPTS_DIR/05_filter_max-dp.sh \
	$VCF_DIR/$SET_ID.populations.snps.04filt-hard.vcf $VCF_DIR/$SET_ID.populations.snps.05filt.vcf 

## Apply final round of filtering for missing data across individuals and genotypes
MAXMISS_GENO=0.9 # One minus maxmimum missingness across genotypes, i.e., maximum missingness = 1-$MAXMISS_GENO
FILTER_INDS=TRUE # Boolean specifying whether to filter for missingness across individuals
MAXMISS_IND=0.5 # Maxmimum missingness across individuals
sbatch --job-name=vcf_filter_pip --dependency=singleton --output=$VCF_DIR/logFiles/06_filter_missing-2.$SET_ID.oe $SCRIPTS_DIR/06_filter_missing-2.sh \
	$VCF_DIR/$SET_ID.populations.snps.05filt.vcf $VCF_DIR/$SET_ID.populations.snps.06filt.vcf $MAXMISS_GENO $FILTER_INDS $MAXMISS_IND

## Apply minor allele count filter
MAC=3
sbatch --job-name=vcf_filter_pip --dependency=singleton --output=$VCF_DIR/logFiles/07_filter_mac.$SET_ID.oe $SCRIPTS_DIR/07_filter_mac.sh \
	$VCF_DIR/$SET_ID.populations.snps.06filt.vcf $VCF_DIR/$SET_ID.populations.snps.07filt.vcf 

## Remove outgroup individuals and those not needed for population structure analyses
DROP_INDS=$VCF_DIR/$SET_ID.samples.drop.txt # List with individuals to remove from VCF (i.e., outgroups and undesired populations), without row or column names
sbatch --job-name=vcf_filter_pip --dependency=singleton --output=$VCF_DIR/logFiles/08_drop_inds.$SET_ID.oe $SCRIPTS_DIR/08_drop_inds.sh \
	$VCF_DIR/$SET_ID.populations.snps.07filt.vcf $VCF_DIR/$SET_ID.populations.snps.08filt.vcf 
# bgzip final file
bgzip -c $VCF_DIR/$SET_ID.populations.snps.08filt.vcf  > $VCF_DIR/$SET_ID.populations.snps.08filt.vcf.gz

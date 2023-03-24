################################################################################
#### POPULATION STRUCTURE ANALYSES ####
################################################################################
SCRIPTS_DIR=/home/nibtve93/scripts/populationStructure

SET_ID=lehilahytsara
POP_DIR=$PWORK/$SET_ID/populationStructure
BEAGLE=$PWORK/$SET_ID/angsd/$SET_ID.beagle.gz # Genotype likelihood file created in genotype_likelihoods_sub.sh

mkdir -p $POP_DIR/logFiles


#################################################################
#### 1 PRINCIPAL COMPONENT ANALYSIS ####
#################################################################

##################################################
#### 1.1 GENOTYPE CALL BASED ####
##################################################

## See Dryad digital repository (https://doi.org/10.5061/dryad.dncjsxkvz).

##################################################
#### 1.2 GENOTYPE LIKELIHOOD BASED ####
##################################################
mkdir -p $POP_DIR/pca/genotypeLikelihoods

IND_FILE=$POP_DIR/$SET_ID.txt # File with individual IDs in first columns and population assignments in second column
NT=20

sbatch --output=$POP_DIR/logFiles/pca.$SET_ID.oe $SCRIPTS_DIR/pca.sh $NC $BEAGLE $POP_DIR/pca/genotypeLikelihoods $SCRIPTS_DIR $IND_FILE $SET_ID


#################################################################
#### 2 ADMIXTURE ####
#################################################################

##################################################
#### 2.1 GENOTYPE CALL BASED ####
##################################################

## See Dryad digital repository (https://doi.org/10.5061/dryad.dncjsxkvz).

##################################################
#### 2.2 GENOTYPE LIKELIHOOD BASED ####
##################################################
mkdir -p $POP_DIR/ngsadmix/genotypeLikelihoods

CLUSTERS=12 # Maximum number of clusters to assume in admixture analysis
REPEATS=10 # Number of independent runs
PERCENTAGE="75/100" # Minimum percentage of represented individuals
MININD=$(( ($(zcat $BEAGLE | head -1 | wc -w)/3-1) * $PERCENTAGE )) # Minimum number of represented individuals
NC=80

## Submit array job to infer individual ancestries for each number of clusters (ranging from 1 to $CLUSTERS), using $REPEATS repetitions 
for K in $(seq 1 $CLUSTERS)
do
	[[ $K < $CLUSTERS ]] && sbatch --array=1-$ITERATIONS --output=$POP_DIR/logFiles/ngsadmix$K.$SET_ID.%A_%a.oe $SCRIPTS_DIR/ngsadmix.sh $NT $K $BEAGLE $POP_DIR/ngsadmix/genotypeLikelihoods $MININD $SET_ID
	[[ $K == $CLUSTERS ]] && sbatch --wait --array=1-$ITERATIONS --output=$POP_DIR/logFiles/ngsadmix$K.$SET_ID.%A_%a.oe $SCRIPTS_DIR/ngsadmix.sh $NT $K $BEAGLE $POP_DIR/ngsadmix/genotypeLikelihoods $MININD $SET_ID
done

## Print likelihood values file
rm $POP_DIR/ngsadmix/genotypeLikelihoods/likevalues.$SET_ID.txt; touch $POP_DIR/ngsadmix/genotypeLikelihoods/likevalues.$SET_ID.txt
for K in $(seq 1 $CLUSTERS); 
do
	for SEED in $(seq 1 $REPEATS)
	do
		grep "best" $POP_DIR/ngsadmix/genotypeLikelihoods/$SET_ID.K$K.seed$SEED.log | awk '{print $K}' | cut -d'=' -f2- | sort -g | sed "s/after/$K/g" | sed "s/iterations/$SEED/g" >> $POP_DIR/ngsadmix/genotypeLikelihoods/likevalues.$SET_ID.txt
		rm $POP_DIR/ngsadmix/genotypeLikelihoods/$SET_ID.K$K.seed$SEED.fopt.gz
	done
done

## Plot results
IND_FILE=$POP_DIR/$SET_ID.txt # File with individual IDs in first columns and population assignments in second column

sbatch --output=$POP_DIR/logFiles/plot_ngsadmix.$SET_ID.oe $SCRIPTS_DIR/plot_ngsadmix.sh $SCRIPTS_DIR $POP_DIR/ngsadmix/genotypeLikelihoods $POP_DIR/ngsadmix/genotypeLikelihoods/likevalues.$SET_ID.txt $IND_FILE $SET_ID


#################################################################
#### 3 ISOLATION-BY-DISTANCE ####
#################################################################

##################################################
#### 3.1 GENOTYPE CALL BASED ####
##################################################
mkdir -p $POP_DIR/ibd/genotypeCalls

VCF_IN=$PWORK/$SET_ID/vcf/$SET_ID.populations.snps.08filt.vcf.gz # Filtered input VCF without outgroups and undesired individuals created in vcf_filtering_sub.sh
SAMPLE_FILE=$POP_DIR/ibd/genotypeCalls/$SET_ID.samples.pops.txt # List with indviduals in first column and associated populations in second column, without row or column names
OUT=$POP_DIR/ibd/genotypeCalls/$SET_ID

## Submit script to infer pairwise F_ST between populations
sbatch --output=$POP_DIR/logFiles/$SET_ID.hierfstat.oe $SCRIPTS_DIR/hierfstat.sh $SCRIPTS_DIR $VCF_IN $SAMPLE_FILE $OUT

## Conduct Mantel tests and plot IBD
GEO_DIST=$POP_DIR/ibd/genotypeCalls/geo_dist.txt # Distance matrix with mean geographic distances between population pairs as estimated with Geographic Distance Matrix Generator v1.2.3 (https://biodiversityinformatics.amnh.org/open_source/gdmg/), with row and column names
GEN_DIST=$POP_DIR/ibd/genotypeCalls/gen_dist.txt # Distance matrix with pairwise F_ST between populations as estimated with hierfstat, with row and column names
sbatch --output=$POP_DIR/logFiles/$SET_ID.genotypeCalls.ibd.oe $SCRIPTS_DIR/ibd.sh $GEO_DIST $GEN_DIST $POP_DIR/ibd/genotypeCalls/$SET_IT

##################################################
#### 3.2 GENOTYPE LIKELIHOOD BASED ####
##################################################
mkdir -p $POP_DIR/ibd/genotypeLikelihoods

SAF_DIR=$PWORK/$SET_ID/angsd # Directory with site allele frequency likelihoods for each population inferred in genotype_likelihoods_sub.sh
POPS="ambavala ambatovy ambohitantely anjanaharibe anjiahely ankafobe marojejy tsinjoarivo"
NT=24

## Get pairwise population comparisons
awk -f $SCRIPTS_DIR/combinations.awk <<< $POPS > $POP_DIR/ibd/genotypeLikelihoods/pop_combinations.txt

## Initialize summary file
echo "pair unweighted weighted" > $POP_DIR/ibd/genotypeLikelihoods/$SET_ID.fst_sumstats.txt

## Get pairwise F_ST between populations
while read combination
do
	FIRST=$(awk '{print $1}' <<< $combination)
	SECOND=$(awk '{print $2}' <<< $combination)
	
	echo -e "#### Processing combination $FIRST $SECOND ...\n"
	# Estimate joint minor allele frequency spectrum
	OUT_FILE=$POP_DIR/ibd/genotypeLikelihoods/$SET_ID.$FIRST.$SECOND.ml
	sbatch --job-name=fst --output=$POP_DIR/logFiles/$SET_ID.realsfs.2d.$FIRST.$SECOND.oe $SCRIPTS_DIR/realsfs.sh $NT "$SAF_DIR/$FIRST.saf.idx $SAF_DIR/$SECOND.saf.idx" $OUT_FILE
	
	# Estimate F_ST values
	OUT=$POP_DIR/ibd/genotypeLikelihoods/$SET_ID.$FIRST.$SECOND
	sbatch --wait --job-name=fst --dependency=singleton --output=$POP_DIR/logFiles/$SET_ID.realsfs_fst.$FIRST.$SECOND.oe $SCRIPTS_DIR/realsfs_fst.sh $NT "$SAF_DIR/$FIRST.saf.idx $SAF_DIR/$SECOND.saf.idx" $OUT
	
	# Write to summary file 
	echo ${FIRST}_$SECOND $(cat $OUT.fst) >> $POP_DIR/ibd/genotypeLikelihoods/$SET_ID.fst_sumstats.txt
done < $POP_DIR/ibd/genotypeLikelihoods/pop_combinations.txt

## Conduct Mantel tests and plot IBD
GEO_DIST=$POP_DIR/ibd/genotypeLikelihoods/geo_dist.txt # Distance matrix with mean geographic distances between population pairs as estimated with Geographic Distance Matrix Generator v1.2.3 (https://biodiversityinformatics.amnh.org/open_source/gdmg/), with row and column names
GEN_DIST=$POP_DIR/ibd/genotypeLikelihoods/gen_dist.txt # Distance matrix with weighted pairwise F_ST between populations as estimated with realSFS, with row and column names
sbatch --output=$POP_DIR/logFiles/$SET_ID.genotypeLikelihoods.ibd.oe $SCRIPTS_DIR/ibd.sh $GEO_DIST $GEN_DIST $POP_DIR/ibd/genotypeLikelihoods/$SET_IT


#################################################################
#### AMOVA ####
#################################################################
mkdir -p $POP_DIR/amova
VCF_IN=$PWORK/$SET_ID/vcf/$SET_ID.populations.snps.08filt.vcf.gz # Filtered input VCF without outgroups and undesired individuals created in vcf_filtering_sub.sh
SAMPLE_FILE=$POP_DIR/amova/$SET_ID.samples.pops.txt # List with indviduals in first column, populations in second column and major clusters (north/south/central) in third column, without row or column names
OUT=$POP_DIR/amova/$SET_ID

## Submit script for AMOVA
sbatch --output=$POP_DIR/logFiles/$SET_ID.amova.oe $SCRIPTS_DIR/amova.sh $SCRIPTS_DIR $VCF_IN $SAMPLE_FILE $OUT

################################################################################
#### POPULATION STRUCTURE ANALYSES ####
################################################################################
SCRIPTS_DIR=/home/nibtve93/scripts/populationStructure

SET_ID=gerpi
POP_DIR=$PWORK/$SET_ID/populationStructure
BEAGLE=$PWORK/$SET_ID/angsd/$SET_ID.beagle.gz # Genotype likelihood file created in genotype_likelihoods_sub.sh

mkdir -p $POP_DIR/logFiles

#################################################################
#### 1 PRINCIPAL COMPONENT ANALYSIS ####
#################################################################
mkdir -p $POP_DIR/pca

IND_FILE=$POP_DIR/$SET_ID.txt # File with individual IDs in first column and population assignments in second column
NT=20

sbatch --account=nib00015 --output=$POP_DIR/logFiles/pca.$SET_ID.oe $SCRIPTS_DIR/pca.sh $NT $BEAGLE $POP_DIR/pca $SCRIPTS_DIR $IND_FILE $SET_ID

#################################################################
#### 2 ADMIXTURE ####
#################################################################
mkdir -p $POP_DIR/ngsadmix

CLUSTERS=10 # Maximum number of clusters to assume in admixture analysis
REPEATS=10 # Number of independent runs
PERCENTAGE="75/100" # Minimum percentage of represented individuals
MININD=$(( ($(zcat $BEAGLE | head -1 | wc -w)/3-1) * $PERCENTAGE )) # Minimum number of represented individuals
NT=80

## Submit array job to infer individual ancestries for each number of clusters (ranging from 1 to $CLUSTERS), using $REPEATS repetitions 
for K in $(seq 1 $CLUSTERS)
do
	JID=$(sbatch --array=1-$REPEATS --output=$POP_DIR/logFiles/ngsadmix$K.$SET_ID.%A_%a.oe $SCRIPTS_DIR/ngsadmix.sh $NT $K $BEAGLE $POP_DIR/ngsadmix $MININD $SET_ID)
	declare RUNID_$K=${JID##* }
done

## Print likelihood values file
LIKE_FILE=$POP_DIR/ngsadmix/likevalues.$SET_ID.txt # File for likelihoods summary
rm $LIKE_FILE; touch $LIKE_FILE
for K in $(seq 1 $CLUSTERS); 
do
	for SEED in $(seq 1 $REPEATS)
	do
		[[ $K == 1 ]] && [[ $SEED == 1 ]] && VARNAME=RUNID_$K && JID=$(sbatch --account=nib00015 --dependency=afterok:${!VARNAME} --output=$POP_DIR/logFiles/print_likes.$SET_ID.oe $SCRIPTS_DIR/print_likes.sh $POP_DIR/ngsadmix/$SET_ID.K$K.seed$SEED.log $LIKE_FILE)
		[[ $K != 1 ]] || [[ $SEED != 1 ]] && VARNAME=RUNID_$K && JID=$(sbatch --account=nib00015 --dependency=afterok:${!VARNAME}:${JID##* } --output=$POP_DIR/logFiles/print_likes.$SET_ID.oe $SCRIPTS_DIR/print_likes.sh $POP_DIR/ngsadmix/$SET_ID.K$K.seed$SEED.log $LIKE_FILE)
	done
done

## Plot results
IND_FILE=$POP_DIR/$SET_ID.txt # File with individual IDs in first columns and population assignments in second column

until [[ $(cat $LIKE_FILE | wc -l) == $(( $CLUSTERS*$REPEATS )) ]]
do
	sleep 5
done

sbatch --account=nib00015 --output=$POP_DIR/logFiles/plot_ngsadmix.$SET_ID.oe $SCRIPTS_DIR/plot_ngsadmix.sh $SCRIPTS_DIR $POP_DIR/ngsadmix $LIKE_FILE $IND_FILE $SET_ID

#################################################################
#### 3 ISOLATION-BY-DISTANCE ####
#################################################################

##################################################
#### 3.1 BASED ON F_ST ####
##################################################
mkdir -p $POP_DIR/ibd/fst

COMB_FILE=$PWORK/$SET_ID/angsd/maf/combinations.txt # File with pairwise population comparisons estimated in genotype_likelihoods_sub.sh

## Initialize summary file
echo "pair unweighted weighted" > $POP_DIR/ibd/fst/$SET_ID.fst_sumstats.txt

## Get pairwise F_ST between populations
while read COMB
do
	FIRST=$(awk '{print $1}' <<< $COMB)
	SECOND=$(awk '{print $2}' <<< $COMB)
	
	echo -e "#### Processing combination $FIRST $SECOND ...\n"	
	# Estimate F_ST values
	OUT=$POP_DIR/ibd/fst/$SET_ID.$FIRST.$SECOND
	sbatch --wait --account=nib00015 --output=$POP_DIR/logFiles/$SET_ID.realsfs_fst.$FIRST.$SECOND.oe $SCRIPTS_DIR/realsfs_fst.sh $NT "$PWORK/$SET_ID/angsd/saf.$FIRST.idx $PWORK/$SET_ID/angsd/saf.$SECOND.idx" $OUT
	
	# Write to summary file 
	echo ${FIRST}_$SECOND $(cat $OUT.fst) >> $POP_DIR/ibd/fst/$SET_ID.fst_sumstats.txt
done < $COMB_FILE

## Conduct Mantel tests and plot IBD
GEO_DIST=$POP_DIR/ibd/geo_dist.txt # Distance matrix with mean geographic distances between population pairs as estimated with Geographic Distance Matrix Generator v1.2.3 (https://biodiversityinformatics.amnh.org/open_source/gdmg/), with row and column names
GEN_DIST=$POP_DIR/ibd/fst/gen_dist.txt # Distance matrix with weighted pairwise F_ST between populations as estimated with realSFS, with row and column names
sbatch --account=nib00015 --output=$POP_DIR/logFiles/$SET_ID.ibd.fst.oe $SCRIPTS_DIR/ibd.sh $GEO_DIST $GEN_DIST $POP_DIR/ibd/fst/$SET_IT "_fst"

##################################################
#### 3.1 BASED ON GENETIC DISTANCES ####
##################################################
mkdir -p $POP_DIR/ibd/geneticDistances

VCF_FILE=$PWORK/$SET_ID/vcf/$SET_ID.allScaffolds.snps.07filt.vcf
## Estimate genetic distance between individuals
sbatch --wait --account=nib00015 --output=$POP_DIR/logFiles/$SET_ID.vcfr.oe $SCRIPTS_DIR/vcfr.sh $SCRIPTS_DIR $VCF_FILE $POP_DIR/ibd/geneticDistances/geneticDistances.csv

## Conduct Mantel tests and plot IBD
GEO_DIST=$POP_DIR/ibd/geo_dist.txt # Distance matrix with mean geographic distances between population pairs as estimated with Geographic Distance Matrix Generator v1.2.3 (https://biodiversityinformatics.amnh.org/open_source/gdmg/), with row and column names
GEN_DIST=$POP_DIR/ibd/geneticDistances/gen_dist.txt # Distance matrix with mean genetic distances between populations as estimated with vcfR, with row and column names
GEN_DIST_SD=$POP_DIR/ibd/geneticDistances/gen_dist_sd.txt # Distance matrix with standard deviations of mean genetic distances between populations as estimated with vcfR, with row and column names
sbatch --account=nib00015 --output=$POP_DIR/logFiles/$SET_ID.ibd.geneticDistances.oe $SCRIPTS_DIR/ibd.sh $GEO_DIST $GEN_DIST $GEN_DIST_SD $POP_DIR/ibd/geneticDistances/$SET_IT "_geneticDistances"


#################################################################
#### 4 ESTIMATED EFFECTIVE MIGRATION SURFACES (EEMS) ####
#################################################################
mkdir -p $POP_DIR/eems/results

VCF_FILE=$PWORK/$SET_ID/vcf/$SET_ID.allScaffolds.snps.07filt.vcf

## Estimate average genetic dissimilarity matrix
NT=4
CHROM_FILE=$POP_DIR/eems/renameChromosomes.txt # File with chromosome names in first column and integers in second column (no header), which is required for formating of VCF file
sbatch --wait --account=nib00015 --output=$POP_DIR/logFiles/$SET_ID.bed2diffs.oe $SCRIPTS_DIR/bed2diffs.sh $NT $POP_DIR/eems $VCF_FILE $CHROM_FILE

## Two files need to be created manually:
# $POP_DIR/eems/$(basename $VCF_FILE .vcf).coord which contains coordinates of each indivdidual
# $POP_DIR/eems/$(basename $VCF_FILE .vcf).outer which contains coordinate boundaries of focal area

## Create configuration file
for NDEMES in 200 500 1000
do
	> $POP_DIR/eems/$(basename $VCF_FILE .vcf).ndemes$NDEMES.ini
	echo "datapath = $POP_DIR/eems/$(basename $VCF_FILE .vcf)" >> $POP_DIR/eems/$(basename $VCF_FILE .vcf).ndemes$NDEMES.ini
	echo "mcmcpath = $POP_DIR/eems/results/$(basename $VCF_FILE .vcf).ndemes$NDEMES" >> $POP_DIR/eems/$(basename $VCF_FILE .vcf).ndemes$NDEMES.ini
	echo "nIndiv = $(bcftools query -l $VCF_FILE | wc -l)" >> $POP_DIR/eems/$(basename $VCF_FILE .vcf).ndemes$NDEMES.ini
	echo "nSites = $(egrep -v "#" $VCF_FILE | wc -l)" >> $POP_DIR/eems/$(basename $VCF_FILE .vcf).ndemes$NDEMES.ini
	echo "nDemes = $NDEMES" >> $POP_DIR/eems/$(basename $VCF_FILE .vcf).ndemes$NDEMES.ini
	echo "diploid = true" >> $POP_DIR/eems/$(basename $VCF_FILE .vcf).ndemes$NDEMES.ini
	echo "numMCMCIter = 4000000" >> $POP_DIR/eems/$(basename $VCF_FILE .vcf).ndemes$NDEMES.ini
	echo "numBurnIter = 1000000" >> $POP_DIR/eems/$(basename $VCF_FILE .vcf).ndemes$NDEMES.ini
	echo "numThinIter = 9999" >> $POP_DIR/eems/$(basename $VCF_FILE .vcf).ndemes$NDEMES.ini
done

## Estimate effective migration surfaces and plot
SEED=123
POP_COORDS=$POP_DIR/eems/results/$(basename $VCF_FILE .vcf).ndemes$NDEMES/pop_coords.txt # File with coordinates to plot populations on EEMS
SHAPE=$POP_DIR/eems/results/$(basename $VCF_FILE .vcf).ndemes$NDEMES/River_Mada_1 # Prefix of river shape file
for NDEMES in 200 500 1000
do
	# Infer EEMS
	JID=$(sbatch --account=nib00015 --output=$POP_DIR/logFiles/eems.nDemes$NDEMES.oe $SCRIPTS_DIR/eems.sh $POP_DIR/eems/$(basename $VCF_FILE .vcf).ndemes$NDEMES.ini $SEED)
	
	# Plot results
	sbatch --account=nib00015 --dependency=afterok:${JID##* } --output=$POP_DIR/logFiles/plot_eems.nDemes$NDEMES.oe plot_eems.sh $SCRIPTS_DIR $POP_DIR/eems/$(basename $VCF_FILE .vcf).ndemes$NDEMES $POP_DIR/eems/results/$(basename $VCF_FILE .vcf).ndemes$NDEMES $POP_COORDS $SHAPE_FILE
done



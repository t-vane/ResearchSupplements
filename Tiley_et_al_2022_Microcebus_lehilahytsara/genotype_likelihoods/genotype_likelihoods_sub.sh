################################################################################
#### GENOTYPE LIKELIHOOD INFERENCE ####
################################################################################
SCRIPTS_DIR=/home/nibtve93/scripts/genotypeLikelihoods

SET_ID=lehilahytsara
BAM_DIR=$PWORK/bamFiles/$SET_ID
REFERENCE_DIR=$PWORK/references/mmur3
REFERENCE=$REFERENCE_DIR/GCF_000165445.2_Mmur_3.0_genomic.fna # Reference genome in fasta format
ANGSD_DIR=$PWORK/$SET_ID/angsd
NT=80

mkdir -p $ANGSD_DIR/logFiles


#################################################################
#### 1 ESTIMATE GENOTYPE LIKELIHOODS ACROSS ALL SAMPLES ####
#################################################################
IN_FILE_ALL=$ANGSD_DIR/angsd_all.txt # List of all samples (without file extensions) that shall be included for global genotype likelihood estimation
NO_INDS_ALL=$(cat $IN_FILE_ALL | wc -l)

mkdir -p $ANGSD_DIR/bamHits

## Create bamHits file with locus coverages for each individual
for i in PE SE
do
	IN_FILE=$ANGSD_DIR/angsd_$i.txt # List of samples (without file extensions) that shall be included for global genotype likelihood estimation, separated by sequencing mode (PE or SE)
	NO_INDS=$(cat $IN_FILE | wc -l)

	[[ $i == PE ]] && sbatch --array=1-$NO_INDS --output=$ANGSD_DIR/logFiles/bamHits.$SET_ID.%A_%a.oe $SCRIPTS_DIR/coverage.sh $i $IN_FILE $BAM_DIR $ANGSD_DIR/bamHits
	[[ $i == SE ]] && sbatch --wait --array=1-$NO_INDS --output=$ANGSD_DIR/logFiles/bamHits.$SET_ID.%A_%a.oe $SCRIPTS_DIR/coverage.sh $i $IN_FILE $BAM_DIR $ANGSD_DIR/bamHits
done

## Estimate and plot coverage distributions for each individual
sbatch --wait --output=$ANGSD_DIR/logFiles/cov_plot.$SET_ID.oe $SCRIPTS_DIR/cov_plot.sh $SCRIPTS_DIR $IN_FILE_ALL $ANGSD_DIR/bamHits $SET_ID

## Set thresholds for angsd
MINDEPTHIND=$(cat $ANGSD_DIR/bamHits/statistics/$SET_ID.minmax.txt | cut -d " " -f2 | sort -n | head -1) # Minimum depth per individual
MAXDEPTHIND=$(cat $ANGSD_DIR/bamHits/statistics/$SET_ID.minmax.txt | cut -d " " -f3 | sort -n | tail -1) # Maximum depth per individual
GMIN=$(cat $ANGSD_DIR/bamHits/statistics/$SET_ID.minmax.txt | cut -d " " -f2 | paste -sd+ | bc) # Minimum depth across all individuals
GMAX=$(cat $ANGSD_DIR/bamHits/statistics/$SET_ID.minmax.txt | cut -d " " -f3 | paste -sd+ | bc) # Maximum depth across all individuals
PERCENTAGE="75/100" # Minimum percentage of represented individuals
MININD=$(($NO_INDS_ALL * $PERCENTAGE )) # Minimum number of represented individuals

## Create BAM list as input to angsd
rm $ANGSD_DIR/$SET_ID.bamlist
while read INDV
do
echo $BAM_DIR/$INDV.auto.bam >> $ANGSD_DIR/$SET_ID.bamlist
done < $IN_FILE_ALL

## Run angsd
FILTERS="-setMinDepth $GMIN -setMaxDepth $GMAX -setMaxDepthInd $MAXDEPTHIND -setMinDepthInd $MINDEPTHIND -minInd $MININD -SNP_pval 1e-5 -minQ 20 -minMapQ 20 -minMaf 0.05 -uniqueOnly 1 -remove_bads 1 -skipTriallelic 1 -only_proper_pairs 1 -baq 1 -C 50"
TODO="-GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 1 -doCounts 1"
sbatch --output=$ANGSD_DIR/logFiles/angsd.$SET_ID.oe $SCRIPTS_DIR/angsd.sh $NT $REFERENCE $ANGSD_DIR/$SET_ID.bamlist "$TODO" "$FILTERS" $ANGSD_DIR/$SET_ID


#################################################################
#### 2 ESTIMATE SITE ALLELE FREQUENCY LIKELIHOODS PER POPULATION ####
#################################################################
POPS="riamalandy ambavala ambatovy ambohitantely anjanaharibe anjiahely ankafobe marojejy tsinjoarivo metapopulation"

## Estimate site allele frequency likelihoods per population
for i in $POPS
do
	echo -e "#### Processing population $i ...\n"
	IN_FILE=$ANGSD_DIR/angsd_$i.txt # List of samples in the population
	NO_INDS=$(cat $IN_FILE | wc -l)

	## Estimate and plot coverage distributions for each individual
	sbatch --wait --output=$ANGSD_DIR/logFiles/cov_plot.$i.oe $SCRIPTS_DIR/cov_plot.sh $SCRIPTS_DIR $IN_FILE $ANGSD_DIR/bamHits $i

	## Set thresholds for angsd
	MINDEPTHIND=$(cat $ANGSD_DIR/bamHits/statistics/$i.minmax.txt | cut -d " " -f2 | sort -n | head -1) # Minimum depth per individual
	MAXDEPTHIND=$(cat $ANGSD_DIR/bamHits/statistics/$i.minmax.txt | cut -d " " -f3 | sort -n | tail -1) # Maximum depth per individual
	GMIN=$(cat $ANGSD_DIR/bamHits/statistics/$i.minmax.txt | cut -d " " -f2 | paste -sd+ | bc) # Minimum depth across all individuals
	GMAX=$(cat $ANGSD_DIR/bamHits/statistics/$i.minmax.txt | cut -d " " -f3 | paste -sd+ | bc) # Maximum depth across all individuals
	PERCENTAGE="75/100" # Minimum percentage of represented individuals
	MININD=$(($NO_INDS * $PERCENTAGE )) # Minimum number of represented individuals

	## Create BAM list as input to angsd
	rm $ANGSD_DIR/$i.bamlist
	while read INDV
	do
	echo $BAM_DIR/$INDV.auto.bam >> $ANGSD_DIR/$i.bamlist
	done < $IN_FILE

	## Run angsd
	TODO="-GL 1 -doSaf 1 -doCounts 1 -anc $REFERENCE"
	FILTERS="-setMinDepth $GMIN -setMaxDepth $GMAX -setMaxDepthInd $MAXDEPTHIND -setMinDepthInd $MINDEPTHIND -minInd $MININD -minMapQ 20 -minQ 20 -uniqueOnly 1 -remove_bads 1 -C 50 -baq 1 -only_proper_pairs 1"
	sbatch --output=$ANGSD_DIR/logFiles/angsd.saf.$i.oe $SCRIPTS_DIR/angsd.sh $NT $REFERENCE $ANGSD_DIR/$i.bamlist "$TODO" "$FILTERS" $ANGSD_DIR/$i
done

################################################################################
#### DEMOGRAPHIC MODELLING ####
################################################################################
SCRIPTS_DIR=/home/nibtve93/scripts/demographicModeling

SET_ID=lehilahytsara
DEMO_DIR=$PWORK/$SET_ID/demographicModeling
SAF_DIR=$PWORK/$SET_ID/angsd # Directory with site allele frequency likelihoods for each population inferred in genotype_likelihoods_sub.sh

mkdir -p $DEMO_DIR/logFiles

#################################################################
#### 1 MODEL DEMOGRAPHIC HISTORY ####
#################################################################
POPS="riamalandy ambavala ambatovy ambohitantely anjanaharibe anjiahely ankafobe marojejy tsinjoarivo metapopulation"
NT=24

## Get pairwise population comparisons
awk -f $SCRIPTS_DIR/combinations.awk <<< $POPS > $DEMO_DIR/pop_combinations.txt

## Estimate joint minor allele frequency (MAF) spectra for selected population comparisons
while read combination
do
	FIRST=$(awk '{print $1}' <<< $combination)
	SECOND=$(awk '{print $2}' <<< $combination)
	
	echo -e "#### Processing combination $FIRST $SECOND ...\n"
	
	# Across all sites for demographic modeling
	OUT_FILE=$DEMO_DIR/$SET_ID.$FIRST.$SECOND.ml
	sbatch --output=$DEMO_DIR/logFiles/$SET_ID.realsfs.2d.$FIRST.$SECOND.oe $SCRIPTS_DIR/realsfs.sh $NT "$SAF_DIR/$FIRST.saf.idx $SAF_DIR/$SECOND.saf.idx" $OUT_FILE
	# Over blocks of 10000 bp (with bootstrapping) to generate confidence intervals in demographic modelling
	OPTIONS="-nSites 10000 -bootstrap 100"
	OUT_FILE=$DEMO_DIR/$SET_ID.$FIRST.$SECOND.block.ml
	sbatch --output=$DEMO_DIR/logFiles/$SET_ID.realsfs.2d.block.$FIRST.$SECOND.oe $SCRIPTS_DIR/realsfs.sh $NT "$SAF_DIR/$FIRST.saf.idx $SAF_DIR/$SECOND.saf.idx" $OUT_FILE "$OPTIONS"
done < $DEMO_DIR/pop_combinations.txt

## See Dryad digital repository (https://doi.org/10.5061/dryad.dncjsxkvz) for scripts to assemble block bootstraps and run fastsimcoal v2.6.0.3 (http://cmpg.unibe.ch/software/fastsimcoal27/).


#################################################################
#### 2 CHARACTERIZE POPULATION SIZE CHANGES THROUGH TIME ####
#################################################################
## Software:
STAIRWAYPLOT=/home/nibtve93/software/stairway_plot_v2.1.1/stairway_plot_es # (v2.0; https://github.com/xiaoming-liu/stairway-plot-v2)

POPS="ambavala ambatovy ankafobe metapopulation"
NT=24
MU="1.2e-8" # Mutation rate for study system
GEN_TIME="2.5" # Generation time for study system
SEED=$RANDOM

for i in $POPS
do
	echo -e "#### Processing population $i ...\n"
	
	## Estimate minor allele frequency (MAF) spectra for selected population
	OUT_FILE=$DEMO_DIR/$i.sfs
	sbatch --wait --output=$DEMO_DIR/logFiles/realsfs.1d.$i.oe $SCRIPTS_DIR/realsfs.sh $NT $SAF_DIR/$i.saf.idx $OUT_FILE
	
	## Print blueprint file for Stairway Plot
	IN_FILE=$PWORK/$SET_ID/angsd/angsd_$i.txt # List of samples in the population
	NO_SEQ=$(( $(cat $IN_FILE | wc -l) * 2 )) # Multiplied by 2 because of diploid genome
	SAF_OUT=$PWORK/$SET_ID/angsd/logFiles/angsd.saf.$i.oe
	
	rm $DEMO_DIR/populationSize/$i.blueprint; touch $DEMO_DIR/populationSize/$i.blueprint
	echo -e "popid: $i" >> $DEMO_DIR/populationSize/$i.blueprint
	echo -e "nseq: $NO_SEQ" >> $DEMO_DIR/populationSize/$i.blueprint
	echo -e "L: $(cat $SAF_OUT | grep nSites | sed -n 1p | awk '{ print $3}')" >> $DEMO_DIR/populationSize/$i.blueprint
	echo -e "whether_folded: true" >> $DEMO_DIR/populationSize/$i.blueprint
	echo -e "SFS: $(cut -d' ' -f2-$(( 1 + $NSEQ/2)) $OUT_FILE)" >> $DEMO_DIR/populationSize/$i.blueprint
	echo -e "smallest_size_of_SFS_bin_used_for_estimation: 1" >> $DEMO_DIR/populationSize/$i.blueprint
	echo -e "largest_size_of_SFS_bin_used_for_estimation: $(cat $IN_FILE | wc -l)" >> $DEMO_DIR/populationSize/$i.blueprint
	echo -e "pct_training: 0.67" >> $DEMO_DIR/populationSize/$i.blueprint
	echo -e "nrand: $(( ($NO_SEQ-2)/4 )) $(( ($NO_SEQ-2)/2 )) $(( ($NO_SEQ-2)*3/4 )) $(( $NO_SEQ-2 ))" >> $DEMO_DIR/populationSize/$i.blueprint
	echo -e "project_dir: $DEMO_DIR/populationSize" >> $DEMO_DIR/populationSize/$i.blueprint
	echo -e "stairway_plot_dir: $STAIRWAYPLOT" >> $DEMO_DIR/populationSize/$i.blueprint
	echo -e "ninput: 200" >> $DEMO_DIR/populationSize/$i.blueprint
	echo -e "random_seed: $SEED" >> $DEMO_DIR/populationSize/$i.blueprint
	echo -e "mu: $MU" >> $DEMO_DIR/populationSize/$i.blueprint
	echo -e "year_per_generation: $GEN_TIME" >> $DEMO_DIR/populationSize/$i.blueprint
	echo -e "plot_title: $i" >> $DEMO_DIR/populationSize/$i.blueprint
	echo -e "xrange: 0,0" >> $DEMO_DIR/populationSize/$i.blueprint
	echo -e "yrange: 0,0" >> $DEMO_DIR/populationSize/$i.blueprint
	echo -e "xspacing: 2" >> $DEMO_DIR/populationSize/$i.blueprint
	echo -e "yspacing: 2" >> $DEMO_DIR/populationSize/$i.blueprint
	echo -e "fontsize: 12" >> $DEMO_DIR/populationSize/$i.blueprint
	
	## Submit script to run Stairway Plot
	sbatch --output=$DEMO_DIR/logFiles/stairwayplot.$i.oe $SCRIPTS_DIR/stairwayplot.sh $DEMO_DIR/populationSize/$i.blueprint
done

## Plot Stairway Plot output
Rscript $SCRIPTS_DIR/plot_stairwayplot.R $DEMO_DIR/populationSize/$(echo $POPS | awk '{print $1}').final.summary $DEMO_DIR/populationSize/$(echo $POPS | awk '{print $2}').final.summary $DEMO_DIR/populationSize/$(echo $POPS | awk '{print $3}').final.summary $DEMO_DIR/populationSize/$(echo $POPS | awk '{print $4}').final.summary


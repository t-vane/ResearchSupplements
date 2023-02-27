################################################################################
#### DIVERGENCE DATING ####
################################################################################
## Software:
# spruceup.py needs to be included in $PATH
AMAS=/global/homes/jg/t_vane02/software/AMAS-master/amas/AMAS.py

HOME=/global/homes/jg/t_vane02
SCRIPTS=$HOME/scripts
WD=$HOME/uce-myrmecocystus
TAXON_SET=dating # A configuration file called taxon-set-$TAXON_SET.conf with a list of samples needs to be present in $WD; a reduced taxon set was used for dating compared to phylogenetic inference

#################################################################
#### 1 SETUP
#################################################################

##################################################
### 1.1 EXTRACT UCES
##################################################
LOCUS_DB=uces.sqlite # Name for the database created in UCE_extraction.sh
ASSEMBLY_DIR=$WD/metaspades

## Extract UCEs for desired taxon set
qsub -N UCE_extraction_$TAXON_SET -o $WD/uces/logs -e $WD/uces/logs $SCRIPTS/UCE_extraction.sh $ASSEMBLY_DIR $WD/uces $LOCUS_DB $TAXON_SET

## Copy loci with 90% taxon-completeness (more than 37 samples) to separate directory
mkdir -p $WD/uces/$TAXON_SET/exploded-fastas-90
grep -c ">" $WD/uces/$TAXON_SET/exploded-fastas-all/* | awk -v myvar=$WD/uces/$TAXON_SET/exploded-fastas-90 -F : '$2 > 37 {print ("cp "$1" myvar/"$1)}' | /bin/bash

##################################################
### 7.1.2 ALIGNMENT
##################################################
mkdir -p $WD/uces/$TAXON_SET/mafft-aligned-90

NT=12

## Create list of loci 
ls $WD/uces/$TAXON_SET/exploded-fastas-90 > $WD/uces/$TAXON_SET/loci_$TAXON_SET.txt

## Align and trim loci
qsub -pe smp $NT -t 1-$(cat $WD/$TAXON_SET/loci_$TAXON_SET.txt | wc -l) -N alignment_trimming_$TAXON_SET -o $WD/uces/logs -e $WD/uces/logs $SCRIPTS/alignment_trimming.sh $NT $WD/uces/$TAXON_SET/loci_$TAXON_SET.txt $WD/uces/$TAXON_SET/mafft-aligned-90

## Concatenate alignments
mkdir -p $WD/uces/$TAXON_SET/mafft-aligned-90/concat
python $AMAS concat -i $WD/uces/$TAXON_SET/mafft-aligned-90/*-mafft-trimal -f fasta -d dna -t $WD/uces/$TAXON_SET/mafft-aligned-90/concat/mafft-trimal-concat.fas -p $WD/uces/$TAXON_SET/mafft-aligned-90/concat/mafft-trimal-concat.fas-part

## Trim with spruceup
SPRUCEUP_CONF=$WD/uces/concat/spruceup.conf # Manually created configuration file
spruceup.py $SPRUCEUP_CONF

## Convert final concatenated alignment to PHYLIP format
cd $WD/uces/$TAXON_SET/mafft-aligned-90/concat
python $AMAS convert -i $WD/uces/$TAXON_SET/mafft-aligned-90/concat/mafft-trimal-concat.fas -f fasta -d dna -u phylip

## Calculate summary statistics 
#Calculate summary statistics on the final alignment using AMAS
python $AMAS summary -i $WD/uces/$TAXON_SET/mafft-aligned-90/concat/mafft-trimal-spruceup-concat.fas -f fasta -d dna -o $WD/uces/$TAXON_SET/mafft-aligned-90/concat/mafft-trimal-spruceup-concat-summary.txt

#################################################################
#### 7 DIVERGENCE TIME ESTIMATION
#################################################################
mkdir -p $WD/divergence_dating/logs

## Create directory structure
for TOPOLOGY in concat astral
do
	for CALIBRATION in lasius paratrechinanylanderia
	do
		mkdir -p $WD/divergence_dating/$TOPOLOGY/$CALIBRATION/usedata0
		mkdir -p $WD/divergence_dating/$TOPOLOGY/$CALIBRATION/usedata3
		for i in 1 2
		do
			mkdir -p $WD/divergence_dating/$TOPOLOGY/$CALIBRATION/usedata2_run$i
		done
	done
done

##################################################
### 7.2.1 RUN BASEML
##################################################
mkdir -p $WD/divergence_dating/concatenated $WD/divergence_dating/astral

CONCAT_CTL=$WD/divergence_dating/concat_baseml.ctl # Manually created control file with instructions for baseml
ASTRAL_CTL=$WD/divergence_dating/astral_baseml.ctl # Manually created control file with instructions for baseml

## Submit baseml with topology of ML inference of concatenated alignment
qsub -N baseml_concat -o $WD/divergence_dating/logs -e $WD/divergence_dating/logs $SCRIPTS/baseml.sh $CONCAT_CTL
## Submit baseml with topology of ASTRAL analysis
qsub -N baseml_astral -o $WD/divergence_dating/logs -e $WD/divergence_dating/logs $SCRIPTS/baseml.sh $ASTRAL_CTL

##################################################
### 7.2.2 RUN MCMCTREE
##################################################

## Run MCMCTree with usedata=3 to get estimate of Gradient and Hessian
## Create control files manually following the pattern ${TOPOLOGY}_${CALIBRATION}_usedata3.ctl
for TOPOLOGY in concat astral
do
	for CALIBRATION in lasius paratrechinanylanderia
	do
		qsub -N mcmctree_${TOPOLOGY}_${CALIBRATION}_usedata3 -o $WD/divergence_dating/logs -e $WD/divergence_dating/logs $SCRIPTS/mcmctree.sh ${TOPOLOGY}_${CALIBRATION}_usedata3.ctl
	done
done

## Copy file with Gradient and Hessian to usedata2 directories and run MCMCTree with usedata=2
for TOPOLOGY in concat astral
do
	for CALIBRATION in lasius paratrechinanylanderia
	do 
		for i in 1 2
		do
			cp $WD/divergence_dating/$TOPOLOGY/$CALIBRATION/usedata3/out.BV $WD/divergence_dating/$TOPOLOGY/$CALIBRATION/usedata2_run$i/in.BV
			qsub -N mcmctree_${TOPOLOGY}_${CALIBRATION}_usedata2 -o $WD/divergence_dating/logs -e $WD/divergence_dating/logs $SCRIPTS/mcmctree.sh ${TOPOLOGY}_${CALIBRATION}_usedata2.ctl
		done
	done
done

## Run MCMCTREE with usedata=0 to get prior estimates
for TOPOLOGY in concat astral
do
	for CALIBRATION in lasius paratrechinanylanderia
	do
		qsub -N mcmctree_${TOPOLOGY}_${CALIBRATION}_usedata0 -o $WD/divergence_dating/logs -e $WD/divergence_dating/logs $SCRIPTS/mcmctree.sh ${TOPOLOGY}_${CALIBRATION}_usedata0.ctl
	done
done
################################################################################
#### ALIGNMENT AND TRIMMING ####
################################################################################
## Software:
# spruceup needs to be included in $PATH (https://github.com/marekborowiec/spruceup)
AMAS=/global/homes/jg/t_vane02/software/AMAS-master/amas/AMAS.py # (https://github.com/marekborowiec/AMAS)

HOME=/global/homes/jg/t_vane02
SCRIPTS=$HOME/scripts
TAXON_SET=genus
WD=$HOME/uce-myrmecocystus/uces/$TAXON_SET/alignments

mkdir -p $WD/logs

#################################################################
#### 1 ALIGN AND TRIM SINGLE LOCI ####
#################################################################
NT=12

## Create list of loci 
ls $HOME/uce-myrmecocystus/uces/$TAXON_SET/exploded-fastas-all > $WD/loci.txt
## Align and trim loci
qsub -sync y -pe smp $NT -t 1-$(cat $WD/loci.txt | wc -l) -N alignment_trimming -o $WD/logs -e $WD/logs $SCRIPTS/alignment_trimming.sh $NT $WD/loci.txt $WD

#################################################################
#### 2 CONCATENATE ALIGNMENTS AND TRIM ####
#################################################################
## Concatenate all alignments
mkdir -p $WD/concat
python $AMAS concat -i $WD/*-mafft-trimal -f fasta -d dna -t $WD/concat/mafft-trimal-concat.fas -p $WD/concat/mafft-trimal-concat.fas-part

## Trim with spruceup
SPRUCEUP_CONF=$WD/concat/spruceup.conf # Manually created configuration file for spruceup
spruceup.py $SPRUCEUP_CONF

## Convert final concatenated alignment to NEXUS and PHYLIP format
cd $WD/concat
python $AMAS convert -i $WD/concat/mafft-trimal-concat.fas -f fasta -d dna -u nexus
python $AMAS convert -i $WD/concat/mafft-trimal-concat.fas -f fasta -d dna -u phylip


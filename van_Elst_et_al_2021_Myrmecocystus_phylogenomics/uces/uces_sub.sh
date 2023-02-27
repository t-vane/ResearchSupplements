################################################################################
#### UCE EXTRACTION ####
################################################################################
HOME=/global/homes/jg/t_vane02
SCRIPTS=$HOME/scripts
WD=$HOME/uce-myrmecocystus/uces
ASSEMBLY_DIR=$HOME/uce-myrmecocystus/metaspades
PROBES=$WD/uces/hymenoptera-v2-ANT-SPECIFIC-uce-baits.fasta # UCE probe sequences

mkdir -p $WD/logs

#################################################################
#### 1 HARVEST UCES FROM PUBLISHED GENOMES
#################################################################
NT=16
GENOMES=$WD/uce-harvesting/genomes.txt # List of genome names without file extension

qsub -pe smp $NT -N uce_harvesting -o $WD/logs -e $WD/logs $SCRIPTS/uce_harvesting.sh $NT $WD/uce-harvesting $GENOMES $PROBES $ASSEMBLY_DIR

#################################################################
#### 2 EXTRACT UCES
#################################################################
LOCUS_DB=uces.sqlite # Name for the database created in UCE_extraction.sh
TAXON_SET=genus # A configuration file called taxon-set-$TAXON_SET.conf with a list of samples needs to be present in $WD/uce-extraction

## Match contigs to probes for all samples
ID=$(qsub -N uce_match_contigs -o $WD/logs -e $WD/logs $SCRIPTS/uce_match_contigs.sh $ASSEMBLY_DIR $PROBES $WD/uce-extraction)
## Extract UCEs for desired taxon set
qsub -N uce_extraction -o $WD/logs -e $WD/logs -W depend=afterok:$ID $SCRIPTS/uce_extraction.sh $ASSEMBLY_DIR $WD/uce-extraction $LOCUS_DB $TAXON_SET


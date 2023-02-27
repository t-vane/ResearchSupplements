################################################################################
#### 5 PARTITIONING ####
################################################################################
HOME=/global/homes/jg/t_vane02
SCRIPTS=$HOME/scripts
TAXON_SET=genus
WD=$HOME/uce-myrmecocystus/uces/$TAXON_SET/partitioning

mkdir -p $WD/logs

#################################################################
#### 1 CREATE TRIPLET PARTITIONS FOR EACH LOCUS
#################################################################
ALIGNMENT=$HOME/uce-myrmecocystus/uces/$TAXON_SET/alignments/concat/mafft-trimal-concat.fas-nex.out # NEXUS-formatted alignment
LOCUS_INFO=$WD/locuspartitions.nex # NEXUS file with UCE locus information; can be taken from mafft-trimal-concat.fas-part (see alignment section)

cat $ALIGNMENT $LOCUS_INFO > $WD/alignment-concatenated-partitions.nex

## Run SWSC-EN
qsub -N SWSC-EN -o $WD/logs -e $WD/logs $SCRIPTS/swscen.sh $WD/alignment-concatenated-partitions.nex

#################################################################
#### 2 RUN PARTITIONFINDER
#################################################################
mkdir -p $WD/locus $WD/triplet

## Two partitioning schemes are run (per triplet and per locus)
LOCUS_CONF=$WD/partitionfinder_locus.cfg # Configuration file for PartitionFinder2 considering loci as partitions
TRIPLET_CONF=$WD/partitionfinder_triplet.cfg #Configuration file for PartitionFinder2 considering triplets as partitions

qsub -N partitionfinder_locus -o $WD/logs -e $WD/logs $SCRIPTS/partitionfinder.sh $LOCUS_CONF
qsub -N partitionfinder_triplet -o $WD/logs -e $WD/logs $SCRIPTS/partitionfinder.sh $TRIPLET_CONF

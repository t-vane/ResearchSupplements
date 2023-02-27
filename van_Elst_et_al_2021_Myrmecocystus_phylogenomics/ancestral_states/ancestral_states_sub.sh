################################################################################
#### ANCESTRAL STATE ESTIMATION ####
################################################################################
HOME=/global/homes/jg/t_vane02
SCRIPTS=$HOME/scripts
WD=/global/homes/jg/t_vane02/uce-myrmecocystus/ancestral-states

mkdir -p $WD/logs

## Estimate ancestral states for two alternative topologies (based on concatenated ML inference and ASTRAL species tree inference) 
TREE_CONCAT=$HOME/uce-myrmecocystus/uce-extraction/genus/phylogenetic-inference/ML/run1/partitions.nex.treefile
TREE_ASTRAL=$HOME/uce-myrmecocystus/uce-extraction/genus/phylogenetic-inference//astral/unbinned/speciestree_unbinned.tre
CHARS=$WD/characters.csv # CSV formatted file with two columns named 'species' and 'char'

# Concatenated topology
qsub -N ancestral_states_concat -o $WD/logs -e $WD/logs $SCRIPTS/ancestral_states.sh $TREE_CONCAT $CHARS $WD/concat
# ASTRAL species tree topology
qsub -N ancestral_states_astral -o $WD/logs -e $WD/logs $SCRIPTS/ancestral_states.sh $TREE_ASTRAL $CHARS $WD/astral

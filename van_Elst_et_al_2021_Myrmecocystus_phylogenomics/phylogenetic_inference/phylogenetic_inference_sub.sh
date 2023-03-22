################################################################################
#### PHYLOGENETIC INFERENCE ####
################################################################################
## Software:
# iqtree needs to be included in $PATH (v1.6.11; http://www.iqtree.org/)
AMAS=/global/homes/jg/t_vane02/software/AMAS-master/amas/AMAS.py # (https://github.com/marekborowiec/AMAS)
NWED=/global/homes/jg/t_vane02/software/newick-utils-1.6/src/nw_ed # (https://github.com/tjunier/newick_utils)

## Variables:
HOME=/global/homes/jg/t_vane02
SCRIPTS=$HOME/scripts
TAXON_SET=genus
WD=$HOME/uce-myrmecocystus/uce-extraction/$TAXON_SET/phylogenetic-inference

ALIGNMENT=$HOME/uce-myrmecocystus/uce-extraction/$TAXON_SET/alignments/concat/mafft-trimal-spruceup-concat.fas-nex.out # Concatenated final alignment in NEXUS format
PARTITION_SCHEME=$WD/partitions.nex # File with paragraph 'Nexus formatted character sets for IQtree' from best_scheme.txt of partitioning analysis

mkdir -p $WD/logs

#################################################################
#### 1 MAXIMUM LIKELIHOOD INFERENCE ####
#################################################################
## Submit four independent runs of IQ-TREE to increase probability to find global likelihood maximum
REPEATS=4
for i in $(seq 1 $REPEATS)
do	
	mkdir -p $WD/ML/run$i
	ln -s $ALIGNMENT $WD/ML/run$i
	ln -s $PARTITION_SCHEME $WD/ML/run$i
	qsub -N ML_inference_$i -o $WD/logs -e $WD/logs $SCRIPTS/ML_inference.sh $WD/ML/run$i/mafft-trimal-spruceup-concat.fasta $WD/ML/run$i/partitions.nex
done

## Estimate Robinson-Foulds distances
cat $WD/ML/run*/partitions.nex.treefile > $WD/ML/ml_trees.txt
iqtree -rf_all $WD/ML/ml_trees.txt

#################################################################
#### 2 SPECIES TREE INFERENCE UNDER THE MULTISPECIES COALESCENT ####
#################################################################

##################################################
### 2.1 GENE TREE INFERENCE ####
##################################################
LOCUS_PARTITIONS=$HOME/uce-myrmecocystus/uce-extraction/$TAXON_SET/partitioning/locuspartitions.nex # Contains location of each UCE in concatenated alignment

## Split concatenated alignment into locus alignments
cd $(dirname $ALIGNMENT)
python $AMAS split -i $ALIGNMENT -u nexus -d dna -f fasta -l $LOCUS_PARTITIONS

## Rename the alignment files to uceXXXX.fas
for i in $(dirname $ALIGNMENT)/*-out.fas
do
	rename.ul uce- uce $i
	rename.ul -out.fas .fas $i
done

## Move alignment files to directory for phylogenetic inference
mkdir -p $WD/ML/single-locus
mv $(dirname $ALIGNMENT)/UCE*.fas $WD/ML/single-locus

## Create alignment list
ls $WD/ML/single-locus/UCE*.fas > $WD/ML/single-locus/alignments.txt

## Maximum likelihood inference per locus
qsub -sync y -t 1-$(cat $WD/ML/single-locus/alignments.txt | wc -l) -N gene_tree_inference -o $WD/logs -e $WD/logs $SCRIPTS/gene_tree_inference.sh $WD/ML/single-locus/alignments.txt

##################################################
### 2.2 STATISTICAL BINNING ####
##################################################
mkdir -p $WD/statistical-binning

## Create directory for each alignment
cat $WD/ML/single-locus/alignments.txt | awk -v myvar=$WD/statistical-binning '{print("mkdir "$1" myvar")}' | sed 's/.fas//g' | /bin/bash

## Copy alignments and trees there
for i in $WD/ML/single-locus/UCE*.fas
do 
	cp $i $(dirname $i)/$(basename$i. fas)
	cp $i.treefile $(dirname $i)/$(basename$i. fas)
done

## Run statistical binning
SUPPORT=95
qsub -sync y -N statistical_binning -o $WD/logs -e $WD/logs $SCRIPTS/statistical_binning.sh $WD/ML/single-locus/ $SUPPORT $WD/ML/run1/partitions.nex.treefile

## Identify bins that are comprised of more than one locus (i.e., two or three loci) and run ML inference on them
find $WD/ML/single-locus/output/supergenes -type f -name "bin*" -print0 | xargs -0 wc -l | grep "2 .\|3 ." | sed -e 's/2 .\///g' | sed -e 's/3 .\///g' > $WD/ML/single-locus/output/supergenes/bin_multiple_loci.txt
qsub -sync y -t 1-$(cat $WD/ML/single-locus/output/supergenes/bin_multiple_loci.txt | wc -l) -N gene_tree_inference_supergenes -o $WD/logs -e $WD/logs $SCRIPTS/gene_tree_inference.sh $WD/ML/single-locus/output/supergenes/bin_multiple_loci.txt

##################################################
### 2.3 SPECIES TREE ESTIMATION WITH ASTRAL FOR BINNED AND UNBINNED LOCI ####
##################################################
MAPPING=$WD/astral/mapping.txt # Assigns specimens to species
mkdir -p $WD/astral/unbinned $WD/astral/binned

## Create input trees file for unbinned analysis in ASTRAL
cat $WD/ML/single-locus/UCE*/*.treefile > $WD/astral/unbinned/gene_trees.txt

## Create weighted input trees file for binned analysis in ASTRAL
# Single-locus bins
find $WD/ML/single-locus/output/ -type f -name "bin*" -print0 | xargs -0 wc -l | grep "1 ." | sed -e 's/1 .\///g' | xargs grep "uce" | sed -e 's/bin.\{1,9\}\+.txt://g' | while read line
do 
	cat $WD/ML/single-locus/$line/$line.treefile > $WD/astral/binned/gene_trees.txt
done 
# Dual-locus bins
find $WD/ML/single-locus/output/ -type f -name "bin*" -print0 | xargs -0 wc -l | grep "2 ." | sed -e 's/2 .\///g' | while read line
do 
	for i in {1..2}
	do
		cat $WD/ML/single-locus/output/supergenes/$line.treefile >> $WD/astral/binned/gene_trees.txt
	done
done 
# Triple-locus bins
find $WD/ML/single-locus/output/ -type f -name "bin*" -print0 | xargs -0 wc -l | grep "3 ." | sed -e 's/3 .\///g' | while read line
do 
	for i in {1..3}
	do
		cat $WD/ML/single-locus/output/supergenes/$line.treefile >> $WD/astral/binned/gene_trees.txt
	done
done

## Collapse branches in input trees files with bootstrap support below 20 and submit analysis
for i in binnend unbinned
do
	$NWED $WD/astral/$i/gene_trees.txt 'i & b <=20' o > $WD/astral/$i/gene_trees_bs20.txt
	qsub -N species_tree_astral_$i -o $WD/logs -e $WD/logs $SCRIPTS/species_tree_astral.sh $WD/astral/$i/gene_trees_bs20.txt $MAPPING $WD/astral/$i/speciestree_$i.tre

done

##################################################
### 2.3 SPECIES TREE ESTIMATION WITH SVDquartets ####
##################################################
mkdir -p $WD/svdq

## Create PAUP block file
echo "BEGIN PAUP;" > $WD/svdq/blockfile.txt
echo -e "\toutgroup SPECIES.Nylanderia_terricola SPECIES.Paratrechina_longicornis;" >> $WD/svdq/blockfile.txt
echo -e "\tset root=outgroup outroot=monophyl;" >> $WD/svdq/blockfile.txt
echo -e "\tsvdq nthreads=8 evalQuartets=all taxpartition=SPECIES loci=LOCI bootstrap=multilocus treeFile=svdq_bootstrap.tre;" >> $WD/svdq/blockfile.txt
echo "END;" >> $WD/svdq/blockfile.txt

## Create NEXUS file to run SVDQ in PAUP
CHAR_PARTITIONS=$WD/svdq/char_partitions.txt # Character partitions
TAX_PARTITIONS=$WD/svdq/tax_partitions.txt # Taxon partitions

cat $ALIGNMENT $CHAR_PARTITIONS $TAX_PARTITIONS $WD/svdq/blockfile.txt > $WD/svdq/concat_alignment_paup.nex

## Run SVDQuartets
qsub -N species_tree_svdq -o $WD/logs -e $WD/logs $SCRIPTS/species_tree_svdq.sh $WD/svdq/concat_alignment_paup.nex $WD/svdq/concat_alignment_paup.log

##################################################
### 2.4 CONCORDANCE FACTOR ANALYSIS ####
##################################################
mkdir -p $WD/concordance-factors
NT=16

qsub -N concordance_factors -o $WD/logs -e $WD/logs $SCRIPTS/concordance_factors.sh $NT $WD/ML/run1/partitions.nex.treefile $WD/astral/unbinned/gene_trees.txt $ALIGNMENT $WD/concordance-factors/concord
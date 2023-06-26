################################################################################
#### PHYLOGENETIC INFERENCE ####
################################################################################
SCRIPTS_DIR=/home/nibtve93/scripts/phylogeneticInference

SET_ID=gerpi
VCF_DIR=$PWORK/$SET_ID/vcf
LOCUS_DIR=$PWORK/$SET_ID/locusExtraction/fasta/${SET_ID}_bylocus_final
ALIGNMENT_DIR=$PWORK/$SET_ID/alignment
PHYL_DIR=$PWORK/$SET_ID/phylogeneticInference

mkdir -p $ALIGNMENT_DIR/logFiles
mkdir -p $PHYL_DIR/logFiles

#################################################################
#### 0 ALIGNMENT ####
#################################################################
## Align loci with muscle and change headers to remove locus IDs (otherwise concatenation will not work)
NT=64 # Number of loci processed in parallel
JID1=$(sbatch --account=nib00015 --output=$ALIGNMENT_DIR/logFiles/alignment.oe $SCRIPTS_DIR/alignment.sh $NT $LOCUS_DIR $ALIGNMENT_DIR)

## Calculate statistics for locus alignments
sbatch --account=nib00015 --dependency=afterok:${JID1##* } --output=$ALIGNMENT_DIR/logFiles/locus_statistics.oe $SCRIPTS_DIR/locus_statistics.sh $ALIGNMENT_DIR $ALIGNMENT_DIR/locus.stats

## Concatenate locus alignments
NT=8
JID2=$(sbatch --account=nib00015 --dependency=afterok:${JID1##* } --output=$ALIGNMENT_DIR/logFiles/concatenate.oe $SCRIPTS_DIR/concatenate.sh $NT $ALIGNMENT_DIR $SET_ID)

## Calculate statistics for concatenated alignment (including per-taxon)
sbatch --account=nib00015 --dependency=afterok:${JID2##* } --output=$ALIGNMENT_DIR/logFiles/concatenated_statistics.oe $SCRIPTS_DIR/concatenated_statistics.sh $ALIGNMENT_DIR/$SET_ID.concatenated.nex

#################################################################
#### 1 MAXIMUM LIKELIHOOD INFERENCE ####
#################################################################
mkdir -p $PHYL_DIR/ml

## Remove invariant sites (high memory consumption)
JID3=$(sbatch --account=nib00015 --dependency=afterok:${JID2##* } --output=$PHYL_DIR_DIR/logFiles/remove_invariants.oe $SCRIPTS_DIR/remove_invariants.sh $ALIGNMENT_DIR/$SET_ID.concatenated.phy $ALIGNMENT_DIR/$SET_ID.concatenated.noinv.phy)

## Run phylogenetic inference with ascertainment bias correction in RAxML-NG
NT=80
BOOTSTRAP=100 # Number of bootstrap replicates
OUTGROUP="Mmur_RMR44,Mmur_RMR45,Mmur_RMR49" # Outgroup individuals in alignment file
sbatch --account=nib00015 --dependency=afterok:${JID3##* } --output=$PHYL_DIR/logFiles/raxml_asc.$SET_ID.oe $SCRIPTS_DIR/raxml_asc.sh $NT $BOOTSTRAP $OUTGROUP $ALIGNMENT_DIR/$SET_ID.concatenated.noinv.phy $PHYL_DIR/ml/$SET_ID.concatenated.noinv.tre

#################################################################
#### 2 QUARTET-BASED INFERENCE FOR INDIVIDUAL AND POPULATION ASSIGNMENT ####
#################################################################
mkdir -p $PHYL_DIR/quartet

## Wait until alignment concatenation has been terminated (i.e., the following job has started)
until [ -f $ALIGNMENT_DIR/logFiles/concatenated_statistics.oe ]
do
	sleep 20
done

## Create locus partitions block file
echo "BEGIN SETS;" > $PHYL_DIR/quartet/$SET_ID.locusPartitions.txt
echo -e "\t CHARPARTITION LOCI =" >> $PHYL_DIR/quartet/$SET_ID.locusPartitions.txt
while read LINE
do
	LOC_NAME=$(cut -f1 -d' ' <<< $LINE)
	LOC_COORD=$(cut -f3 -d' ' <<< $LINE)
	echo -e "#### Processing locus $LOC_NAME ..."
	echo -e "\t\t$LOC_NAME:$LOC_COORD," >> $PHYL_DIR/quartet/$SET_ID.locusPartitions.txt
done < $ALIGNMENT_DIR/$SET_ID.partitions.txt
sed -i '$ s/,$//' $PHYL_DIR/quartet/$SET_ID.locusPartitions.txt # Remove last comma
echo -e "\t\t;" >> $PHYL_DIR/quartet/$SET_ID.locusPartitions.txt
echo "END;" >> $PHYL_DIR/quartet/$SET_ID.locusPartitions.txt
echo "" >> $PHYL_DIR/quartet/$SET_ID.locusPartitions.txt

## Create taxon partitions block files
# For population assignment, the file has to be created manually
# For individual assignment:
echo "BEGIN SETS;" > $PHYL_DIR/quartet/$SET_IT.taxPartitions.individual.nex
echo -e "\t TAXPARTITION SPECIES =" >> echo "BEGIN SETS;" > $PHYL_DIR/quartet/$SET_IT.taxPartitions.individual.nex
for IND in $(awk '{ print $1 }' $ALIGNMENT_DIR/$SET_ID.concatenated.noinv.phy | tail -n+2)
do
	echo -e "\t\t$IND:$IND," >> echo "BEGIN SETS;" > $PHYL_DIR/quartet/$SET_IT.taxPartitions.individual.nex
done
sed -i '$ s/,$//' echo "BEGIN SETS;" > $PHYL_DIR/quartet/$SET_IT.taxPartitions.individual.nex # Remove last comma
echo -e "\t\t;" >> echo "BEGIN SETS;" > $PHYL_DIR/quartet/$SET_IT.taxPartitions.individual.nex
echo "END;" >> echo "BEGIN SETS;" > $PHYL_DIR/quartet/$SET_IT.taxPartitions.individual.nex
echo "" >> echo "BEGIN SETS;" > $PHYL_DIR/quartet/$SET_IT.taxPartitions.individual.nex

## Create PAUP block files
NT=80
SEED=$RANDOM
# For population assignment:
echo "BEGIN PAUP;" > $PHYL_DIR/quartet/$SET_IT.paup.population.nex
echo -e "\toutgroup SPECIES.Mmurinus;" >> $PHYL_DIR/quartet/$SET_IT.paup.population.nex
echo -e "\tset root=outgroup outroot=monophyl;" >> $PHYL_DIR/quartet/$SET_IT.paup.population.nex
echo -e "\tsvdq nthreads=$NT evalQuartets=all taxpartition=SPECIES loci=LOCI bootstrap=standard seed=${SEED};" >> $PHYL_DIR/quartet/$SET_IT.paup.population.nex
echo -e "\tsavetrees format=Newick file=$PHYL_DIR/quartet/$SET_ID.concatenated.noinv.population.svdq.tre savebootp=nodelabels;" >> $PHYL_DIR/quartet/$SET_IT.paup.population.nex
echo -e "\tquit;" >> $PHYL_DIR/quartet/$SET_IT.paup.population.nex
echo "END;" >> $PHYL_DIR/quartet/$SET_IT.paup.population.nex
# For individual assignment
echo "BEGIN PAUP;" > $PHYL_DIR/quartet/$SET_IT.paup.individual.nex
echo -e "\toutgroup SPECIES.Mmur_RMR44 SPECIES.Mmur_RMR45 SPECIES.Mmur_RMR49;" >> $PHYL_DIR/quartet/$SET_IT.paup.individual.nex
echo -e "\tset root=outgroup outroot=monophyl;" >> $PHYL_DIR/quartet/$SET_IT.paup.individual.nex
echo -e "\tsvdq nthreads=$NT evalQuartets=all taxpartition=SPECIES loci=LOCI bootstrap=standard seed=${SEED};" >> $PHYL_DIR/quartet/$SET_IT.paup.individual.nex
echo -e "\tsavetrees format=Newick file=$PHYL_DIR/quartet/$SET_ID.concatenated.noinv.individual.svdq.tre savebootp=nodelabels;" >> $PHYL_DIR/quartet/$SET_IT.paup.individual.nex
echo -e "\tquit;" >> $PHYL_DIR/quartet/$SET_IT.paup.individual.nex
echo "END;" >> $PHYL_DIR/quartet/$SET_IT.paup.individual.nex

## Concatenate files and submit SVDquartets job
for i in population individual
do
	cat $ALIGNMENT_DIR/$SET_ID.concatenated.nex $PHYL_DIR/quartet/$SET_IT.taxPartitions.$i.nex $PHYL_DIR/quartet/$SET_IT.paup.$i.nex > $PHYL_DIR/quartet/$SET_IT.paup.$i.concat.nex
	sbatch --account=nib00015 --output=$PHYL_DIR/logFiles/svdq.$SET_ID.oe $SCRIPTS_DIR/svdq.sh $PHYL_DIR/quartet/$SET_IT.paup.$i.concat.nex $PHYL_DIR/quartet/$SET_IT.paup.$i.concat.nex.log
done

#################################################################
#### 3 NEIGHBOR-NET INFERENCE IN SPLITSTREE ####
#################################################################
mkdir -p $PHYL_DIR/neighbor-net

## Prepare NEXUS file which can be opened in the GUI version of SplitsTree to obtain the Neighbor-Net network
sbatch --account=nib00015 --output=$PHYL_DIR/logFiles/splitstree.oe $SCRIPTS_DIR/splitstree.sh $ALIGNMENT_DIR/$SET_ID.concatenated.nex $PHYL_DIR/neighbor-net/splitstree.out.nex

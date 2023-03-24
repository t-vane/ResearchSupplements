################################################################################
#### PHYLOGENETIC INFERENCE ####
################################################################################

SCRIPTS_DIR=/home/nibtve93/scripts/phylogeneticInference

SET_ID=lehilahytsara
PHYL_DIR=$PWORK/$SET_ID/phylogeneticInference
VCF_IN=$PWORK/$SET_ID/vcf/populations.snps.07filt.vcf 

mkdir -p $PHYL_DIR/logFiles

#################################################################
#### 1 MAXIMUM LIKELIHOOD INFERENCE WITH ASCERTAINMENT BIAS CORRECTION ####
#################################################################
mkdir -p $PHYL_DIR/ml

## Convert VCF file to PHYLIP format and calculate basic alignment statistics
FORMAT=phylip # Output format of alignment (phylip or nexus)
sbatch --job-name=ml_inference --output=$PHYL_DIR/logFiles/vcf_convert.$FORMAT.$SET_ID.oe $SCRIPTS_DIR/vcf_convert.sh $SCRIPTS_DIR $VCF_IN $PHYL_DIR/ml $FORMAT

## Run phylogenetic inference with ascertainment bias correction in RAxML-NG
NT=80
BOOTSTRAP=100 # Number of bootstrap replicates
OUTGROUP="Mmur_RMR44,Mmur_RMR45,Mmur_RMR49" # Outgroup individuals in alignment file
sbatch --job-name=ml_inference --dependency=singleton --output=$PHYL_DIR/logFiles/raxml_asc.$SET_ID.oe $SCRIPTS_DIR/raxml_asc.sh $NT $BOOTSTRAP $OUTGROUP $PHYL_DIR/ml/populations.snps.07filt.noinv.phy $PHYL_DIR/ml/populations.snps.07filt.noinv.tre

#################################################################
#### 2 QUARTET-BASED INFERENCE FOR INDIVIDUAL AND POPULATION ASSIGNMENT ####
#################################################################
mkdir -p $PHYL_DIR/quartet

## Thin VCF file by 10,000 bp intervals to ensure independence of SNPs
vcftools --vcf $VCF_IN --thin 10000 --recode --recode-INFO-all --stdout > $(dirname VCF_IN)/$(basename $VCF_IN .vcf).thin10k.vcf

## Convert VCF file to NEXUS format and calculate basic alignment statistics
FORMAT=nexus # Output format of alignment (phylip or nexus)
sbatch --wait --output=$PHYL_DIR/logFiles/vcf_convert.$FORMAT.$SET_ID.oe $SCRIPTS_DIR/vcf_convert.sh $SCRIPTS_DIR $VCF_IN $PHYL_DIR/quartet $FORMAT

## Create taxon partitions block files
# For population assignment, the file has to be created manually
# For individual assignment:
echo "BEGIN SETS;" > $PHYL_DIR/quartet/$SET_IT.taxPartitions.individual.nex
echo -e "\t TAXPARTITION SPECIES =" >> echo "BEGIN SETS;" > $PHYL_DIR/quartet/$SET_IT.taxPartitions.individual.nex
for i in $(awk '{ print $1 }' $PHYL_DIR/ml/populations.snps.07filt.noinv.phy | tail -n+2)
do
	echo -e "\t\t${i}:${i}," >> echo "BEGIN SETS;" > $PHYL_DIR/quartet/$SET_IT.taxPartitions.individual.nex
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
echo -e "\tsvdq nthreads=$NT evalQuartets=all taxpartition=SPECIES bootstrap=standard seed=${SEED};" >> $PHYL_DIR/quartet/$SET_IT.paup.population.nex
echo -e "\tsavetrees format=Newick file=$PHYL_DIR/quartet/$SET_ID.populations.snps.07filt.population.svdq.tre savebootp=nodelabels;" >> $PHYL_DIR/quartet/$SET_IT.paup.population.nex
echo -e "\tquit;" >> $PHYL_DIR/quartet/$SET_IT.paup.population.nex
echo "END;" >> $PHYL_DIR/quartet/$SET_IT.paup.population.nex
# For individual assignment
echo "BEGIN PAUP;" > $PHYL_DIR/quartet/$SET_IT.paup.individual.nex
echo -e "\toutgroup SPECIES.Mmur_RMR44 SPECIES.Mmur_RMR45 SPECIES.Mmur_RMR49;" >> $PHYL_DIR/quartet/$SET_IT.paup.individual.nex
echo -e "\tset root=outgroup outroot=monophyl;" >> $PHYL_DIR/quartet/$SET_IT.paup.individual.nex
echo -e "\tsvdq nthreads=$NT evalQuartets=all taxpartition=SPECIES bootstrap=standard seed=${SEED};" >> $PHYL_DIR/quartet/$SET_IT.paup.individual.nex
echo -e "\tsavetrees format=Newick file=$PHYL_DIR/quartet/$SET_ID.populations.snps.07filt.individual.svdq.tre savebootp=nodelabels;" >> $PHYL_DIR/quartet/$SET_IT.paup.individual.nex
echo -e "\tquit;" >> $PHYL_DIR/quartet/$SET_IT.paup.individual.nex
echo "END;" >> $PHYL_DIR/quartet/$SET_IT.paup.individual.nex

## Concatenate files and submit SVDquartets job
for i in population individual
do
	cat $PHYL_DIR/quartet/populations.snps.07filt.thin10k.noinv.nex $PHYL_DIR/quartet/$SET_IT.taxPartitions.$i.nex $PHYL_DIR/quartet/$SET_IT.paup.$i.nex > $PHYL_DIR/quartet/$SET_IT.paup.$i.concat.nex
	sbatch --job-name=quartet_inference --output=$PHYL_DIR/logFiles/svdq.$SET_ID.oe $SCRIPTS_DIR/svdq.sh $PHYL_DIR/quartet/$SET_IT.paup.$i.concat.nex $PHYL_DIR/quartet/$SET_IT.paup.$i.concat.nex.log
done

#################################################################
#### 3 APPROXIMATELY UNBIASED (AU) TEST ####
#################################################################
mkdir -p $PHYL_DIR/au_test

## Infer phylogenies under specified constraints
NT=80
BOOTSTRAP=0 # Number of bootstrap replicates
OUTGROUP="Mmur_RMR44,Mmur_RMR45,Mmur_RMR49" # Outgroup individuals in alignment file
CONSTRAINTS="ambatovy.anjz.monophyl ambatovy.monophyl east.west north.south snapp snapp.polytomy svdq svdq.polytomy" # String with prefixes of files containing constrained trees for AU test (format $PREFIX.constraint.nwk)
for i in $CONSTRAINTS
do
	CONSTRAINT_FILE=$i.constraint.nwk
	COMMAND="-g $CONSTRAINT_FILE"
	sbatch --dependency=afterok:$(squeue --noheader --format %i --name ml_inference):$(squeue --noheader --format %i --name quartet_inference) --output=$PHYL_DIR/logFiles/raxml_asc.$SET_ID.constraint.$i.oe $SCRIPTS_DIR/raxml_asc.sh \
		$NT $BOOTSTRAP $OUTGROUP $PHYL_DIR/ml/populations.snps.07filt.noinv.phy $PHYL_DIR/au_test/populations.snps.07filt.noinv.$i.tre "$COMMAND"
done

## Concatenate tree constrained and unconstrained tree files
cat $PHYL_DIR/au_test/populations.snps.07filt.noinv.*.tree.*support $PHYL_DIR/ml/populations.snps.07filt.noinv.tre.*support $PHYL_DIR/quartet/$SET_ID.populations.snps.07filt.individual.svdq.tre > $PHYL_DIR/au_test/au_test.trees

## Perform AU test
NT=40
sbatch --output=$PHYL_DIR/logFiles/au_test.$SET_ID.oe $SCRIPTS_DIR/au_test.sh $NT $PHYL_DIR/ml/populations.snps.07filt.noinv.phy $PHYL_DIR/au_test/au_test.trees $PHYL_DIR/au_test/$SET_ID.au


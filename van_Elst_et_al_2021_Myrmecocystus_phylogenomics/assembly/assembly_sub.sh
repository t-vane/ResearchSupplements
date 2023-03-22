################################################################################
#### ASSEMBLY ####
################################################################################
## Software:
# phyluce needs to be included in $PATH (v1.6.7; https://phyluce.readthedocs.io/en/latest/)

HOME=/global/homes/jg/t_vane02
SCRIPTS=$HOME/scripts
WD=$HOME/uce-myrmecocystus/metaspades

mkdir -p $WD/logs

NT=8
MEM=30
INDS=$WD/individuals.txt # List of individuals for assembly
INDS_DOWNSAMPLE=$WD/individuals_downsample.txt # List of individuals for which reads need to be downsampled before assembly

## Assembly for individuals for which downsampling has to be done
qsub -pe smp $NT -l h_vmem=${MEM}G -t 1-$(cat $INDS_DOWNSAMPLE | wc -l) -N assembly_downsample -o $WD/logs -e $WD/logs $SCRIPTS/assembly.sh $NT $MEM $INDS $PROJECT/reads/clean-fastq $WD TRUE 1500000
## Assembly for individuals without downsampling
qsub -sync y -pe smp $NT -l h_vmem=${MEM}G -t 1-$(cat $INDS | wc -l) -N assembly -o $WD/logs -e $WD/logs $SCRIPTS/assembly.sh $NT $MEM $INDS $PROJECT/reads/clean-fastq $WD FALSE

## Calculate basic summary statistics
for FILE in $WD/*.fasta
do
	echo -e "#### processing file $FILE ... \n"
	phyluce_assembly_get_fasta_lengths --input $FILE --csv
done > $WD/fasta-lengths-metaspades.txt

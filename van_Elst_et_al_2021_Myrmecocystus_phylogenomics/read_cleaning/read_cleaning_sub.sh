################################################################################
#### READ CLEANING ####
################################################################################
HOME=/global/homes/jg/t_vane02
SCRIPTS=$HOME/scripts
WD=$HOME/uce-myrmecocystus/reads

mkdir -p $WD/logs

NT=16
CONF=$WD/illumiprocessor.conf # Configuration/adapter file for newly generated reads
CONF_SRA=$WD/illumiprocessor_SRA.conf # Configuration/adapter file for reads downloaded from SRA
R1_PATTERN=_1
R2_PATTERN=_2

## Read cleaning for newly generated reads
qsub -pe smp $NT -N read_cleaning -o $WD/logs -e $WD/logs  $SCRIPTS/read_cleaning.sh $NT $WD/raw-fastq $WD/clean-fastq $CONF $R1_PATTERN $R2_PATTERN
## Read cleaning for SRA reads
qsub -pe smp $NT -N read_cleaning_SRA -o $WD/logs -e $WD/logs $SCRIPTS/read_cleaning.sh $NT $WD/raw-fastq-SRA $WD/clean-fastq $CONF_SRA $R1_PATTERN $R2_PATTERN



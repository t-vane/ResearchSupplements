#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
# gatk needs to be included in $PATH (v4.1.9.0; https://gatk.broadinstitute.org/hc/en-us)

## Command-line args:
NT=$1
REFERENCE=$2
IND_FILE=$3
REGION_FILE=$4
GVCF_DIR=$5
DB_DIR=$6
VCF_SCAFFOLD_DIR=$7
TMP_DIR=$8

SCAFFOLD_NAME=$(sed -n "$SLURM_ARRAY_TASK_ID"p $REGION_FILE | cut -f 1)
SCAFFOLD_START=$(sed -n "$SLURM_ARRAY_TASK_ID"p $REGION_FILE | cut -f 2)
SCAFFOLD_END=$(sed -n "$SLURM_ARRAY_TASK_ID"p $REGION_FILE | cut -f 3)

INDS_COMMAND=$(for IND in `cat $IND_FILE`; do printf " --variant $GVCF_DIR/$IND.rawvariants.g.vcf.gz"; done)

## Activate conda environment
conda activate java

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### joint_genotyping.sh: Starting script."
echo -e "#### joint_genotyping.sh: Number of threads: $NT"
echo -e "#### joint_genotyping.sh: Reference genome: $REFERENCE"
echo -e "#### joint_genotyping.sh: List with individuals: $IND_FILE"
echo -e "#### joint_genotyping.sh: List with regions and coordinates: $REGION_FILE"
echo -e "#### joint_genotyping.sh: Directory with GVCF files: $GVCF_DIR"
echo -e "#### joint_genotyping.sh: Database directory: $DB_DIR"
echo -e "#### joint_genotyping.sh: Output directory for per-scaffold VCF files: $VCF_SCAFFOLD_DIR"
echo -e "#### joint_genotyping.sh: Temporary directory for running GenotypeGVCFs: $TMP_DIR"
echo -e "#### joint_genotyping.sh: Scaffold name: $SCAFFOLD_NAME"
echo -e "#### joint_genotyping.sh: Scaffold start coordinate: $SCAFFOLD_START"
echo -e "#### joint_genotyping.sh: Scaffold end coordinate: $SCAFFOLD_END \n\n"

################################################################################
#### CONDUCT JOINT GENOTYPING ####
################################################################################
echo -e "#### joint_genotyping.sh: Creating GenomicsDB data storage for scaffold $SCAFFOLD_NAME ...\n"
gatk GenomicsDBImport $INDS_COMMAND --genomicsdb-shared-posixfs-optimizations --genomicsdb-workspace-path $DB_DIR/$SCAFFOLD_NAME --batch-size 0 -L $SCAFFOLD_NAME:$SCAFFOLD_START-$SCAFFOLD_END --reader-threads $NT --interval-padding 100

echo -e "#### joint_genotyping.sh: Conducting joint genotyping for scaffold $SCAFFOLD_NAME ...\n"
gatk GenotypeGVCFs -R $REFERENCE -V "gendb://$DB_DIR/$SCAFFOLD_NAME" -O $VCF_SCAFFOLD_DIR/$SCAFFOLD_NAME.vcf.gz --tmp-dir $TMP_DIR --genomicsdb-shared-posixfs-optimizations

## Report:
echo -e "\n#### joint_genotyping.sh: Done with script."
date


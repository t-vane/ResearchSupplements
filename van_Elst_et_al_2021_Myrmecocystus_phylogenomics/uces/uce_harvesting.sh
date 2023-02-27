#!/bin/bash

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
# ucsc-fatotwobit needs to be included in $PATH
# ucsc-twoBitInfo needs to be included in $PATH
# phyluce scripts (https://phyluce.readthedocs.io/en/latest/) need to be included in $PATH

## Command-line args:
NT=$1
IN_DIR=$2
GENOMES=$3
PROBES=$4
CONTIG_DIR=$5

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### UCE_harvesting.sh: Starting script."
echo -e "#### UCE_harvesting.sh: Number of threads: $NT"
echo -e "#### UCE_harvesting.sh: Input directory: $IN_DIR"
echo -e "#### UCE_harvesting.sh: File with list of genomes: $GENOMES"
echo -e "#### UCE_harvesting.sh: Probe file: $PROBES"
echo -e "#### UCE_harvesting.sh: Contig directory: $CONTIG_DIR \n\n"

################################################################################
#### CONVERT AND CALCULATE SUMMARY STATISTICS####
################################################################################
echo -e "#### UCE_harvesting.sh: Unzipping, converting to 2bit format and calculating summary statistics ... \n"
for i in $(cat $GENOMES)
do
	gunzip $IN_DIR/$i/${i}_genomic.fna.gz
	faToTwoBit $IN_DIR/$i/${i}_genomic.fna. $IN_DIR/$i/$i.2bit
	twoBitInfo $IN_DIR/$i/$i.2bit $IN_DIR/$i/$i.sizes.tab
done

################################################################################
#### HARVEST UCES FROM GENOMES####
################################################################################
cd $IN_DIR

echo -e "#### UCE_harvesting.sh: Creating database ... \n"
GENOMES_STRING=$(awk '$1=$1' ORS=' ' $GENOMES)
phyluce_probe_run_multiple_lastzs_sqlite --cores $NT --identity 75 --coverage 77.5 --db harvesting.sqlite --output harvesting-lastz --scaffoldlist $GENOMES_STRING --genome-base-path ./ --probefile $PROBES

echo -e "#### UCE_harvesting.sh: Slice sequences from genomes ... \n"
for i in $(cat $GENOMES)
do
	echo -e "## UCE_harvesting.sh: Processing genome $i... \n"
	printf "[scaffolds]\n$i:$IN_DIR/$i/$i.2bit" > $i.conf
	phyluce_probe_slice_sequence_from_genomes --lastz outgroups-lastz --conf $i.conf --flank 400 --name-pattern $(basename $PROBES)_v_{}.lastz.clean --output $i-genome-fasta
	ln -s $i-genome-fasta/$i.fasta $CONTIG_DIR/$i.contigs.fasta
done

## Report:
echo -e "\n#### UCE_harvesting.sh: Done with script."
date

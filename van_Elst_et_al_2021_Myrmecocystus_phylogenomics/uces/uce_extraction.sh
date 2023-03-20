#!/bin/bash

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
# phyluce needs to be included in $PATH (v1.6.7; https://phyluce.readthedocs.io/en/latest/)

## Command-line args:
CONTIG_DIR=$1
OUT_DIR=$2
LOCUS_DB=$3
TAXON_SET=$4

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### UCE_extractions.sh: Starting script."
echo -e "#### UCE_extractions.sh: Contig directory: $CONTIG_DIR"
echo -e "#### UCE_extractions.sh: Output directory: $OUT_DIR"
echo -e "#### UCE_extractions.sh: Locus database: $LOCUS_DB"
echo -e "#### UCE_extractions.sh: Taxon set: $TAXON_SET \n\n"

################################################################################
#### EXTRACT UCES FROM CONTIGS####
################################################################################
mkdir -p $OUT_DIR/$TAXON_SET
cd $OUT_DIR/$TAXON_SET

echo -e "#### UCE_extractions.sh: Creating monolithic FASTA file with all loci from all taxa ... \n"
phyluce_assembly_get_match_counts --locus-db $LOCUS_DB --taxon-list-config $OUT_DIR/taxon-set-$TAXON_SET.conf --taxon-group $TAXON_SET --incomplete-matrix --output $OUT_DIR/$TAXON_SET/$TAXON_SET-taxa-incomplete.conf
phyluce_assembly_get_fastas_from_match_counts --contigs $CONTIG_DIR --locus-db $LOCUS_DB --match-count-output $OUT_DIR/$TAXON_SET/$TAXON_SET-taxa-incomplete.conf --output $OUT_DIR/$TAXON_SET/$TAXON_SET-taxa-incomplete.fasta --incomplete-matrix $OUT_DIR/$TAXON_SET/$TAXON_SET-taxa-incomplete.incomplete

echo -e "#### UCE_extractions.sh: Exploding monolithic FASTA file into one file per taxon ... \n"
phyluce_assembly_explode_get_fastas_file --input $OUT_DIR/$TAXON_SET/$TAXON_SET-taxa-incomplete.fasta --output $OUT_DIR/$TAXON_SET/exploded-fastas-all
cd $OUT_DIR/$TAXON_SET/exploded-fasta-all
for i in $OUT_DIR/$TAXON_SET/exploded-fastas-all/*
do 
	sed -r -i 's/>uce-[0-9]+_/>/g;s/ \|uce-[0-9]+//g' $i
done

## Report:
echo -e "\n#### UCE_extractions.sh: Done with script."
date

#!/usr/bin/env Rscript

cat("\n#### hierfstat.R: Starting script.\n")

## Command-line arguments
options(warn = 1)
args <- commandArgs(trailingOnly = TRUE)
vcf_file <- args[1]
sample_file <- args[2]
out <- args[3]

## Report:
cat("#### hierfstat.R: Geographic distance matrix:", geo_dist, "\n")
cat("#### hierfstat.R: Genetic distance matrix:", gen_dist, "\n")
cat("#### hierfstat.R: Output prefix:", out, "\n\n")

## Packages
library(pegas)
library(adegenet)
library(hierfstat)

## Process command-line arguments
cat("#### hierfstat.R: Reading files ... \n")
vcf <- read.vcf(vcf_file, which.loci=1:1000000)
samples <- read.table(sample_file)
pops <- samples[,2]

## Convert
cat("#### hierfstat.R: Converting to hierfstat format ... \n")
vcf_conv <- genind2hierfstat(loci2genind(vcf), pop=pops)

## Estimate pairwise F_ST between populations
cat("#### hierfstat.R: Estimating pairwise F_ST between populations ... \n")
results <- pairwise.WCfst(vcf_conv, diploid=TRUE)

## Write output
cat("#### hierfstat.R: Writing output ... \n")
write.table(results, file=paste0(out, ".fst.txt"))

## Report:
cat("\n#### hierfstat.R: Done with script.\n")


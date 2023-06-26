#!/usr/bin/env Rscript

cat("\n#### vcfr.R: Starting script.\n")

## Command-line arguments
options(warn = 1)
args <- commandArgs(trailingOnly = TRUE)
vcf <- args[1]
out <- args[2]

## Report:
cat("#### vcfr.R: VCF file:", vcf, "\n")
cat("#### vcfr.R: Output file:", out, "\n\n")

## Packages
library(vcfR)

## Process command-line arguments
cat("#### vcfr.R: Processing and transforming input ... \n")
vcf <- read.vcfR(vcf)
x <- vcfR2genlight(vcf)

## Estimate genetic distances
cat("#### vcfr.R: Estimating genetic distances ... \n")
x.dist <- dist(x)

## Write output
cat("#### vcfr.R: Writing output ... \n")
mydf <- as.data.frame(as.matrix(x.dist))
write.csv(mydf, file = out)

## Report:
cat("\n#### vcfr.R: Done with script.\n")


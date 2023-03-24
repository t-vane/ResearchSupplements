#!/usr/bin/env Rscript

cat("\n#### amova.R: Starting script.\n")

## Command-line arguments
options(warn = 1)
args <- commandArgs(trailingOnly = TRUE)
vcf_file <- args[1]
sample_file <- args[2]
out <- args[3]

## Report:
cat("#### amova.R: Working directory:", WD, "\n")
cat("#### amova.R: Geographic distance matrix:", geo_dist, "\n")
cat("#### amova.R: Genetic distance matrix:", gen_dist, "\n")
cat("#### amova.R: Output prefix:", out, "\n\n")

## Packages
library(pegas)
library(adegenet)
library(hierfstat)
library(poppr)

## Process command-line arguments
cat("#### amova.R: Reading VCF file ... \n")
vcf <- read.vcf(vcf_file, which.loci=1:1000000)
samples <- read.table(sample_file)
samples[,2] <- as.factor(samples[,2])
samples[,3] <- as.factor(samples[,3])

## Convert
cat("#### amova.R: Converting ... \n")
genind <- loci2genind(vcf)
pop(genind) <- samples[,2]
strata(genind) <- samples[,-1]
genclone <- as.genclone(genind)

## Run AMOVA and write output
cat("#### amova.R: Conducting AMOVA ... \n")
amova.results <- poppr.amova(genclone, ~V3/V2, cutoff=0.3)
capture.output(amova.results, file = paste0(out, ".amova.txt"))

amova.test.results <- randtest(amova.results)
capture.output(amova.test.results, file = paste0(out, ".amova.test.txt"))

pdf(paste0(out, ".amova.plot.pdf"))
plot(amova.test.results)
dev.off()

## Report:
cat("\n#### amova.R: Done with script.\n")














#read in hierarchy file
#header: ind pop cluster
#1st column: sample names, 2nd column: individual numbers; 3rd column: population numbers; 4th column cluster numbers)
hierarchy=read.table(samplefile, header = T)



#assign populations to genind object
#convert before to numbers with as.numeric(as.factor())
pop(x2) <- hierarchy[,2]
#assign analysis strata to genind object
strata(x2) <- hierarchy[,-1]
#convert to genclone object
x3 <- as.genclone(x2, strata=hierarchy)

#run actual AMOVA
amova.results <- poppr.amova(x3, ~V3/V2, cutoff=0.3)
capture.output(amova.results, file = paste0(out, ".amova.txt"))

amova.test.results <- randtest(amova.results)
capture.output(amova.test.results, file = paste0(out, ".amova.test.txt"))

pdf(paste0(out, ".plot.pdf"))
plot(amova.test.results)
dev.off()

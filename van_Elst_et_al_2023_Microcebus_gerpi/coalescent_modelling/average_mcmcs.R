#!/usr/bin/env Rscript

cat("\n#### average_mcmcs.R: Starting script.\n")

## Command-line arguments
options(warn = 1)
args <- commandArgs(trailingOnly = TRUE)
mcmc1 <- args[1]
mcmc2 <- args[2]
mcmc3 <- args[3]
mcmc4 <- args[4]
out_file <- args[5]

## Report
cat("\n#### average_mcmcs.R: First MCMC:", mcmc1, "\n")
cat("#### average_mcmcs.R: Second MCMC:", mcmc2, "\n")
cat("#### average_mcmcs.R: Third MCMC:", mcmc3, "\n")
cat("#### average_mcmcs.R: Fourth MCMC:", mcmc4, "\n\n")

## Process command-line arguments
cat("#### average_mcmcs.R: Reading MCMC files ... \n")
x1 <- read.table(mcmc1,header=TRUE,sep="\t")
x2 <- read.table(mcmc2,header=TRUE,sep="\t")
x3 <- read.table(mcmc3,header=TRUE,sep="\t")
x4 <- read.table(mcmc4,header=TRUE,sep="\t")

## Estimate average and write
cat("#### average_mcmcs.R: Estimating and outputting average ... \n")
write.table((x1+x2+x3+x4)/4, file = out_file, row.names=FALSE, sep="\t", quote = FALSE)

## Report:
cat("\n#### average_mcmcs.R: Done with script.\n")


#!/usr/bin/env Rscript

cat("\n#### ancestral_states.R: Starting script.\n")

## Command-line arguments
options(warn = 1)
args <- commandArgs(trailingOnly = TRUE)
phylogeny = args[1]
characters = args[2]
WD = args[3]

## Report
cat("\n#### ancestral_states.R: Phylogeny:", phylogeny, "\n")
cat("#### ancestral_states.R: Characters:", characters, "\n")
cat("#### ancestral_states.R: Working directory:", WD, "\n\n")

## Packages
library(phytools)
library(geiger)

## Process command-line arguments
tree = read.tree(phylogeny)
chars = read.csv(characters)
setwd(WD)

## Estimate model fit
cat("#### ancestral_states.R: Estimating model fit ... \n")
foraging = factor(chars$char)
names(foraging) = chars$species

species.to.exclude = tree$tip.label[!(tree$tip.label %in% chars$species)]
tree = drop.tip(tree, species.to.exclude)

results.anc = data.frame(model=c("ER","SYM","ARD"), lnL=numeric(3), AICc=numeric(3), params=numeric(3))

ER_fit <- fitDiscrete(tree, foraging, model="ER")
SYM_fit <- fitDiscrete(tree, foraging, model="SYM")
ARD_fit <- fitDiscrete(tree, foraging, model="ARD")
results.anc[1,-1] <- c(lnL=ER_fit$opt$lnL, AICc=ER_fit$opt$aicc, ER_fit$opt$k)
results.anc[2,-1] <- c(lnL=SYM_fit$opt$lnL, AICc=SYM_fit$opt$aicc, SYM_fit$opt$k)
results.anc[3,-1] <- c(lnL=ARD_fit$opt$lnL, AICc=ARD_fit$opt$aicc, ARD_fit$opt$k)

write.table(results.anc, "modelfit.txt")

# Ancestral state estimation and plotting
cat("#### ancestral_states.R: Estimating ancestral states ... \n")
ER <- make.simmap(tree, matrix, nsim=500, model="ER")
SYM <- make.simmap(tree, matrix, nsim=500, model="SYM")
ARD <- make.simmap(tree, matrix, nsim=500, model="ARD")

pdf(paste0("anc_states.pdf")
describe.simmap(ER, plot=TRUE)
describe.simmap(SYM, plot=TRUE)
describe.simmap(ARD, plot=TRUE)
dev.off()

## Report:
cat("\n#### ancestral_states.R: Done with script.\n")
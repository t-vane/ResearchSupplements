#!/usr/bin/env Rscript

#Packages
library(phytools)
library(geiger)

args = commandArgs(trailingOnly=TRUE)
phylogeny = args[1]
characters = args[2]
directory = args[3]

# Read the data
tree = read.tree(phylogeny)
chars = read.csv(characters)

# Estimate model fit
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

write.table(results.anc, paste0(directory, "/modelfit.txt"))

# Ancestral state estimation and plotting
ER <- make.simmap(tree, matrix, nsim=500, model="ER")
SYM <- make.simmap(tree, matrix, nsim=500, model="SYM")
ARD <- make.simmap(tree, matrix, nsim=500, model="ARD")

pdf(paste0(directory, "/anc_states.pdf"))
describe.simmap(ER, plot=TRUE)
describe.simmap(SYM, plot=TRUE)
describe.simmap(ARD, plot=TRUE)
dev.off()
#!/usr/bin/env Rscript

cat("\n#### pca.R: Starting script.\n")

## Command-line arguments
options(warn = 1)
args <- commandArgs(trailingOnly = TRUE)
WD <- args[1]
cov_matrix <- args[2]
ind_list <- args[3]
set_id <- args[4]

## Packages
library(ggplot2)
library(grid)
library(gridExtra)

## Report
cat("\n#### pca.R: Working directory:", WD, "\n")
cat("#### pca.R: Covariance matrix:", cov_matrix, "\n")
cat("#### pca.R: List with individuals mapped to populations:", ind_list, "\n")
cat("#### pca.R: Set ID:", set_id, "\n\n")

## Process command-line arguments
setwd(WD)
tab <- read.table(cov_matrix)
popdat <- read.table(ind_list)

## Calculate eigenvalues
cat("#### pca.R: Estimating eigenvalues ... \n")
m <- as.matrix(tab)
e <- eigen(m)

## Calculate Euclidean distances and save them
cat("#### pca.R: Estimating Euclidean distances ... \n")
cov_euc_dd <- dist(m)
save(cov_euc_dd, file=paste0(cov_matrix,".euc_dist.rda"))

## Create vector that encodes colors and shapes for subsequent plotting
populationColors <- c()
populationShapes <- c()
for (i in popdat[,2]) {
  if (i == "Anjahamana") {
    populationColors <- c(populationColors, "#7f3e9b")
    populationShapes <- c(populationShapes, 0)
  }
  else if (i == "Sahamamy") {
    populationColors <- c(populationColors, "#003c8d")
    populationShapes <- c(populationShapes, 1)
  }
  else if (i == "Andobo") {
    populationColors <- c(populationColors, "#cc3a91")
    populationShapes <- c(populationShapes, 2)
  }
  else if (i == "Sahafina") {
    populationColors <- c(populationColors, "#f50035")
    populationShapes <- c(populationShapes, 3)
  }
  else if (i == "Vohiposa") {
    populationColors <- c(populationColors, "#ff5a1f")
    populationShapes <- c(populationShapes, 4)
  }
  else if (i == "Ambodisakoana") {
    populationColors <- c(populationColors, "#e09d00")
    populationShapes <- c(populationShapes, 0)
  }
  else if (i == "Antanambao") {
    populationColors <- c(populationColors, "#f5f100")
    populationShapes <- c(populationShapes, 1)
  }
}

## Plot
cat("#### pca.R: Plotting ... \n")
# Distribution of explained variance
svg(paste0(set_id,".svg"), 11.7, 11.7)
par(mfrow=c(2,2))
plot(1:dim(tab)[1], e$values/sum(e$values)*100, pch=1, col="black", xlab="PCs", ylab="% of variance explained", ylim=c(0,60), cex.axis=1.5, cex.lab=1.5)
# PCA
for(i in c(1,3,5)){ j=i+1;
	plot(e$vectors[,i:j], xlab=paste0("PC", i, "; ", round(e$values[i]/sum(e$values)*100,2), "%"), ylab=paste0("PC", j, "; ",round(e$values[j]/sum(e$values)*100,2), "%"), lwd=2, cex=2.5, cex.lab=1.5, cex.axis=1.5, col=populationColors, pch=populationShapes)
	legend("topright", pch=c(0,1,2,3,4,5,6),
       col=c("#7F3E9B","#003C8D","#CC3A91","#F50035","#FF5A1F","#E09D00","#F5F100"),
       legend=c("Anjahamana","Sahamamy","Andobo","Sahafina","Vohiposa","Ambodisakoana","Antanambao"), cex=1.25)
}
dev.off()

## Report:
cat("\n#### pca.R: Done with script.\n")




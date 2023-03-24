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
  if (i == "Marojejy") {
    populationColors <- c(populationColors, "#0A253B")
    populationShapes <- c(populationShapes, 0)
  }
  else if (i == "Anjanaharibe") {
    populationColors <- c(populationColors, "#235179")
    populationShapes <- c(populationShapes, 1)
  }
  else if (i == "Anjiahely") {
    populationColors <- c(populationColors, "#3E7DBB")
    populationShapes <- c(populationShapes, 2)
  }
  else if (i == "Ambavala") {
    populationColors <- c(populationColors, "#67488B")
    populationShapes <- c(populationShapes, 3)
  }
  else if (i == "Riamalandy") {
    populationColors <- c(populationColors, "#8F3559")
    populationShapes <- c(populationShapes, 4)
  }
  else if (i == "Ambatovy") {
    populationColors <- c(populationColors, "#C90000")
    populationShapes <- c(populationShapes, 0)
  }
  else if (i == "Anjozorobe") {
    populationColors <- c(populationColors, "#EB7C05")
    populationShapes <- c(populationShapes, 1)
  }
  else if (i == "Tsinjoarivo") {
    populationColors <- c(populationColors, "#E09406")
    populationShapes <- c(populationShapes, 2)
  }
  else if (i == "Ambohitantely") {
    populationColors <- c(populationColors, "#BA7839")
    populationShapes <- c(populationShapes, 3)
  }
  else if (i == "Ankafobe") {
    populationColors <- c(populationColors, "#92583F")
    populationShapes <- c(populationShapes, 4)
    
  }
}

## Plot
cat("#### pca.R: Plotting ... \n")
# Distribution of explained variance
svg(paste0(set_id,".svg"), 8.3, 11.7)
par(mfrow=c(3,2))
plot(1:dim(tab)[1], e$values/sum(e$values)*100, pch=1, col="black", xlab="PCs", ylab="% of variance explained", ylim=c(0,15))
# PCA
for(i in c(1,3,5,7,9)){ j=i+1;
	plot(e$vectors[,i:j], xlab=paste0("PC", i, "; ", round(e$values[i]/sum(e$values)*100,2), "%"), ylab=paste0("PC", j, "; ",round(e$values[j]/sum(e$values)*100,2), "%"), lwd=1, col=populationColors, pch=populationShapes)
	legend("bottomright",pch=c(0,1,2,3,4,0,1,2,3,4),
       col=c("#0A253B","#235179","#3E7DBB","#67488B","#8F3559","#C90000","#EB7C05","#E09406","#BA7839","#92583F"),
       legend=c("Marojejy","Anjanaharibe","Anjiahely","Ambavala","Riamalandy","Ambatovy","Anjozorobe","Tsinjoarivo","Ambohitantely","Ankafobe"), cex=0.8)
}
dev.off()

## Report:
cat("\n#### pca.R: Done with script.\n")




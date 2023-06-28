#!/usr/bin/env Rscript

cat("\n#### ibd_geneticDistances.R: Starting script.\n")

## Command-line arguments
options(warn = 1)
args <- commandArgs(trailingOnly = TRUE)
WD <- args[1]
geo_dist <- args[2]
gen_dist <- args[3]
out <- args[4]
gen_dist_sd <- args[5]

## Report
cat("\n#### ibd_geneticDistances.R: Working directory:", WD, "\n")
cat("#### ibd_geneticDistances.R: Geographic distance matrix:", geo_dist, "\n")
cat("#### ibd_geneticDistances.R: Genetic distance matrix:", gen_dist, "\n")
cat("#### ibd_geneticDistances.R: Genetic distance standard deviation matrix:", gen_dist_sd, "\n")
cat("#### ibd_geneticDistances.R: Output prefix:", out, "\n\n")

## Packages
library(vegan)
library(ggplot2)
library(ggrepel)

## Process command-line arguments
cat("#### ibd_geneticDistances.R: Processing and transforming input ... \n")
setwd(WD)
geo <- read.table(geo_dist, row.names=1 header=TRUE)
gen <- read.table(gen_dist, row.names=1 header=TRUE)
gen_sd <- read.table(gen_dist_sd, row.names=1 header=TRUE)
# Remove NA values
geo[is.na(geo)]<-0
gen[is.na(gen)]<-0
gen_sd[is.na(gen_sd)]<-0
# Transform
geo <- as.data.frame(sapply(geo+0.1, log))
# Convert to distance matrix
dist_geo <- as.dist(as(geo, "matrix"))
dist_gen <- as.dist(as(gen, "matrix"))

## Conduct and write Mantel test
cat("#### ibd_geneticDistances.R: Conducting Mantel tests ... \n")
mantel_all <- mantel(dist_geo, dist_gen, method = "spearman", permutations = 9999)
capture.output(mantel_all, file = paste0(out, ".all.mantel.txt"))

## Reformat
populationPairs <- c()
geneticDistances<- c()
geneticDistances_sd <- c()
geographicDistances <- c()
for(row in 1:nrow(gen)) {
  for(col in 1:ncol(gen)) {
    if (row < col) {
      populationPairs <- c(populationPairs,paste0(rownames(gen)[row],"-",colnames(gen)[col]))
      geneticDistances <- c(geneticDistances, gen[row, col])
	  geneticDistances_sd <- c(geneticDistances_sd, gen_sd[row, col])
      geographicDistances <- c(geographicDistances, geo[row, col])
    }
  }
}
comb <- data.frame(geographicDistances, geneticDistances, populationPairs, geneticDistances_sd)
colnames(comb) <- c("geographicDistances", "geneticDistances", "populationPair", "geneticDistances_sd")

## Annotate which is north-north/south-south/north-south comparison
comb$category <- rep("all", nrow(comb))
comb[7,4] <- "north"
comb[10,4] <- "north"
comb[14,4] <- "north"
comb[3,4] <- "south"
comb[4,4] <- "south"
comb[6,4] <- "south"
comb[16,4] <- "south"
comb[18,4] <- "south"
comb[20,4] <- "south"

## Plot isolation-by-distance
cat("#### ibd_geneticDistances.R: Plotting ... \n")
svg(paste0(out, ".plot.svg"))
ggplot(comb, aes(x = geographicDistances, y = geneticDistances, fontface=1)) + 
  geom_errorbar(aes(ymin=geneticDistances-geneticDistances_sd, ymax=geneticDistances+geneticDistances_sd, color=category), width=.075, position="identity") +
  geom_point(aes(color=category), size=1.8) + 
  scale_color_manual(labels = c("North-South", "North", "South"), values=c("#000000", "#6e3c93", "#f27a15")) + 
  xlab("Geographic distance (log km)") + ylab("Mean genetic distance") + 
  geom_text_repel(aes(label = populationPair), size=3) +
  theme_bw() + 
  theme(legend.position = c(0.85, 0.13), legend.direction = "vertical", legend.title = element_blank(), legend.spacing.y = unit(0, "mm"), 
	panel.border = element_rect(colour = "black", fill=NA), legend.background = element_blank(), legend.box.background = element_rect(colour = "black").
	axis.text = element_text(size=10), axis.title=element_text(size=12)) 
dev.off()

## Report:
cat("\n#### ibd_geneticDistances.R: Done with script.\n")


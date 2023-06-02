#!/usr/bin/env Rscript

cat("\n#### ibd.R: Starting script.\n")

## Command-line arguments
options(warn = 1)
args <- commandArgs(trailingOnly = TRUE)
geo_dist <- args[1]
gen_dist <- args[2]
out <- args[3]

## Report:
cat("#### ibd.R: Geographic distance matrix:", geo_dist, "\n")
cat("#### ibd.R: Genetic distance matrix:", gen_dist, "\n")
cat("#### ibd.R: Output prefix:", out, "\n\n")

## Packages
library(vegan)
library(ggplot2)
library(ggrepel)

## Process command-line arguments
cat("#### ibd.R: Processing and transforming input ... \n")
geo <- read.table(geo_dist, row.names=1 header=TRUE)
gen <- read.table(gen_dist, row.names=1 header=TRUE)
# Remove NA values
geo[is.na(geo)]<-0
gen[is.na(gen)]<-0
# Transform
slatkinLinearized <- function(val) {
  transformed <- val/(1-val)
  return(transformed)
}
gen <- as.data.frame(sapply(gen, slatkinLinearized))
geo <- as.data.frame(sapply(geo+0.1, log))
# Convert to distance matrix
dist_geo <- as.dist(as(geo, "matrix"))
dist_gen <- as.dist(as(gen, "matrix"))

## Conduct and write Mantel tests
cat("#### ibd.R: Conducting Mantel tests ... \n")
# Complete data set
cat("## ibd.R: For all populations ... \n")
mantel_all <- mantel(dist_geo, dist_gen, method = "spearman", permutations = 9999)
capture.output(mantel_all, file = paste0(out, ".all.mantel.txt"))
# Only northern populations
cat("## ibd.R: For northern populations ... \n")
pops_n <- c("ambavala", "anjanaharibe", "anjiahely", "marojejy")
dist_geo_n <- dist_geo[rownames(dist_geo) %in% pops_n, colnames(dist_geo) %in% pops_n]
dist_gen_n <- dist_gen[rownames(dist_gen) %in% pops_n, colnames(dist_gen) %in% pops_n]
mantel_n <- mantel(dist_geo_n, dist_gen_n, method = "spearman", permutations = 9999)
capture.output(mantel_n, file = paste0(out, ".north.mantel.txt"))
# Only southern and central populations
cat("## ibd.R: For southern and central populations ... \n")
pop_sc <- c("ambatovy", "ambohitantely", "ankafobe", "tsinjoarivo")
dist_geo_sc <- dist_geo[rownames(dist_geo) %in% pops_sc, colnames(dist_geo) %in% pops_sc]
dist_gen_sc <- dist_gen[rownames(dist_gen) %in% pops_sc, colnames(dist_gen) %in% pops_sc]
mantel_sc <- mantel(dist_geo_sc, dist_gen_sc, method = "spearman", permutations = 9999)
capture.output(mantel_sc, file = paste0(out, ".south.mantel.txt"))

## Reformat
populationPairs <- c()
geneticDistances<- c()
geographicDistances <- c()
for(row in 1:nrow(gen)) {
  for(col in 1:ncol(gen)) {
    if (row < col) {
      populationPairs <- c(populationPairs,paste0(rownames(gen)[row],"-",colnames(gen)[col]))
      geneticDistances <- c(geneticDistances, gen[row, col])
      geographicDistances <- c(geographicDistances, geo[row, col])
    }
  }
}
comb <- data.frame(geographicDistances, geneticDistances, populationPairs)
colnames(comb) <- c("geographicDistances", "geneticDistances", "populationPair")

## Annotate which is north-north/south-south/north-south comparison
comb$category <- rep("all", nrow(comb))
comb[9,4] <- "north"
comb[10,4] <- "north"
comb[12,4] <- "north"
comb[19,4] <- "north"
comb[21,4] <- "north"
comb[24,4] <- "north"
comb[2,4] <- "south"
comb[5,4] <- "south"
comb[7,4] <- "south"
comb[16,4] <- "south"
comb[18,4] <- "south"
comb[27,4] <- "south"

## Plot isolation-by-distance
cat("#### ibd.R: Plotting ... \n")
svg(paste0(out, ".plot.svg"))
ggplot(comb, aes(x = geographicDistances, y = geneticDistances, fontface=1)) + 
  geom_smooth(data=subset(comb), color="black", show.legend = FALSE, method='lm', se=FALSE, linetype = "dashed") +
  geom_point(aes(color=category), size=2) + 
  scale_color_manual(labels = c("North - South/Central Highlands", "North", "South/Central Highlands"), values=c("#000000", "#3E7DBB", "#C90000")) + 
  xlab("Geographic distance (log km)") + ylab(expression(Genetic~distance~(italic(F[ST])/(1-italic(F[ST]))))) + 
  geom_text_repel(aes(label = populationPair), size=2.75) +
  theme_bw() + 
  theme(legend.position = c(0.14, 0.897), legend.direction = "vertical", legend.title = element_blank(), legend.spacing.y = unit(0, "mm"), 
	panel.border = element_rect(colour = "black", fill=NA), legend.background = element_blank(), legend.box.background = element_rect(colour = "black")) 
dev.off()

## Report:
cat("\n#### ibd.R: Done with script.\n")


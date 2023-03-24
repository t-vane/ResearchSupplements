#!/usr/bin/env Rscript

#### Script adapted and modified from J. Salmona

cat("\n#### plot_ngsadmix.R: Starting script.\n")

## Command-line arguments
options(warn = 1)
args <- commandArgs(trailingOnly = TRUE)
WD <- args[1]
like_values <- args[2]
ind_list <- args[3]
set_id <- args[4]

## Report
cat("\n#### plot_ngsadmix.R: Working directory:", WD, "\n")
cat("#### plot_ngsadmix.R: File with likelihood values:", like_values, "\n")
cat("#### plot_ngsadmix.R: List with individuals mapped to populations:", ind_list, "\n")
cat("#### plot_ngsadmix.R: Set ID:", set_id, "\n\n")

## Process command-line arguments
setwd(WD)
tab <- read.table(like_values)
popdat <- read.table(ind_list)

## Estimate deltaK following Evanno et al. (2005), Mol. Ecol.
cat("#### plot_ngsadmix.R: Estimating deltaK ... \n")
simsum <- tab
colnames(simsum) <-c ("LnPr","K","it","rep")
maxK <- max(simsum[,2])
Kmeans <- tapply(simsum[,1],simsum[,2],mean)
Kvars <- tapply(simsum[,1],simsum[,2],var)
Ksds <- Kvars^0.5
deltaK <- rep(NA,maxK)
for(x in 2:maxK){
	deltaK[x] <- abs(Kmeans[x+1] - 2*Kmeans[x] + Kmeans[x-1])/Ksds[x]
}

## Plot deltaK
cat("#### plot_ngsadmix.R: Plotting deltaK ... \n")
pdf(paste(set_id, ".plot_ngsadmix.pdf",sep=""), 16, 12)
par(mfrow=c(1,2))
bp <- boxplot(simsum[,1]~simsum[,2], xlab="K", ylab='Log likelihood')
plot(deltaK, type = 'b', xlab='K', ylab=expression(paste(Delta, "K")))
axis(labels=c(1,2,3,4,5,6,7,8,9,10), at=c(1,2,3,4,5,6,7,8,9,10), cex.axis=1, side=1)

## Plot ancestry for best likelihood run of each K
par(mfrow=c(3,1), mar=c(9,6,4,2))
for(i in 2:max(tab[,2])) {
	cat("#### plot_ngsadmix.R: Plotting ancestry for K=", i, " ... \n")
	tab_tmp <- tab[tab[,2]==i,]
	bestlike <- tab_tmp[tab_tmp[,1] == max(tab_tmp[,1], na.rm = T), 4]
	admix <- t(as.matrix(read.table(paste(set_id,".K",i,".seed",bestlike,".qopt",sep=""))))
	admix <- admix[,order(popdat[,2])]; 
	popdat_tmp <- popdat[order(popdat[,2]),]
  
	# Plot with individual labels
	barplot(admix, col=rainbow(i), space=0, border=NA, ylab="", yaxt="n", horiz=F, names.arg=popdat_tmp[,1], las=3, cex.names=.8, oma=c(12, 2, 2, 2), mar=c(20, 2, 20, 2))
	axis(2, pos=-1, cex.axis=1)
	title(ylab="Admixture", line=1, cex.lab=1)
	title(paste("K =",i), adj = 0, line = 3, cex.main=1.5)
  
	# Plot with population labels
	level <- levels(as.factor(popdat_tmp[,2]))
	level <- gsub("[ABCDEFGHIJ]_", "", level)
	pop <- popdat_tmp[order(match(popdat_tmp[,2], level)),]
	tempo <- tapply(1:nrow(pop),pop[,2],mean);
	tempo <- tempo[order(tempo)]
	tempomax <- tapply(1:nrow(pop),pop[,2],max);
	tempomax <- tempomax[order(tempomax)]
	barplot(admix, col=rainbow(i), space=0, border=NA, ylab="", yaxt="n", , xaxt='n', oma=c(12, 2, 2, 2))
	axis(2, pos=-1, cex.axis=1)
	title(ylab="Admixture", line=1, cex.lab=1)
	text(tempo-0.5, 1.08, labels=level, xpd=T, srt=0, cex.lab=.8)
	abline(v=head(tempomax,-1), lty=2, lwd=1.2, col="white")
}
dev.off()

## Report:
cat("\n#### plot_ngsadmix.R: Done with script.\n")


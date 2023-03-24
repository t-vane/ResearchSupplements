#!/usr/bin/env Rscript

cat("\n#### plot_stairwayplot.R: Starting script.\n")

## Command-line arguments
options(warn = 1)
args <- commandArgs(trailingOnly = TRUE)
WD <- args[1]
in1 <- args[2]
in2 <- args[3]
in3 <- args[4]
in4 <- args[5]

## Report
cat("\n#### plot_stairwayplot.R: Working directory:", WD, "\n")
cat("#### plot_stairwayplot.R: Input file for population 1:", in1, "\n")
cat("#### plot_stairwayplot.R: Input file for population 2:", in2, "\n")
cat("#### plot_stairwayplot.R: Input file for population 3:", in3, "\n")
cat("#### plot_stairwayplot.R: Input file for population 4:", in4, "\n\n")

## Process command-line arguments
setwd(WD)
pop1 = read.table(in1, header=TRUE)
pop2 = read.table(in2, header=TRUE)
pop3 = read.table(in3, header=TRUE)
pop4 = read.table(in4, header=TRUE)
pops <- c(pop1, pop2, pop3, pop4)

names <- c("ankafobe", "ambatovy", "ambavala", "metapopulation")
cols <- c("#92583F", "#C90000","#67488B", "#4C7DBF")

## Plot all populations together
cat("#### plot_stairwayplot.R: Plotting all populations together ... \n")
pdf("all.plot_stairwayplot.pdf", 16, 12)
plot(pop1$year, pop1$Ne_median, log="x", xlab="Time before present (years)", ylab=expression(Effective~population~size~(x10^'4')), xlim=c(100,1000000), ylim=c(0,25), type="l", lwd = 4, col=cols[1])
lines(pop2$year, pop2$Ne_median, lwd = 4, col=cols[2])
lines(pop3$year, pop3$Ne_median, lwd = 4, col=cols[3])
lines(pop4$year, pop4$Ne_median, lwd = 4, col=cols[4])
abline(v=c(1000, 2000, 19000, 26500), col="black", lwd = 1, lty=2) 
legend(100, 350000, legend=c(paste(names[1], "(N=7)"), paste(names[2], "(N=7)"), paste(names[3], "(N=7)"), paste(names[4], "(N=7)"), col=cols, lty=1:1, cex=0.9)
dev.off()

## Plot single populations with confidence intervals
for (i in 1:4) {
	cat("#### plot_stairwayplot.R: Plotting population ", names[i], " ... \n")
	pdf(paste0(names[i], ".plot_stairwayplot.pdf",sep=""), 16, 12)
	plot(pops[i]$year, pops[i]$Ne_median, log="x", xlab="Time before present (years)", ylab=expression(Effective~population~size~(x10^'4')), xlim=c(100,1000000), ylim=c(0,25), type="l", lwd = 4, col=cols[i])
	lines(pops[i]$year, pops[i]$Ne_2.5., lwd = 1, col=cols[i])
	lines(pops[i]$year, pops[i]$Ne_97.5., lwd = 1, col=cols[i])
	dev.off()
}

## Report:
cat("\n#### plot_stairwayplot.R: Done with script.\n")


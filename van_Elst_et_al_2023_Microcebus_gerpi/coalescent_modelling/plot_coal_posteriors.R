#!/usr/bin/env Rscript

cat("\n#### plot_coal_posteriors.R: Starting script.\n")

## Command-line arguments
options(warn = 1)
args <- commandArgs(trailingOnly = TRUE)
mcmc_nm1 <- args[1]
mcmc_nm2 <- args[2]
mcmc_nm3 <- args[3]
mcmc_nm4 <- args[4]
mcmc_m1 <- args[5]
mcmc_m2 <- args[6]
mcmc_m3 <- args[7]
mcmc_m4 <- args[8]
m_scale <- as.numeric(args[9])
t_scale <- as.numeric(args[10])
out_dir <- args[11]

## Report
cat("\n#### plot_coal_posteriors.R: First MCMC (no migration):", mcmc_nm1, "\n")
cat("#### plot_coal_posteriors.R: Second MCMC (no migration):", mcmc_nm2, "\n")
cat("#### plot_coal_posteriors.R: Third MCMC (no migration):", mcmc_nm3, "\n")
cat("#### plot_coal_posteriors.R: Fourth MCMC (no migration):", mcmc_nm4, "\n")
cat("#### plot_coal_posteriors.R: First MCMC (migration):", mcmc_m1, "\n")
cat("#### plot_coal_posteriors.R: Second MCMC (migration):", mcmc_m2, "\n")
cat("#### plot_coal_posteriors.R: Third MCMC (migration):", mcmc_m3, "\n")
cat("#### plot_coal_posteriors.R: Fourth MCMC (migration):", mcmc_m4, "\n")
cat("#### plot_coal_posteriors.R: Inverse scaling factor used in the G-PhoCS configuration file for migration parameter:", m_scale, "\n")
cat("#### plot_coal_posteriors.R: Inverse scaling factor used in the G-PhoCS configuration file for tau and theta:", t_scale, "\n")
cat("#### plot_coal_posteriors.R: Output directory:", out_dir, "\n\n")

## Packages
library(ggplot2)
library(gridExtra)

## Process command-line arguments
cat("#### plot_coal_posteriors.R: Reading MCMC files ... \n")
no_mig1 <- read.table(mcmc_nm1,header=TRUE,sep="\t")
no_mig2 <- read.table(mcmc_nm2,header=TRUE,sep="\t")
no_mig3 <- read.table(mcmc_nm3,header=TRUE,sep="\t")
no_mig4 <- read.table(mcmc_nm4,header=TRUE,sep="\t")
mig1 <- read.table(mcmc_m1,header=TRUE,sep="\t")
mig2 <- read.table(mcmc_m2,header=TRUE,sep="\t")
mig3 <- read.table(mcmc_m3,header=TRUE,sep="\t")
mig4 <- read.table(mcmc_m4,header=TRUE,sep="\t")

## Convert MCMCs with scaling factors
cat("#### plot_coal_posteriors.R: Converting MCMCs with scaling factors ... \n")
mcmcs <- c("no_mig1", "no_mig2", "no_mig3", "no_mig4", "mig1", "mig2", "mig3", "mig4")
for (mcmc_name in mcmcs) {
  mcmc <- get(mcmc_name)
  mcmc[,2:20] <- mcmc[,2:20]/t_scale
  assign(mcmc_name, mcmc)
}
mcmcs <- c("mig1", "mig2", "mig3", "mig4")
for (mcmc_name in mcmcs) {
  mcmc <- get(mcmc_name)
  mcmc[,21:32] <- mcmc[,21:32]/m_scale
  assign(mcmc_name, mcmc)
}

## Assign names
paramTracker <- c("Dummy",
                  expression("\u03B8"[italic("M. jollyae")]),expression("\u03B8"[italic("M. marohita")]),expression("\u03B8"[Sahamamy]),expression("\u03B8"[Andobo]),expression("\u03B8"[VohiposaSahafina]),expression("\u03B8"[Antanambao]),expression("\u03B8"[Ambodisakoana]),expression("\u03B8"[italic("M. jollyae-M. marohita")]),expression("\u03B8"[Sahamamy-Andobo]),expression("\u03B8"[Ambodisakoana-VohiposaSahafina]),expression("\u03B8"[(Ambodisakoana-VohiposaSahafina)-Antanambao]),expression("\u03B8"[paste(italic("M. gerpi"), " root")]),expression("\u03B8"[Root]),
                  expression("\u03C4"[italic("M. jollyae-M. marohita")]), expression("\u03C4"[Sahamamy-Andobo]), expression("\u03C4"[Ambodisakoana-VohiposaSahafina]), expression("\u03C4"[(Ambodisakoana-VohiposaSahafina)-Antanambao]), expression("\u03C4"[paste(italic("M. gerpi"), " root")]),expression("\u03C4"[Root]),
                  expression("m"[italic("M. jollyae") *symbol('\256')* italic("M. marohita")]), expression("m"[italic("M. marohita") *symbol('\256')* italic("M. jollyae")]), expression("m"[Sahamamy *symbol('\256')* Andobo]), expression("m"[Ambodisakoana *symbol('\256')* Antanambao]), expression("m"[Antanambao *symbol('\256')* Ambodisakoana]),expression("m"[Ambodisakoana *symbol('\256')* italic("M. marohita")]), expression("m"[VohiposaSahafina *symbol('\256')* Andobo]), expression("m"[Andobo *symbol('\256')* VohiposaSahafina]), expression("m"[VohiposaSahafina *symbol('\256')* Sahamamy]), 
                  expression("m"[Sahamamy *symbol('\256')* VohiposaSahafina]), expression("m"[VohiposaSahafina *symbol('\256')* Ambodisakoana]), expression("m"[Antanambao *symbol('\256')* italic("M. marohita")]),"Average data log-likelihood","Full log-likelihood")

## Step 1: Plot theta
cat("#### plot_coal_posteriors.R: Plotting theta ... \n")
mnh12 <- c()
mnh34 <- c()
mnhA12 <- c()
mnhA34 <- c()
mnhCounter <- 1

for (i in 2:14){
	# Select specific theta and convert to data frame
	no_mig_d1 <- data.frame(value = no_mig1[,i])
	no_mig_d2 <- data.frame(value = no_mig2[,i])
	no_mig_d3 <- data.frame(value = no_mig3[,i])
	no_mig_d4 <- data.frame(value = no_mig4[,i])
	mig_d1 <- data.frame(value = mig1[,i])
	mig_d2 <- data.frame(value = mig2[,i])
	mig_d3 <- data.frame(value = mig3[,i])
	mig_d4 <- data.frame(value = mig4[,i])

	# Concatenate data frames
	dat <- rbind(no_mig_d1,no_mig_d2,no_mig_d3,no_mig_d4,mig_d1,mig_d2,mig_d3,mig_d4)
	dat$variable <- factor(dat$variable)
		
	# Plot and rename
	p <- ggplot(dat, aes(x = variable, y = value)) + geom_violin(scale = "width", adjust = 1, width = 0.5,aes(fill=variable)) + 
		scale_fill_manual(values=c("#5ab4ac", "#5ab4ac", "#5ab4ac", "#5ab4ac","#C99700", "#C99700", "#C99700", "#C99700")) + 
		theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),legend.position="none") + labs(title=paramTracker[i])
	
	assign(paste("p",i, sep=""),p)
		
	# Create list of medians for node heigh plot
	temp1 <- c(no_mig1[,i],no_mig2[,i])
	temp2 <- c(no_mig3[,i],no_mig4[,i])
	mnh12[mnhCounter] <- median(temp1)
	mnh34[mnhCounter] <- median(temp2)
	tempA1 <- c(mig1[,i],mig2[,i])
	tempA2 <- c(mig3[,i],mig4[,i]) 
	mnhA12[mnhCounter] <- median(tempA1)
	mnhA34[mnhCounter] <- median(tempA2)
	mnhCounter <- mnhCounter + 1
}

# Concatenate list of medians and plot
mnhData <- as.data.frame(cbind(mnh12,mnh34))
mnhAData <- as.data.frame(cbind(mnhA12,mnhA34))
p_mnh <- ggplot(data=mnhData, aes(x=mnh12,y=mnh34,alpha=0.7)) + geom_point(size=3) + geom_segment(aes(x=0.0,y=0.0,xend=0.005,yend=0.005), linetype=2) + theme(legend.position="none") + labs(title="Median Node Heights (no migration)",x="chains 1 + 2",y="chains 3 + 4") + xlim(0.0,0.005) + ylim(0.0,0.005)
p_mnhA <- ggplot(data=mnhAData, aes(x=mnhA12,y=mnhA34,alpha=0.7)) + geom_point(size=3) + geom_segment(aes(x=0.0,y=0.0,xend=0.005,yend=0.005), linetype=2) + theme(legend.position="none") + labs(title="Median Node Heights (migration)",x="chains 1 + 2",y="chains 3 + 4") + xlim(0.0,0.005) + ylim(0.0,0.005)

# Write output
png(paste0(out_dir, "/theta_posteriors.png"),res=300, height=10*300, width=8*300)
grid.arrange(p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,empty,p_mnh,p_mnhA, ncol=2, nrow=8)
dev.off()

## Step 2: Plot tau
cat("#### plot_coal_posteriors.R: Plotting tau ... \n")
mnh12 <- c()
mnh34 <- c()
mnhA12 <- c()
mnhA34 <- c()
mnhCounter <- 1

for (i in 15:20){
	# Select specific tau and convert to data frame
	no_mig_d1 <- data.frame(value = no_mig1[,i])
	no_mig_d2 <- data.frame(value = no_mig2[,i])
	no_mig_d3 <- data.frame(value = no_mig3[,i])
	no_mig_d4 <- data.frame(value = no_mig4[,i])
	mig_d1 <- data.frame(value = mig1[,i])
	mig_d2 <- data.frame(value = mig2[,i])
	mig_d3 <- data.frame(value = mig3[,i])
	mig_d4 <- data.frame(value = mig4[,i])

	# Concatenate data frames
	dat <- rbind(no_mig_d1,no_mig_d2,no_mig_d3,no_mig_d4,mig_d1,mig_d2,mig_d3,mig_d4)
	dat$variable <- factor(dat$variable)
		
	# Plot and rename
	p <- ggplot(dat, aes(x = variable, y = value)) + geom_violin(scale = "width", adjust = 1, width = 0.5,aes(fill=variable)) + 
		scale_fill_manual(values=c("#5ab4ac", "#5ab4ac", "#5ab4ac", "#5ab4ac","#C99700", "#C99700", "#C99700", "#C99700")) + 
		theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),legend.position="none") + labs(title=paramTracker[i])
	
	assign(paste("p",i, sep=""),p)
		
	# Create list of medians for node heigh plot
	temp1 <- c(no_mig1[,i],no_mig2[,i])
	temp2 <- c(no_mig3[,i],no_mig4[,i])
	mnh12[mnhCounter] <- median(temp1)
	mnh34[mnhCounter] <- median(temp2)
	tempA1 <- c(mig1[,i],mig2[,i])
	tempA2 <- c(mig3[,i],mig4[,i]) 
	mnhA12[mnhCounter] <- median(tempA1)
	mnhA34[mnhCounter] <- median(tempA2)
	mnhCounter <- mnhCounter + 1
}

# Concatenate list of medians and plot
mnhData <- as.data.frame(cbind(mnh12,mnh34))
mnhAData <- as.data.frame(cbind(mnhA12,mnhA34))
p_mnh <- ggplot(data=mnhData, aes(x=mnh12,y=mnh34,alpha=0.7)) + geom_point(size=3) + geom_segment(aes(x=0.0,y=0.0,xend=0.005,yend=0.005), linetype=2) + theme(legend.position="none") + labs(title="Median Node Heights (no migration)",x="chains 1 + 2",y="chains 3 + 4") + xlim(0.0,0.005) + ylim(0.0,0.005)
p_mnhA <- ggplot(data=mnhAData, aes(x=mnhA12,y=mnhA34,alpha=0.7)) + geom_point(size=3) + geom_segment(aes(x=0.0,y=0.0,xend=0.005,yend=0.005), linetype=2) + theme(legend.position="none") + labs(title="Median Node Heights (migration)",x="chains 1 + 2",y="chains 3 + 4") + xlim(0.0,0.005) + ylim(0.0,0.005)

# Write output
png(paste0(out_dir, "/tau_posteriors.png"),res=300, height=10*300, width=8*300)
grid.arrange(p15,p16,p17,p18,p19,p20,p_mnh,p_mnhA, ncol=2, nrow=4)
dev.off()

## Step 3: Plot migration
cat("#### plot_coal_posteriors.R: Plotting migration parameters ... \n")
mnh12 <- c()
mnh34 <- c()
mnhCounter <- 1

for (i in 21:32){
	# Select specific migration parameters and convert to data frame
	mig_d1 <- data.frame(value = mig1[,i])
	mig_d2 <- data.frame(value = mig2[,i])
	mig_d3 <- data.frame(value = mig3[,i])
	mig_d4 <- data.frame(value = mig4[,i])

	# Concatenate data frames
	dat <- rbind(mig_d1,mig_d2,mig_d3,mig_d4)
	dat$variable <- factor(dat$variable)
		
	# Plot and rename
	p <- ggplot(dat, aes(x = variable, y = value)) + geom_violin(scale = "width", adjust = 1, width = 0.5,aes(fill=variable)) + 
		scale_fill_manual(values=c("#C99700", "#C99700", "#C99700", "#C99700")) + 
		theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),legend.position="none") + labs(title=paramTracker[i])
 	
	assign(paste("p",i, sep=""),p)
		
	# Create list of medians for node heigh plot
	temp1 <- c(no_mig1[,i],no_mig2[,i])
	temp2 <- c(no_mig3[,i],no_mig4[,i])
	mnh12[mnhCounter] <- median(temp1)
	mnh34[mnhCounter] <- median(temp2)
	mnhCounter <- mnhCounter + 1
}

# Concatenate list of medians and plot
mnhData <- as.data.frame(cbind(mnh12,mnh34))
p_mnh <- ggplot(data=mnhData, aes(x=mnh12,y=mnh34,alpha=0.7)) + geom_point(size=3) + geom_segment(aes(x=0.0,y=0.0,xend=0.005,yend=0.005), linetype=2) + theme(legend.position="none") + labs(title="Median Node Heights (no migration)",x="chains 1 + 2",y="chains 3 + 4") + xlim(0.0,0.005) + ylim(0.0,0.005)

# Write output
png(paste0(out_dir, "/migration_posteriors.png"),res=300,height=10*300,width=8*300)
grid.arrange(p21,p22,p23,p24,p25,p26,p27,p28,p29,p30,p31,p32,p_mnh, ncol=3, nrow=5)
dev.off()

## Step 4: Plot likelihoods (full and average data likelihood)
cat("#### plot_coal_posteriors.R: Plotting likelihoods ... \n")
for (i in 0:1){
	# Select specific tau and convert to data frame
	no_mig_d1 <- data.frame(value = no_mig1[,ncol(no_mig1)-i])
	no_mig_d2 <- data.frame(value = no_mig2[,ncol(no_mig2)-i])
	no_mig_d3 <- data.frame(value = no_mig3[,ncol(no_mig3)-i])
	no_mig_d4 <- data.frame(value = no_mig4[,ncol(no_mig4)-i])
	mig_d1 <- data.frame(value = mig1[,ncol(mig1)-i])
	mig_d2 <- data.frame(value = mig2[,ncol(mig2)-i])
	mig_d3 <- data.frame(value = mig3[,ncol(mig3)-i])
	mig_d4 <- data.frame(value = mig4[,ncol(mig4)-i])

	# Concatenate data frames
	dat <- rbind(no_mig_d1,no_mig_d2,no_mig_d3,no_mig_d4)
	datm <- rbind(mig_d1,mig_d2,mig_d3,mig_d4)
	dat$variable <- factor(dat$variable)
	datm$variable <- factor(datm$variable)
		
	# Plot and rename
	p <- ggplot(dat, aes(x = variable, y = value)) + geom_violin(scale = "width", adjust = 1, width = 0.5,aes(fill=variable)) + 
		scale_fill_manual(values=c("#5ab4ac", "#5ab4ac", "#5ab4ac", "#5ab4ac")) + 
		theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),legend.position="none") + labs(title=paste(paramTracker[ncol(mig1)-1-i], " (no migration)", sep=""))
	pm <- ggplot(datm, aes(x = variable, y = value)) + geom_violin(scale = "width", adjust = 1, width = 0.5,aes(fill=variable)) + 
		scale_fill_manual(values=c("#C99700", "#C99700", "#C99700", "#C99700")) + 
		theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),legend.position="none") + labs(title=paste(paramTracker[ncol(mig1)-1-i], " (migration)", sep=""))
  
	
	assign(paste("p",ncol(x5)-i, sep=""), p)
	assign(paste("pm",ncol(x5)-i, sep=""), pm)
}

# Write output
png(paste0(out_dir, "/likelihoods.png"),res=300, height=10*300, width=8*300)
grid.arrange(p34,pm34,p35,pm35, ncol=2, nrow=2)
dev.off()

## Report:
cat("\n#### plot_coal_posteriors.R: Done with script.\n")


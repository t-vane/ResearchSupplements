#!/usr/bin/env Rscript

cat("\n#### plot_coal_models.R: Starting script.\n")

## Command-line arguments
options(warn = 1)
args <- commandArgs(trailingOnly = TRUE)
mcmc_nm <- args[1]
mcmc_m <- args[2]
out_dir <- args[3]

## Report
cat("\n#### plot_coal_models.R: Mean MCMC (no migration):", mcmc_nm, "\n")
cat("#### plot_coal_models.R: Mean MCMC (migration):", mcmc_m, "\n")
cat("#### plot_coal_models.R: Output directory:", out_dir, "\n\n")

## Packages
library(ggplot2)
library(pacman)
library(dplyr)
library(tidyr)
library(bayestestR)
library(grid)

## Process command-line arguments
cat("#### plot_coal_models.R: Reading MCMC files ... \n")
no_mig <- read.table(mcmc_nm, header=TRUE)
mig <- read.table(mcmc_m, header=TRUE)

## Information on parent and child populations
pop <- c("Mjollyae", "Mmarohita", "Sahamamy", "Andobo", "Ambodisakoana", "VohiposaSahafina", "sako_vohiFina", "Antanambao", "saha_dobo", "bao_sakoVohiFina", "maro_joll", "gerpiRoot", "root")
parentPop <- c("maro_joll", "maro_joll", "saha_dobo", "saha_dobo", "sako_vohiFina", "sako_vohiFina", "bao_sakoVohiFina", "bao_sakoVohiFina", "gerpiRoot", "gerpiRoot", "root", "root", NA)
poplist <- data.frame(pop, parentPop)

## Plot models
for (model in c("no_mig", "mig")) {
	
	cat("#### plot_coal_models.R: Plotting model", model, "... \n")
	# Extract means
	summaryTable <- data.frame()
	for (population in poplist$pop) {
		print(population)
		meanPopSize <- mean(eval(parse(text = model))[[paste("Ne_", population,sep = "")]])
		meanDivTime <- mean(eval(parse(text = model))[[paste("div.time_", population,sep = "")]])
  
		summaryTable <- rbind(summaryTable, data.frame(population, meanPopSize,meanDivTime))
	}
	
	# Set margins for rectangles for each population 
	ylimit <- 550000
	summaryTable$xmin <- NA
	summaryTable$xmax <- NA
	summaryTable$ymin <- NA
	summaryTable$ymax <- NA

	for (population in poplist$pop) {
		# Get parent and daughter populations from poplist
		parent <- poplist[which(poplist[,1]==population),2]
		child <- poplist[which(poplist[,2]==population),1]
		
		# Set ymax
		if (is.na(parent)) {
			ymax <- ylimit
		} else {
			ymax <- summaryTable[which(summaryTable$population==parent),3] # Parent divergence time
		}
		
		# Set ymin
		if (length(child) == 0) {
			ymin <- 0
		} else {
			ymin <- summaryTable[which(summaryTable$population==population),3] # Own divergence time
		}
  
		summaryTable[which(summaryTable$population==population),7] <- ymax
		summaryTable[which(summaryTable$population==population),6] <- ymin
	}

	# Sahamamy
	population = "Sahamamy"
	xmin_Sahamamy <- 25000
	xmax_Sahamamy <- xmin_Sahamamy + summaryTable[which(summaryTable$population==population),2]
	summaryTable[which(summaryTable$population==population),5] <- xmax_Sahamamy
	summaryTable[which(summaryTable$population==population),4] <- xmin_Sahamamy
	
	# saha_dobo
	population = "saha_dobo"
	xmin_saha_dobo <- xmax_Sahamamy
	xmax_saha_dobo <-  xmax_Sahamamy + summaryTable[which(summaryTable$population==population),2]
	summaryTable[which(summaryTable$population==population),5] <- xmax_saha_dobo
	summaryTable[which(summaryTable$population==population),4] <- xmin_saha_dobo

	# Andobo
	population = "Andobo"
	xmin_Andobo <- xmax_saha_dobo
	xmax_Andobo <- xmax_saha_dobo + summaryTable[which(summaryTable$population==population),2]
	summaryTable[which(summaryTable$population==population),5] <- xmax_Andobo
	summaryTable[which(summaryTable$population==population),4] <- xmin_Andobo

	# Insert horizontal line beneath M. gerpi root to stretch phylogeny
	ymin_line1 <- summaryTable[which(summaryTable$population=="gerpiRoot"),3]
	ymax_line1 <- summaryTable[which(summaryTable$population=="gerpiRoot"),3]
	xmin_line1 <- xmax_saha_dobo
	xmax_line1 <- xmax_saha_dobo + 50000
	summaryTable <- rbind(summaryTable, c(NA, NA, NA, xmin_line1, xmax_line1, ymin_line1, ymax_line1))

	# gerpiRoot
	population = "gerpiRoot"
	xmin_gerpiRoot <- xmin_line1 + (xmax_line1 - xmin_line1 - summaryTable[which(summaryTable$population=="gerpiRoot"),2])/2
	xmax_gerpiRoot <- xmax_line1 - (xmax_line1 - xmin_line1 - summaryTable[which(summaryTable$population=="gerpiRoot"),2])/2
	summaryTable[which(summaryTable$population==population),5] <- xmax_gerpiRoot
	summaryTable[which(summaryTable$population==population),4] <- xmin_gerpiRoot
	
	# bao_sakoVohiFina
	population="bao_sakoVohiFina"
	xmin_bao_sakoVohiFina <- xmax_line1
	xmax_bao_sakoVohiFina <- xmax_line1 + summaryTable[which(summaryTable$population==population),2]
	summaryTable[which(summaryTable$population==population),5] <- xmax_bao_sakoVohiFina
	summaryTable[which(summaryTable$population==population),4] <- xmin_bao_sakoVohiFina

	# Insert horizontal line beneath root to stretch phylogeny
	ymin_line2 <- summaryTable[which(summaryTable$population=="root"),3]
	ymax_line2 <- summaryTable[which(summaryTable$population=="root"),3]
	xmin_line2 <- xmax_gerpiRoot
	xmax_line2 <- xmax_gerpiRoot + 90000
	summaryTable <- rbind(summaryTable, c(NA, NA, NA, xmin_line2, xmax_line2, ymin_line2, ymax_line2))

	# root
	population="root"
	xmin_root <- xmin_line2 + (xmax_line2 - xmin_line2 - summaryTable[which(summaryTable$population=="root"),2])/2
	xmax_root <- xmax_line2 - (xmax_line2 - xmin_line2 - summaryTable[which(summaryTable$population=="root"),2])/2
	summaryTable[which(summaryTable$population==population),5] <- xmax_root
	summaryTable[which(summaryTable$population==population),4] <- xmin_root
	
	# Antanambao
	population="Antanambao"
	xmin_Antanambao <- xmax_bao_sakoVohiFina
	xmax_Antanambao <- xmax_bao_sakoVohiFina + summaryTable[which(summaryTable$population==population),2]
	summaryTable[which(summaryTable$population==population),5] <- xmax_Antanambao
	summaryTable[which(summaryTable$population==population),4] <- xmin_Antanambao
	
	# sako_vohiFina
	population="sako_vohiFina"
	xmin_sako_vohiFina <- xmin_bao_sakoVohiFina - summaryTable[which(summaryTable$population==population),2]
	xmax_sako_vohiFina <- xmin_bao_sakoVohiFina
	summaryTable[which(summaryTable$population==population),5] <- xmax_sako_vohiFina
	summaryTable[which(summaryTable$population==population),4] <- xmin_sako_vohiFina

	# Ambodisakoana
	population="Ambodisakoana"
	xmin_Ambodisakoana <- xmax_sako_vohiFina
	xmax_Ambodisakoana <- xmax_sako_vohiFina + summaryTable[which(summaryTable$population==population),2]
	summaryTable[which(summaryTable$population==population),5] <- xmax_Ambodisakoana
	summaryTable[which(summaryTable$population==population),4] <- xmin_Ambodisakoana

	# VohiposaSahafina
	population="VohiposaSahafina"
	xmin_VohiposaSahafina <- xmin_sako_vohiFina - summaryTable[which(summaryTable$population==population),2]
	xmax_VohiposaSahafina <- xmin_sako_vohiFina
	summaryTable[which(summaryTable$population==population),5] <- xmax_VohiposaSahafina
	summaryTable[which(summaryTable$population==population),4] <- xmin_VohiposaSahafina
	
	# maro_joll
	population="maro_joll"
	xmin_maro_joll <- xmax_line2
	xmax_maro_joll <- xmax_line2 + summaryTable[which(summaryTable$population==population),2]
	summaryTable[which(summaryTable$population==population),5] <- xmax_maro_joll
	summaryTable[which(summaryTable$population==population),4] <- xmin_maro_joll
	
	# Mmarohita
	population <- "Mmarohita"
	xmin_Mmarohita <- xmin_maro_joll - summaryTable[which(summaryTable$population==population),2]
	xmax_Mmarohita <- xmin_maro_joll
	summaryTable[which(summaryTable$population==population),5] <- xmax_Mmarohita
	summaryTable[which(summaryTable$population==population),4] <- xmin_Mmarohita

	# Mjollyae
	population <- "Mjollyae"
	xmin_Mjollyae <- xmax_maro_joll
	xmax_Mjollyae <- xmax_maro_joll + summaryTable[which(summaryTable$population==population),2]
	summaryTable[which(summaryTable$population==population),5] <- xmax_Mjollyae
	summaryTable[which(summaryTable$population==population),4] <- xmin_Mjollyae
	
	# Divide dataframe by factor 1000 to facilitate formatting
	summaryTableAdapted <- summaryTable
	summaryTableAdapted[,2:7] <- summaryTableAdapted[,2:7]/1000
	
	# Plot and assign name
	pop <- c("Mjollyae", "Mmarohita", "Sahamamy", "Andobo", "Ambodisakoana", "VohiposaSahafina", "sako_vohiFina", "Antanambao", "saha_dobo", "bao_sakoVohiFina", "maro_joll", "gerpiRoot", "root")
	colors <- c("#000000", "#000000", "#003C8D", "#CC3A91", "#E09D00", "#FA2D2A", "#ED6515" , "#F5F100", "#663B8F", "#F1AB0B", "#000000", "#AC734D", "grey90", "#000000", "#000000")

	p <- ggplot() + geom_rect(data = summaryTableAdapted, aes(xmin = ymin, xmax = ymax, ymin = 200 - xmin, ymax = 200 - xmax), fill= colors, color = "black") +
		theme_bw() + #general theme
		scale_y_continuous(expand = c(0, 1), breaks = seq(from = -60, to = 200, by = 10)) +
		scale_x_continuous(expand = c(0, 1), breaks = seq(from = 0, to = 600, by = 100)) +
		labs(x = "ka ago") + coord_cartesian(clip="off") +
		labs(y = expression(paste("N"[e], "(tick mark = 10k)"))) +
		theme(
			axis.text.y = element_blank(),
			axis.text.x = element_blank(),
			axis.title.x = element_blank(),
			axis.title.y = element_text(margin = margin(0, 0, 0, 0, 'cm'), colour ='black', size = 12, vjust=3))
			
	assign(paste0("p_",model), p)
}

## Write output
cat("#### plot_coal_models.R: Writing output ... \n")
svg(paste0(out_dir, "/coal_models.svg"), height=12, width=5)
grid.arrange(p_no_mig, p_mig, ncol=1, nrow=2)
dev.off()

## Report:
cat("\n#### plot_coal_models.R: Done with script.\n")


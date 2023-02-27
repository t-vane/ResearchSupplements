#!/usr/bin/env Rscript

#The following R code was adapted from http://www.robertlanfear.com/blog/files/concordance_factors.html)

#Packages
library(dplyr)
library(viridis)
library(ggplot2)
library(ggrepel)
library(GGally)
library(entropy)

args = commandArgs(trailingOnly=TRUE)
directory <- args[1]

# Read the data
d = read.delim(paste0(directory, "/concord.cf.stat"), header = T, comment.char='#')           
names(d)[10] = "bootstrap"
names(d)[11] = "branchlength"

# Test ILS assumptions
# First we use a slightly modified chisq function which behaves nicely when you feed it zeros
chisq = function(DF1, DF2, N){
    tryCatch({
        # Converts percentages to counts, runs chisq, gets pvalue
        chisq.test(c(round(DF1*N)/100, round(DF2*N)/100))$p.value
    },
    error = function(err) {
        # eErrors come if you give chisq two zeros
        # but here we're sure that there's no difference
        return(1.0)
    })
}

e = d %>% 
    group_by(ID) %>%
    mutate(gEF_p = chisq(gDF1, gDF2, gN)) %>%
    mutate(sEF_p = chisq(sDF1, sDF2, sN))
    
subset(data.frame(e), (gEF_p < 0.05 | sEF_p < 0.05))

write.table(e, paste0(directory, "/chisq.txt"))

# Calculate internode certainty
IC = function(CF, DF1, DF2, N){
    
    # Converts to counts
    X = CF * N / 100
    Y = max(DF1, DF2) * N / 100
        
    pX = X/(X+Y)
    pY = Y/(X+Y)
    
    IC = 1 + pX * log2(pX) +
             pY * log2(pY)

    return(IC)
}

e = e %>% 
    group_by(ID) %>%
    mutate(gIC = IC(gCF, gDF1, gDF2, gN)) %>%
    mutate(sIC = IC(sCF, sDF1, sDF2, sN))
	
write.table(e, paste0(directory, "/internode_certainty.txt"))
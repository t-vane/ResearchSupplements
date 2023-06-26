#!/usr/bin/env Rscript

cat("\n#### plot_coal_params.R: Starting script.\n")
options(scipen=10000)

## Command-line arguments
options(warn = 1)
args <- commandArgs(trailingOnly = TRUE)
sum_table <- args[1]
out_dir <- args[2]

## Report
cat("\n#### plot_coal_params.R: Summary table with divergence times, migration rates and genealogical divergence indices:", sum_table, "\n")
cat("#### plot_coal_params.R: Output directory:", out_dir, "\n\n")

## Packages
library(ggplot2)
library(gridExtra)
library(readxl)

## Step 1: Plot divergence time
# Read table
cat("#### plot_coal_params.R: Reading summary table for divergence time... \n")
sum_div <- as.data.frame(read_excel(sum_table, sheet = "div"))
sum_div$Model <- as.factor(sum_div$Model)

# Create rectangle data frame for shading
rect_left <- c(seq(from = 0.5, to = 12, by = 2))
rectangles <- data.frame(x.min = rect_left, x.max = rect_left + 2)

# Plot
cat("#### plot_coal_params.R: Plotting divergence time... \n")
div <- ggplot() + 
  geom_rect(data = rectangles, aes(xmin = x.min, xmax = x.max, ymin = -Inf, ymax = Inf), 
    fill = c("white", "grey90", "white", "grey90", "white", "grey90"), alpha = 0.5) +
  geom_point(data=sum_div, aes(x=x, y=Mean, fill=Model), size=3, 
    color=c("#5ab4ac", "#C99700", "#5ab4ac", "#C99700", "#5ab4ac", "#C99700", "#5ab4ac", "#C99700", "#5ab4ac", "#C99700", "#5ab4ac", "#C99700")) +
  geom_errorbar(data=sum_div, aes(x=x, y=Mean, ymin=lower95, ymax=upper95), width=.2,position=position_dodge(0.05), size=1, 
    color=c("#5ab4ac", "#C99700", "#5ab4ac", "#C99700", "#5ab4ac", "#C99700", "#5ab4ac", "#C99700", "#5ab4ac", "#C99700", "#5ab4ac", "#C99700")) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 12, colour = 'black', face = 'plain'),
    axis.text.x = element_text(size = 12, hjust = 1, colour = 'black', face = 'plain'),
    axis.title.x = element_blank(),
    axis.title.y = element_text(margin = margin(0, 0, 0, 0, 'cm'), colour ='black', size = 12, vjust=3),
    legend.text = element_text(size = 12),
    panel.border = element_rect(colour = "grey20", fill = NA, size = 1),
    plot.margin = margin(0.5, 0.6, 0, 0.2, 'cm')) + 
  labs(y = "Divergence time [ka] ") + 
  coord_cartesian(xlim = c(1, 12)) +
  scale_x_continuous(breaks = c(1.5,3.5,5.5,7.5,9.5,11.5), labels=c("1.5" = "A", "3.5" = "B", "5.5" = "C", "7.5" = "D", "9.5" = "E", "11.5" ="F")) +
  theme(
    legend.title = element_text(size = 12, colour = 'black', face = 'plain'),
    legend.text = element_text(size = 12, colour = 'black', face = 'plain'),
    legend.key.height = unit(0.5, "cm"),
    legend.key.width = unit(0.5, "cm"),
    legend.background = element_rect(fill = "grey90", colour = "grey30"),
    legend.key = element_rect(fill = "grey90"),
    legend.position = "top")


## Step 2: Plot migration
# Read table
cat("#### plot_coal_params.R: Reading summary table for migration ... \n")
sum_mig <- as.data.frame(read_excel(sum_table, sheet = "mig"))
sum_mig$Model <- as.factor(sum_mig$Model)

# Create rectangle data frame for shading
rect_left <- c(seq(from = 0.5, to = 12, by = 1))
rectangles <- data.frame(x.min = rect_left, x.max = rect_left + 1)

# Plot
cat("#### plot_coal_params.R: Plotting migration ... \n")
mig <- ggplot() + 
  geom_point(data=sum_mig, aes(x=x, y=Mean, fill=Model), size=3, color="#C99700") +
  geom_errorbar(data=sum_mig, aes(x=x, y=Mean, ymin=lower95, ymax=upper95), width=.2,position=position_dodge(0.05), size=1, color="#C99700") +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 12, colour = 'black', face = 'plain'),
    axis.text.x = element_text(size = 12, angle = 40, hjust = 1, colour = 'black', face = 'plain'),
    axis.title.x = element_blank(),
    axis.title.y = element_text(margin = margin(0, 0, 0, 0, 'cm'), colour ='black', size = 12, vjust=3),
    legend.text = element_text(size = 14),
    panel.border = element_rect(colour = "grey20", fill = NA, size = 1),
    plot.margin = margin(0.5, 0.6, 0, 0.2, 'cm')) + 
  labs(y = expression(paste("Population migration rate (", italic("2Nm"), ")"))) +
  coord_cartesian(xlim = c(1, 12)) +
  theme(legend.position = "none")


## Step 3: Plot genealogical divergence
# Read table
cat("#### plot_coal_params.R: Reading summary table for genealogical divergence ... \n")
sum_gdi <- as.data.frame(read_excel(sum_table, sheet = "gdi"))
sum_gdi$Model <- as.factor(sum_gdi$Model)

# Create rectangle data frame for shading
rect_left <- c(seq(from = 0.5, to = 10, by = 2))
rectangles <- data.frame(x.min = rect_left, x.max = rect_left + 2)

# Plot
cat("#### plot_coal_params.R: Plotting genealogical divergence ... \n")
gdi <- ggplot() + 
  geom_rect(data = rectangles, aes(xmin = x.min, xmax = x.max, ymin = -Inf, ymax = Inf), 
    fill = c("white", "grey90", "white", "grey90", "white"), alpha = 0.5) +
  geom_point(data=sum_gdi, aes(x=x, y=Mean, fill=Model), size=3, color="#5ab4ac") +
  geom_errorbar(data=sum_gdi, aes(x=Category, y=Mean, ymin=lower95, ymax=upper95), width=.2,position=position_dodge(0.05), size=1, color="#5ab4ac") +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 12, colour = 'black', face = 'plain'),
    axis.text.x = element_text(size = 12, , colour = 'black', face = 'plain'),
    axis.title.x = element_blank(),
    axis.title.y = element_text(margin = margin(0, 0, 0, 0, 'cm'), colour ='black', size = 12, vjust=3),
    legend.text = element_text(size = 12),
    panel.border = element_rect(colour = "grey20", fill = NA, size = 1),
    plot.margin = margin(0.5, 0.6, 0, 0.2, 'cm')) + 
  labs(y = expression(paste("Genealogical divergence index (", italic("gdi"), ")"))) + 
  coord_cartesian(xlim = c(1, 10), ylim = c(0,1)) +
   scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
  theme(
    legend.title = element_text(size = 12, colour = 'black', face = 'plain'),
    legend.text = element_text(size = 12, colour = 'black', face = 'plain'),
    legend.key.height = unit(0.5, "cm"),
    legend.key.width = unit(0.5, "cm"),
    legend.background = element_rect(fill = "grey90", colour = "grey30"),
    legend.key = element_rect(fill = "grey90"),
    legend.position = "none") +
  geom_hline(yintercept=0.2, linetype="dashed", 
             color = "black", size=0.8) +
  geom_hline(yintercept=0.7, linetype="dashed", 
             color = "black", size=0.8)

## Write output
cat("#### plot_coal_params.R: Writing output ... \n")
svg(paste0(out_dir, "/coal_params.svg"), width=14, height=5)
grid.arrange(div,mig,gdi, ncol=3, nrow=1)
dev.off()

## Report:
cat("\n#### plot_coal_params.R: Done with script.\n")


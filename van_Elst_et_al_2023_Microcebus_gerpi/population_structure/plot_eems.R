#!/usr/bin/env Rscript

cat("\n#### plot_eems.R: Starting script.\n")

## Command-line arguments
options(warn = 1)
args <- commandArgs(trailingOnly = TRUE)
input <- args[1]
out_dir <- args[2]
pop_coords <- args[3]
shape <- args[4]

## Report
cat("\n#### plot_eems.R: Number of demes:", ndemes, "\n")
cat("#### plot_eems.R: Input prefix:", input, "\n")
cat("#### plot_eems.R: Output directory:", out_dir, "\n")
cat("#### plot_eems.R: File with population coordinates directory:", pop_coords, "\n")
cat("#### plot_eems.R: Shape file prefix:", shape, "\n\n")

## Packages
library("rgdal")
library("rworldmap")
library("rworldxtra")
library("maps")
library("mapdata")
library("sp")
library("maptools")
library("raster")
library("rgeos")
library("rgdal")
library("scales")
library("ggplot2")
library("mapplots")
library("rEEMSplots")

## Process command-line arguments
cat("#### plot_eems.R: Reading population coordinates ... \n")
coords <- as.matrix(read.table(pop_coords, row.names=1))
coords <- coords[,c(2, 1)]
cat("#### plot_eems.R: Loading shape file ... \n")
madarivers_major_shapes <- readOGR(dirname(shape), basename(shape))
madarivers_major_shapes_t = spTransform(madarivers_major_shapes, "+proj=merc +datum=WGS84")

## Prepare map
cat("#### plot_eems.R: Preparing map ... \n")
projection_none <- "+proj=longlat +datum=WGS84"
projection_mercator <- "+proj=merc +datum=WGS84"
map_world <- getMap()
map_africa <- map_world[which(map_world@data$continent == "Africa"), ]
map_africa <- spTransform(map_africa, CRSobj = CRS(projection_mercator))
coords_merc <- sp::spTransform(SpatialPoints(coords, CRS(projection_none)),CRS(projection_mercator))
coords_merc <- coords_merc@coords

## Plot deltaK
cat("#### plot_eems.R: Plotting estimated effective migration surface ... \n")
# Assign population colors
colors <- c("#E09D00", "#CC3A91", "#7F3e9B", "#F5F100", "black", "#F50035", "#003C8D", "#FF5A1F")
# Plot
pdf(paste0(out_dir, "/plot_eems.pdf"), 16, 14)
eems.plots(mcmcpath = input, plotpath = out_dir, longlat = FALSE,
           plot.height = 8, plot.width = 7, res = 600, add.outline= TRUE, col.outline = "black", lwd.outline = 3,
           projection.in=projection_none, projection.out=projection_mercator,
           add.map = TRUE, col.map = "black", lwd.map = 3, out.png = FALSE, 
           m.plot.xy = { plot(madarivers_major_shapes_t, col=alpha("royalblue", 1), border=FALSE, add=TRUE, lwd=2);
             points(coords_merc, pch = 21, cex= 3.5, col = "black", bg = colors, lwd = 1); },
           q.plot.xy = { plot(madarivers_major_shapes_t, col=alpha("royalblue", 1), border=FALSE, add=TRUE, lwd=2);
             points(coords_merc, pch = 21, cex= 3., col = "black", bg = colors, lwd = 1); }
)
dev.off()

## Report:
cat("\n#### plot_eems.R: Done with script.\n")


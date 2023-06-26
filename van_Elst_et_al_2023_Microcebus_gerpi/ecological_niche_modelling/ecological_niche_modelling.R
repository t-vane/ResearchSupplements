#### Script by D. Schüßler


#### packages ####
library(raster)     # handling of raster data
library(rJava)      # to work with maxent in R
library(ENMTools)   # ENM using maxent algorithm
library(ENMeval)    # ENM validation and fine tuning
library(RStoolbox)  # calculation of PCA on raster data
library(dplyr)      # data handling
library(readxl)     # loading data from excel table
library(ecospat)    # CBI calculation in 
library(biomod2)    # ENM calculation using random forest algorithm


#### data preparation ####
setwd("E://136-Distribution_of_M_gerpi") # set working directory
set.seed(21) # set seed for random operation in the script
options(max.print=100000) # increase print out for maxent models
data <- read_excel("C://Users//domin//OneDrive//Dokumente//019-Distribution_of_M_gerpi//vanElstetal_EcologyandEvolution_SupplementaryTables.xlsx", 
                   sheet = "S4 - Occurrence records", col_names = TRUE) # load in occurrence data of Microcebus spp.

PO_data <- filter(data, Species == "M. gerpi")     # Presence-only data for M. gerpi
PO_data <- dplyr::select(PO_data, Species, Longitude, Latitude)       # downsize dataframe
PA_data <- dplyr::select(data, Species, Longitude, Latitude) 

# preparation of environmental data; if no bioclimatic variables are available, load and stack provided PC1-3 below to continue
env_variables <- stack(c(bio01 = "E://999_Bioclim_Mada//Chelsa_2.1//CH_bio01.tif",
                         bio02 = "E://999_Bioclim_Mada//Chelsa_2.1//CH_bio02.tif",
                         bio03 = "E://999_Bioclim_Mada//Chelsa_2.1//CH_bio03.tif",
                         bio04 = "E://999_Bioclim_Mada//Chelsa_2.1//CH_bio04.tif",
                         bio05 = "E://999_Bioclim_Mada//Chelsa_2.1//CH_bio05.tif",
                         bio06 = "E://999_Bioclim_Mada//Chelsa_2.1//CH_bio06.tif",
                         bio07 = "E://999_Bioclim_Mada//Chelsa_2.1//CH_bio07.tif",
                         bio08 = "E://999_Bioclim_Mada//Chelsa_2.1//CH_bio08.tif",
                         bio09 = "E://999_Bioclim_Mada//Chelsa_2.1//CH_bio09.tif",
                         bio10 = "E://999_Bioclim_Mada//Chelsa_2.1//CH_bio10.tif",
                         bio11 = "E://999_Bioclim_Mada//Chelsa_2.1//CH_bio11.tif",
                         bio12 = "E://999_Bioclim_Mada//Chelsa_2.1//CH_bio12.tif",
                         bio13 = "E://999_Bioclim_Mada//Chelsa_2.1//CH_bio13.tif",
                         bio14 = "E://999_Bioclim_Mada//Chelsa_2.1//CH_bio14.tif",
                         bio15 = "E://999_Bioclim_Mada//Chelsa_2.1//CH_bio15.tif",
                         bio16 = "E://999_Bioclim_Mada//Chelsa_2.1//CH_bio16.tif",
                         bio17 = "E://999_Bioclim_Mada//Chelsa_2.1//CH_bio17.tif",
                         bio18 = "E://999_Bioclim_Mada//Chelsa_2.1//CH_bio18.tif",
                         bio19 = "E://999_Bioclim_Mada//Chelsa_2.1//CH_bio19.tif")) # loading of bioclimatic layers

gerpi_study_region <- shapefile("gerpi_study_region.shp") # load study region as shape file
env_variables <- mask(env_variables, gerpi_study_region)  # clip environmental variables to study region

pca1 <- rasterPCA(env_variables, spca = TRUE) # calculate PCS
summary(pca1$model) # display summary of PCA, select first three PCs

PC1 <- pca1$map$PC1 # convert PCs to single raster layers
PC2 <- pca1$map$PC2
PC3 <- pca1$map$PC3

# => PC1-3 can be loaded here from folder to be stacked in next line (e.g., PC1 <- raster("C://Users//...//PC1.tif"))
climate_variables <- stack(PC1, PC2, PC3) # stack first 3 PCs together


#### Maxent model with presence-only data ####

# parameter fine-tuning for maxent models, AIC-based
validation <- ENMevaluate(PO_data[,2:3], envs = climate_variables,
                          tune.args = list(fc = c("L", "P", "Q", "LP", "LQ"), # try different feature classes
                                           rm = 1:6),                         # try different rendomization multipliers
                          n.bg = 1000, partitions = "jackknife",              # 1,000 background points, n-1 jackknife cross-validation
                          algorithm = 'maxnet', overlap = FALSE, 
                          bin.output = FALSE, clamp = FALSE, 
                          parallel = TRUE, progbar = TRUE, numCores = 4)
validation@results # show results, choose model with deltaAIC = 0, use validation metrics for reporting

enm.gerp <- enmtools.species()   # creating an enmtools.species object for later modeling
enm.gerp$species.name <- "M.gerpi"
enm.gerp$presence.points <- PO_data[,2:3]
enm.gerp$range <- background.raster.buffer(enm.gerp$presence.points, 100000, mask = climate_variables)
enm.gerp$background.points <- background.points.buffer(points = enm.gerp$presence.points,
                                                       radius = 100000, n = 1000, mask = climate_variables[[1]])

args.enm.gerp = c("betamultiplier=2", "linear=TRUE", "quadratic=FALSE",
                  "product=FALSE", "hinge=FALSE", "threshold=FALSE", "autofeature=FALSE") # define maxent parameters
gerp.mx <- enmtools.maxent(enm.gerp, climate_variables, test.prop = 0.0, args = args.enm.gerp, clamp = FALSE) # maxent model
plot(gerp.mx) # see model on map (geographic space)

visualize.enm(gerp.mx, climate_variables, layers = c("PC1", "PC2"), plot.points = TRUE) # plot model in environmental space, for checking plausibility

CBI <- ecospat.boyce(gerp.mx$suitability, gerp.mx$analysis.df[gerp.mx$analysis.df$presence == 1,1:2]) # calculate CBI metric for validation
CBI             # CBI 
gerp.mx         # AUC
gerp.mx$model   # for creating binary maps: 10% and minimum training presence (10pctTP, MTP)

# Raes & ter Steege test
rave.mx_rts <- enmtools.maxent(enm.gerp, climate_variables, rts.reps = 10, test.prop = 0.3, 
                               args = args.enm.gerp, clamp = FALSE) # args defined above
rave.mx_rts$rts.test$rts.pvalues # show results

#### Save results as raster for GIS
suit.map.gerp <- gerp.mx$suitability
writeRaster(suit.map.gerp, filename= "gerpi_maxent", format="GTiff", overwrite = TRUE) # give species name here


#### Random forest with presence-absence data ####
PA_data$species <- recode(PA_data$species, gerpi = "1", lehilahytsara = "0") # recode for modeling
PA_data$species <- as.numeric(PA_data$species)
RF_data <- BIOMOD_FormatingData(resp.var = PA_data$species, expl.var = climate_variables,
                                  resp.xy = PA_data[, c("long", "lat")], resp.name = "PA_gerpi",
                                  PA.nb.rep = 0, na.rm = TRUE)
RF_data # check output of formatting

# set model option beforehand
opt <- BIOMOD_ModelingOptions(GLM = list(type = "quadratic", interaction.level = 1, family = "binomial")) # more options can be predifened, also if other models will be used!

unregister_dopar <- function() {    # to un-register parallel computing; otherwise leads to errors
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}
unregister_dopar()

# modeling
RF_model <- BIOMOD_Modeling(RF_data,
                          models = c("GLM", "ANN", "RF"), # Generalized linear model, artificial neural network, random forest
                          bm.options = opt,
                          nb.rep = 2,
                          data.split.perc = 70,
                          weights = NULL,
                          prevalence = 0.5,
                          var.import = 99,
                          metric.eval = c("KAPPA", "TSS", "ROC"),
                          save.output = TRUE,
                          scale.models = FALSE,
                          do.progress = TRUE)

# save/load calculated models
save(RF_model, file="RF_model.RData")
load("RF_model.RData")

model_scores <- get_evaluations(RF_model)
model_scores # shows evaluation scores

VarImp <- get_variables_importance(RF_model) # extract variable importance from model
VarImp # show it
apply(VarImp, c(1,2), mean) # show simplified (best use if several models were calculated)

GF_RF <- BIOMOD_LoadModels(RF_model, models = "RF")

eval_strip <- bm_PlotResponseCurves(bm.out = RF_model, new.env = get_formal_data(RF_model, "expl.var"),
                                 show.variables = get_formal_data(RF_model, "expl.var.names"),
                                 do.bivariate = FALSE, fixed.var = "median", legend = FALSE,
                                 data_species = get_formal_data(RF_model, "resp.var"))

gc()

get_built_models(RF_model)
RF_proj <- BIOMOD_Projection(RF_model, new.env = climate_variables, proj.name = "PA_gerpi",
                             output.format = ".img",
                             models.chosen = "PA.gerpi_AllData_Full_RF")
plot(RF_proj)

# => browse for proj_PA_gerpi_PA.gerpi and select band 9 for RF model in GIS





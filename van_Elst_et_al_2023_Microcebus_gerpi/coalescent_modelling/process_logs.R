#!/usr/bin/env Rscript

## Script adapted and modified from Poelstra et al. 2021, Systematic Biology (https://doi.org/10.1093/sysbio/syaa053)

cat("\n#### process_logs.R: Starting script.\n")

## Command-line arguments
options(warn = 1)
args <- commandArgs(trailingOnly = TRUE)
model <- args[1]
run_id <- args[2]
mcmc <- args[3]
burnin <- as.double(args[4])
last_sample <- as.double(args[5])
mutrate <- as.double(args[6])
mutrate_var <- as.double(args[7])
gentime <- as.double(args[8])
gentime_sd <- as.double(args[9])
m_scale <- as.double(args[10])
t_scale <- as.double(args[11])
poplist_file <- args[12]

## Report
cat("\n#### process_logs.R: Model:", model, "\n")
cat("#### process_logs.R: Run ID:", run_id, "\n")
cat("#### process_logs.R: MCMC file:", mcmc, "\n")
cat("#### process_logs.R: Burn-in:", burnin, "\n")
cat("#### process_logs.R: Last sample to keep:", last_sample, "\n")
cat("#### process_logs.R: Gamma distribution of mutation rate will have mean:", mutrate, "* 10 e-8 \n")
cat("#### process_logs.R: Gamma distribution of mutation rate will have variance:", mutrate_var, "* 10 e-8 \n")
cat("#### process_logs.R: Lognormal distribution of generation time will have mean: ln(", gentime, ")\n")
cat("#### process_logs.R: Lognormal distribution of generation time will have standard deviation: ln(", gentime_sd, ")\n")
cat("#### process_logs.R: Inverse scaling factor used in the G-PhoCS configuration file for migration parameter:", m_scale, "\n")
cat("#### process_logs.R: Inverse scaling factor used in the G-PhoCS configuration file for tau and theta:", t_scale, "\n")
cat("#### process_logs.R: File with information on parent and child populations:", poplist_file, "\n\n")

## Packages
if(!'pacman' %in% rownames(installed.packages())) install.packages('pacman')
packages <- c('data.table', 'TeachingDemos', 'tidyverse')
pacman::p_load(char = packages, install = TRUE)

## Process command-line arguments
logdir <- dirname(mcmc)
logfile <- basename(mcmc)
poplist_prelim <- read.table(poplist_file, header = TRUE)
poplist <- as.list(setNames(poplist_prelim$parent, poplist_prelim$child))

## Step 1: Define function to open one logfile and cut off burn-in and/or last samples
cut_log <- function(mcmc, burnin, last_sample, subsample, return_log = TRUE, write_log = TRUE) {

  Log <- read.table(mcmc, header = TRUE)
  cat('\n\n## cut_log():', logfile, "\tnrows:", nrow(Log), '\n')

  # Remove burn-in
  if(burnin > 0) Log <- Log[-which(Log$Sample < burnin), ]
  if(nrow(Log) < 1000) {
    warning('\n## cut_log(): Skipping file: less than 1000 rows left aber removing burn-in \n')
    return(NULL)
  }

  # Cut after last sample
  if(!is.null(last_sample) && any(Log$Sample > last_sample)) {
    Log <- Log[-which(Log$Sample > last_sample), ]
  }
  
  # Subsample from the Log (to retain only every nth sample, where n=subsample
  Log <- Log[seq(from = 1, to = nrow(Log), by = subsample), ]
  
  # Write and return log
  if(write_log == TRUE) {
    if(!dir.exists(paste0(logdir, '/cut/'))) dir.create(paste0(logdir, '/cut/'))
    write.table(Log, paste0(logdir, '/cut/cut.', logfile),
                sep = '\t', quote = FALSE, row.names = FALSE)
  }

  if(return_log == TRUE) return(Log)
}

## Step 2: Define function to melt dataframe and prepare variables
prep_log <- function(Log, logfile, model, run_id, poplist, lookup, mutrate, mutrate_var, gentime, gentime_sd, m_scale, t_scale, gdi = TRUE, 
  rename_pops_before = FALSE, # Rename populations (using "lookup" dataframe) before converting to long format, i.e. in column names
  rename_pops_after = FALSE # Rename populations after converting to long format, i.e. in dataframe itself
  ) {

  # Rename columns (second row renames last three columns):
  colnames(Log) <- gsub('\\.\\.', 2, colnames(Log))
  colnames(Log)[(ncol(Log)-2):ncol(Log)] <- c('mut_NA', 'dataLd_NA', 'FullLd_NA')

  # Rename populations (before)
  if(rename_pops_before == TRUE) {
    varcols <- grep('theta_|tau_', colnames(Log))
    for(col_idx in varcols) for(pop_idx in 1:nrow(lookup)) {
      pop <- gsub('.*_(.*)', '\\1', colnames(Log)[col_idx])
      long <- lookup$popname.long[pop_idx]
      short <- lookup$popname.short[pop_idx]
      if(pop == long) colnames(Log)[col_idx] <- sub(long, short, colnames(Log)[col_idx])
    }
    cat("\n## prep_log(): Renamed colnames:", colnames(Log), '\n')
  }

  # MIGRATION:
  # Identify migration columns
  migcols <- grep('m_', colnames(Log))
  
  for(migcol in migcols) {
    # Extract migration column
	m <- Log[, migcol]

    # Get migration pattern (e.g., VohiposaSahafina2Antanambao)
	migpattern <- unlist(strsplit(colnames(Log)[migcol], split = '_'))[2]
	# Get destination population (e.g., Antanambao)
    migto <- unlist(strsplit(migpattern, split = '2'))[2]
	# Rename destination population if specified
    if(rename_pops_before == TRUE) migto <- poprename(migto, lookup)
	# Get theta of destination population
    migto.column <- grep(paste0('theta_', migto, '$'), colnames(Log))
    th <- Log[, migto.column]

    # Estimate population migration rate
    # Get name for population migration rate
	colname.mig2 <- paste0('2Nm_', migpattern)
	# Calculate population migration rate
	Log$newcolumn1 <- (m / m_scale) * (th / t_scale) / 4
    colnames(Log)[grep('newcolumn1', colnames(Log))] <- colname.mig2

    # Estimate proportion of migrants
	# Define shape parameter for gamma function
    my_shape <- mutrate^2 / mutrate_var
	# Define rate parameter for gamma function
    my_rate <- mutrate / mutrate_var
	# Create gamma distribution
    mutrate_dist <- rgamma(nrow(Log), shape = my_shape, rate = my_rate) * 1e-8

	# Get name for proportion of migrants
    colname.mig3 <- paste0('m.prop_', migpattern)
	# Calculate proportion of migrants in percent
    Log$newcolumn2 <- (m / m_scale) * mutrate_dist * 100
	# Rename column
    colnames(Log)[grep('newcolumn2', colnames(Log))] <- colname.mig3
  }

  # GENEALOGICAL DIVERGENCE:
  if(gdi == TRUE) {
	# Get list of child population names
    childpops <- names(poplist)
	# Get list of parent population names
    parentpops <- as.character(poplist)
	
	#define function to get parent population for a population (level 2 is for grandparent, etc.)
    get_parent <- function(pop, level = 1) {
      if(level == 1) parentpop <- poplist[[pop]]
      if(level == 2) parentpop <- poplist[[poplist[[pop]]]]
	  if(level == 3) parentpop <- poplist[[poplist[[poplist[[pop]]]]]]
      return(parentpop)
    }
	
	# Define function to calculate genealogical divergence index (gdi)
    get_gdi <- function(pop, level) {
      parentpop <- get_parent(pop, level = level)
      theta <- Log[, paste0('theta_', pop)]
      tau <- Log[, paste0('tau_', parentpop)]
      gdi <- 1 - exp((-2 * tau) / theta)
    }
	
	# Run get_gdi() on level 1 to get gdis to sister populations
    for(pop in childpops) {
      Log$gdi_new <- get_gdi(pop, level = 1)
      colnames(Log)[grep('gdi_new', colnames(Log))] <- paste0('gdi_', pop)
    }
	
	# Run get_gdi() on level 2 to get gdis to more distantly related populations
    for(pop in childpops) {
      if(!is.null(get_parent(pop, level = 2))) {
        Log$gdi_new <- get_gdi(pop, level = 2)
        colnames(Log)[grep('gdi_new', colnames(Log))] <- paste0('gdi2_', pop)
      }
    }
	
	# Run get_gdi on level 3 to get gdis to even more distantly related populations
	for (pop in childpops) {
		if (get_parent(pop, 1) != "root") {
			if (get_parent(pop, 2) != "root") {
				if (get_parent(pop, 3) != "root") {
					Log$gdi_new <- get_gdi(pop, level = 3)
					colnames(Log)[grep('gdi_new', colnames(Log))] <- paste0('gdi3_', pop)
				}
			}
		}
	}
  }
	
	# EFFECTIVE POPULATION SIZE:
	thetacols <- grep('theta_', colnames(Log))
	
	for(thetacol in thetacols) {
	
	th <- Log[, thetacol]
	
    # Get name for effective population size
	popname <- unlist(str_split(colnames(Log)[thetacol], '_', n=2))[2]
	colname.popsize <- paste0('Ne_', popname)
	
	# Define shape parameter for gamma function
	my_shape <- mutrate^2 / mutrate_var
	# Define rate parameter for gamma function
    my_rate <- mutrate / mutrate_var
	# Create gamma distribution
    mutrate_dist <- rgamma(nrow(Log), shape = my_shape, rate = my_rate) * 1e-8
	
	# Estimate effective population size
    Log$newcolumn3 <- (th / t_scale) / (4 * mutrate_dist)
    colnames(Log)[grep('newcolumn3', colnames(Log))] <- colname.popsize
	}
	
	# DIVERGENCE TIME:
	taucols <- grep('tau_', colnames(Log))
	
	for(taucol in taucols) {
	
	ta <- Log[, taucol]
	
    # Get name for divergence time
	popname <- unlist(str_split(colnames(Log)[taucol], '_', n=2))[2]
	colname.divTime <- paste0('div.time_', popname)
	
	# Define generation time lognormal distribution
	gentime_dist <- rlnorm(nrow(Log), meanlog = log(gentime), sdlog = log(gentime_sd))
	# Define shape parameter for gamma function
	my_shape <- mutrate^2 / mutrate_var
	# Define rate parameter for gamma function
    my_rate <- mutrate / mutrate_var
	# Create gamma distribution
    mutrate_dist <- rgamma(nrow(Log), shape = my_shape, rate = my_rate) * 1e-8
	
	# Estimate divergence time
    Log$newcolumn4 <- (ta / t_scale) * gentime_dist / mutrate_dist
    colnames(Log)[grep('newcolumn4', colnames(Log))] <- colname.divTime
	}
	
  # Write output
	if(!dir.exists(paste0(logdir, '/cut/'))) dir.create(paste0(logdir, '/cut/'))
		write.table(Log, paste0(logdir, '/cut/prep.', logfile),
			sep = '\t', quote = FALSE, row.names = FALSE)

  # Rename populations(after)
  if(rename_pops_after == TRUE) {
    if(rename_pops_before == FALSE) {
      warning('## prep_log(): These pops were not found in the lookup:',
              unique(mlog$pop)[!unique(mlog$pop) %in% lookup$popname.long], '\n')
      mlog$pop <- poprename(mlog$pop, lookup)
    }
    mlog$migfrom <- poprename(mlog$migfrom, lookup)
    mlog$migto <- poprename(mlog$migto, lookup)
  }
}

## Step 3: Define single-log wrapper function: cut_log(), then prep_log()
cutprep_log <- function(logfile, logdir, model, run_id, poplist, lookup,
  burnin, last_sample, subsample, cut = TRUE, gdi = TRUE, mutrate, mutrate_var, gentime, gentime_sd, m_scale, t_scale,
  rename_pops_before = FALSE, rename_pops_after = FALSE, return_log = TRUE, write_log = TRUE
  ) {

  if(cut == TRUE) {
    Log <- cut_log(
      logfile = logfile, logdir = logdir, burnin = burnin,
      last_sample = last_sample, subsample = subsample,
      return_log = return_log, write_log = write_log
      )
  }

  if(cut == FALSE) {
    cat('\n\n## cutprep_log(): Not cutting log, just reading ...')
    cat('## cutprep_log(): ', logfile, "\tnrows:", nrow(Log))
    Log <- read.table(paste0(logdir, '/cut/', logfile), header = TRUE)
  }

  if(!is.null(Log)) {
    cat('## cutprep_log(): Column names of Log:\n', colnames(Log), '\n')

    Log <- prep_log(
      Log = Log, logfile = logfile,
      model = model, run_id = run_id, poplist = poplist, lookup = lookup,
      rename_pops_before = rename_pops_before, rename_pops_after = rename_pops_after,
      mutrate = mutrate, mutrate_var = mutrate_var, gentime = gentime, gentime_sd = gentime_sd, m_scale = m_scale, t_scale = t_scale
      )
    return(Log)
  }
}

## Step 4: Run single-log wrapper function
cat("#### process_logs.R: Processing logs and converting to demographic values ... \n")
cutprep_log(logfile=logfile, logdir=logdir, burnin=burnin, last_sample=last_sample, model=model, run_id=run_id,
			poplist=poplist, mutrate=mutrate, mutrate_var=mutrate_var, gentime=gentime, gentime_sd=gentime_sd,
			m_scale=m_scale, t_scale=t_scale)

## Report:
cat("\n#### process_logs.R: Done with script.\n")


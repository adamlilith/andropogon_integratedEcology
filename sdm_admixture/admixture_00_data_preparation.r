### MODELING ANDROPOGON GERARDI DISTRIBUTION, PHENOTYPE, PHYSIOLOGY, GENOTYPE, and ASSOCIATED MICROBIAL COMMUNITIES
### Erica Newman | Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2023-12
###
### This script compiles response and environmental data for an integrated model combining occurrences of Andropogon gerardi and results from an ADMIXTURE analysis.
###
### source('C:/Ecology/R/andropogon_integratedEcology/sdm_admixture/admixture_00_data_preparation.r')
### source('C:/Subarashi/R/andropogon_integratedEcology/sdm_admixture/admixture_00_data_preparation.r')
###
### CONTENTS ###
### setup ###
### match ADMIXTURE results to coordinates and environment ###

#############
### setup ###
#############

	rm(list = ls())

	drive <- 'C:/Ecology/'
	# drive <- 'C:/Subarashi/'

	library(data.table)
	library(enmSdmX)
	library(omnibus)
	library(readxl)
	library(terra)

	setwd(paste0(drive, '/Research/Andropogon/Andropogon'))

say('##############################################################')
say('### match ADMIXTURE results to coordinates and environment ###')
say('##############################################################')

	### ADMIXTURE results
	admix <- read_xlsx('./data_from_loretta/jack_sytsma_2024_05_16_admixture_analysis/Admixture Spatial Sampling 2023.xlsx', sheet = 'Sheet1')

	admix <- as.data.table(admix)
	names(admix) <- tolower(names(admix))

	### add coordinates
	# will use the pre-prepared microbial data to get coordinates since these have been checked by Erica
	xy <- fread('./data_from_sonny/erica_harmonized_taxa_between_files_2024_07_26/environment_rhz_26JUL2024.csv')
	xy <- xy[ , c('site', 'latitude', 'longitude')]
	xy <- aggregate(xy, by = list(xy$site), FUN = mean)
	xy$site <- NULL
	names(xy)[1] <- 'site'

	### rename ADMIXTURE sites so they correspond to data frame with coordinates
	admix$population_proper <- NA_character_
	for (i in 1:nrow(admix)) {

		improper <- admix$population[i]
		this <- grepl(xy$site, pattern = improper)
		if (sum(this) != 1) stop('More than one site match.')
		admix$population_proper[i] <- xy$site[this]
	
	}

	### add coordinates
	matches <- match(admix$population_proper, xy$site)
	admix$longitude <- xy$longitude[matches]
	admix$latitude <- xy$latitude[matches]

	### match with soil and climate data
	cols <- c('aridity', 'bio7', 'bio12', 'bio15', 'elevation_m', 'sampling_ppt_mm', 'sampling_tmean_c', 'sand_field', 'silt_field', 'clay_field', 'ph_field', 'nitrogen_field_perc', 'soc_field_perc', 'ph_soilgrids', 'sand_soilgrids', 'silt_soilgrids', 'clay_soilgrids', 'soc_soilgrids_perc', 'nitrogen_soilgrids_perc')

	for (j in seq_along(cols)) {
		admix[ , DUMMY := NA_real_]
		names(admix)[ncol(admix)] <- cols[j]
	}

	env <- fread('./data_from_sonny/collated_for_hmsc/environment_combined.csv')
	for (i in 1:nrow(admix)) {
	
		site <- admix$population_proper[i]
		index <- which(env$location == site)
		index <- index[1]
		site_env <- env[index]

		for (j in seq_along(cols)) {
			col <- cols[j]
			admix[i, (col) := site_env[[col]]]
		}
	
	}

	dirCreate('./data_from_loretta/!admixture_data_collated')
	write.csv(admix, './data_from_loretta/!admixture_data_collated/admixture_with_env.csv', row.names = FALSE)

say('FINIS', deco = '~', level = 1)



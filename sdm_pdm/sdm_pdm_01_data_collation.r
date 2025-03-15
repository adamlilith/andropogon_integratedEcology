### MODELING ANDROPOGON GERARDI DISTRIBUTION, PHENOTYPE, PHYSIOLOGY, GENOTYPE, and ASSOCIATED MICROBIAL COMMUNITIES
### Erica Newman | Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2023-12
###
### This script constructs a non-integrated model for AG geographic distribution.
###
### source('C:/Ecology/R/andropogon_integratedEcology/sdm_pdm/sdm_pdm_01_data_collation.r')
### source('C:/Subarashi/R/andropogon_integratedEcology/sdm_pdm/sdm_pdm_01_data_collation.r')
###
### CONTENTS ###
### setup ###
### collate morphology/physiology/community data with climate ###

#############
### setup ###
#############

	rm(list = ls())

	# drive <- 'C:/Ecology/'
	drive <- 'C:/Subarashi/'

	setwd(paste0(drive, '/Research/Andropogon/Andropogon'))

	library(data.table) # fast data tables
	library(enmSdmX) # GIS & SDMing
	library(omnibus) # utilities
	library(predicts) # GIS & SDMing
	library(readxl) # Excel
	library(terra) # spatial objects

say('#################################################################')
say('### collate morphology/physiology/community data with climate ###')
say('#################################################################')

	### AG data

	# sites: coordinates, mean values
	sites <- read_excel('./data_from_loretta/!plant_sitelevel_data_11NOV2024 - Google Sheets [aggregated by Erica].xlsx', sheet = 'plantmaster_bysite_19NOV2024')

	# biomass
	biomass <- read_excel('./data_from_loretta/!plant_sitelevel_data_11NOV2024 - Google Sheets [aggregated by Erica].xlsx', sheet = 'biomass_and_size')

	# morphology & physiology
	morpho_phys <- read_excel('./data_from_loretta/!plant_sitelevel_data_11NOV2024 - Google Sheets [aggregated by Erica].xlsx', sheet = 'morphology_and_phys')

	sites <- as.data.table(sites)
	biomass <- as.data.table(biomass)
	morpho_phys <- as.data.table(morpho_phys)

	### match sites with environments
	ppt <- rast(paste0(drive, '/Research Data/ClimateNA/v 7.3 AdaptWest/1991-2020/', paste0('Normal_1991_2020_PPT', prefix(1:12, 2), '.tif')))
	tmin <- rast(paste0(drive, '/Research Data/ClimateNA/v 7.3 AdaptWest/1991-2020/', paste0('Normal_1991_2020_Tmin', prefix(1:12, 2), '.tif')))
	tmax <- rast(paste0(drive, '/Research Data/ClimateNA/v 7.3 AdaptWest/1991-2020/', paste0('Normal_1991_2020_Tmax', prefix(1:12, 2), '.tif')))

	### extract environment at sites
	sites_spatial <- vect(sites, geom = c('LONGITUDE', 'LATITUDE'), crs = getCRS('NAD83'), keepgeom = TRUE)
	sites_spatial <- project(sites_spatial, ppt)

	ppt_at_sites <- extract(ppt, sites_spatial, ID = FALSE)
	tmax_at_sites <- extract(tmax, sites_spatial, ID = FALSE)
	tmin_at_sites <- extract(tmin, sites_spatial, ID = FALSE)

	ppt_at_sites <- as.matrix(ppt_at_sites)
	tmax_at_sites <- as.matrix(tmax_at_sites)
	tmin_at_sites <- as.matrix(tmin_at_sites)

	bioclims <- bcvars(ppt_at_sites, tmin_at_sites, tmax_at_sites)
	bioclims <- as.data.frame(bioclims)
	bioclims$aridity <- (bioclims$bio1 + 10) / (bioclims$bio12 / 1000) # add aridity
	
	sites <- cbind(sites, bioclims)

	### create slightly different data frames for each type of data
	# each will have coordinates and BIOCLIM predictors that match to site
	# need to do this because different number of plants were sampled for different metrics

	bioclims <- as.data.table(bioclims)
	bioclims[ , SITE := sites$site_id]


	biomass <- merge(biomass, bioclims, by = 'SITE', suffixes = '')
	morpho_phys <- merge(morpho_phys, bioclims, by = 'SITE', suffixes = '')

	dirCreate('./data_from_loretta/sdm_pdm_00_merged_site_data_with_climate')
	saveRDS(sites, './data_from_loretta/sdm_pdm_00_merged_site_data_with_climate/sites.rds')
	saveRDS(biomass, './data_from_loretta/sdm_pdm_00_merged_site_data_with_climate/biomass.rds')
	saveRDS(morpho_phys, './data_from_loretta/sdm_pdm_00_merged_site_data_with_climate/morpho_phys.rds')

say(date())
say('FINIS!', deco = '+', level = 1)

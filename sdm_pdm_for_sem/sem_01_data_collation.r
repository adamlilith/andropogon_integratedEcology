### MODELING ANDROPOGON GERARDI DISTRIBUTION, PHENOTYPE, PHYSIOLOGY, GENOTYPE, and ASSOCIATED MICROBIAL COMMUNITIES
### Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2023-12
###
### This script constructs a non-integrated model for AG geographic distribution.
###
### source('C:/Kaji/R/andropogon_integratedEcology/sdm_pdm_for_sem/sem_01_data_collation.r')
###
### CONTENTS ###
### setup ###
### collate morphology/physiology/community data with climate ###

#############
### setup ###
#############

	rm(list = ls())

	drive <- 'C:/Kaji/'
	setwd(paste0(drive, '/Research/Andropogon/Andropogon'))
	.libPaths(paste0(getwd(), '/libraries'))

	library(data.table) # fast data tables
	library(enmSdmX) # GIS & SDMing
	library(ggplot2) # graphics
	library(ggspatial) # GIS graphics
	library(omnibus) # utilities
	library(predicts) # GIS & SDMing
	library(readxl) # Excel
	library(terra) # spatial objects

say('#################################################################')
say('### collate morphology/physiology/community data with climate ###')
say('#################################################################')

	### data

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

	# pisr <- rast('E:/ecology/Potential Annual Insolation (SAGA)/Based on ClimateNA 7.03 from SAGA 9.3.0/Annual_Insolation_1990_kW_hr_per_m2.tif')
	# pisr <- extract(pisr, sites_spatial, ID = FALSE)
	# pisr <- rowSums(pisr)

	sites <- cbind(sites, bioclims)
	# sites$solar_rad_kW_hr_per_m2 <- pisr

	### create slightly different data frames for each type of data
	# each will have coordinates and BIOCLIM predictors that match to site
	# need to do this because different number of plants were sampled for different metrics

	bioclims <- as.data.table(bioclims)
	bioclims[ , SITE := sites$site_id]

	biomass <- merge(biomass, bioclims, by = 'SITE', suffixes = '')
	morpho_phys <- merge(morpho_phys, bioclims, by = 'SITE', suffixes = '')

	### define data folds
	#####################

	# We define geo-folds of phenotypically sampled sites and of AG occurrence using the locations of the sampled sites

	# geo-folds of sample sites
	sites_vect <- vect(sites, geom = c('LONGITUDE', 'LATITUDE'), crs = 'wgs84')
	sites_gfolds <- geoFold(sites_vect, k = 4, minIn = 6, method = 'single')
	sites_vect$geofold <- sites_gfolds
	sites$geofold <- sites_gfolds

	# geo-folds of biomass values
	biomass$geofold <- sites_gfolds[match(biomass$SITE, sites$site_id)]

	# geo-folds of occurrences
	occs <- vect('./data_from_adam_and_loretta/andropogon_gerardi_occurrences_with_environment_1961_2020_climatena.gpkg')

	occs$focal_region <- !is.na(occs$n_andropogon_gerardi)
	occs_cents <- centroids(occs)
	sites_vect <- project(sites_vect, occs)
	occs_gfolds <- geoFoldContrast(contrast = occs_cents, pres = sites_vect, presFolds = sites_gfolds)
	occs$geofold <- occs_gfolds

	writeVector(occs, './outputs_loretta/sdm_pdm_for_sem/andropogon_gerardi_occurrences_with_environment_1961_2020_for_integration.gpkg', overwrite = TRUE)

	# map
	sites_vect$geofold <- as.factor(sites_vect$geofold)
	occs$geofold <- as.factor(occs$geofold)
	map <- ggplot() +
		layer_spatial(occs, aes(fill = geofold), ) +
		layer_spatial(sites_vect, aes(fill = geofold), pch = 21, size = 4) +
		ggtitle('Geo-folds')
	ggsave(map, filename = './outputs_loretta/sdm_pdm_for_sem/map_geofolds.png', width = 12, height = 12, dpi = 600)

	dirCreate('./data_from_loretta/sem_00_merged_site_data_with_climate')
	saveRDS(sites, './data_from_loretta/sem_00_merged_site_data_with_climate/sites.rds')
	saveRDS(biomass, './data_from_loretta/sem_00_merged_site_data_with_climate/biomass.rds')
	saveRDS(morpho_phys, './data_from_loretta/sem_00_merged_site_data_with_climate/morpho_phys.rds')

say(date())
say('FINIS!', deco = '+', level = 1)

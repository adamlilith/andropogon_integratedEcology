### MODELING ANDROPOGON GERARDI DISTRIBUTION, MORPHOLOGY, PHYSIOLOGY, and ASSOCIATED MICROBIAL COMMUNITIES
### Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2023-11
###
### This script compiles occurrence data and environmental data for modeling the biogeographic distribution of Andropogon gerardi (Poaceae). We will use county-level occurrence data compiled by Smith et al. (2017 Global Change Biology), coupled with climate data from AdaptWest.
###
### source('C:/Kaji/R/andropogon_integratedEcology/shared_00_compile_occurrence_and_environmental_data.r')
###
### CONTENTS ###
### setup ###
### download BIEN Andropogon gerardi and Poaceae data ###
### port occurrence data from Smith et al. (2017 GCB) and add soil and present-day climate values for each county ###
### extract future climate data for each county ###
### extract 1930s, 1950s, and 2010s climate data from PRISM for each county ###
### extract 1950s PRISM climate and SoilGrids variables to McMillan 1964 AG sites ###

#############
### setup ###
#############

	rm(list = ls())

	drive <- 'C:/Kaji/'

	setwd(paste0(drive, '/Research/Andropogon/Andropogon'))

	library(BIEN) # BIEN plant data
	library(enmSdmX) # SDMing & GIS
	library(ggplot2) # graphics
	library(lubridate) # dates
	library(omnibus) # helper functions
	library(patchwork) # combining ggplots
	library(fasterRaster) # GIS
	library(terra) # GIS

# say('#########################################################')
# say('### download BIEN Andropogon gerardi and Poaceae data ###')
# say('#########################################################')

# 	ag <- BIEN_occurrence_species(species = 'Andropogon gerardii')
# 	saveRDS(ag, './data_from_BIEN/andropogon_gerardi.rds')

# 	poaceae <- BIEN_occurrence_family(family = 'Poaceae')
# 	saveRDS(poaceae, './data_from_BIEN/poaceae.rds')

# say('#####################################################################################################################')
# say('### port occurrence data from Smith et al. (2017 GCB) and add soil and present-day climate values for each county ###')
# say('#####################################################################################################################')

# 	# As SDM predictors, Smith et al. used:
# 	# * Diurnal temperature range (BIO 2)
# 	# * Maximum temperature of the warmest month (BIO 5)
# 	# * Temperature annual range (BIO 7)
# 	# * Total annual precipitation (BIO 12)
# 	# * Potential solar radiation (from SAGA)

# 	### port data from Smith et al. 2017
# 	####################################

# 	say('port data from Smith et al. 2017')

# 	# from Smith et al. 2017
# 	load(paste0(drive, '/Research/Andropogon/Analysis - Phenotype Modeling/Species Records V3/!13d GADM Ver 2 - Multipart - North America - WORLDCLIM Ver 2 Rel June 1 2016 & AG Pheno and Geno Records & Removed Lake Counties.Rdata'))
# 	occs <- vect(gadm)

# 	# select columns
# 	occs <- occs[ , c('NAME_0', 'NAME_1', 'NAME_2', 'area_km2', 'numCrd1', 'numCrd2', 'numCrd3', 'anyAg1to3', 'poaRec', 'agDensity', 'poaDensity')]

# 	# aggregate--some counties were separate-parts polygons and so models can consider them separate observations
# 	occs_df <- as.data.frame(occs)
# 	occs$state_prov_county <- apply(occs_df[ , c('NAME_1', 'NAME_2')], 1, paste, collapse = ' ')
# 	occs <- aggregate(occs, by = 'state_prov_county', fun = 'mean', dissolve = TRUE, count = FALSE)

# 	occs_meta <- data.frame(
# 		country = occs$NAME_0,
# 		state_province = occs$NAME_1,
# 		county = occs$NAME_2,
# 		n_andropogon_gerardi = occs$mean_anyAg1to3,
# 		n_poaceae = occs$mean_poaRec
# 	)

# 	occs_meta$area_km2 <- expanse(occs, unit = 'km')
# 	for (i in ncol(occs):1) occs[ , i] <- NULL # remove columns
# 	for (i in 1:ncol(occs_meta)) {
# 		occs$DUMMY <- occs_meta[ , i]
# 		names(occs)[i] <- names(occs_meta)[i]
# 	}

# 	### clean BIEN data and crop to North America
# 	#############################################

# 	say('clean BIEN data')

# 	# outline of Mexico, US, and Canada from GADM version 4.1 pre-downloaded
# 	nam1 <- vect(paste0(drive, '/Research Data/GADM/Version 4.1/High Res North America Level 1 sans Great Lakes SpatVector WGS84.gpkg'))
# 	nam1 <- nam1[!(nam1$NAME_1 %in% c('Alaska', 'Hawaii'))]

# 	# BIEN records
# 	ag_bien <- readRDS('./data_from_BIEN/andropogon_gerardi.rds')
# 	poa_bien <- readRDS('./data_from_BIEN/poaceae.rds')

# 	ag_bien$date_collected <- as.Date(ag_bien$date_collected)
# 	poa_bien$date_collected <- as.Date(poa_bien$date_collected)

# 	# remove records collected before 2014 (last date of AG records in Smith et al. was 2013) and after 2020 (last year of climate data)
# 	ag_bien <- ag_bien[!is.na(ag_bien$date_collected), ]
# 	poa_bien <- poa_bien[!is.na(poa_bien$date_collected), ]

# 	ag_bien <- ag_bien[year(ag_bien$date_collected) >= 2014 & year(ag_bien$date_collected) <= 2020, ]
# 	poa_bien <- poa_bien[year(poa_bien$date_collected) >= 2014 & year(poa_bien$date_collected) <= 2020, ]

# 	# remove records outside North America
# 	ag_bien <- vect(ag_bien, geom = c('longitude', 'latitude'), crs = getCRS('WGS84'))
# 	in_nam <- extract(nam1, ag_bien)
# 	in_nam <- !is.na(in_nam$COUNTRY)
# 	ag_bien <- ag_bien[in_nam]

# 	poa_bien <- vect(poa_bien, geom = c('longitude', 'latitude'), crs = getCRS('WGS84'))
# 	in_nam <- extract(nam1, poa_bien)
# 	in_nam <- !is.na(in_nam$COUNTRY)
# 	poa_bien <- poa_bien[in_nam]

# 	### add BIEN records
# 	####################

# 	say('add BIEN records')

# 	# We're only adding records to counties in the existing data set.
# 	ag_bien_county <- extract(occs, ag_bien)
# 	poa_bien_county <- extract(occs, poa_bien)

# 	for (i in 1:nrow(ag_bien_county)) {
		
# 		state <- ag_bien_county$state_province[i]
# 		county <- ag_bien_county$county[i]

# 		if (!is.na(county)) {
			
# 			index <- which(occs$state_province == state & occs$county == county)
# 			occs$n_andropogon_gerardi[index] <- occs$n_andropogon_gerardi[index] + 1

# 		}
		
# 	}

# 	for (i in 1:nrow(poa_bien_county)) {
		
# 		state <- poa_bien_county$state_province[i]
# 		county <- poa_bien_county$county[i]
		
# 		if (!is.na(county)) {
			
# 			index <- which(occs$state_province == state & occs$county == county)
# 			occs$n_poaceae[index] <- occs$n_poaceae[index] + 1

# 		}
	
# 	}
	
# 	# MAºzquiz county in Coahuila, Mexico has no AG or Poaceae records, so force it to 0 for both
# 	occs$n_andropogon_gerardi[occs$county == 'MAºzquiz'] <- 0
# 	occs$n_poaceae[occs$county == 'MAºzquiz'] <- 0

# 	# # # ### remove largest parts of Ontario and Manitoba occurrences on basis that they are just too big to indicate species-environment relationships
# 	# # # redacts <- c(
# 	# # # 	which(occs_meta$state_province == 'Ontario' & occs_meta$county == 'Kenora'),
# 	# # # 	which(occs_meta$state_province == 'Ontario' & occs_meta$county == 'Cochrane'),
# 	# # # 	which(occs_meta$state_province == 'Manitoba' & occs_meta$county == 'Division No. 19')
# 	# # # )

# 	# # # for (redact in redacts) {
# 	# # # 	occs$n_andropogon_gerardi[redact] <- 0
# 	# # # }

# 	### extract solar insolation
# 	############################
# 	say('extract solar GDD, BIOCLIMs, insolation, elevation', level = 2)

# 	# insolation
# 	target_dates <- c('1990-03-01', '1990-09-30')
# 	target_dates <- as.Date(target_dates)
# 	target_dates <- seq(target_dates[1], target_dates[2], by = '1 day')

# 	insol <- rast('E:/Ecology/Potential Annual Insolation (SAGA)/Based on ClimateNA 7.03 from SAGA 9.3.0/Annual_Insolation_1990_kW_hr_per_m2.tif')
# 	target_rast_names <- paste0('Annual Insolation.', target_dates)
# 	insol <- insol[[target_rast_names]]
# 	insol <- sum(insol)
# 	names(insol) <- 'insolation_1990_growing_season_kWh_per_m2'

# 	# BIOCLIMs
# 	bc <- rast(paste0(drive, '/Research Data/ClimateNA/v 7.3 AdaptWest/1961-2020/bioclim_variables_1961_2020.tif'))
	
# 	# GDD
# 	gdd5 <- rast(paste0(drive, '/Research Data/ClimateNA/v 7.3 AdaptWest/1961-2020/gdd5.tif'))
# 	names(gdd5) <- 'gdd_5_deg'

# 	# CMI
# 	cmi <- rast(paste0(drive, '/Research Data/ClimateNA/v 7.3 AdaptWest/1961-2020/climaticMoistureIndex.tif'))
# 	names(cmi) <- 'climatic_moisture_index'

# 	# PET
# 	pet <- rast(paste0(drive, '/Research Data/ClimateNA/v 7.3 AdaptWest/1961-2020/petExtremes.tif'))
# 	pet <- pet[[c('PETWarmestQuarter')]]
# 	names(pet) <- 'pet_warmest_quarter_mm'

# 	# elevation
# 	elevation <- rast(paste0(drive, '/Research Data/ClimateNA/v 7.3 AdaptWest/elevation.tif'))
# 	names(elevation) <- 'elevation_m'

# 	env <- c(bc, gdd5, cmi, pet, insol, elevation)

# 	# extract
# 	occs <- project(occs, env)
# 	env_at_occs_by_cell <- extract(env, occs, exact = TRUE)

# 	# aridity
# 	env_at_occs_by_cell$aridity <- (env_at_occs_by_cell$bio1 + 10) / ((env_at_occs_by_cell$bio12 + 1) / 1000)

# 	# calculate weighted average values
# 	# weights are proportion of each cell covered by the polygon
# 	env_at_occs <- data.frame()

# 	vars <- names(env_at_occs_by_cell)
# 	vars <- vars[!(vars %in% c('ID', 'fraction'))]
# 	IDs <- unique(env_at_occs_by_cell$ID)

# 	for (ID in IDs) {

# 		vals <- rep(NA_real_, length(vars))
# 		names(vals) <- vars
		
# 		fraction <- env_at_occs_by_cell$fraction[env_at_occs_by_cell$ID == ID]
# 		fraction_sum <- sum(fraction, na.rm = TRUE)
		
# 		for (var in vars) {
		
# 			var_vals <- env_at_occs_by_cell[env_at_occs_by_cell$ID == ID, var]
# 			val <- sum(var_vals * fraction, na.rm = TRUE) / fraction_sum

# 			vals[[var]] <- val
		
# 		}
		
# 		vals <- round(vals, 2)
# 		vals <- rbind(vals)
# 		env_at_occs <- rbind(env_at_occs, vals, make.row.names = FALSE)

# 	}

# 	occs <- cbind(occs, env_at_occs)

# 	### extract soil variables
# 	##########################
# 	say('extract soil variables', level = 2)

# 	# Rasters are too fine resolution to extract across a county then average without running into memory issues. To fix this, we'll crop the raster to the county extent (plus a buffer), then extract from there.

# 	# Average values across counties, using weighted means.
# 	ph <- rast(paste0(drive, '/Research Data/SoilGrids/SoilGrids 2.0/phh2o_0-5cm_mean_northAmerica.tif'))
# 	cec <- rast(paste0(drive, '/Research Data/SoilGrids/SoilGrids 2.0/cec_0-5cm_mean_northAmerica.tif'))
# 	clay <- rast(paste0(drive, '/Research Data/SoilGrids/SoilGrids 2.0/clay_0-5cm_mean_northAmerica.tif'))
# 	silt <- rast(paste0(drive, '/Research Data/SoilGrids/SoilGrids 2.0/silt_0-5cm_mean_northAmerica.tif'))
# 	sand <- rast(paste0(drive, '/Research Data/SoilGrids/SoilGrids 2.0/sand_0-5cm_mean_northAmerica.tif'))
# 	soc <- rast(paste0(drive, '/Research Data/SoilGrids/SoilGrids 2.0/soc_0-5cm_mean_northAmerica.tif'))

# 	names(ph) <- 'ph'
# 	names(cec) <- 'cec'
# 	names(clay) <- 'clay'
# 	names(silt) <- 'silt'
# 	names(sand) <- 'sand'
# 	names(soc) <- 'soc'

# 	ph <- ph / 10
# 	cec <- cec / 1000
# 	clay <- clay / 1000
# 	silt <- silt / 1000
# 	sand <- sand / 1000
# 	soc <- soc / 100

# 	vars <- c('ph', 'cec', 'clay', 'silt', 'sand', 'soc')
# 	soil <- c(ph, cec, clay, silt, sand, soc)
# 	names(soil) <- vars
# 	occs <- project(occs, soil)
	
# 	soil <- aggregate(soil, 4, mean, na.rm = TRUE) # doing this to speed up extraction (raw cells are 250 m resolution)

# 	occs$soc <- occs$sand <- occs$silt <- occs$clay <- occs$cec <- occs$ph <- NA_real_

# 	for (i in 1:nrow(occs)) {

# 		county <- occs[i]
# 		county <- buffer(county, 500)
# 		county <- ext(county)
# 		county <- as.polygons(county, crs = crs(occs))
		
# 		county_soils <- crop(soil, county)
# 		county_soils <- extract(county_soils, county, exact = TRUE, ID = FALSE)

# 		for (var in vars) {
# 			county_soils[ , var] <- county_soils[ , var] * county_soils$fraction
# 		}
		
# 		county_soils <- colSums(county_soils, na.rm = TRUE)
# 		for (var in vars) {
# 			county_soils[[var]] <- county_soils[[var]] / county_soils[['fraction']]
# 		}
		
# 		occs$soc[i] <- county_soils[['soc']]
# 		occs$sand[i] <- county_soils[['sand']]
# 		occs$silt[i] <- county_soils[['silt']]
# 		occs$clay[i] <- county_soils[['clay']]
# 		occs$cec[i] <- county_soils[['cec']]
# 		occs$ph[i] <- county_soils[['ph']]

# 	}

# 	occs <- project(occs, insol)
# 	writeVector(occs, './data_from_adam_and_loretta/andropogon_gerardi_occurrences_with_environment_1961_2020_climatena.gpkg', overwrite = TRUE)

# 	### meta-data
# 	#############

# 	sink('./data_from_adam_and_loretta/andropogon_gerardi_occurrences_with_environment_README.txt', split = TRUE)

# 		say('ANDROPOGON GERARDI COUNTY-LEVEL OCCURRENCES AND ASSOCIATED ENVIRONMENTAL DATA')
# 		say('Compiled by Adam B. Smith ', date(), ' from Smith et al. (2017 Global Change Bology)', post = 2)

# 		say('These files are SpatVectors (terra package) representing all the counties of Mexico, the US (excluding Hawaii and Alaska), and Canada. The spatial representation is from GADM 3.6. Occurrence data was acquired by Loretta Johnson and Adam Smith and collated by Smith circa 2016. Occurrences were retained if they were collected since 1950. The associated data tables have the following fields:', breaks = 80)
		
# 		say('SPATIAL DATA', pre = 1)
# 		say('country .................. country')
# 		say('state_province ........... state/province name (from GADM 3.6)')
# 		say('county ................... name of county or equivalent (from GADM 3.6)')
# 		say('area_km2 ................. county area in km2')
		
# 		say('OCCURRENCE DATA: Values are usually number of specimens. NA values indicate data were not collected for this unit.', pre = 1)
# 		say('n_andropogon_gerardi ..... number of specimens of Andropogon gerardi')
# 		say('n_poaceae ................ number of specimens of Poaceae')
		
# 		say('CLIMATE AND ELEVATION DATA from or derived from ClimateNA Version 7.3 (https://adaptwest.databasin.org/pages/adaptwest-climatena/)', pre = 1)
# 		say('Environmental data from ClimateNA at the county level has been aggregated using the mean value across all 1-km2 ClimateNA cells overlapping the county, with cell weight given by proportion of the cell overlapped by the county. Values represent averages across 1961-2020, which were generated using average of the ClimateNA 1961-1990 and 1991-2020 normals.', breaks = 80)
# 		say('bio01 through bio19 ...... BIOCLIM variables 1 through 19, units of deg C, mm, or unit-less (for definitions see https://www.worldclim.org/data/bioclim.html)')
# 		say('gdd5 ..................... total annual growing degree days above 5 deg C, in deg C')
# 		say('climatic_moisture_index .. climatic moisture index, unit-less: [-1, 1]')
# 		say('pet_warmest_quarter_mm ... potential evapotranspiration of the warmest 3 month run, in mm')
# 		say('elevation_m .............. elevation, in m')	
# 		say('aridity .................. (mean annual temperature + 10) / ((total annual precipitation + 1) / 1000)')
		
# 		say('SOIL DATA from SoilGrids Version 2.0 (https://www.isric.org/explore/soilgrids/soilgrids-access)', pre = 1)
# 		say('All values are for depth 0 to 5 cm.')
# 		say('ph ....................... pH, measured in water')
# 		say('cec ...................... cation exchange capacity, meq/100 g of soil')
# 		say('sand, silt, clay ......... unit-less (proportion: [0, 1])')
# 		say('soc ...................... soil organic carbon (kg of C / m2 ???, or % dry weight ???)')

# 	sink()

# say('###################################################')
# say('### extract future climate data for each county ###')
# say('###################################################')

# 	### present-day vector
# 	cna <- rast(paste0(drive, '/Research Data/ClimateNA/v 7.3 AdaptWest/1961-2020/climaticMoistureIndex.tif'))

# 	ag_sq <- vect('./data_from_adam_and_loretta/andropogon_gerardi_occurrences_with_environment_1961_2020_climatena.gpkg')
# 	ag_sq <- project(ag_sq, cna)

# 	# This chunk extracts future climate data to the spatial vector used to store AG occurrence data. Future climates are from ClimateNA 7.3 (AdaptWest versions: https://adaptwest.databasin.org/pages/adaptwest-climatena/)

# 	faster(grassDir = 'C:/Program Files/GRASS GIS 8.4/', verbose = TRUE, useDataTable = TRUE)

# 	futs <- c('ensemble_8GCMs_ssp245_2041_2070', 'ensemble_8GCMs_ssp245_2071_2100', 'ensemble_8GCMs_ssp370_2041_2070', 'ensemble_8GCMs_ssp370_2071_2100')

# 	keep_cols <- c('country', 'state_province', 'county', 'elevation_m', 'insolation_1990_growing_season_kWh_per_m2', 'ph', 'cec', 'clay', 'silt', 'sand', 'soc')

# 	for (fut in futs) {
	
# 		say(fut)

# 		ag_fut <- ag_sq

# 		# remove present-day climate
# 		removes <- names(ag_fut)[!(names(ag_fut) %in% keep_cols)]
# 		for (remove in removes) ag_fut[ , remove] <- NULL

# 		# extract future
# 		clim_dir <- paste0(drive, '/Research Data/ClimateNA/v 7.3 AdaptWest/', fut, '_monthly')

# 		ppt <- listFiles(clim_dir, pattern = 'PPT')
# 		tmin <- listFiles(clim_dir, pattern = 'Tmin')
# 		tmax <- listFiles(clim_dir, pattern = 'Tmax')
# 		tmean <- listFiles(clim_dir, pattern = 'Tave')

# 		ppt <- fast(ppt)
# 		tmin <- fast(tmin)
# 		tmax <- fast(tmax)
# 		tmean <- fast(tmean)

# 		bc <- bioclims(ppt = ppt, tmin = tmin, tmax = tmax, tmean = tmean, bios = 1:19, verbose = TRUE)
# 		names <- names(bc)
		
# 		bc <- rast(bc)
# 		names(bc) <- names

# 		env_at_occs_by_cell <- extract(bc, ag_fut, exact = TRUE)

# 		# calculate weighted average values
# 		# weights are proportion of each cell covered by the polygon

# 		vars <- names(env_at_occs_by_cell)
# 		vars <- vars[!(vars %in% c('ID', 'fraction'))]
# 		IDs <- unique(env_at_occs_by_cell$ID)

# 		env_at_occs <- data.frame()
# 		for (ID in IDs) {

# 			vals <- rep(NA_real_, length(vars))
# 			names(vals) <- vars
			
# 			fraction <- env_at_occs_by_cell$fraction[env_at_occs_by_cell$ID == ID]
# 			fraction_sum <- sum(fraction, na.rm = TRUE)
			
# 			for (var in vars) {
			
# 				var_vals <- env_at_occs_by_cell[env_at_occs_by_cell$ID == ID, var]
# 				val <- sum(var_vals * fraction, na.rm = TRUE) / fraction_sum

# 				vals[[var]] <- val
			
# 			}
			
# 			vals <- round(vals, 2)
# 			vals <- rbind(vals)
# 			env_at_occs <- rbind(env_at_occs, vals, make.row.names = FALSE)

# 		}
		
# 		env_at_occs$aridity <- (env_at_occs$bio1 + 10) / ((env_at_occs$bio12 + 1) / 1000)

# 		ag_fut <- cbind(ag_fut, env_at_occs)
# 		writeVector(ag_fut, paste0('./data_from_adam_and_loretta/andropogon_gerardi_occurrences_with_environment_', fut, '.gpkg'), overwrite = TRUE)

# 		mow(ask = FALSE)

# 	}

# say('####################################################################')
# say('### calculate time trajectory of climate from 21 Kybp to present ###')
# say('####################################################################')

# 	# select centroids of counties with AG
# 	occs <- vect('./data_from_adam_and_loretta/andropogon_gerardi_occurrences_with_environment_1961_2020_climatena.gpkg')
# 	occs <- occs[occs$n_andropogon_gerardi > 0]
# 	occs  <- centroids(occs)

# 	# centuries
# 	cents_bce  <- 1:200
# 	cents_bce <- paste0('-', ifelse(nchar(cents_bce) == 2, '0', ''), ifelse(nchar(cents_bce) == 1, '00', ''), cents_bce)
# 	cents_bce <- rev(cents_bce)

# 	cents_ce  <- 0:20
# 	cents_ce <- paste0('+', prefix(cents_ce, 3))

# 	cents <- c(cents_bce, cents_ce)

# 	env <- data.table()
# 	for (cent in cents) {
	
# 		say(cent)

# 		bio1 <- rast(paste0('E:/ecology/CHELSA/North America/', cent, ' Century/bio01.tif'))
# 		bio12 <- rast(paste0('E:/ecology/CHELSA/North America/', cent, ' Century/bio12.tif'))
	
# 		this_env <- c(bio1, bio12)
	
# 		vals <- extract(this_env, occs, ID = FALSE)
# 		vals$bio01 <- vals$bio01 / 10
# 		vals$bio12 <- vals$bio12 / 10
# 		means <- colMeans(vals, na.rm = TRUE)
# 		medians <- apply(vals, 2, median, na.rm = TRUE)
# 		lowers <- apply(vals, 2, quantile, 0.05, na.rm = TRUE)
# 		uppers <- apply(vals, 2, quantile, 0.95, na.rm = TRUE)

# 		env <- rbind(
# 			env,
# 			data.table(
# 				century = as.numeric(cent),
# 				variable = c('bio1', 'bio12'),
# 				lower = lowers,
# 				mean = means,
# 				median = medians,
# 				upper = uppers
# 			)
# 		)

# 	}

# 	fwrite(env, './outputs_loretta/climate_21_to_0_Kybp_at_centroids_of_counties_with_AG.csv', row.names = FALSE)

# 	bio1 <- ggplot(env[variable == 'bio1'], aes(x = century, y = mean)) +
# 		geom_line() +
# 		ggtitle('BIO1') +
# 		scale_x_continuous(breaks = scales::pretty_breaks(n = 20)) +
# 		theme(
# 			panel.grid.major.x = element_line(color = "grey80"),
# 			panel.grid.minor.x = element_line(color = "grey90")
# 		)

# 	bio12 <- ggplot(env[variable == 'bio12'], aes(x = century, y = mean)) +
# 		geom_line() +
# 		ggtitle('BIO12') +
# 		scale_x_continuous(breaks = scales::pretty_breaks(n = 20)) +
# 		theme(
# 			panel.grid.major.x = element_line(color = "grey80"),
# 			panel.grid.minor.x = element_line(color = "grey90")
# 		)

# 	ggsave(bio1, filename = './outputs_loretta/climate_21_to_0_Kybp_at_centroids_of_counties_with_AG_bio01.png', width = 12, height = 8)
# 	ggsave(bio12, filename = './outputs_loretta/climate_21_to_0_Kybp_at_centroids_of_counties_with_AG_bio12.png', width = 12, height = 8)

say('###############################################################################')
say('### extract 1930s, 1950s, and 2010s climate data from PRISM for each county ###')
say('###############################################################################')

	### present-day vector
	occs <- vect('./data_from_adam_and_loretta/andropogon_gerardi_occurrences_with_environment_1961_2020_climatena.gpkg')
	occs <- occs[ , c('country', 'state_province', 'county', 'ph', 'sand', 'silt', 'clay')]
	occs <- occs[occs$country == 'United States']

	# periods <- c('1931_1940', '1952_1961', '2013_2022')
	periods <- c('2013_2022')

	for (period in periods) {

		say(period)

		# BIOCLIMs
		bcs <- rast(paste0('C:/Kaji/Research Data/PRISM/lt81m/bioclims_', period, '.tif'))

		# extract
		this_occs <- project(occs, bcs)
		this_occs <- crop(this_occs, bcs)
		env_at_occs_by_cell <- extract(bcs, this_occs, exact = TRUE)

		# aridity
		env_at_occs_by_cell$aridity <- (env_at_occs_by_cell$bio1 + 10) / ((env_at_occs_by_cell$bio12 + 1) / 1000)

		# calculate weighted average values
		# weights are proportion of each cell covered by the polygon
		env_at_occs <- data.frame()

		vars <- names(env_at_occs_by_cell)
		vars <- vars[!(vars %in% c('ID', 'fraction'))]
		IDs <- unique(env_at_occs_by_cell$ID)

		for (ID in IDs) {

			vals <- rep(NA_real_, length(vars))
			names(vals) <- vars
			
			fraction <- env_at_occs_by_cell$fraction[env_at_occs_by_cell$ID == ID]
			fraction_sum <- sum(fraction, na.rm = TRUE)
			
			for (var in vars) {
			
				var_vals <- env_at_occs_by_cell[env_at_occs_by_cell$ID == ID, var]
				val <- sum(var_vals * fraction, na.rm = TRUE) / fraction_sum

				vals[[var]] <- val
			
			}
			
			vals <- round(vals, 2)
			vals <- rbind(vals)
			env_at_occs <- rbind(env_at_occs, vals, make.row.names = FALSE)

		}

		this_occs <- cbind(this_occs, env_at_occs)

		writeVector(this_occs, paste0('./data_from_adam_and_loretta/andropogon_gerardi_occurrences_with_environment_', period, '_prism.gpkg')	, overwrite = TRUE)

	} # next period

# say('#####################################################################################')
# say('### extract 1950s PRISM climate and SoilGrids variables to McMillan 1964 AG sites ###')
# say('#####################################################################################')

# 	mcm <- vect('./data_from_mcmillan/mcmillan_1964_fig_2_andropogon_gerardi.shp')

# 	### 1950s climate
# 	fifties <- rast('C:/Kaji/Research Data/PRISM/lt81m/bioclims_1952_1961.tif')
# 	clim <- extract(fifties, mcm, ID = FALSE)
# 	clim$aridity <- (clim$bio1 + 10) / ((clim$bio12 + 1) / 1000)

# 	### soil
# 	ph <- rast('C:/Kaji/Research Data/SoilGrids/SoilGrids 2.0/phh2o_0-5cm_mean_northAmerica.tif')
# 	sand <- rast('C:/Kaji/Research Data/SoilGrids/SoilGrids 2.0/sand_0-5cm_mean_northAmerica.tif')
# 	silt <- rast('C:/Kaji/Research Data/SoilGrids/SoilGrids 2.0/silt_0-5cm_mean_northAmerica.tif')
# 	clay <- rast('C:/Kaji/Research Data/SoilGrids/SoilGrids 2.0/clay_0-5cm_mean_northAmerica.tif')

# 	soil <- c(ph, sand, silt, clay)
# 	names(soil) <- c('ph', 'sand', 'silt', 'clay')

# 	# fill NA cells
# 	soil <- focal(soil, w = 3, fun = 'mean', na.policy = 'only', na.rm = TRUE)
# 	soil_env <- extract(soil, mcm, ID = FALSE)

# 	mcm <- cbind(mcm, clim, soil_env)
# 	writeVector(mcm, './data_from_mcmillan/mcmillan_1964_fig_2_ag_prism_1952_1961_soilgrids_top_5_cm.gpkg')

say('DONE!', level = 1, deco = '!')

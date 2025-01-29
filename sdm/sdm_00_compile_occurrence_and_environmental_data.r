### MODELING ANDROPOGON GERARDI DISTRIBUTION, MORPHOLOGY, PHYSIOLOGY, and ASSOCIATED MICROBIAL COMMUNITIES
### Erica Newman | Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2023-11
###
### This script compiles occurrence data and environmental data for modeling the biogeographic distribution of Andropogon gerardi (Poaceae). We will use county-level occurrence data compiled by Smith et al. (2017 Global Change Biology), coupled with climate data from AdaptWest.
###
### source('C:/Ecology/R/andropogon_integratedEcology/sdm/sdm_00_compile_occurrence_and_environmental_data.r')
### source('C:/Subarashi/R/andropogon_integratedEcology/sdm/sdm_00_compile_occurrence_and_environmental_data.r')
###
### CONTENTS ###
### setup ###
### download BIEN Andropogon gerardi and Poaceae data ###
### port occurrence data from Smith et al. (2017 GCB) and add soil and present-day climate values for each county ###
### extract future climate data for each county ###

#############
### setup ###
#############

	rm(list = ls())

	drive <- 'C:/Ecology/'
	# drive <- 'C:/Subarashi/'

	setwd(paste0(drive, '/Research/Andropogon/Andropogon'))

	library(BIEN) # BIEN plant data
	library(enmSdmX) # SDMing & GIS
	library(lubridate) # dates
	library(omnibus) # helper functions
	library(terra) # GIS
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

	# # As SDM predictors, Smith et al. used:
	# # * Diurnal temperature range (BIO 2)
	# # * Maximum temperature of the warmest month (BIO 5)
	# # * Temperature annual range (BIO 7)
	# # * Total annual precipitation (BIO 12)
	# # * Potential solar radiation (from SAGA)

	# ### clean BIEN data and crop to North America
	# #############################################

	# say('clean BIEN data')

	# # outline of Mexico, US, and Canada from GADM version 4.1 pre-downloaded
	# nam1 <- vect(paste0(drive, '/Research Data/GADM/Version 4.1/High Res North America Level 1 sans Great Lakes SpatVector WGS84.gpkg'))
	# nam1 <- nam1[!(nam1$NAME_1 %in% c('Alaska', 'Hawaii'))]

	# # BIEN records
	# ag_bien <- readRDS('./data_from_BIEN/andropogon_gerardi.rds')
	# poa_bien <- readRDS('./data_from_BIEN/poaceae.rds')

	# ag_bien$date_collected <- as.Date(ag_bien$date_collected)
	# poa_bien$date_collected <- as.Date(poa_bien$date_collected)

	# # remove records collected before 2014 (last date of AG records in Smith et al. was 2013) and 2020 (last year of climate data)
	# ag_bien <- ag_bien[!is.na(ag_bien$date_collected), ]
	# poa_bien <- poa_bien[!is.na(poa_bien$date_collected), ]

	# ag_bien <- ag_bien[year(ag_bien$date_collected) >= 2014 & year(ag_bien$date_collected) <= 2020, ]
	# poa_bien <- poa_bien[year(poa_bien$date_collected) >= 2014 & year(poa_bien$date_collected) <= 2020, ]

	# # remove records outside North America
	# ag_bien <- vect(ag_bien, geom = c('longitude', 'latitude'), crs = getCRS('WGS84'))
	# in_nam <- extract(nam1, ag_bien)
	# in_nam <- !is.na(in_nam$COUNTRY)
	# ag_bien <- ag_bien[in_nam]

	# poa_bien <- vect(poa_bien, geom = c('longitude', 'latitude'), crs = getCRS('WGS84'))
	# in_nam <- extract(nam1, poa_bien)
	# in_nam <- !is.na(in_nam$COUNTRY)
	# poa_bien <- poa_bien[in_nam]

	# ### port data from Smith et al. 2017
	# ####################################

	# say('port data from Smith et al. 2017')

	# # from Smith et al. 2017
	# load(paste0(drive, '/Research/Andropogon/Analysis - Phenotype Modeling/Species Records V3/!13d GADM Ver 2 - Multipart - North America - WORLDCLIM Ver 2 Rel June 1 2016 & AG Pheno and Geno Records & Removed Lake Counties.Rdata'))
	# occs <- vect(gadm)

	# # select columns
	# occs <- occs[ , c('NAME_0', 'NAME_1', 'NAME_2', 'area_km2', 'numCrd1', 'numCrd2', 'numCrd3', 'anyAg1to3', 'poaRec', 'agDensity', 'poaDensity')]

	# # rename columns
	# names(occs)[names(occs) == 'NAME_0'] <- 'country'
	# names(occs)[names(occs) == 'NAME_1'] <- 'state_province'
	# names(occs)[names(occs) == 'NAME_2'] <- 'county'
	# names(occs)[names(occs) == 'numCrd1'] <- 'num_ag_quality_1'
	# names(occs)[names(occs) == 'numCrd2'] <- 'num_ag_quality_2'
	# names(occs)[names(occs) == 'numCrd3'] <- 'num_ag_quality_3'
	# names(occs)[names(occs) == 'anyAg1to3'] <- 'any_ag_quality_1_to_3'
	# names(occs)[names(occs) == 'poaRec'] <- 'num_poaceae_records'
	# names(occs)[names(occs) == 'agDensity'] <- 'ag_density'
	# names(occs)[names(occs) == 'poaDensity'] <- 'poa_density'

	# ### add BIEN records
	# ####################

	# say('add BIEN records')

	# # We're only adding records to counties in the existing data set.

	# # outline of Mexico, US, and Canada from GADM version 4.1 pre-downloaded
	# nam2 <- vect(paste0(drive, '/Research Data/GADM/Version 4.1/High Res North America Level 2 sans Great Lakes SpatVector WGS84.gpkg'))
	# nam2 <- nam2[!(nam2$NAME_1 %in% c('Alaska', 'Hawaii'))]

	# # get state/county of each record
	# ag_bien_county <- extract(nam2, ag_bien)
	# poa_bien_county <- extract(nam2, poa_bien)

	# for (i in 1:nrow(ag_bien_county)) {
		
		# state <- ag_bien_county$NAME_1[i]
		# county <- ag_bien_county$NAME_2[i]

		# if (!is.na(county)) {
			
			# index <- which(occs$state_province == state & occs$county == county)
			# occs$any_ag_quality_1_to_3[index] <- occs$any_ag_quality_1_to_3[index] + 1

		# }
		
	# }

	# for (i in 1:nrow(poa_bien_county)) {
		
		# state <- poa_bien_county$NAME_1[i]
		# county <- poa_bien_county$NAME_2[i]
		
		# if (!is.na(county)) {
			
			# index <- which(occs$state_province == state & occs$county == county)
			# occs$num_poaceae_records[index] <- occs$num_poaceae_records[index] + 1

		# }
	
	# }
	
	# occs$ag_density <- occs$any_ag_quality_1_to_3 / occs$area_km2
	# occs$poa_density <- occs$num_poaceae_records / occs$area_km2

	# ### extract solar insolation
	# ############################
	# say('extract solar GDD, BIOCLIMs, insolation, elevation', level = 2)

	# # insolation
	# target_dates <- c('1990-03-01', '1990-09-30')
	# target_dates <- as.Date(target_dates)
	# target_dates <- seq(target_dates[1], target_dates[2], by = '1 day')

	# insol <- rast('E:/Ecology/Potential Annual Insolation (SAGA)/Based on ClimateNA 7.03 from SAGA 9.3.0/Annual_Insolation_1990_kW_hr_per_m2.tif')
	# target_rast_names <- paste0('Annual Insolation.', target_dates)
	# insol <- insol[[target_rast_names]]
	# insol <- sum(insol)
	# names(insol) <- 'insolation_1990_growing_season_kWh_per_m2'

	# # BIOCLIMs
	# bc <- rast(paste0(drive, '/Research Data/ClimateNA/v 7.3 AdaptWest/1961-2020/bioclim_variables_1961_2020.tif'))
	
	# # GDD
	# gdd5 <- rast(paste0(drive, '/Research Data/ClimateNA/v 7.3 AdaptWest/1961-2020/gdd5.tif'))
	# names(gdd5) <- 'gdd_5_deg'

	# # CMI
	# cmi <- rast(paste0(drive, '/Research Data/ClimateNA/v 7.3 AdaptWest/1961-2020/climaticMoistureIndex.tif'))
	# names(cmi) <- 'climatic_moisture_index'

	# # PET
	# pet <- rast(paste0(drive, '/Research Data/ClimateNA/v 7.3 AdaptWest/1961-2020/petExtremes.tif'))
	# pet <- pet[[c('PETWarmestQuarter')]]
	# names(pet) <- 'pet_warmest_quarter_mm'

	# # elevation
	# elevation <- rast(paste0(drive, '/Research Data/ClimateNA/v 7.3 AdaptWest/elevation.tif'))
	# names(elevation) <- 'elevation_m'

	# env <- c(bc, gdd5, cmi, pet, insol, elevation)

	# # extract
	# occs <- project(occs, env)
	# env_at_occs_by_cell <- extract(env, occs, exact = TRUE)

	# # aridity
	# env_at_occs_by_cell$aridity <- (env_at_occs_by_cell$bio1 + 10) / (env_at_occs_by_cell$bio12 / 1000)

	# # calculate weighted average values
	# # weights are proportion of each cell covered by the polygon
	# env_at_occs <- data.frame()

	# vars <- names(env_at_occs_by_cell)
	# vars <- vars[!(vars %in% c('ID', 'fraction'))]
	# IDs <- unique(env_at_occs_by_cell$ID)

	# for (ID in IDs) {

		# vals <- rep(NA_real_, length(vars))
		# names(vals) <- vars
		
		# fraction <- env_at_occs_by_cell$fraction[env_at_occs_by_cell$ID == ID]
		# fraction_sum <- sum(fraction, na.rm = TRUE)
		
		# for (var in vars) {
		
			# var_vals <- env_at_occs_by_cell[env_at_occs_by_cell$ID == ID, var]
			# val <- sum(var_vals * fraction, na.rm = TRUE) / fraction_sum

			# vals[[var]] <- val
		
		# }
		
		# vals <- round(vals, 2)
		# vals <- rbind(vals)
		# env_at_occs <- rbind(env_at_occs, vals, make.row.names = FALSE)

	# }

	# occs <- cbind(occs, env_at_occs)

	# ### extract soil variables
	# ##########################
	# say('extract soil variables', level = 2)

	# # Rasters are too fine resolution to extract across a county then average without running into memory issues. To fix this, we'll crop the raster to the county extent (plus a buffer), then extract from there.

	# # Average values across counties, using weighted means.

	# ph <- rast(paste0(drive, '/Research Data/SoilGrids/SoilGrids 2.0/phh2o_0-5cm_mean_northAmerica.tif'))
	# cec <- rast(paste0(drive, '/Research Data/SoilGrids/SoilGrids 2.0/cec_0-5cm_mean_northAmerica.tif'))
	# clay <- rast(paste0(drive, '/Research Data/SoilGrids/SoilGrids 2.0/clay_0-5cm_mean_northAmerica.tif'))
	# silt <- rast(paste0(drive, '/Research Data/SoilGrids/SoilGrids 2.0/silt_0-5cm_mean_northAmerica.tif'))
	# sand <- rast(paste0(drive, '/Research Data/SoilGrids/SoilGrids 2.0/sand_0-5cm_mean_northAmerica.tif'))
	# soc <- rast(paste0(drive, '/Research Data/SoilGrids/SoilGrids 2.0/soc_0-5cm_mean_northAmerica.tif'))

	# names(ph) <- 'ph'
	# names(cec) <- 'cec'
	# names(clay) <- 'clay'
	# names(silt) <- 'silt'
	# names(sand) <- 'sand'
	# names(soc) <- 'soc'

	# ph <- ph / 10
	# cec <- cec / 1000
	# clay <- clay / 1000
	# silt <- silt / 1000
	# sand <- sand / 1000
	# soc <- soc / 100

	# vars <- c('ph', 'cec', 'clay', 'silt', 'sand', 'soc')
	# soil <- c(ph, cec, clay, silt, sand, soc)
	# names(soil) <- vars
	# occs <- project(occs, soil)
	
	# soil <- aggregate(soil, 4, mean, na.rm = TRUE)

	# occs$soc <- occs$sand <- occs$silt <- occs$clay <- occs$cec <- occs$ph <- NA_real_

	# for (i in 1:nrow(occs)) {

		# county <- occs[i]
		# county <- buffer(county, 500)
		# county <- ext(county)
		# county <- as.polygons(county, crs = crs(occs))
		
		# county_soils <- crop(soil, county)
		# county_soils <- extract(county_soils, county, exact = TRUE, ID = FALSE)

		# for (var in vars) {
			# county_soils[ , var] <- county_soils[ , var] * county_soils$fraction
		# }
		
		# county_soils <- colSums(county_soils, na.rm = TRUE)
		# for (var in vars) {
			# county_soils[[var]] <- county_soils[[var]] / county_soils[['fraction']]
		# }
		
		# occs$soc[i] <- county_soils[['soc']]
		# occs$sand[i] <- county_soils[['sand']]
		# occs$silt[i] <- county_soils[['silt']]
		# occs$clay[i] <- county_soils[['clay']]
		# occs$cec[i] <- county_soils[['cec']]
		# occs$ph[i] <- county_soils[['ph']]

	# }

	# occs <- project(occs, insol)
	# writeVector(occs, './data_from_adam_and_loretta/andropogon_gerardi_occurrences_with_environment_1961_2020.gpkg', overwrite = TRUE)

	# ### meta-data
	# #############

	# sink('./data_from_adam_and_loretta/andropogon_gerardi_occurrences_with_environment_README.txt', split = TRUE)

		# say('ANDROPOGON GERARDI COUNTY-LEVEL OCCURRENCES AND ASSOCIATED ENVIRONMENTAL DATA')
		# say('Compiled by Adam B. Smith ', date(), ' from Smith et al. (2017 Global Change Bology)', post = 2)

		# say('These files are SpatVectors (terra package) representing all the counties of Mexico, the US (excluding Hawaii and Alaska), and Canada. The spatial representation is from GADM 3.6. Occurrence data was acquired by Loretta Johnson and Adam Smith and collated by Smith circa 2016. Occurrences were retained if they were collected since 1950. The associated data tables have the following fields:', breaks = 80)
		
		# say('SPATIAL DATA', pre = 1)
		# say('country .................. country')
		# say('state_province ........... state/province name (from GADM 3.6)')
		# say('county ................... name of county or equivalent (from GADM 3.6)')
		# say('area_km2 ................. county area in km2')
		
		# say('OCCURRENCE DATA: Values are usually number of specimens. NA values indicate data were not collected for this unit.', pre = 1)
		# say('num_ag_quality_1 ......... number of specimens of Andropogon gerardi of high spatial quality (ie, coordinates with small uncertainty)')
		# say('num_ag_quality_2 ......... number of specimens of Andropogon gerardi of moderate spatial quality (ie, coordinates only)')
		# say('num_ag_quality_3 ......... number of specimens of Andropogon gerardi of high spatial quality (ie, political unit only)')
		# say('any_ag_quality_1_to_3 .... number of specimens of Andropogon gerardi of any quality... NA means not scored for AG')
		# say('num_poaceae_records ...... number of specimens of Poaceae')
		# say('ag_density ............... Andropogon gerardi specimens / km2')
		# say('poa_density .............. Poaceae specimens / km2')
		
		# say('CLIMATE AND ELEVATION DATA from or derived from ClimateNA Version 7.3 (https://adaptwest.databasin.org/pages/adaptwest-climatena/)', pre = 1)
		# say('Environmental data from ClimateNA at the county level has been aggregated using the mean value across all 1-km2 ClimateNA cells overlapping the county, with cell weight given by proportion of the cell overlapped by the county. Values represent averages across 1961-2020, which were generated using average of the ClimateNA 1961-1990 and 1991-2020 normals.', breaks = 80)
		# say('bio01 through bio19 ...... BIOCLIM variables 1 through 19, units of deg C, mm, or unit-less (for definitions see https://www.worldclim.org/data/bioclim.html)')
		# say('gdd5 ..................... total annual growing degree days above 5 deg C, in deg C')
		# say('climatic_moisture_index .. climatic moisture index, unit-less: [-1, 1]')
		# say('pet_warmest_quarter_mm ... potential evapotranspiration of the warmest 3 month run, in mm')
		# say('elevation_m .............. elevation, in m')	
		# say('aridity .................. (mean annual temperature + 10) / (total annual precipitation / 1000)')
		
		# say('SOIL DATA from SoilGrids Version 2.0 (https://www.isric.org/explore/soilgrids/soilgrids-access)', pre = 1)
		# say('All values are for depth 0 to 5 cm.')
		# say('ph ....................... pH, measured in water')
		# say('cec ...................... cation exchange capacity, meq/100 g of soil')
		# say('sand, silt, clay ......... unit-less (proportion: [0, 1])')
		# say('soc ...................... soil organic carbon (kg of C / m2 ???, or % dry weight ???)')

	# sink()

say('###################################################')
say('### extract future climate data for each county ###')
say('###################################################')

	### user-defined
	################
		
	load(paste0(drive, '/Research/Andropogon/Analysis - Phenotype Modeling/Species Records V3/!13d GADM Ver 2 - Multipart - North America - WORLDCLIM Ver 2 Rel June 1 2016 & AG Pheno and Geno Records & Removed Lake Counties.Rdata'))
	occs <- vect(gadm)

	# select columns
	occs <- occs[ , c('NAME_0', 'NAME_1', 'NAME_2', 'area_km2')]

	# rename columns
	names(occs)[names(occs) == 'NAME_0'] <- 'country'
	names(occs)[names(occs) == 'NAME_1'] <- 'stateProvince'
	names(occs)[names(occs) == 'NAME_2'] <- 'county'

	sq <- vect('./data_from_adam_and_loretta/andropogon_gerardi_occurrences_with_environment_1961_2020.gpkg')
	sq <- as.data.frame(sq)
	sq <- sq[ , c('insolation_1990_growing_season_kWh_per_m2', 'ph', 'cec', 'clay', 'silt', 'sand', 'soc')]
	
	occs <- cbind(occs, sq)

	cna <- rast(paste0(drive, '/Research Data/ClimateNA/v 7.3 AdaptWest/1961-2020/climaticMoistureIndex.tif'))

	ag_nam <- occs
	ag_nam <- project(ag_nam, cna)

	# This chunk extracts future climate data to the spatial vector used to store AG occurrence data. Future climates are from ClimateNA 7.3 (AdaptWest versions: https://adaptwest.databasin.org/pages/adaptwest-climatena/)

	faster(grassDir = 'C:/Program Files/GRASS GIS 8.4/', verbose = TRUE, useDataTable = TRUE)

	futs <- c('ensemble_8GCMs_ssp245_2041_2070', 'ensemble_8GCMs_ssp245_2071_2100', 'ensemble_8GCMs_ssp370_2041_2070', 'ensemble_8GCMs_ssp370_2071_2100')

	for (fut in futs) {
	
		say(fut)

		ag_fut <- ag_nam
		ag_fut[ , paste0('bio', 1:19)] <- NULL
		ag_fut$pet_warmest_quarter_mm <- NULL
		ag_fut$climatic_moisture_index <- NULL
		ag_fut$aridity <- NULL
		ag_fut$gdd_5_deg <- NULL
		ag_fut$ag_lambda <- NULL

		clim_dir <- paste0(drive, '/Research Data/ClimateNA/v 7.3 AdaptWest/', fut, '_monthly')

		ppt <- listFiles(clim_dir, pattern = 'PPT')
		tmin <- listFiles(clim_dir, pattern = 'Tmin')
		tmax <- listFiles(clim_dir, pattern = 'Tmax')
		tmean <- listFiles(clim_dir, pattern = 'Tave')

		ppt <- fast(ppt)
		tmin <- fast(tmin)
		tmax <- fast(tmax)
		tmean <- fast(tmean)

		bc <- bioclims(ppt = ppt, tmin = tmin, tmax = tmax, tmean = tmean, bios = c(1, 7, 12, 15), verbose = TRUE)
		names <- names(bc)
		
		bc <- rast(bc)
		names(bc) <- names

		env_at_occs_by_cell <- extract(bc, ag_fut, exact = TRUE)

		# calculate weighted average values
		# weights are proportion of each cell covered by the polygon

		vars <- names(env_at_occs_by_cell)
		vars <- vars[!(vars %in% c('ID', 'fraction'))]
		IDs <- unique(env_at_occs_by_cell$ID)

		env_at_occs <- data.frame()
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
		
		env_at_occs$aridity <- (env_at_occs$bio1 + 10) / (env_at_occs$bio12 / 1000)

		ag_fut <- cbind(ag_fut, env_at_occs)
		writeVector(ag_fut, paste0('./data_from_adam_and_loretta/andropogon_gerardi_occurrences_with_environment_', fut, '.gpkg'), overwrite = TRUE)

		mow(ask = FALSE)

	}

say('DONE!', level = 1, deco = '!')

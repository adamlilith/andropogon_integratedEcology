### MODELING ANDROPOGON GERARDI DISTRIBUTION, MORPHOLOGY, PHYSIOLOGY, and ASSOCIATED MICROBIAL COMMUNITIES
### Erica Newman | Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2023-11
###
### This script compiles occurrence data and environmental data for modeling the biogeographic distribution of Andropogon gerardi (Poaceae). We will use county-level occurrence data compiled by Smith et al. (2017 Global Change Biology), coupled with climate data from AdaptWest.
###
### source('C:/Adam/R/andropogon_integratedEcology/sdm_00_compile_occurrence_and_environmental_data.r')
### source('E:/Adam/R/andropogon_integratedEcology/sdm_00_compile_occurrence_and_environmental_data.r')
###
### CONTENTS ###
### setup ###
### port occurrence data from Smith et al. (2017 GCB) ###

#############
### setup ###
#############

	rm(list = ls())

	# drive <- 'C:/Adam/'
	drive <- 'E:/Adam/'

	setwd(paste0(drive, '/Research/Andropogon/Andropogon'))

	library(omnibus)
	library(terra)
	# library(fasterRaster) # from www.github.com/adamlilith/fasterRaster NB... using "intuitive" branch as of 2023-11

say('#########################################################')
say('### port occurrence data from Smith et al. (2017 GCB) ###')
say('#########################################################')

	# As SDM predictors, Smith et al. used:
	# * Diurnal temperature range (BIO 2)
	# * Maximum temperature of the warmest month (BIO 5)
	# * Temperature annual range (BIO 7)
	# * Total annual precipitation (BIO 12)
	# * Potential solar radiation (from SAGA)

	### port data from Smith et al. 2017
	####################################

	# from Smith et al. 2017
	load(paste0(drive, '/Research/Andropogon/Analysis - Phenotype Modeling/Species Records V3/!13d GADM Ver 2 - Multipart - North America - WORLDCLIM Ver 2 Rel June 1 2016 & AG Pheno and Geno Records & Removed Lake Counties.Rdata'))
	occs <- vect(gadm)

	# select columns
	occs <- occs[ , c('NAME_0', 'NAME_1', 'NAME_2', 'area_km2', 'numCrd1', 'numCrd2', 'numCrd3', 'anyAg1to3', 'poaRec', 'agDensity', 'poaDensity')]

	# rename columns
	names(occs)[names(occs) == 'NAME_0'] <- 'country'
	names(occs)[names(occs) == 'NAME_1'] <- 'stateProvince'
	names(occs)[names(occs) == 'NAME_2'] <- 'county'
	names(occs)[names(occs) == 'numCrd1'] <- 'num_ag_quality1'
	names(occs)[names(occs) == 'numCrd2'] <- 'num_ag_quality2'
	names(occs)[names(occs) == 'numCrd3'] <- 'num_ag_quality3'
	names(occs)[names(occs) == 'anyAg1to3'] <- 'any_ag_quality1to3'
	names(occs)[names(occs) == 'poaRec'] <- 'num_poaceae_records'
	names(occs)[names(occs) == 'agDensity'] <- 'ag_density'
	names(occs)[names(occs) == 'poaDensity'] <- 'poa_density'

	### extract solar insolation
	############################
	say('extract solar GDD, BIOCLIMs, insolation, elevation', level = 2)

	# insolation
	target_dates <- c('1990-03-01', '1990-09-30')
	target_dates <- as.Date(target_dates)
	target_dates <- seq(target_dates[1], target_dates[2], by = '1 day')

	insol <- rast('F:/Ecology/Potential Annual Insolation (SAGA)/Based on ClimateNA 7.03 from SAGA 9.3.0/Annual_Insolation_1990_kW_hr_per_m2.tif')
	target_rast_names <- paste0('Annual Insolation.', target_dates)
	insol <- insol[[target_rast_names]]
	insol <- sum(insol)
	names(insol) <- 'insolation_1990_growing_season_kWh_per_m2'

	# BIOCLIMs
	bc <- rast(paste0(drive, '/Research Data/ClimateNA/ClimateNA v7.3/1961-2020/bioclim_variables_1961_2020.tif'))

	# GDD
	gdd5 <- rast(paste0(drive, '/Research Data/ClimateNA/ClimateNA v7.3/1961-2020/gdd5.tif'))
	names(gdd5) <- 'gdd_5_deg'

	# CMI
	cmi <- rast(paste0(drive, '/Research Data/ClimateNA/ClimateNA v7.3/1961-2020/climaticMoistureIndex.tif'))
	names(cmi) <- 'climatic_moisture_index'

	# PET
	pet <- rast(paste0(drive, '/Research Data/ClimateNA/ClimateNA v7.3/1961-2020/petExtremes.tif'))
	pet <- pet[[c('PETWarmestQuarter')]]
	names(pet) <- 'pet_warmest_quarter_mm'

	# elevation
	elevation <- rast(paste0(drive, '/Research Data/ClimateNA/ClimateNA v7.3/elevation.tif'))
	names(elevation) <- 'elevation_m'

	env <- c(bc, gdd5, cmi, pet, insol, elevation)

	occs <- project(occs, env)
	env_at_occs_by_cell <- extract(env, occs, weights = TRUE)

	# calculate weighted average values
	# weights are proportion of each cell covered by the polygon
	env_at_occs <- data.frame()

	vars <- names(env_at_occs_by_cell)
	vars <- vars[vars %notin% c('ID', 'weight')]
	IDs <- unique(env_at_occs_by_cell$ID)

	for (ID in IDs) {

		vals <- rep(NA_real_, length(vars))
		names(vals) <- vars
		
		weight <- env_at_occs_by_cell$weight[env_at_occs_by_cell$ID == ID]
		weight_sum <- sum(weight, na.rm = TRUE)
		
		for (var in vars) {
		
			var_vals <- env_at_occs_by_cell[env_at_occs_by_cell$ID == ID, var]
			val <- sum(var_vals * weight, na.rm = TRUE) / weight_sum

			vals[[var]] <- val
		
		}
		
		vals <- round(vals, 2)
		vals <- rbind(vals)
		env_at_occs <- rbind(env_at_occs, vals, make.row.names = FALSE)

	}


	occs <- cbind(occs, env_at_occs)

	### calculate aridity
	#####################
	say('calculate aridity')

	occs$aridity <- occs$bio1 / occs$bio12

	### extract soil variables
	##########################
	say('extract soil variables', level = 2)

	# Rasters are too fine resolution to extract across a county then average without running into memory issues. To fix this, we'll crop the raster to the county extent (plus a buffer), then extract from there.

	# Average values across counties, using weighted means.

	ph <- rast(paste0(drive, '/Research Data/SoilGrids/SoilGrids 2.0/phh2o_0-5cm_mean_northAmerica.tif'))
	cec <- rast(paste0(drive, '/Research Data/SoilGrids/SoilGrids 2.0/cec_0-5cm_mean_northAmerica.tif'))
	clay <- rast(paste0(drive, '/Research Data/SoilGrids/SoilGrids 2.0/clay_0-5cm_mean_northAmerica.tif'))
	silt <- rast(paste0(drive, '/Research Data/SoilGrids/SoilGrids 2.0/silt_0-5cm_mean_northAmerica.tif'))
	sand <- rast(paste0(drive, '/Research Data/SoilGrids/SoilGrids 2.0/sand_0-5cm_mean_northAmerica.tif'))
	soc <- rast(paste0(drive, '/Research Data/SoilGrids/SoilGrids 2.0/soc_0-5cm_mean_northAmerica.tif'))

	names(ph) <- 'pH'
	names(cec) <- 'CEC'
	names(clay) <- 'clay'
	names(silt) <- 'silt'
	names(sand) <- 'sand'
	names(soc) <- 'SOC'

	ph <- ph / 10
	cec <- cec / 1000
	clay <- clay / 1000
	silt <- silt / 1000
	sand <- sand / 1000
	soc <- soc / 100

	vars <- c('pH', 'CEC', 'clay', 'silt', 'sand', 'SOC')
	soil <- c(ph, cec, clay, silt, sand, soc)
	names(soil) <- vars
	occs <- project(occs, soil)

	occs$soc <- occs$sand <- occs$silt <- occs$clay <- occs$cec <- occs$ph <- NA_real_

	for (i in 1:nrow(occs)) {

		say(i, ' ', nrow(occs))

		county <- occs[i]
		county <- buffer(county, 500)
		county <- ext(county)
		county <- as.polygons(county, crs = crs(occs))
		
		county_soils <- crop(soil, county)
		county_soils <- extract(county_soils, county, weights = TRUE, ID = FALSE)

		for (var in vars) {
		
			county_soils[ , var] <- county_soils[ , var] * county_soils$weight
		
		}
		
		county_soils <- colSums(county_soils, na.rm = TRUE)
		for (var in vars) {
		
			county_soils[[var]] <- county_soils[[var]] / county_soils[['weight']]
		
		}
		
		occs$soc[i] <- county_soils[['SOC']]
		occs$sand[i] <- county_soils[['sand']]
		occs$silt[i] <- county_soils[['silt']]
		occs$clay[i] <- county_soils[['clay']]
		occs$cec[i] <- county_soils[['CEC']]
		occs$ph[i] <- county_soils[['pH']]

	}

	occs <- project(occs, insol)
	writeVector(occs, './data/occurrence_data/andropogon_gerardi_occurrences_with_environment2.gpkg', overwrite = TRUE)

	### meta-data
	#############

	sink('./data/occurrence_data/andropogon_gerardi_occurrences_with_environment_README.txt', split = TRUE)

		say('ANDROPOGON GERARDI COUNTY-LEVEL OCCURRENCES AND ASSOCIATED ENVIRONMENTAL DATA')
		say('Compiled by Adam B. Smith ', date(), ' from Smith et al. (2017 Global Change Bology)', post = 2)

		say('This file is a SpatVector (terra package) representing all the counties of Mexico, the US (excluding Hawaii and Alaska), and Canada. The spatial representation is from GADM 3.6. Occurrence data was acquired by Loretta Johnson and collated by Adam Smith circa 2016, and were retained if they were collected since 1950. The associated data table has the following fields:', breaks = 80)
		
		say('SPATIAL DATA', pre = 1)
		say('country .................. country')
		say('stateProvince ............ state/province name (from GADM 3.6)')
		say('county ................... name of county or equivalent (from GADM 3.6)')
		say('area_km2 ................. county area in km2')
		
		say('OCCURRENCE DATA: Values are usually number of specimens. NA values indicate data were not collected for this unit.', pre = 1)
		say('num_ag_quality1 .......... number of specimens of Andropogon gerardi of high spatial quality (ie, coordinates with small uncertainty)')
		say('num_ag_quality2 .......... number of specimens of Andropogon gerardi of moderate spatial quality (ie, coordinates only)')
		say('num_ag_quality3 .......... number of specimens of Andropogon gerardi of high spatial quality (ie, political unit only)')
		say('any_ag_quality1to3 ....... number of specimens of Andropogon gerardi of any quality... NA means not scored for AG')
		say('num_poaceae_records ...... number of specimens of Poaceae')
		say('ag_density ............... Andropogon gerardi specimens / km2')
		say('poa_density .............. Poaceae specimens / km2')
		
		say('CLIMATE AND ELEVATION DATA from or derived from ClimateNA Version 7.3 (https://adaptwest.databasin.org/pages/adaptwest-climatena/)', pre = 1)
		say('Environmental data from ClimateNA at the county level has been aggregated using the mean value across all 1-km2 ClimateNA cells overlapping the county, with cell weight given by proportion of the cell overlapped by the county. Values represent averages across 1961-2020, which were generated using average of the ClimateNA 1961-1990 and 1991-2020 normals.', breaks = 80)
		say('bio01 through bio19 ...... BIOCLIM variables 1 through 19, units of deg C, mm, or unit-less (for definitions see https://www.worldclim.org/data/bioclim.html)')
		say('gdd5 ..................... total annual growing degree days above 5 deg C, in deg C')
		say('climatic_moisture_index .. climatic moisture index, unit-less: [-1, 1]')
		say('pet_warmest_quarter_mm ... potential evapotranspiration of the warmest 3 month run, in mm')
		say('elevation_m .............. elevation, in m')	
		say('aridity .................. meat annual temperature / total annual precipitation (deg C / mm)')
		
		say('SOIL DATA from SoilGrids Version 2.0 (https://www.isric.org/explore/soilgrids/soilgrids-access)', pre = 1)
		say('All values are for depth 0 to 5 cm.')
		say('ph ....................... pH, measured in water')
		say('CEC ...................... cation exchange capacity, meq/100 g of soil')
		say('sand, silt, clay ......... unit-less (proportion: [0, 1])')
		say('SOC ...................... soil organic carbon (kg of C / m2 ???)')
		

	sink()

say('DONE!', level = 1, deco = '!')

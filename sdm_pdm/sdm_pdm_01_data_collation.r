### MODELING ANDROPOGON GERARDI DISTRIBUTION, PHENOTYPE, PHYSIOLOGY, GENOTYPE, and ASSOCIATED MICROBIAL COMMUNITIES
### Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2023-12
###
### This script constructs a non-integrated model for AG geographic distribution.
###
### source('C:/Kaji/R/andropogon_integratedEcology/sdm_pdm/sdm_pdm_01_data_collation.r')
###
### CONTENTS ###
### setup ###
### collate morphology/physiology/community data with climate ###
### plot geo-folds in environmental and geographic space ###

#############
### setup ###
#############

	rm(list = ls())

	drive <- 'C:/Kaji/'

	setwd(paste0(drive, '/Research/Andropogon/Andropogon'))
	source(paste0(drive, '/R/andropogon_integratedEcology/sdm_pdm/sdm_pdm_00_shared_functions_and_variables.r'))

	library(readxl)

say('#################################################################')
say('### collate morphology/physiology/community data with climate ###')
say('#################################################################')

	### data
	occs <- vect('./data_from_adam_and_loretta/andropogon_gerardi_occurrences_with_environment_1961_2020_climatena.gpkg')

	# # adjust aridity so values of 0 precipitation do not cause issued
	# occs$aridity <- (occs$bio1 + 10) / ((occs$bio12 + 1)/1000)

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

	ppt_at_sites <- terra::extract(ppt, sites_spatial, ID = FALSE)
	tmax_at_sites <- terra::extract(tmax, sites_spatial, ID = FALSE)
	tmin_at_sites <- terra::extract(tmin, sites_spatial, ID = FALSE)

	ppt_at_sites <- as.matrix(ppt_at_sites)
	tmax_at_sites <- as.matrix(tmax_at_sites)
	tmin_at_sites <- as.matrix(tmin_at_sites)

	bioclims <- bcvars(ppt_at_sites, tmin_at_sites, tmax_at_sites)
	bioclims <- as.data.frame(bioclims)
	bioclims$aridity <- (bioclims$bio1 + 10) / ((1 + bioclims$bio12) / 1000) # add aridity

	# # pisr <- rast('E:/ecology/Potential Annual Insolation (SAGA)/Based on ClimateNA 7.03 from SAGA 9.3.0/Annual_Insolation_1990_kW_hr_per_m2.tif')
	# # pisr <- terra::extract(pisr, sites_spatial, ID = FALSE)
	# # pisr <- rowSums(pisr)

	sites <- cbind(sites, bioclims)
	# # sites$solar_rad_kW_hr_per_m2 <- pisr

	### create  different data frames for each type of data
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
		# created these manually from map so that as much as possible, sites close to one another are in the same fold, folds are spatially distinct, and each fold as nearly the same number of sites
		sites$geofold <- NA_integer_
		sites$geofold[sites$site_id %in% c('MO_1', 'MT_1', 'SD_1', 'KS_2', 'NC_1', 'SC_1')] <- 1
		sites$geofold[sites$site_id %in% c('IA_1', 'MN_1', 'ND_1', 'TX_1', 'TX_3', 'MS_1', 'AL_1')] <- 2
		sites$geofold[sites$site_id %in% c('WI_1', 'CO_1', 'NE_2', 'OK_1', 'AR_1', 'IL_1')] <- 3
		sites$geofold[sites$site_id %in% c('IN_1', 'MI_1', 'KS_1', 'NE_1', 'NM_1', 'LA_1', 'TX_2')] <- 4

		sites_vect <- vect(sites, geom = c('LONGITUDE', 'LATITUDE'), crs = 'wgs84')
		# sites_gfolds <- geoFold(sites_vect, k = 4, minIn = 6, method = 'single')

		# geo-folds of biomass values
		biomass$geofold <- sites$geofold[match(biomass$SITE, sites$site_id)]
		morpho_phys$geofold <- sites$geofold[match(morpho_phys$SITE, sites$site_id)]

		# geo-folds of occurrences
		occs$focal_region <- !is.na(occs$n_andropogon_gerardi)
		occs_cents <- centroids(occs)
		sites_vect <- project(sites_vect, occs)
		occs_gfolds <- geoFoldContrast(contrast = occs_cents, pres = sites_vect, presFolds = sites$geofold)
		occs$geofold <- occs_gfolds

		writeVector(occs, './outputs_loretta/integrated_sdm_pdm/andropogon_gerardi_occurrences_with_environment_1961_2020_for_integration.gpkg', overwrite = TRUE)

		dirCreate('./data_from_loretta/sdm_pdm_00_merged_site_data_with_climate')
		saveRDS(sites, './data_from_loretta/sdm_pdm_00_merged_site_data_with_climate/sites.rds')
		saveRDS(biomass, './data_from_loretta/sdm_pdm_00_merged_site_data_with_climate/biomass.rds')
		saveRDS(morpho_phys, './data_from_loretta/sdm_pdm_00_merged_site_data_with_climate/morpho_phys.rds')

# say('############################################################')
# say('### plot geo-folds in environmental and geographic space ###')
# say('############################################################')

# 	# Create a map of the geofolds and plot them in environmental space

# 	# North America, occurrence,s and biomass
# 	nam <- vect('./data_from_gadm/gadm_4pt1_level_1_north_america_sans_alaska_lambert.gpkg')
# 	ag_vect_sq <- vect('./outputs_loretta/integrated_sdm_pdm/andropogon_gerardi_occurrences_with_environment_1961_2020_for_integration.gpkg')
# 	sites <- readRDS('./data_from_loretta/sdm_pdm_00_merged_site_data_with_climate/sites.rds')

# 	sites$id <- sub(sites$site_id, pattern = '_', replacement = '')

# 	### map
# 	#######

# 		sites_vect <- vect(sites, geom = c('LONGITUDE', 'LATITUDE'), crs = getCRS('WGS84'))

# 		# extent
# 		sites_vect <- project(sites_vect, ag_vect_sq)
# 		extent <- ext(sites_vect)
# 		extent <- as.polygons(extent, crs = ag_vect_sq)
# 		extent <- buffer(extent, width = 200 * 1000) # nominal plot extent
# 		extent_display <- buffer(extent, width = 300 * 1000) # larger than plot extent
# 		extent <- ext(extent)
# 		extent <- as.vector(extent)

# 		ag_vect_display <- crop(ag_vect_sq, extent_display)

# 		counties_with_ag <- ag_vect_display[ag_vect_display$n_andropogon_gerardi > 0]
# 		counties_with_ag <- centroids(counties_with_ag)

# 		map <- ggplot() +
# 			layer_spatial(ag_vect_display, aes(fill = factor(geofold)), color = alpha('gray', 0.5)) +
# 			layer_spatial(nam, color = 'gray40', fill = NA, linewidth = 0.3) +
# 			layer_spatial(counties_with_ag, pch = 16, color = alpha('gray20', 0.5), size = 0.35) +
# 			layer_spatial(sites_vect, mapping = aes(size = mBiomass, fill = factor(geofold)), pch = 21, color = 'red') +
# 			scale_fill_viridis_d() +
# 			geom_sf_text(data = st_as_sf(sites_vect), aes(label = id), color = 'white', size = 3, fontface = 'bold') +
# 			scale_size_continuous(name = 'Mean\nbiomass (g)') +
# 			ggtitle('A) Map of geo-folds') +
# 			xlim(extent[1], extent[2]) + ylim(extent[3], extent[4]) +
# 			labs(fill = 'Fold')

# 	### environmental space
# 	#######################

# 		ag_env <- as.data.table(ag_vect_sq)
# 		ag_env_occs <- ag_env[n_andropogon_gerardi > 0]
		
# 		env <- ggplot() +
# 			geom_point(data = ag_env, aes(x = bio1, y = bio12, color = factor(geofold)), alpha = 0.1, pch = 3) +
# 			geom_point(data = ag_env_occs, aes(x = bio1, y = bio12, color = factor(geofold)), alpha = 0.5) +
# 			geom_point(data = sites, aes(x = bio1, y = bio12, fill = factor(geofold), size = mBiomass), color = 'red', pch = 21) +
# 			geom_text(data = sites, aes(x = bio1, y = bio12, label = id), color = 'black', size = 3) +
# 			scale_color_viridis_d(name = 'Fold') +
# 			scale_fill_viridis_d(name = 'Fold') +
# 			scale_size_continuous(name = 'Mean\nbiomass (g)') +
# 			xlab('Mean annual temperature (°C)') + ylab('Total annual precipitation (mm)') +
# 			ggtitle('B) Geo-folds in environmental space')


# 	combo <- map + env
# 	ggsave(combo, filename = './outputs_loretta/integrated_sdm_pdm/map_geofolds.png', width = 17, height = 8, dpi = 300)

say(date())
say('FINIS!', deco = '+', level = 1)

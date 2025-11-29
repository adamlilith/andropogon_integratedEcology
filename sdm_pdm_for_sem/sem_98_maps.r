### MODELING ANDROPOGON GERARDI DISTRIBUTION, PHENOTYPE, PHYSIOLOGY, GENOTYPE, and ASSOCIATED MICROBIAL COMMUNITIES
### Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2023-12
###
### This script constructs select maps of model inputs and outputs.
###
### source('C:/Kaji/R/andropogon_integratedEcology/sdm_pdm_for_sem/sem_98_maps.r')
### 
### CONTENTS
###
### setup ###
### map of occurrence model input and output plus sample sites ###
### map of present-day biomass and probability of zero biomass ###

#############
### setup ###
#############

	rm(list = ls())

	drive <- 'C:/Kaji/'

	setwd(paste0(drive, '/Research/Andropogon/Andropogon'))
	source(paste0(drive, '/R/andropogon_integratedEcology/sdm_pdm_for_sem/sem_00_shared_functions_and_variables.r'))

	library(viridis)

###########################
### user-defined values ###
###########################


# say('##################################################################')
# say('### map of occurrence model input and output plus sample sites ###')
# say('##################################################################')

# 	occs <- vect('./data_from_adam_and_loretta/andropogon_gerardi_occurrences_with_environment_1961_2020_climatena.gpkg')
# 	preds <- vect('./outputs_loretta/sdm_pdm_for_sem/models_occurrence/[occs_poisson~normal_homoscedastic_bio1^2_bio12^2_bio15^2_[bias~1]]/prediction_vector_nam.gpkg')
# 	sites <- readRDS('./data_from_loretta/sem_00_merged_site_data_with_climate/sites.rds')
# 	sites <- vect(sites, geom = c('LONGITUDE', 'LATITUDE'), crs = getCRS('WGS84'))
# 	sites <- project(sites, occs)

# 	nam <- vect('./data_from_gadm/gadm_4pt1_level_1_north_america_sans_alaska_lambert.gpkg')

# 	# extent
# 	extent <- ext(sites)
# 	extent <- as.polygons(extent, crs = occs)
# 	extent <- buffer(extent, width = 200 * 1000) # nominal plot extent
# 	extent_display <- buffer(extent, width = 300 * 1000) # larger than plot extent
# 	extent <- ext(extent)
# 	extent <- as.vector(extent)

# 	pred_vect_display <- crop(preds, extent_display)

# 	pred_vect_display$color <- NA
# 	quants <- quantile(pred_vect_display$N_ag_county_sq, c(0.25, 0.5, 0.75, 0.95))
# 	pred_vect_display$color[pred_vect_display$N_ag_county_sq < quants[1]] <- 0
# 	pred_vect_display$color[pred_vect_display$N_ag_county_sq >= quants[1] & pred_vect_display$N_ag_county_sq < quants[2]] <- 1
# 	pred_vect_display$color[pred_vect_display$N_ag_county_sq >= quants[2] & pred_vect_display$N_ag_county_sq < quants[3]] <- 2
# 	pred_vect_display$color[pred_vect_display$N_ag_county_sq >= quants[3] & pred_vect_display$N_ag_county_sq < quants[4]] <- 3
# 	pred_vect_display$color[pred_vect_display$N_ag_county_sq >= quants[4]] <- 4
# 	pred_vect_display$color <- factor(pred_vect_display$color)

# 	viridis_palette <- viridis_pal(option = 'viridis')(256)
# 	occ_map <- ggplot() +
# 		layer_spatial(pred_vect_display, aes(fill = n_andropogon_gerardi), color = NA) +
# 		scale_fill_gradient(trans = 'log10', na.value = 'gray80', name = 'No. records', low = viridis_palette[1], high = viridis_palette[256]) +
# 		layer_spatial(nam, color = 'gray40', fill = NA, linewidth = 0.3) +
# 		layer_spatial(sites, pch = 21, size = 4, fill = 'white') +
# 		xlim(extent[1], extent[2]) + ylim(extent[3], extent[4]) +
# 		ggtitle('A) Occurrences') +
# 		theme(
# 			plot.title = element_text(size = 22),
# 			plot.subtitle = element_text(size = 14),
# 			legend.title = element_text(size = 20),
# 			legend.text = element_text(size = 20)
# 		)

# 	pred_map <- ggplot() +
# 		layer_spatial(pred_vect_display, aes(fill = color), color = NA) +
# 		scale_fill_manual(
# 			values = c('0' = '#f7fcf5', '1' = '#c7e9c0', '2' = '#74c476', '3' = '#238b45', '4' = '#00441b'),
# 			labels = c('≥0 to <0.25', '≥0.25 to <0.50', '≥0.50 to <0.75', '≥0.75 to <0.95', '≥0.95 to 1'),
# 			name = 'Relative\nabundance\nquantile'
# 		) +
# 		guides(fill = guide_legend(reverse = TRUE)) +
# 		layer_spatial(nam, color = 'gray40', fill = NA, linewidth = 0.3) +
# 		layer_spatial(sites, pch = 21, size = 4, fill = 'white') +
# 		xlim(extent[1], extent[2]) + ylim(extent[3], extent[4]) +
# 		ggtitle('B) Distribution model predictions') +
# 		theme(
# 			plot.title = element_text(size = 22),
# 			plot.subtitle = element_text(size = 14),
# 			legend.title = element_text(size = 20),
# 			legend.text = element_text(size = 20)
# 		)
		
# 	maps <- occ_map + pred_map

# 	filename <- paste0('./outputs_loretta/sem/occs_plus_sdm_from_sdm_pdm_bio1^2_bio12^2_bio15^2_bias~1.png')
# 	ggsave(plot = maps, filename = filename, width = 22, height = 10, dpi = 600)
	
say('##################################################################')
say('### map of present-day biomass and probability of zero biomass ###')
say('##################################################################')

	out_dir <- './outputs_loretta/sdm_pdm_for_sem/model_biomass~zig~bio1^2_bio12'
	bioclims <- c(1, 12)
	zero_inflated <- TRUE

	nam <- vect('./data_from_gadm/gadm_4pt1_level_1_north_america_sans_alaska_lambert.gpkg')
	pred_vect_nam <- vect(paste0(out_dir, '/prediction_vector_nam.gpkg'))

	data_biomass <- prepare_biomass(formula_biomass = ~ 1 + bio1 + bio12 + I(bio1^2), n_response_curve_values = n_response_curve_values, calib = FALSE)

	# extent
	site_vect <- data_biomass$site_vect_biomass
	site_vect <- project(site_vect, pred_vect_nam)
	extent <- ext(site_vect)
	extent <- as.polygons(extent, crs = pred_vect_nam)
	extent <- buffer(extent, width = 300 * 1000) # nominal plot extent
	extent_display <- buffer(extent, width = 300 * 1000) # larger than plot extent
	extent <- ext(extent)

	### BIOCLIMs
	############

	ppt <- rast(paste0('C:/Kaji/Research Data/ClimateNA/v 7.3 AdaptWest/1991-2020/Normal_1991_2020_PPT', prefix(1:12, 2), '.tif'))
	tmin <- rast(paste0('C:/Kaji/Research Data/ClimateNA/v 7.3 AdaptWest/1991-2020/Normal_1991_2020_Tmin', prefix(1:12, 2), '.tif'))
	tmax <- rast(paste0('C:/Kaji/Research Data/ClimateNA/v 7.3 AdaptWest/1991-2020/Normal_1991_2020_Tmax', prefix(1:12, 2), '.tif'))
	tmean <- rast(paste0('C:/Kaji/Research Data/ClimateNA/v 7.3 AdaptWest/1991-2020/Normal_1991_2020_Tave', prefix(1:12, 2), '.tif'))

	ppt <- crop(ppt, extent_display)
	tmin <- crop(tmin, extent_display)
	tmax <- crop(tmax, extent_display)
	tmean <- crop(tmean, extent_display)

	fact <- 8
	ppt <- aggregate(ppt, fact = fact, fun = 'mean', na.rm = TRUE)
	tmin <- aggregate(tmin, fact = fact, fun = 'mean', na.rm = TRUE)
	tmax <- aggregate(tmax, fact = fact, fun = 'mean', na.rm = TRUE)
	tmean <- aggregate(tmean, fact = fact, fun = 'mean', na.rm = TRUE)

	bios_rasts <- fasterRaster::bioclims(ppt = ppt, tmin = tmin, tmax = tmax, tmean = tmean, bios = bioclims, verbose = TRUE)

	extent <- as.vector(extent)
	pred_vect_display <- crop(pred_vect_nam, extent_display)

	### predict biomass
	chains <- readRDS(paste0(out_dir, '/chains.rds'))
	meta <- readRDS(paste0(out_dir, '/!meta_generic.rds'))
	form <- meta$formulae$formula_biomass_mu

	bios_cells <- as.data.table(bios_rasts, cell = TRUE)
	bioclim_names <- paste0('bio', bioclims)
	bios <- bios_cells[ , ..bioclim_names]
	bios <- scale(bios, center = data_biomass$x_centers_biomass, scale = data_biomass$x_scales_biomass)
	x <- model.matrix(form, as.data.table(bios))
	mu <- predict_biomass(chains = chains, x = x, zero_inflated = zero_inflated, type = 'mu', pre_averaged = TRUE)

	mu <- setValueByCell(ppt[[1]], val = mu, cell = bios_cells$cell)
	names(mu) <- 'mu'
	core_quant <- global(mu, fun = quantile, prob = ag_core_quant, na.rm = TRUE)

	### predict pzero
	form <- meta$formulae$formula_pzero

	bios_cells <- as.data.table(bios_rasts, cell = TRUE)
	bioclim_names <- paste0('bio', bioclims)
	bios <- bios_cells[ , ..bioclim_names]
	bios <- scale(bios, center = data_biomass$x_centers_biomass, scale = data_biomass$x_scales_biomass)
	x <- model.matrix(form, as.data.table(bios))
	pzero <- predict_biomass(chains = chains, x = x, zero_inflated = zero_inflated, type = 'pzero', pre_averaged = TRUE)
	# pzero <- colMeans(preds)

	pzero <- setValueByCell(ppt[[1]], val = pzero, cell = bios_cells$cell)
	names(pzero) <- 'pzero'

	### BIOMASS
	###########

	range_core <- mu >= core_quant[1, 1]
	range_core[range_core == 0] <- NA
	range_core <- as.polygons(range_core)

	# get range of values for plotting
	biomass_stats <- data_biomass$raw_data_biomass[ , .(biomass_mean = mean(Biomass), biomass_sd = sd(Biomass)), by = SITE]
	max_val <- max(biomass_stats$biomass_mean, globalx(mu, 'max'))
	max_val <- 1.01 * max_val
	resp_limits <- c(1, max_val)

	counties_with_ag <- pred_vect_display[pred_vect_display$n_andropogon_gerardi > 0]
	counties_with_ag <- centroids(counties_with_ag)

	biomass_map <- ggplot() +
		layer_spatial(mu) +
		scale_fill_gradientn(
			name = 'Biomass (g)',
			colors = c('#edf8e9', '#74c476', '#005a32'),
			limits = resp_limits,
			na.value = 'transparent',
			trans = 'log10'
		) +
		layer_spatial(nam, color = 'black', fill = NA, linewidth = 0.3) +
		layer_spatial(counties_with_ag, pch = 16, color = 'black', size = 0.5) +
		layer_spatial(range_core, color = 'cyan', fill = NA, linewidth = 0.8) +
		layer_spatial(data_biomass$site_vect_biomass, aes(fill = biomass_mean), pch = 21, size = 4) +
		xlim(extent[1], extent[2]) + ylim(extent[3], extent[4]) +
		coord_sf(expand = FALSE) +
		ggtitle('A) Individual aboveground biomass') +
		theme(
			legend.position = 'bottom',
			legend.key.width = unit(1, 'in'),
			plot.title = element_text(size = 30),
			legend.title = element_text(size = 26),
			legend.text = element_text(size = 16)
		)

	### PROBABILITY OF ZERO BIOMASS
	###############################

	say('pzero')

	core_quant <- global(pzero, fun = quantile, prob = ag_core_quant, na.rm = TRUE)

	range_core <- pzero >= core_quant[1, 1]
	range_core[range_core == 0] <- NA
	range_core <- as.polygons(range_core)

	pzero_map <- ggplot() +
		layer_spatial(pzero) +
		scale_fill_gradientn(
			name = 'Probability',
			colors = c('#253494', '#41b6c4', '#ffffcc'),
			na.value = 'transparent'
		) +
		layer_spatial(nam, color = 'black', fill = NA, linewidth = 0.3) +
		layer_spatial(counties_with_ag, pch = 16, color = 'gray', size = 0.5) +
		layer_spatial(range_core, color = 'cyan', fill = NA, linewidth = 0.8) +
		layer_spatial(site_vect, pch = 3, size = 4, color = 'orange') +
		xlim(extent[1], extent[2]) + ylim(extent[3], extent[4]) +
		coord_sf(expand = FALSE) +
		ggtitle('B) Probability of zero biomass (absence)') +
		theme(
			legend.position = 'bottom',
			legend.key.width = unit(1, 'in'),
			plot.title = element_text(size = 30),
			legend.title = element_text(size = 26),
			legend.text = element_text(size = 22)
		)

	maps <- biomass_map + pzero_map
	ggsave(plot = maps, filename = paste0(out_dir, '/map_biomass_pzero.png'), width = 20, height = 10, dpi = 600)

say(date())
say('FINIS!', deco = '+', level = 1)

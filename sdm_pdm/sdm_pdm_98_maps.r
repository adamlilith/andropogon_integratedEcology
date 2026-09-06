### MODELING ANDROPOGON GERARDI DISTRIBUTION, PHENOTYPE, PHYSIOLOGY, GENOTYPE, and ASSOCIATED MICROBIAL COMMUNITIES
### Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2023-12
###
### This script constructs select maps of model inputs and outputs.
###
### source('C:/Kaji/R/andropogon_integratedEcology/sdm_pdm/sdm_pdm_98_maps.r')
###
### CONTENTS ###
### setup ###
### map of occurrence model input and output plus sample sites ###
### predict integrated models to North America and Dust Bowl ###
### predict abundance-only model to North America while enforcing presence ###
### predict biomass and non-biomass facet models to North America while enforcing presence ###
###
### maps of abundance and change in abundance for main text ###
### maps of biomass and change in biomass for main text ###
###
### maps of abundance and change in abundance for supplement ###
### maps of biomass and change in biomass for supplement ###
### maps of non-biomass traits and change in non-biomass traits for supplement ###
### 
### maps of abundance in Dust Bowl region in 1930s for supplement ###
### maps of biomass in Dust Bowl region in 1930s for supplement ###
### maps of non-biomass in Dust Bowl region in 1930s for supplement ###
### maps of change in Dust Bowl region in 1930s for main text ###
###
#############
### setup ###
#############

	rm(list = ls())

	drive <- 'C:/Kaji/'

	setwd(paste0(drive, '/Research/Andropogon/Andropogon'))
	source(paste0(drive, '/R/andropogon_integratedEcology/sdm_pdm/sdm_pdm_00_shared_functions_and_variables.r'))

	library(viridis)

###########################
### user-defined values ###
###########################

# say('##################################################################')
# say('### map of occurrence model input and output plus sample sites ###')
# say('##################################################################')

# 	occs <- vect('./data_from_adam_and_loretta/andropogon_gerardi_occurrences_with_environment_1961_2020.gpkg')
# 	preds <- vect('./outputs_loretta/integrated_sdm_pdm/models_occurrence/[occs_poisson~normal_homoscedastic_bio1^2_bio12^2_bio15^2_[bias~1]]/prediction_vector_nam.gpkg')
# 	sites <- readRDS('./data_from_loretta/sdm_pdm_00_merged_site_data_with_climate/sites.rds')
# 	sites <- vect(sites, geom = c('LONGITUDE', 'LATITUDE'), crs = getCRS('WGS84'))
# 	sites <- project(sites, occs)

# 	nam <- vect('./data_from_gadm/gadm_4pt1_level_1_north_america_sans_alaska_lambert.gpkg')
#	nam <- simplifyGeom(nam, tolerance = 1000)

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

# say('################################################################')
# say('### predict integrated models to North America and Dust Bowl ###')
# say('################################################################')

# 	# model_dir <- './outputs_loretta/integrated_sdm_pdm/models_integrated/non_centered_MVN([occs~hurdlepoisson_offset]_[biomass~hurdleln]_[bla_can_hei])_eta=2'
# 	model_dir <- './outputs_loretta/integrated_sdm_pdm/models_integrated/non_centered_MVN([occs~hurdlepoisson_offset]_[biomass~hurdleln]_[cnr_pho_spa_sto_tra])_eta=2'

# 	say('Predicting to:\n', model_dir)

# 	formulae <- readRDS(paste0(model_dir, '/formulae.rds'))

# 	### chains
# 	##########
# 	# load chains, subset to fewest number of iterations done, then thin to 1000 each
# 	chains <- readRDS(paste0(model_dir, '/chains.rds'))
# 	niter <- Inf
# 	for (i in 1:4) niter <- min(niter, nrow(chains$samples[[i]]))
# 	keeps <- round(seq((niter / 2) + 1, niter, length.out = 1000))
# 	for (i in 1:4) chains$samples[[i]] <- chains$samples[[i]][keeps, ]

# 	### prepare data: North America
# 	###############################

# 	# data for OCCURRENCES at counties
# 	data_occs_counties <- prepare_occurrence_data(
# 		formula_occs = formulae$formula_occs,
# 		formula_bias = ~ 1,
# 		n_response_curve_values = n_response_curve_values,
# 		psa_quant = psa_quant,
# 		calib = FALSE
# 	)

# 	# data for ZERO-INFLATION at counties
# 	data_psi_counties <- prepare_occurrence_data(
# 		formula_occs = formulae$formula_psi,
# 		formula_bias = ~ 1,
# 		n_response_curve_values = n_response_curve_values,
# 		psa_quant = psa_quant,
# 		calib = FALSE
# 	)

# 	# data for BIOMASS using county-level environment
# 	data_biomass_counties <- prepare_occurrence_data(
# 		formula_occs = formulae$formula_biomass,
# 		formula_bias = ~ 1,
# 		n_response_curve_values = n_response_curve_values,
# 		psa_quant = psa_quant,
# 		calib = FALSE
# 	)

# 	# data for NON-BIOMASS FACETS at sites and counties
# 	nonbiomass_facets <- formulae$nonbiomass_facets
# 	data_nonbiomass_sites <- data_nonbiomass_counties <- list()
# 	for (f in seq_along(nonbiomass_facets)) {

# 		facet <- names(nonbiomass_facets)[f]
# 		say('preparing data for ', facet, ' ', date())
		
# 		formula_facet <- nonbiomass_facets[[facet]]$formula

# 		data_nonbiomass_counties[[f]] <- prepare_occurrence_data(
# 			formula_occs = formula_facet,
# 			formula_bias = ~ 1,
# 			n_response_curve_values = n_response_curve_values,
# 			psa_quant = psa_quant,
# 			calib = FALSE
# 		)

# 	}
# 	names(data_nonbiomass_counties) <- names(nonbiomass_facets)

# 	### predict
# 	pred_vect_nam <- burn_fully_integrated_into_vector(
# 		demesne = 'nam',
# 		chains = chains,
# 		formula_occs = formulae$formula_occs,
# 		formula_bias = formulae$formula_bias,
# 		formula_psi = formulae$formula_psi,
# 		formula_biomass = formulae$formula_biomass,
# 		resp_distrib_biomass = formulae$meta_biomass$resp_distrib_biomass,
# 		transform_biomass = formulae$meta_biomass$transform_biomass,
# 		log_precip_biomass = formulae$meta_biomass$log_precip_biomass,
# 		nonbiomass_facets = formulae$nonbiomass_facets
# 	)
# 	writeVector(pred_vect_nam, paste0(model_dir, '/prediction_vector_nam.gpkg'), overwrite = TRUE)

# 	pred_vect_nam <- burn_fully_integrated_into_vector(
# 		demesne = 'nam',
# 		chains = chains,
# 		formula_occs = formulae$formula_occs,
# 		formula_bias = formulae$formula_bias,
# 		formula_psi = formulae$formula_psi,
# 		formula_biomass = formulae$formula_biomass,
# 		resp_distrib_biomass = formulae$meta_biomass$resp_distrib_biomass,
# 		transform_biomass = formulae$meta_biomass$transform_biomass,
# 		log_precip_biomass = formulae$meta_biomass$log_precip_biomass,
# 		nonbiomass_facets = formulae$nonbiomass_facets,
# 		force_presence = TRUE
# 	)
# 	writeVector(pred_vect_nam, paste0(model_dir, '/prediction_vector_nam_force_presence.gpkg'), overwrite = TRUE)

# 	pred_vect_1930s <- burn_fully_integrated_into_vector(
# 		demesne = '1930s',
# 		chains = chains,
# 		formula_occs = formulae$formula_occs,
# 		formula_bias = formulae$formula_bias,
# 		formula_psi = formulae$formula_psi,
# 		formula_biomass = formulae$formula_biomass,
# 		resp_distrib_biomass = formulae$meta_biomass$resp_distrib_biomass,
# 		transform_biomass = formulae$meta_biomass$transform_biomass,
# 		log_precip_biomass = formulae$meta_biomass$log_precip_biomass,
# 		nonbiomass_facets = formulae$nonbiomass_facets
# 	)
# 	writeVector(pred_vect_1930s, paste0(model_dir, '/prediction_vector_1930s.gpkg'), overwrite = TRUE)

# 	pred_vect_1930s <- burn_fully_integrated_into_vector(
# 		demesne = '1930s',
# 		chains = chains,
# 		formula_occs = formulae$formula_occs,
# 		formula_bias = formulae$formula_bias,
# 		formula_psi = formulae$formula_psi,
# 		formula_biomass = formulae$formula_biomass,
# 		resp_distrib_biomass = formulae$meta_biomass$resp_distrib_biomass,
# 		transform_biomass = formulae$meta_biomass$transform_biomass,
# 		log_precip_biomass = formulae$meta_biomass$log_precip_biomass,
# 		nonbiomass_facets = formulae$nonbiomass_facets,
# 		force_presence = TRUE
# 	)
# 	writeVector(pred_vect_1930s, paste0(model_dir, '/prediction_vector_1930s_force_presence.gpkg'), overwrite = TRUE)

# say('##############################################################################')
# say('### predict abundance-only model to North America while enforcing presence ###')
# say('##############################################################################')

# 		dir_model <- './outputs_loretta/integrated_sdm_pdm/models_occurrence/[occs~hurdlepoisson(bio1^2_log(bio12)^2_bio15^2)]_[bias~poaceae_offset_plus_exp_neg_1_WITH_zero_correction_covariate]'

# 		chains <- readRDS(paste0(dir_model, '/chains.rds'))
# 		formula <- readRDS(paste0(dir_model, '/formulae.rds'))
# 		pred_vect_nam_force_presence <- burn_occs_into_vector(
# 			demesne = 'nam',
# 			chains = chains,
# 			formula_occs = formula$formula_occs,
# 			formula_psi = formula$formula_psi,
# 			overdispersed = FALSE,
# 			force_presence = TRUE
# 		)

# 		writeVector(pred_vect_nam_force_presence, paste0(dir_model, '/prediction_vector_nam_force_presence.gpkg'), overwrite = TRUE)

# say('##############################################################################################')
# say('### predict biomass and non-biomass facet models to North America while enforcing presence ###')
# say('##############################################################################################')

# 	# The workflows run after each model was calibrated did not force presence when predictions were made to a vector. For comparisons, we want to have vectors that force presence fromm these single-facet models.

# 	### predict biomass while NOT enforcing presence
# 	################################################

# 		dir_model <- './outputs_loretta/integrated_sdm_pdm/models_biomass/[biomass~hurdleln_log(bio12)]'

# 		chains <- readRDS(paste0(dir_model, '/chains.rds'))
# 		pred_vect_nam <- burn_biomass_into_vector(
# 			demesne = 'nam',
# 			chains = chains,
# 			formula_biomass = ~ 1 + bio12_log10p1,
# 			formula_psi = ~ 1 + bio12_log10p1,
# 			resp_distrib = 'hurdleLN',
# 			transform = 'identity',
# 			force_presence = TRUE
# 		)

# 		writeVector(pred_vect_nam, paste0(dir_model, '/prediction_vector_nam_force_presence.gpkg'), overwrite = TRUE)

# 	### predict non-biomass facets while NOT enforcing presence
# 	###########################################################

# 		for (f in seq_along(nonbiomass_facets)) {

# 			facet <- names(nonbiomass_facets)[f]
# 			say(facet, level = 2)

# 			resp_distrib <- nonbiomass_facets[[f]]$resp_distrib
# 			filename <- nonbiomass_facets[[f]]$filename
# 			formula <- nonbiomass_facets[[f]]$formula
# 			transform <- nonbiomass_facets[[f]]$transform

# 			dir_model <- paste0('./outputs_loretta/integrated_sdm_pdm/models_', facet, '/[', facet, '~', tolower(resp_distrib), '(', filename, ')]')

# 			chains <- readRDS(paste0(dir_model, '/chains.rds'))
# 			pred_vect_nam <- burn_nonbiomass_into_vector(
# 				facet = facet,
# 				demesne = 'nam',
# 				chains = chains,
# 				formula_facet = formula,
# 				formula_psi = formula,
# 				resp_distrib = resp_distrib,
# 				transform = transform,
# 				force_presence = TRUE
# 			)

# 			writeVector(pred_vect_nam , paste0(dir_model, '/prediction_vector_nam_force_presence.gpkg'), overwrite = TRUE)

# 		}

# say('###############################################################')
# say('### maps of abundance and change in abundance for main text ###')
# say('###############################################################')

# 	### maps of abundance in 3-row grid, 1st two rows have 3 columns, last has two:
# 		# abund standalone SQ   abund standalne fut		abund standalone future delta
# 		# abund integrated SQ   abund integrated fut	abund integrated future delta
# 		# delta integrated and standalone sq	delta integrated and standalone sq

# 	nam <- vect('./data_from_gadm/gadm_4pt1_level_1_north_america_sans_alaska_lambert.gpkg')
# 	nam <- simplifyGeom(nam, tolerance = 1000)

# 	multi <- vect('./outputs_loretta/integrated_sdm_pdm/models_integrated/non_centered_MVN([occs~hurdlepoisson_offset]_[biomass~hurdleln]_[bla_can_hei])_eta=2/prediction_vector_nam.gpkg')
# 	multi_force_presence <- vect('./outputs_loretta/integrated_sdm_pdm/models_integrated/non_centered_MVN([occs~hurdlepoisson_offset]_[biomass~hurdleln]_[bla_can_hei])_eta=2/prediction_vector_nam_force_presence.gpkg')

# 	uni <- vect('./outputs_loretta/integrated_sdm_pdm/models_occurrence/[occs~hurdlepoisson(bio1^2_log(bio12)^2_bio15^2)]_[bias~poaceae_offset_plus_exp_neg_1_WITH_zero_correction_covariate]/prediction_vector_nam.gpkg')
# 	uni_force_presence <- vect('./outputs_loretta/integrated_sdm_pdm/models_occurrence/[occs~hurdlepoisson(bio1^2_log(bio12)^2_bio15^2)]_[bias~poaceae_offset_plus_exp_neg_1_WITH_zero_correction_covariate]/prediction_vector_nam_force_presence.gpkg')

# 	sites <- readRDS('./outputs_loretta/integrated_sdm_pdm/models_integrated/data_biomass_sites_bla_can_hei.rds')
# 	site_vect <- vect(sites$site_data_raw, geom = c('LONGITUDE', 'LATITUDE'), crs = getCRS('WGS84'))
# 	site_vect <- project(site_vect, nam)

# 	# extent
# 	extent <- ext(site_vect)
# 	extent <- as.polygons(extent, crs = nam)
# 	extent <- buffer(extent, width = 320 * 1000) # nominal plot extent
# 	extent <- ext(extent)
# 	extent_coords <- as.vector(extent)

# 	nam <- crop(nam, extent)
# 	multi <- crop(multi, extent)
# 	multi_force_presence <- crop(multi_force_presence, extent)
# 	uni <- crop(uni, extent)

# 	# delineate current range "core"
# 	range_core_uni_sq <- delineate_range_core(uni, column = 'N_ag_mean_sq', core_quant = core_quant)
# 	range_core_uni_fut <- delineate_range_core(uni, column = 'N_ag_mean_ssp245_2071_2100', core_quant = core_quant)

# 	range_anticore_uni_sq <- delineate_range_core(uni, column = 'N_ag_mean_sq', core_quant = anticore_quant, rule = '<')
# 	range_anticore_uni_fut <- delineate_range_core(uni, column = 'N_ag_mean_ssp245_2071_2100', core_quant = anticore_quant, rule = '<')

# 	range_core_multi_sq <- delineate_range_core(multi, column = 'N_ag_mean_sq', core_quant = core_quant)
# 	range_core_multi_fut <- delineate_range_core(multi, column = 'N_ag_mean_ssp370_2071_2100', core_quant = core_quant)

# 	range_anticore_multi_sq <- delineate_range_core(multi, column = 'N_ag_mean_sq', core_quant = anticore_quant, rule = '<')
# 	range_anticore_multi_fut <- delineate_range_core(multi, column = 'N_ag_mean_ssp370_2071_2100', core_quant = anticore_quant, rule = '<')

# 	range_core_multi_sq_force_presence <- delineate_range_core(multi_force_presence, column = 'N_ag_mean_sq', core_quant = core_quant)
# 	range_core_multi_fut_force_presence <- delineate_range_core(multi_force_presence, column = 'N_ag_mean_ssp370_2071_2100', core_quant = core_quant)

# 	range_anticore_multi_sq_force_presence <- delineate_range_core(multi_force_presence, column = 'N_ag_mean_sq', core_quant = anticore_quant, rule = '<')
# 	range_anticore_multi_fut_force_presence <- delineate_range_core(multi_force_presence, column = 'N_ag_mean_ssp370_2071_2100', core_quant = anticore_quant, rule = '<')

# 	counties_with_ag <- uni[uni$n_andropogon_gerardi > 0]
# 	counties_with_ag <- centroids(counties_with_ag)

# 	# across both models, get max value of predictions for scaling
# 	vars <- c('N_ag_mean_sq', 'N_ag_mean_ssp245_2041_2070', 'N_ag_mean_ssp245_2071_2100', 'N_ag_mean_ssp370_2041_2070', 'N_ag_mean_ssp370_2071_2100')
# 	# vars <- c('N_ag_median_sq', 'N_ag_median_ssp245_2041_2070', 'N_ag_median_ssp245_2071_2100', 'N_ag_median_ssp370_2041_2070', 'N_ag_median_ssp370_2071_2100')

# 	max_val <- -Inf
# 	for (this_vect in c('multi', 'uni')) {
# 		x <- get(this_vect)
# 		for (var in vars) {
# 			x_var <- x[[var]]
# 			x_var <- unlist(x_var)
# 			max_val <- max(max_val, max(x_var, na.rm = TRUE))
# 			# max_val <- max(max_val, quantile(x_var, 0.99, na.rm = TRUE))
# 		}
# 	}
# 	for (this_vect in c('multi_force_presence', 'uni_force_presence')) {
# 		x <- get(this_vect)
# 		for (var in vars) {
# 			x_var <- x[[var]]
# 			x_var <- unlist(x_var)
# 			max_val <- max(max_val, max(x_var, na.rm = TRUE))
# 			# max_val <- max(max_val, quantile(x_var, 0.99, na.rm = TRUE))
# 		}
# 	}

# 	resp_limits <- c(0, max_val)

# 	title <- 'a) Abundance-only: Present\nAbsence allowed'
# 	legend_title <- 'Abundance'
# 	sq_uni <- ggplot() +
# 		layer_spatial(uni, aes(fill = N_ag_mean_sq), color = NA) +
# 		layer_spatial(nam, color = 'gray40', fill = NA, linewidth = 0.3) +
# 		layer_spatial(counties_with_ag, pch = 16, color = alpha('gray20', 0.3), size = 0.15) +
# 		layer_spatial(range_core_uni_sq, color = 'orange2', fill = NA, linewidth = 0.4) +
# 		layer_spatial(range_anticore_uni_sq, color = 'red', fill = NA, linewidth = 0.4) +
# 		layer_spatial(site_vect, pch = 3, size = 1) +
# 		scale_fill_gradientn(
# 			name = legend_title,
# 			colors = c('#edf8e9', '#74c476', '#005a32'),
# 			limits = resp_limits
# 		) +
# 		ggtitle(title)

# 	title <- 'g) Abundance-only: Present\nPresence forced'
# 	legend_title <- 'Abundance'
# 	sq_uni_force_presence <- ggplot() +
# 		layer_spatial(uni_force_presence, aes(fill = N_ag_mean_sq), color = NA) +
# 		layer_spatial(nam, color = 'gray40', fill = NA, linewidth = 0.3) +
# 		layer_spatial(counties_with_ag, pch = 16, color = alpha('gray20', 0.3), size = 0.15) +
# 		layer_spatial(range_core_uni_sq, color = 'orange2', fill = NA, linewidth = 0.4) +
# 		layer_spatial(range_anticore_uni_sq, color = 'red', fill = NA, linewidth = 0.4) +
# 		layer_spatial(site_vect, pch = 3, size = 1) +
# 		scale_fill_gradientn(
# 			name = legend_title,
# 			colors = c('#edf8e9', '#74c476', '#005a32'),
# 			limits = resp_limits
# 		) +
# 		ggtitle(title)

# 	title <- 'd) Integrated: Present\nAbsence allowed'
# 	legend_title <- 'Abundance'
# 	sq_multi <- ggplot() +
# 		layer_spatial(multi, aes(fill = N_ag_mean_sq), color = NA) +
# 		layer_spatial(nam, color = 'gray40', fill = NA, linewidth = 0.3) +
# 		layer_spatial(counties_with_ag, pch = 16, color = alpha('gray20', 0.3), size = 0.15) +
# 		layer_spatial(range_core_multi_sq, color = 'orange2', fill = NA, linewidth = 0.4) +
# 		layer_spatial(range_anticore_multi_sq, color = 'red', fill = NA, linewidth = 0.4) +
# 		layer_spatial(site_vect, pch = 3, size = 1) +
# 		scale_fill_gradientn(
# 			name = legend_title,
# 			colors = c('#edf8e9', '#74c476', '#005a32'),
# 			limits = resp_limits
# 		) +
# 		ggtitle(title)

# 	title <- 'j) Integrated: Present\nPresence forced'
# 	legend_title <- 'Abundance'
# 	sq_multi_force_presence <- ggplot() +
# 		layer_spatial(multi_force_presence, aes(fill = N_ag_mean_sq), color = NA) +
# 		layer_spatial(nam, color = 'gray40', fill = NA, linewidth = 0.3) +
# 		layer_spatial(counties_with_ag, pch = 16, color = alpha('gray20', 0.3), size = 0.15) +
# 		layer_spatial(range_core_multi_sq_force_presence, color = 'orange2', fill = NA, linewidth = 0.4) +
# 		layer_spatial(range_anticore_multi_sq_force_presence, color = 'red', fill = NA, linewidth = 0.4) +
# 		layer_spatial(site_vect, pch = 3, size = 1) +
# 		scale_fill_gradientn(
# 			name = legend_title,
# 			colors = c('#edf8e9', '#74c476', '#005a32'),
# 			limits = resp_limits
# 		) +
# 		ggtitle(title)

# 	title <- 'b) Abundance-only: Future\nAbsence allowed'
# 	legend_title <- 'Abundance'
# 	fut_uni <- ggplot() +
# 		layer_spatial(uni, aes(fill = N_ag_mean_ssp370_2071_2100), color = NA) +
# 		layer_spatial(nam, color = 'gray40', fill = NA, linewidth = 0.3) +
# 		layer_spatial(counties_with_ag, pch = 16, color = alpha('gray20', 0.3), size = 0.15) +
# 		layer_spatial(range_core_uni_fut, color = 'orange2', fill = NA, linewidth = 0.4) +
# 		layer_spatial(range_anticore_uni_fut, color = 'red', fill = NA, linewidth = 0.4) +
# 		layer_spatial(site_vect, pch = 3, size = 1) +
# 		scale_fill_gradientn(
# 			name = legend_title,
# 			colors = c('#edf8e9', '#74c476', '#005a32'),
# 			limits = resp_limits
# 		) +
# 		ggtitle(title)

# 	title <- 'h) Abundance-only: Future\nPresence forced'
# 	legend_title <- 'Abundance'
# 	fut_uni_force_presence <- ggplot() +
# 		layer_spatial(uni_force_presence, aes(fill = N_ag_mean_ssp370_2071_2100), color = NA) +
# 		layer_spatial(nam, color = 'gray40', fill = NA, linewidth = 0.3) +
# 		layer_spatial(counties_with_ag, pch = 16, color = alpha('gray20', 0.3), size = 0.15) +
# 		layer_spatial(range_core_uni_sq, color = 'orange2', fill = NA, linewidth = 0.4) +
# 		layer_spatial(range_anticore_uni_sq, color = 'red', fill = NA, linewidth = 0.4) +
# 		layer_spatial(site_vect, pch = 3, size = 1) +
# 		scale_fill_gradientn(
# 			name = legend_title,
# 			colors = c('#edf8e9', '#74c476', '#005a32'),
# 			limits = resp_limits
# 		) +
# 		ggtitle(title)

# 	title <- 'e) Integrated: Future\nAbsence allowed'
# 	legend_title <- 'Abundance'
# 	fut_multi <- ggplot() +
# 		layer_spatial(multi, aes(fill = N_ag_mean_ssp370_2071_2100), color = NA) +
# 		layer_spatial(nam, color = 'gray40', fill = NA, linewidth = 0.3) +
# 		layer_spatial(counties_with_ag, pch = 16, color = alpha('gray20', 0.3), size = 0.15) +
# 		layer_spatial(range_core_multi_fut, color = 'orange2', fill = NA, linewidth = 0.4) +
# 		layer_spatial(range_anticore_multi_fut, color = 'red', fill = NA, linewidth = 0.4) +
# 		layer_spatial(site_vect, pch = 3, size = 1) +
# 		scale_fill_gradientn(
# 			name = legend_title,
# 			colors = c('#edf8e9', '#74c476', '#005a32'),
# 			limits = resp_limits
# 		) +
# 		ggtitle(title)

# 	title <- 'k) Integrated: Future\nPresence forced'
# 	legend_title <- 'Abundance'
# 	fut_multi_force_presence <- ggplot() +
# 		layer_spatial(multi_force_presence, aes(fill = N_ag_mean_ssp370_2071_2100), color = NA) +
# 		layer_spatial(nam, color = 'gray40', fill = NA, linewidth = 0.3) +
# 		layer_spatial(counties_with_ag, pch = 16, color = alpha('gray20', 0.3), size = 0.15) +
# 		layer_spatial(range_core_multi_fut_force_presence, color = 'orange2', fill = NA, linewidth = 0.4) +
# 		layer_spatial(range_anticore_multi_fut_force_presence, color = 'red', fill = NA, linewidth = 0.4) +
# 		layer_spatial(site_vect, pch = 3, size = 1) +
# 		scale_fill_gradientn(
# 			name = legend_title,
# 			colors = c('#edf8e9', '#74c476', '#005a32'),
# 			limits = resp_limits
# 		) +
# 		ggtitle(title)

# 	uni$delta_sq_vs_fut <- uni$N_ag_mean_ssp370_2071_2100 - uni$N_ag_mean_sq
# 	multi$delta_sq_vs_fut <- multi$N_ag_mean_ssp370_2071_2100 - multi$N_ag_mean_sq

# 	uni_force_presence$delta_sq_vs_fut <- uni_force_presence$N_ag_mean_ssp370_2071_2100 - uni_force_presence$N_ag_mean_sq
# 	multi_force_presence$delta_sq_vs_fut <- multi_force_presence$N_ag_mean_ssp370_2071_2100 - multi_force_presence$N_ag_mean_sq

# 	delta_limits <- range(c(uni$delta_sq_vs_fut, uni_force_presence$delta_sq_vs_fut, multi$delta_sq_vs_fut, multi_force_presence$delta_sq_vs_fut))
	
# 	rescale_delta_limits <- c(-max(abs(delta_limits)), max(abs(delta_limits)))

# 	title <- 'c) Abundance-only: Change\nAbsence allowed'
# 	legend_title <- 'Change'
# 	pm <- c(rep('', 10), rep('+', 11))
# 	fut_delta_uni <- ggplot() +
# 		layer_spatial(uni, aes(fill = delta_sq_vs_fut), color = NA) +
# 		layer_spatial(nam, color = 'gray40', fill = NA, linewidth = 0.3) +
# 		layer_spatial(counties_with_ag, pch = 16, color = alpha('gray20', 0.3), size = 0.15) +
# 		layer_spatial(site_vect, pch = 3, size = 1) +
# 		scale_fill_gradientn(
# 			name = legend_title,
# 			colors = c('#8c510a', '#d8b365', '#f6e8c3', '#f5f5f5', '#c7eae5', '#5ab4ac', '#01665e'),
# 			limits = rescale_delta_limits
# 		) +
# 		ggtitle(title)

# 	title <- 'f) Integrated: Change\nAbsence allowed'
# 	legend_title <- 'Change'
# 	fut_delta_multi <- ggplot() +
# 		layer_spatial(multi, aes(fill = delta_sq_vs_fut), color = NA) +
# 		layer_spatial(nam, color = 'gray40', fill = NA, linewidth = 0.3) +
# 		layer_spatial(counties_with_ag, pch = 16, color = alpha('gray20', 0.3), size = 0.15) +
# 		layer_spatial(site_vect, pch = 3, size = 1) +
# 		scale_fill_gradientn(
# 			name = legend_title,
# 			colors = c('#8c510a', '#d8b365', '#f6e8c3', '#f5f5f5', '#c7eae5', '#5ab4ac', '#01665e'),
# 			limits = rescale_delta_limits
# 		) +
# 		ggtitle(title)

# 	title <- 'i) Abundance-only: Change\nPresence forced'
# 	legend_title <- 'Change'
# 	pm <- c(rep('', 10), rep('+', 11))
# 	fut_delta_uni_force_presence <- ggplot() +
# 		layer_spatial(uni_force_presence, aes(fill = delta_sq_vs_fut), color = NA) +
# 		layer_spatial(nam, color = 'gray40', fill = NA, linewidth = 0.3) +
# 		layer_spatial(counties_with_ag, pch = 16, color = alpha('gray20', 0.3), size = 0.15) +
# 		layer_spatial(site_vect, pch = 3, size = 1) +
# 		scale_fill_gradientn(
# 			name = legend_title,
# 			colors = c('#8c510a', '#d8b365', '#f6e8c3', '#f5f5f5', '#c7eae5', '#5ab4ac', '#01665e'),
# 			limits = rescale_delta_limits
# 		) +
# 		ggtitle(title)

# 	title <- 'l) Integrated: Change\nPresence forced'
# 	legend_title <- 'Change'
# 	fut_delta_multi_force_presence <- ggplot() +
# 		layer_spatial(multi_force_presence, aes(fill = delta_sq_vs_fut), color = NA) +
# 		layer_spatial(nam, color = 'gray40', fill = NA, linewidth = 0.3) +
# 		layer_spatial(counties_with_ag, pch = 16, color = alpha('gray20', 0.3), size = 0.15) +
# 		layer_spatial(site_vect, pch = 3, size = 1) +
# 		scale_fill_gradientn(
# 			name = legend_title,
# 			colors = c('#8c510a', '#d8b365', '#f6e8c3', '#f5f5f5', '#c7eae5', '#5ab4ac', '#01665e'),
# 			limits = rescale_delta_limits
# 		) +
# 		ggtitle(title)

# 	maps <- list(
# 		sq_uni,
# 		fut_uni,
# 		fut_delta_uni,
# 		sq_multi,
# 		fut_multi,
# 		fut_delta_multi,
# 		sq_uni_force_presence,
# 		fut_uni_force_presence,
# 		fut_delta_uni_force_presence,
# 		sq_multi_force_presence,
# 		fut_multi_force_presence,
# 		fut_delta_multi_force_presence
# 	)

# 	for (i in seq_along(maps)) {

# 		maps[[i]] <- maps[[i]] +
# 			coord_sf(xlim = c(extent_coords[1], extent_coords[2]), ylim = c(extent_coords[3], extent_coords[4]), expand = FALSE) +
# 			theme(
# 				legend.position = 'none',
# 				plot.title = element_text(size = 12),
# 				axis.text = element_blank(),
# 				axis.ticks = element_blank(),
# 				panel.background = element_rect(fill = NA),
# 				plot.margin = margin(6, 0, 0, 0)
# 			)

# 	}

# 	maps_occs <- plot_grid(plotlist = maps, ncol = 3, align = 'hv')

# 	ggsave(maps_occs, filename = './outputs_loretta/integrated_sdm_pdm/maps_abundance_integrated_vs_non_integrated.png', width = 7.4, height = 10.5, dpi = 600, bg = 'white')

say('###########################################################')
say('### maps of biomass and change in biomass for main text ###')
say('###########################################################')

	### maps of biomass in 3-row grid, 1st two rows have 3 columns, last has two:
		# abund standalone SQ   abund standalne fut		abund standalone future delta
		# abund integrated SQ   abund integrated fut	abund integrated future delta
		# delta integrated and standalone sq	delta integrated and standalone sq

	nam <- vect('./data_from_gadm/gadm_4pt1_level_1_north_america_sans_alaska_lambert.gpkg')
	nam <- simplifyGeom(nam, tolerance = 1000)

	multi <- vect('./outputs_loretta/integrated_sdm_pdm/models_integrated/non_centered_MVN([occs~hurdlepoisson_offset]_[biomass~hurdleln]_[bla_can_hei])_eta=2/prediction_vector_nam.gpkg')
	multi_force_presence <- vect('./outputs_loretta/integrated_sdm_pdm/models_integrated/non_centered_MVN([occs~hurdlepoisson_offset]_[biomass~hurdleln]_[bla_can_hei])_eta=2/prediction_vector_nam_force_presence.gpkg')

	uni <- vect('./outputs_loretta/integrated_sdm_pdm/models_biomass/[biomass~hurdleln_log(bio12)]/prediction_vector_nam.gpkg')
	uni_force_presence <- vect('./outputs_loretta/integrated_sdm_pdm/models_biomass/[biomass~hurdleln_log(bio12)]/prediction_vector_nam_force_presence.gpkg')

	sites <- readRDS('./outputs_loretta/integrated_sdm_pdm/models_integrated/data_biomass_sites_bla_can_hei.rds')
	site_vect <- vect(sites$site_data_raw, geom = c('LONGITUDE', 'LATITUDE'), crs = getCRS('WGS84'))
	site_vect <- project(site_vect, nam)

	# extent
	extent <- ext(site_vect)
	extent <- as.polygons(extent, crs = nam)
	extent <- buffer(extent, width = 320 * 1000) # nominal plot extent
	extent <- ext(extent)
	extent_coords <- as.vector(extent)

	nam <- crop(nam, extent)
	multi <- crop(multi, extent)
	multi_force_presence <- crop(multi_force_presence, extent)
	uni <- crop(uni, extent)
	uni_force_presence <- crop(uni_force_presence, extent)

	# delineate current range "core"
	range_core_uni_sq <- delineate_range_core(uni, column = 'biomass_mean_sq', core_quant = core_quant)
	range_core_uni_fut <- delineate_range_core(uni, column = 'biomass_mean_ssp245_2071_2100', core_quant = core_quant)

	range_anticore_uni_sq <- delineate_range_core(uni, column = 'biomass_mean_sq', core_quant = anticore_quant, rule = '<')
	range_anticore_uni_fut <- delineate_range_core(uni, column = 'biomass_mean_ssp245_2071_2100', core_quant = anticore_quant, rule = '<')

	range_core_multi_sq <- delineate_range_core(multi, column = 'biomass_mean_sq', core_quant = core_quant)
	range_core_multi_fut <- delineate_range_core(multi, column = 'biomass_mean_ssp370_2071_2100', core_quant = core_quant)

	range_anticore_multi_sq <- delineate_range_core(multi, column = 'biomass_mean_sq', core_quant = anticore_quant, rule = '<')
	range_anticore_multi_fut <- delineate_range_core(multi, column = 'biomass_mean_ssp370_2071_2100', core_quant = anticore_quant, rule = '<')

	range_core_uni_sq_force_presence <- delineate_range_core(uni_force_presence, column = 'biomass_mean_sq', core_quant = core_quant)
	range_core_uni_fut_force_presence <- delineate_range_core(uni_force_presence, column = 'biomass_mean_ssp245_2071_2100', core_quant = core_quant)

	range_anticore_uni_sq_force_presence <- delineate_range_core(uni_force_presence, column = 'biomass_mean_sq', core_quant = anticore_quant, rule = '<')
	range_anticore_uni_fut_force_presence <- delineate_range_core(uni_force_presence, column = 'biomass_mean_ssp245_2071_2100', core_quant = anticore_quant, rule = '<')

	range_core_multi_sq_force_presence <- delineate_range_core(multi_force_presence, column = 'biomass_mean_sq', core_quant = core_quant)
	range_core_multi_fut_force_presence <- delineate_range_core(multi_force_presence, column = 'biomass_mean_ssp370_2071_2100', core_quant = core_quant)

	range_anticore_multi_sq_force_presence <- delineate_range_core(multi_force_presence, column = 'biomass_mean_sq', core_quant = anticore_quant, rule = '<')
	range_anticore_multi_fut_force_presence <- delineate_range_core(multi_force_presence, column = 'biomass_mean_ssp370_2071_2100', core_quant = anticore_quant, rule = '<')

	counties_with_ag <- uni[uni$n_andropogon_gerardi > 0]
	counties_with_ag <- centroids(counties_with_ag)

	# across both models, get max value of predictions for scaling
	vars <- c('biomass_mean_sq', 'biomass_mean_ssp245_2041_2070', 'biomass_mean_ssp245_2071_2100', 'biomass_mean_ssp370_2041_2070', 'biomass_mean_ssp370_2071_2100')
	# vars <- c('biomass_median_sq', 'biomass_median_ssp245_2041_2070', 'biomass_median_ssp245_2071_2100', 'biomass_median_ssp370_2041_2070', 'biomass_median_ssp370_2071_2100')

	max_val <- -Inf
	for (this_vect in c('multi', 'uni')) {
		x <- get(this_vect)
		for (var in vars) {
			x_var <- x[[var]]
			x_var <- unlist(x_var)
			max_val <- max(max_val, max(x_var, na.rm = TRUE))
			# max_val <- max(max_val, quantile(x_var, 0.99, na.rm = TRUE))
		}
	}
	for (this_vect in c('multi_force_presence', 'uni_force_presence')) {
		x <- get(this_vect)
		for (var in vars) {
			x_var <- x[[var]]
			x_var <- unlist(x_var)
			max_val <- max(max_val, max(x_var, na.rm = TRUE))
			# max_val <- max(max_val, quantile(x_var, 0.99, na.rm = TRUE))
		}
	}

	max_val <- max(max_val, sites$site_data_raw$mBiomass)
	resp_limits <- c(0, max_val)

	box_cox_trans <- function(p = 2/3) {
		trans_new(
			name = 'box_cox',
			transform = function(x) { x^p },
			inverse = function(x) { x^(1 / p) }
		)
	}

	title <- 'a) Biomass-only: Present\nAbsence allowed'
	legend_title <- 'Biomass (g)'
	sq_uni <- ggplot() +
		layer_spatial(uni, aes(fill = biomass_mean_sq), color = NA) +
		layer_spatial(nam, color = 'gray40', fill = NA, linewidth = 0.3) +
		layer_spatial(counties_with_ag, pch = 16, color = alpha('gray20', 0.3), size = 0.15) +
		layer_spatial(range_core_uni_sq, color = 'orange2', fill = NA, linewidth = 0.4) +
		layer_spatial(range_anticore_uni_sq, color = 'red', fill = NA, linewidth = 0.4) +
		layer_spatial(site_vect, aes(fill = mBiomass), pch = 21, size = 2) +
		scale_fill_gradientn(
			name = legend_title,
			colors = c('#edf8e9', '#74c476', '#005a32'),
			limits = resp_limits,
			transform = box_cox_trans(), breaks = scales::breaks_extended(n = 9),
		) +
		ggtitle(title)

	title <- 'b) Biomass-only: Future\nAbsence allowed'
	legend_title <- 'Biomass (g)'
	fut_uni <- ggplot() +
		layer_spatial(uni, aes(fill = biomass_mean_ssp370_2071_2100), color = NA) +
		layer_spatial(nam, color = 'gray40', fill = NA, linewidth = 0.3) +
		layer_spatial(counties_with_ag, pch = 16, color = alpha('gray20', 0.3), size = 0.15) +
		layer_spatial(range_core_uni_fut, color = 'orange2', fill = NA, linewidth = 0.4) +
		layer_spatial(range_anticore_uni_fut, color = 'red', fill = NA, linewidth = 0.4) +
		layer_spatial(site_vect, pch = 3, size = 1) +
		scale_fill_gradientn(
			name = legend_title,
			colors = c('#edf8e9', '#74c476', '#005a32'),
			limits = resp_limits,
			transform = box_cox_trans(), breaks = scales::breaks_extended(n = 9),
		) +
		ggtitle(title)

	title <- 'd) Integrated: Present\nAbsence allowed'
	legend_title <- 'Biomass (g)'
	sq_multi <- ggplot() +
		layer_spatial(multi, aes(fill = biomass_mean_sq), color = NA) +
		layer_spatial(nam, color = 'gray40', fill = NA, linewidth = 0.3) +
		layer_spatial(counties_with_ag, pch = 16, color = alpha('gray20', 0.3), size = 0.15) +
		layer_spatial(range_core_multi_sq, color = 'orange2', fill = NA, linewidth = 0.4) +
		layer_spatial(range_anticore_multi_sq, color = 'red', fill = NA, linewidth = 0.4) +
		layer_spatial(site_vect, aes(fill = mBiomass), pch = 21, size = 2) +
		scale_fill_gradientn(
			name = legend_title,
			colors = c('#edf8e9', '#74c476', '#005a32'),
			limits = resp_limits,
			transform = box_cox_trans(), breaks = scales::breaks_extended(n = 9),
		) +
		ggtitle(title)

	title <- 'e) Integrated: Future\nAbsence allowed'
	legend_title <- 'Biomass (g)'
	fut_multi <- ggplot() +
		layer_spatial(multi, aes(fill = biomass_mean_ssp370_2071_2100), color = NA) +
		layer_spatial(nam, color = 'gray40', fill = NA, linewidth = 0.3) +
		layer_spatial(counties_with_ag, pch = 16, color = alpha('gray20', 0.3), size = 0.15) +
		layer_spatial(range_core_multi_fut, color = 'orange2', fill = NA, linewidth = 0.4) +
		layer_spatial(range_anticore_multi_fut, color = 'red', fill = NA, linewidth = 0.4) +
		layer_spatial(site_vect, pch = 3, size = 1) +
		scale_fill_gradientn(
			name = legend_title,
			colors = c('#edf8e9', '#74c476', '#005a32'),
			limits = resp_limits,
			transform = box_cox_trans(), breaks = scales::breaks_extended(n = 9),
		) +
		ggtitle(title)

	title <- 'g) Biomass-only: Present\nPresence forced'
	legend_title <- 'Biomass (g)'
	sq_uni_force_presence <- ggplot() +
		layer_spatial(uni_force_presence, aes(fill = biomass_mean_sq), color = NA) +
		layer_spatial(nam, color = 'gray40', fill = NA, linewidth = 0.3) +
		layer_spatial(counties_with_ag, pch = 16, color = alpha('gray20', 0.3), size = 0.15) +
		layer_spatial(range_core_uni_sq_force_presence, color = 'orange2', fill = NA, linewidth = 0.4) +
		layer_spatial(range_anticore_uni_sq_force_presence, color = 'red', fill = NA, linewidth = 0.4) +
		layer_spatial(site_vect, aes(fill = mBiomass), pch = 21, size = 2) +
		scale_fill_gradientn(
			name = legend_title,
			colors = c('#edf8e9', '#74c476', '#005a32'),
			limits = resp_limits,
			transform = box_cox_trans(), breaks = scales::breaks_extended(n = 9),
		) +
		ggtitle(title)

	title <- 'h) Biomass-only: Future\nPresence forced'
	legend_title <- 'Biomass (g)'
	fut_uni_force_presence <- ggplot() +
		layer_spatial(uni_force_presence, aes(fill = biomass_mean_ssp370_2071_2100), color = NA) +
		layer_spatial(nam, color = 'gray40', fill = NA, linewidth = 0.3) +
		layer_spatial(counties_with_ag, pch = 16, color = alpha('gray20', 0.3), size = 0.15) +
		layer_spatial(range_core_uni_fut_force_presence, color = 'orange2', fill = NA, linewidth = 0.4) +
		layer_spatial(range_anticore_uni_fut_force_presence, color = 'red', fill = NA, linewidth = 0.4) +
		layer_spatial(site_vect, pch = 3, size = 1) +
		scale_fill_gradientn(
			name = legend_title,
			colors = c('#edf8e9', '#74c476', '#005a32'),
			limits = resp_limits,
			transform = box_cox_trans(), breaks = scales::breaks_extended(n = 9),
		) +
		ggtitle(title)

	title <- 'j) Integrated: Present\nPresence forced'
	legend_title <- 'Biomass (g)'
	sq_multi_force_presence <- ggplot() +
		layer_spatial(multi_force_presence, aes(fill = biomass_mean_sq), color = NA) +
		layer_spatial(nam, color = 'gray40', fill = NA, linewidth = 0.3) +
		layer_spatial(counties_with_ag, pch = 16, color = alpha('gray20', 0.3), size = 0.15) +
		layer_spatial(range_core_multi_sq_force_presence, color = 'orange2', fill = NA, linewidth = 0.4) +
		layer_spatial(range_anticore_multi_sq_force_presence, color = 'red', fill = NA, linewidth = 0.4) +
		layer_spatial(site_vect, aes(fill = mBiomass), pch = 21, size = 2) +
		scale_fill_gradientn(
			name = legend_title,
			colors = c('#edf8e9', '#74c476', '#005a32'),
			limits = resp_limits,
			transform = box_cox_trans(), breaks = scales::breaks_extended(n = 9),
		) +
		ggtitle(title)

	title <- 'k) Integrated: Future\nPresence forced'
	legend_title <- 'Biomass (g)'
	fut_multi_force_presence <- ggplot() +
		layer_spatial(multi_force_presence, aes(fill = biomass_mean_ssp370_2071_2100), color = NA) +
		layer_spatial(nam, color = 'gray40', fill = NA, linewidth = 0.3) +
		layer_spatial(counties_with_ag, pch = 16, color = alpha('gray20', 0.3), size = 0.15) +
		layer_spatial(range_core_multi_fut_force_presence, color = 'orange2', fill = NA, linewidth = 0.4) +
		layer_spatial(range_anticore_multi_fut_force_presence, color = 'red', fill = NA, linewidth = 0.4) +
		layer_spatial(site_vect, pch = 3, size = 1) +
		scale_fill_gradientn(
			name = legend_title,
			colors = c('#edf8e9', '#74c476', '#005a32'),
			limits = resp_limits,
			transform = box_cox_trans(), breaks = scales::breaks_extended(n = 9),
		) +
		ggtitle(title)

	uni$delta_sq_vs_fut <- uni$biomass_mean_ssp370_2071_2100 - uni$biomass_mean_sq
	multi$delta_sq_vs_fut <- multi$biomass_mean_ssp370_2071_2100 - multi$biomass_mean_sq

	uni_force_presence$delta_sq_vs_fut <- uni_force_presence$biomass_mean_ssp370_2071_2100 - uni_force_presence$biomass_mean_sq
	multi_force_presence$delta_sq_vs_fut <- multi_force_presence$biomass_mean_ssp370_2071_2100 - multi_force_presence$biomass_mean_sq

	delta_limits <- range(c(uni$delta_sq_vs_fut, uni_force_presence$delta_sq_vs_fut, multi$delta_sq_vs_fut, multi_force_presence$delta_sq_vs_fut))
	
	rescale_delta_limits <- c(-max(abs(delta_limits)), max(abs(delta_limits)))

	title <- 'c) Biomass-only: Change\nAbsence allowed'
	legend_title <- 'Change (g)'
	pm <- c(rep('', 10), rep('+', 11))
	fut_delta_uni <- ggplot() +
		layer_spatial(uni, aes(fill = delta_sq_vs_fut), color = NA) +
		layer_spatial(nam, color = 'gray40', fill = NA, linewidth = 0.3) +
		layer_spatial(counties_with_ag, pch = 16, color = alpha('gray20', 0.3), size = 0.15) +
		layer_spatial(site_vect, pch = 3, size = 1) +
		scale_fill_gradientn(
			name = legend_title,
			colors = c('#8c510a', '#d8b365', '#f6e8c3', '#f5f5f5', '#c7eae5', '#5ab4ac', '#01665e'),
			limits = rescale_delta_limits,
			breaks = scales::breaks_extended(n = 9),
			# transform = scales::asinh_trans()
		) +
		ggtitle(title)

	title <- 'f) Integrated: Change\nAbsence allowed'
	legend_title <- 'Change (g)'
	fut_delta_multi <- ggplot() +
		layer_spatial(multi, aes(fill = delta_sq_vs_fut), color = NA) +
		layer_spatial(nam, color = 'gray40', fill = NA, linewidth = 0.3) +
		layer_spatial(counties_with_ag, pch = 16, color = alpha('gray20', 0.3), size = 0.15) +
		layer_spatial(site_vect, pch = 3, size = 1) +
		scale_fill_gradientn(
			name = legend_title,
			colors = c('#8c510a', '#d8b365', '#f6e8c3', '#f5f5f5', '#c7eae5', '#5ab4ac', '#01665e'),
			limits = rescale_delta_limits,
			# transform = scales::asinh_trans()
		) +
		ggtitle(title)

	title <- 'i) Biomass-only: Change\nPresence forced'
	legend_title <- 'Change (g)'
	pm <- c(rep('', 10), rep('+', 11))
	fut_delta_uni_force_presence <- ggplot() +
		layer_spatial(uni_force_presence, aes(fill = delta_sq_vs_fut), color = NA) +
		layer_spatial(nam, color = 'gray40', fill = NA, linewidth = 0.3) +
		layer_spatial(counties_with_ag, pch = 16, color = alpha('gray20', 0.3), size = 0.15) +
		layer_spatial(site_vect, pch = 3, size = 1) +
		scale_fill_gradientn(
			name = legend_title,
			colors = c('#8c510a', '#d8b365', '#f6e8c3', '#f5f5f5', '#c7eae5', '#5ab4ac', '#01665e'),
			limits = rescale_delta_limits,
			# transform = scales::asinh_trans()
		) +
		ggtitle(title)

	title <- 'l) Integrated: Change\nPresence forced'
	legend_title <- 'Change (g)'
	fut_delta_multi_force_presence <- ggplot() +
		layer_spatial(multi_force_presence, aes(fill = delta_sq_vs_fut), color = NA) +
		layer_spatial(nam, color = 'gray40', fill = NA, linewidth = 0.3) +
		layer_spatial(counties_with_ag, pch = 16, color = alpha('gray20', 0.3), size = 0.15) +
		layer_spatial(site_vect, pch = 3, size = 1) +
		scale_fill_gradientn(
			name = legend_title,
			colors = c('#8c510a', '#d8b365', '#f6e8c3', '#f5f5f5', '#c7eae5', '#5ab4ac', '#01665e'),
			limits = rescale_delta_limits,
			transform = scales::asinh_trans()
		) +
		ggtitle(title)

	maps <- list(
		sq_uni,
		fut_uni,
		fut_delta_uni,
		sq_multi,
		fut_multi,
		fut_delta_multi,
		sq_uni_force_presence,
		fut_uni_force_presence,
		fut_delta_uni_force_presence,
		sq_multi_force_presence,
		fut_multi_force_presence,
		fut_delta_multi_force_presence
	)

	for (i in seq_along(maps)) {

		maps[[i]] <- maps[[i]] +
			coord_sf(xlim = c(extent_coords[1], extent_coords[2]), ylim = c(extent_coords[3], extent_coords[4]), expand = FALSE) +
			theme(
				legend.position = 'none',
				plot.title = element_text(size = 12),
				axis.text = element_blank(),
				axis.ticks = element_blank(),
				panel.background = element_rect(fill = NA),
				plot.margin = margin(6, 0, 0, 0)
			)

	}

	maps_composite <- plot_grid(plotlist = maps, ncol = 3, align = 'hv')
	ggsave(maps_composite, filename = './outputs_loretta/integrated_sdm_pdm/maps_biomass_integrated_vs_non_integrated.png', width = 7.4, height = 10.5, dpi = 600, bg = 'white')

# say('################################################################')
# say('### maps of abundance and change in abundance for supplement ###')
# say('################################################################')

# 	# ### layout
# 	# 	standalone sq (just one map)
# 	# 	fut1	delta
# 	# 	fut2	delta
# 	# 	fut3	delta
# 	# 	fut4	delta

# 	# 	(repeat on next page for integrated)

# 	# response variables
# 	var_sq <- 'N_ag_mean_sq'
# 	var_fut1 <- 'N_ag_mean_ssp245_2041_2070'
# 	var_fut2 <- 'N_ag_mean_ssp245_2071_2100'
# 	var_fut3 <- 'N_ag_mean_ssp370_2041_2070'
# 	var_fut4 <- 'N_ag_mean_ssp370_2071_2100'

# 	nam <- vect('./data_from_gadm/gadm_4pt1_level_1_north_america_sans_alaska_lambert.gpkg')
# 	nam <- simplifyGeom(nam, tolerance = 1000)

# 	multi <- vect('./outputs_loretta/integrated_sdm_pdm/models_integrated/non_centered_MVN([occs~hurdlepoisson_offset]_[biomass~hurdleln]_[bla_can_hei])_eta=2/prediction_vector_nam.gpkg')

# 	uni <- vect('./outputs_loretta/integrated_sdm_pdm/models_occurrence/[occs~overdispersed_hurdlepoisson(bio1^2_log(bio12)^2_bio15^2)]_[bias~poaceae_offset_plus_exp_neg_1_WITH_zero_correction_covariate]/prediction_vector_nam.gpkg')

# 	sites <- readRDS('./outputs_loretta/integrated_sdm_pdm/models_integrated/data_biomass_sites_bla_can_hei.rds')
# 	site_vect <- vect(sites$site_data_raw, geom = c('LONGITUDE', 'LATITUDE'), crs = getCRS('WGS84'))
# 	site_vect <- project(site_vect, nam)

# 	# extent
# 	extent <- ext(site_vect)
# 	extent <- as.polygons(extent, crs = nam)
# 	extent <- buffer(extent, width = 320 * 1000) # nominal plot extent
# 	extent <- ext(extent)
# 	extent_coords <- as.vector(extent)

# 	nam <- crop(nam, extent)
# 	multi <- crop(multi, extent)
# 	uni <- crop(uni, extent)

# 	# delineate current range "core"
# 	range_core_uni_sq <- delineate_range_core(uni, column = var_sq, core_quant = core_quant)
# 	range_core_multi_sq <- delineate_range_core(multi, column = var_sq, core_quant = core_quant)

# 	range_anticore_uni_sq <- delineate_range_core(uni, column = var_sq, core_quant = anticore_quant, rule = '<')
# 	range_anticore_multi_sq <- delineate_range_core(multi, column = var_sq, core_quant = anticore_quant, rule = '<')

# 	range_core_uni_fut <- range_core_multi_fut <- list()
# 	range_anticore_uni_fut <- range_anticore_multi_fut <- list()
# 	for (fut in 1:4) {
	
# 		range_core_uni_fut[[fut]] <- delineate_range_core(uni, column = get(paste0('var_fut', fut)), core_quant = core_quant)
# 		range_core_multi_fut[[fut]] <- delineate_range_core(multi, column = get(paste0('var_fut', fut)), core_quant = core_quant)

# 		range_anticore_uni_fut[[fut]] <- delineate_range_core(uni, column = get(paste0('var_fut', fut)), core_quant = anticore_quant, rule = '<')
# 		range_anticore_multi_fut[[fut]] <- delineate_range_core(multi, column = get(paste0('var_fut', fut)), core_quant = anticore_quant, rule = '<')

# 	}

# 	counties_with_ag <- uni[uni$n_andropogon_gerardi > 0]
# 	counties_with_ag <- centroids(counties_with_ag)

# 	# across both models, get max value of predictions for scaling
# 	vars <- c(var_sq, var_fut1, var_fut2, var_fut3, var_fut4)
# 	max_val <- -Inf
# 	for (this_vect in c('multi', 'uni')) {
# 		x <- get(this_vect)
# 		for (var in vars) {
# 			x_var <- x[[var]]
# 			x_var <- unlist(x_var)
# 			max_val <- max(max_val, max(x_var, na.rm = TRUE))
# 		}
# 	}

# 	resp_limits <- c(0, max_val)

# 	# across both models, get min/max value of change in predictions for scaling
# 	vars <- c(var_fut1, var_fut2, var_fut3, var_fut4)
# 	min_val <- Inf
# 	max_val <- -Inf
# 	for (this_vect in c('multi', 'uni')) {
# 		x <- get(this_vect)
# 		sq <- x[[var_sq]]
# 		sq <- unlist(sq)
# 		for (var in vars) {
# 			x_var <- x[[var]]
# 			x_var <- unlist(x_var)
# 			delta <- x_var - sq
# 			min_val <- min(min_val, min(delta, na.rm = TRUE))
# 			max_val <- max(max_val, max(delta, na.rm = TRUE))
# 		}
# 	}

# 	max_abs_val <- max(abs(c(min_val, max_val)))
# 	delta_resp_limits <- c(-max_abs_val, max_abs_val)

# 	title <- 'Abundance-only: Present'
# 	legend_title <- 'Abundance'
# 	sq_uni <- ggplot() +
# 		layer_spatial(uni, aes(fill = .data[[var_sq]]), color = NA) +
# 		layer_spatial(nam, color = 'gray40', fill = NA, linewidth = 0.3) +
# 		layer_spatial(counties_with_ag, pch = 16, color = alpha('gray20', 0.3), size = 0.15) +
# 		layer_spatial(range_core_uni_sq, color = 'orange2', fill = NA, linewidth = 0.4) +
# 		layer_spatial(range_anticore_uni_sq, color = 'red', fill = NA, linewidth = 0.4) +
# 		layer_spatial(site_vect, pch = 3, size = 1) +
# 		scale_fill_gradientn(
# 			name = legend_title,
# 			colors = c('#edf8e9', '#74c476', '#005a32'),
# 			limits = resp_limits
# 		) +
# 		ggtitle(title) +
# 		coord_sf(xlim = c(extent_coords[1], extent_coords[2]), ylim = c(extent_coords[3], extent_coords[4]), expand = FALSE) +
# 		theme(
# 			plot.title = element_text(size = 10),
# 			axis.text = element_blank(),
# 			axis.ticks = element_blank(),
# 			legend.title = element_text(size = 9),
# 			legend.text = element_text(size = 6),
# 			legend.key.width = grid::unit(0.24, 'cm'),
# 			legend.key.height = grid::unit(0.57, 'cm'),
# 			legend.background = element_rect(fill = alpha('white', 0))
# 		)

# 	title <- 'Integrated: Abundance-only: Present'
# 	legend_title <- 'Abundance'
# 	sq_multi <- ggplot() +
# 		layer_spatial(multi, aes(fill = .data[[var_sq]]), color = NA) +
# 		layer_spatial(nam, color = 'gray40', fill = NA, linewidth = 0.3) +
# 		layer_spatial(counties_with_ag, pch = 16, color = alpha('gray20', 0.3), size = 0.15) +
# 		layer_spatial(range_core_multi_sq, color = 'orange2', fill = NA, linewidth = 0.4) +
# 		layer_spatial(range_anticore_multi_sq, color = 'red', fill = NA, linewidth = 0.4) +
# 		layer_spatial(site_vect, pch = 3, size = 1) +
# 		scale_fill_gradientn(
# 			name = legend_title,
# 			colors = c('#edf8e9', '#74c476', '#005a32'),
# 			limits = resp_limits
# 		) +
# 		ggtitle(title) +
# 		coord_sf(xlim = c(extent_coords[1], extent_coords[2]), ylim = c(extent_coords[3], extent_coords[4]), expand = FALSE) +
# 		theme(
# 			plot.title = element_text(size = 10),
# 			axis.text = element_blank(),
# 			axis.ticks = element_blank(),
# 			legend.title = element_text(size = 9),
# 			legend.text = element_text(size = 6),
# 			legend.key.width = grid::unit(0.24, 'cm'),
# 			legend.key.height = grid::unit(0.57, 'cm'),
# 			legend.background = element_rect(fill = alpha('white', 0))
# 		)

# 	### standalone models: future predictions and deltas
# 	####################################################

# 	futs <- deltas <- list()
# 	for (fut in 1:4) {

# 		if (fut == 1) {
# 			var_fut <- var_fut1
# 			title_fut <- 'Abundance-only: SSP245 2041-2070'
# 			title_delta <- 'Change: 2041-2070 - Present'
# 		} else if (fut == 2) {
# 			var_fut <- var_fut2
# 			title_fut <- 'Abundance-only: SSP245 2071-2100'
# 			title_delta <- 'Change: 2071-2100 - Present'
# 		} else if (fut == 3) {
# 			var_fut <- var_fut3
# 			title_fut <- 'Abundance-only: SSP370 2041-2070'
# 			title_delta <- 'Change: 2041-2070 - Present'
# 		} else if (fut == 4) {
# 			var_fut <- var_fut4
# 			title_fut <- 'Abundance-only: SSP370 2071-2100'
# 			title_delta <- 'Change: 2071-2100 - Present'
# 		}

# 		legend_title <- 'Abund.'
# 		futs[[fut]] <- ggplot() +
# 			layer_spatial(uni, aes(fill = .data[[var_fut]]), color = NA) +
# 			layer_spatial(nam, color = 'gray40', fill = NA, linewidth = 0.3) +
# 			layer_spatial(counties_with_ag, pch = 16, color = alpha('gray20', 0.3), size = 0.15) +
# 			layer_spatial(range_core_uni_fut[[fut]], color = 'orange2', fill = NA, linewidth = 0.4) +
# 			layer_spatial(range_anticore_uni_fut[[fut]], color = 'red', fill = NA, linewidth = 0.4) +
# 			layer_spatial(site_vect, pch = 3, size = 1) +
# 			scale_fill_gradientn(
# 				name = legend_title,
# 				colors = c('#edf8e9', '#74c476', '#005a32'),
# 				limits = resp_limits
# 			) +
# 			ggtitle(title_fut) +
# 			coord_sf(xlim = c(extent_coords[1], extent_coords[2]), ylim = c(extent_coords[3], extent_coords[4]), expand = FALSE) +
# 			theme(
# 				plot.title = element_text(size = 10),
# 				axis.text = element_blank(),
# 				axis.ticks = element_blank(),
# 				legend.title = element_text(size = 9),
# 				legend.text = element_text(size = 6),
# 				legend.key.width = grid::unit(0.24, 'cm'),
# 				legend.key.height = grid::unit(0.57, 'cm'),
# 				legend.background = element_rect(fill = alpha('white', 0))
# 			)

# 		legend_title <- 'Change'
# 		deltas[[fut]] <- ggplot() +
# 			layer_spatial(uni, aes(fill = .data[[var_fut]] - .data[[var_sq]]), color = NA) +
# 			layer_spatial(nam, color = 'gray40', fill = NA, linewidth = 0.3) +
# 			layer_spatial(counties_with_ag, pch = 16, color = alpha('gray20', 0.3), size = 0.15) +
# 			layer_spatial(range_core_uni_fut[[fut]], color = 'orange2', fill = NA, linewidth = 0.4) +
# 			layer_spatial(range_anticore_uni_fut[[fut]], color = 'red', fill = NA, linewidth = 0.4) +
# 			layer_spatial(site_vect, pch = 3, size = 1) +
# 			scale_fill_gradientn(
# 				name = legend_title,
# 				colors = c('#8c510a', '#d8b365', '#f6e8c3', '#f5f5f5', '#c7eae5', '#5ab4ac', '#01665e'),
# 				limits = delta_resp_limits
# 			) +
# 			ggtitle(title_delta) +
# 			coord_sf(xlim = c(extent_coords[1], extent_coords[2]), ylim = c(extent_coords[3], extent_coords[4]), expand = FALSE) +
# 			theme(
# 				plot.title = element_text(size = 10),
# 				axis.text = element_blank(),
# 				axis.ticks = element_blank(),
# 				legend.title = element_text(size = 9),
# 				legend.text = element_text(size = 6),
# 				legend.key.width = grid::unit(0.24, 'cm'),
# 				legend.key.height = grid::unit(0.57, 'cm'),
# 				legend.background = element_rect(fill = alpha('white', 0))
# 			)

# 	}

# 	maps_futs <- plot_grid(plotlist = futs, ncol = 1, align = 'hv')
# 	maps_deltas <- plot_grid(plotlist = deltas, ncol = 1, align = 'hv')
# 	maps_futs_deltas <- plot_grid(maps_futs, maps_deltas, ncol = 2, align = 'hv')
# 	maps <- plot_grid(sq_uni, maps_futs_deltas, ncol = 1, align = 'hv', rel_heights = c(1.2, 4))

# 	ggsave(maps, filename = './outputs_loretta/integrated_sdm_pdm/maps_abundance_all_maps_non_integrated.png', width = 6, height = 9, dpi = 600, bg = 'white')

# 	### integrated models: future predictions and deltas
# 	####################################################

# 	futs <- deltas <- list()
# 	for (fut in 1:4) {

# 		if (fut == 1) {
# 			var_fut <- var_fut1
# 			title_fut <- 'Integrated: SSP245 2041-2070'
# 			title_delta <- 'Change: 2041-2070 - Present'
# 		} else if (fut == 2) {
# 			var_fut <- var_fut2
# 			title_fut <- 'Integrated: SSP245 2071-2100'
# 			title_delta <- 'Change: 2071-2100 - Present'
# 		} else if (fut == 3) {
# 			var_fut <- var_fut3
# 			title_fut <- 'Integrated: SSP370 2041-2070'
# 			title_delta <- 'Change: 2041-2070 - Present'
# 		} else if (fut == 4) {
# 			var_fut <- var_fut4
# 			title_fut <- 'Integrated: SSP370 2071-2100'
# 			title_delta <- 'Change: 2071-2100 - Present'
# 		}

# 		legend_title <- 'Abund.'
# 		futs[[fut]] <- ggplot() +
# 			layer_spatial(multi, aes(fill = .data[[var_fut]]), color = NA) +
# 			layer_spatial(nam, color = 'gray40', fill = NA, linewidth = 0.3) +
# 			layer_spatial(counties_with_ag, pch = 16, color = alpha('gray20', 0.3), size = 0.15) +
# 			layer_spatial(range_core_multi_fut[[fut]], color = 'orange2', fill = NA, linewidth = 0.4) +
# 			layer_spatial(range_anticore_multi_fut[[fut]], color = 'red', fill = NA, linewidth = 0.4) +
# 			layer_spatial(site_vect, pch = 3, size = 1) +
# 			scale_fill_gradientn(
# 				name = legend_title,
# 				colors = c('#edf8e9', '#74c476', '#005a32'),
# 				limits = resp_limits
# 			) +
# 			ggtitle(title_fut) +
# 			coord_sf(xlim = c(extent_coords[1], extent_coords[2]), ylim = c(extent_coords[3], extent_coords[4]), expand = FALSE) +
# 			theme(
# 				plot.title = element_text(size = 10),
# 				axis.text = element_blank(),
# 				axis.ticks = element_blank(),
# 				legend.title = element_text(size = 9),
# 				legend.text = element_text(size = 6),
# 				legend.key.width = grid::unit(0.24, 'cm'),
# 				legend.key.height = grid::unit(0.57, 'cm'),
# 				legend.background = element_rect(fill = alpha('white', 0))
# 			)

# 		legend_title <- 'Change'
# 		deltas[[fut]] <- ggplot() +
# 			layer_spatial(multi, aes(fill = .data[[var_fut]] - .data[[var_sq]]), color = NA) +
# 			layer_spatial(nam, color = 'gray40', fill = NA, linewidth = 0.3) +
# 			layer_spatial(counties_with_ag, pch = 16, color = alpha('gray20', 0.3), size = 0.15) +
# 			layer_spatial(range_core_multi_fut[[fut]], color = 'orange2', fill = NA, linewidth = 0.4) +
# 			layer_spatial(range_anticore_multi_fut[[fut]], color = 'red', fill = NA, linewidth = 0.4) +
# 			layer_spatial(site_vect, pch = 3, size = 1) +
# 			scale_fill_gradientn(
# 				name = legend_title,
# 				colors = c('#8c510a', '#d8b365', '#f6e8c3', '#f5f5f5', '#c7eae5', '#5ab4ac', '#01665e'),
# 				limits = delta_resp_limits
# 			) +
# 			ggtitle(title_delta) +
# 			coord_sf(xlim = c(extent_coords[1], extent_coords[2]), ylim = c(extent_coords[3], extent_coords[4]), expand = FALSE) +
# 			theme(
# 				plot.title = element_text(size = 10),
# 				axis.text = element_blank(),
# 				axis.ticks = element_blank(),
# 				legend.title = element_text(size = 9),
# 				legend.text = element_text(size = 6),
# 				legend.key.width = grid::unit(0.24, 'cm'),
# 				legend.key.height = grid::unit(0.57, 'cm'),
# 				legend.background = element_rect(fill = alpha('white', 0))
# 			)

# 	}

# 	maps_futs <- plot_grid(plotlist = futs, ncol = 1, align = 'hv')
# 	maps_deltas <- plot_grid(plotlist = deltas, ncol = 1, align = 'hv')
# 	maps_futs_deltas <- plot_grid(maps_futs, maps_deltas, ncol = 2, align = 'hv')
# 	maps <- plot_grid(sq_multi, maps_futs_deltas, ncol = 1, align = 'hv', rel_heights = c(1.2, 4))

# 	ggsave(maps, filename = './outputs_loretta/integrated_sdm_pdm/maps_abundance_all_maps_integrated.png', width = 6, height = 9, dpi = 600, bg = 'white')

# say('############################################################')
# say('### maps of biomass and change in biomass for supplement ###')
# say('############################################################')

# 	# ### layout
# 	# 	standalone sq (just one map)
# 	# 	fut1	delta
# 	# 	fut2	delta
# 	# 	fut3	delta
# 	# 	fut4	delta

# 	# 	(repeat on next page for integrated)

# 	# response variables
# 	var_sq <- 'biomass_mean_sq'
# 	var_fut1 <- 'biomass_mean_ssp245_2041_2070'
# 	var_fut2 <- 'biomass_mean_ssp245_2071_2100'
# 	var_fut3 <- 'biomass_mean_ssp370_2041_2070'
# 	var_fut4 <- 'biomass_mean_ssp370_2071_2100'

# 	nam <- vect('./data_from_gadm/gadm_4pt1_level_1_north_america_sans_alaska_lambert.gpkg')
# 	nam <- simplifyGeom(nam, tolerance = 1000)

# 	multi <- vect('./outputs_loretta/integrated_sdm_pdm/models_integrated/non_centered_MVN([occs~hurdlepoisson_offset]_[biomass~hurdleln]_[bla_can_hei])_eta=2/prediction_vector_nam_force_presence.gpkg')

# 	uni <- vect('./outputs_loretta/integrated_sdm_pdm/models_biomass/[biomass~hurdleln_log(bio12)]/prediction_vector_nam.gpkg')

# 	sites <- readRDS('./outputs_loretta/integrated_sdm_pdm/models_integrated/data_biomass_sites_bla_can_hei.rds')
# 	site_vect <- vect(sites$site_data_raw, geom = c('LONGITUDE', 'LATITUDE'), crs = getCRS('WGS84'))
# 	site_vect <- project(site_vect, nam)

# 	# extent
# 	extent <- ext(site_vect)
# 	extent <- as.polygons(extent, crs = nam)
# 	extent <- buffer(extent, width = 320 * 1000) # nominal plot extent
# 	extent <- ext(extent)
# 	extent_coords <- as.vector(extent)

# 	nam <- crop(nam, extent)
# 	multi <- crop(multi, extent)
# 	uni <- crop(uni, extent)

# 	# delineate current range "core"
# 	range_core_uni_sq <- delineate_range_core(uni, column = var_sq, core_quant = core_quant)
# 	range_core_multi_sq <- delineate_range_core(multi, column = var_sq, core_quant = core_quant)

# 	range_anticore_uni_sq <- delineate_range_core(uni, column = var_sq, core_quant = anticore_quant, rule = '<')
# 	range_anticore_multi_sq <- delineate_range_core(multi, column = var_sq, core_quant = anticore_quant, rule = '<')

# 	range_core_uni_fut <- range_core_multi_fut <- list()
# 	range_anticore_uni_fut <- range_anticore_multi_fut <- list()
# 	for (fut in 1:4) {
	
# 		range_core_uni_fut[[fut]] <- delineate_range_core(uni, column = get(paste0('var_fut', fut)), core_quant = core_quant)
# 		range_core_multi_fut[[fut]] <- delineate_range_core(multi, column = get(paste0('var_fut', fut)), core_quant = core_quant)

# 		range_anticore_uni_fut[[fut]] <- delineate_range_core(uni, column = get(paste0('var_fut', fut)), core_quant = anticore_quant, rule = '<')
# 		range_anticore_multi_fut[[fut]] <- delineate_range_core(multi, column = get(paste0('var_fut', fut)), core_quant = anticore_quant, rule = '<')

# 	}

# 	counties_with_ag <- uni[uni$n_andropogon_gerardi > 0]
# 	counties_with_ag <- centroids(counties_with_ag)

# 	# across both models, get max value of predictions for scaling
# 	vars <- c(var_sq, var_fut1, var_fut2, var_fut3, var_fut4)
# 	max_val <- -Inf
# 	for (this_vect in c('multi', 'uni')) {
# 		x <- get(this_vect)
# 		for (var in vars) {
# 			x_var <- x[[var]]
# 			x_var <- unlist(x_var)
# 			max_val <- max(max_val, max(x_var, na.rm = TRUE))
# 		}
# 	}

# 	resp_limits <- c(0, max_val)

# 	# across both models, get min/max value of change in predictions for scaling
# 	vars <- c(var_fut1, var_fut2, var_fut3, var_fut4)
# 	min_val <- Inf
# 	max_val <- -Inf
# 	for (this_vect in c('multi', 'uni')) {
# 		x <- get(this_vect)
# 		sq <- x[[var_sq]]
# 		sq <- unlist(sq)
# 		for (var in vars) {
# 			x_var <- x[[var]]
# 			x_var <- unlist(x_var)
# 			delta <- x_var - sq
# 			min_val <- min(min_val, min(delta, na.rm = TRUE))
# 			max_val <- max(max_val, max(delta, na.rm = TRUE))
# 		}
# 	}

# 	max_abs_val <- max(abs(c(min_val, max_val)))
# 	delta_resp_limits <- c(-max_abs_val, max_abs_val)

# 	title <- 'Biomass-only: Present'
# 	legend_title <- 'Biomass (g)'
# 	sq_uni <- ggplot() +
# 		layer_spatial(uni, aes(fill = .data[[var_sq]]), color = NA) +
# 		layer_spatial(nam, color = 'gray40', fill = NA, linewidth = 0.3) +
# 		layer_spatial(counties_with_ag, pch = 16, color = alpha('gray20', 0.3), size = 0.15) +
# 		layer_spatial(range_core_uni_sq, color = 'orange2', fill = NA, linewidth = 0.4) +
# 		layer_spatial(range_anticore_uni_sq, color = 'red', fill = NA, linewidth = 0.4) +
# 		layer_spatial(site_vect, pch = 3, size = 1) +
# 		scale_fill_gradientn(
# 			name = legend_title,
# 			colors = c('#edf8e9', '#74c476', '#005a32'),
# 			limits = resp_limits
# 		) +
# 		ggtitle(title) +
# 		coord_sf(xlim = c(extent_coords[1], extent_coords[2]), ylim = c(extent_coords[3], extent_coords[4]), expand = FALSE) +
# 		theme(
# 			plot.title = element_text(size = 10),
# 			axis.text = element_blank(),
# 			axis.ticks = element_blank(),
# 			legend.title = element_text(size = 9),
# 			legend.text = element_text(size = 6),
# 			legend.key.width = grid::unit(0.24, 'cm'),
# 			legend.key.height = grid::unit(0.57, 'cm'),
# 			legend.background = element_rect(fill = alpha('white', 0))
# 		)

# 	title <- 'Integrated: Present'
# 	legend_title <- 'Biomass (g)'
# 	sq_multi <- ggplot() +
# 		layer_spatial(multi, aes(fill = .data[[var_sq]]), color = NA) +
# 		layer_spatial(nam, color = 'gray40', fill = NA, linewidth = 0.3) +
# 		layer_spatial(counties_with_ag, pch = 16, color = alpha('gray20', 0.3), size = 0.15) +
# 		layer_spatial(range_core_multi_sq, color = 'orange2', fill = NA, linewidth = 0.4) +
# 		layer_spatial(range_anticore_multi_sq, color = 'red', fill = NA, linewidth = 0.4) +
# 		layer_spatial(site_vect, pch = 3, size = 1) +
# 		scale_fill_gradientn(
# 			name = legend_title,
# 			colors = c('#edf8e9', '#74c476', '#005a32'),
# 			limits = resp_limits
# 		) +
# 		ggtitle(title) +
# 		coord_sf(xlim = c(extent_coords[1], extent_coords[2]), ylim = c(extent_coords[3], extent_coords[4]), expand = FALSE) +
# 		theme(
# 			plot.title = element_text(size = 10),
# 			axis.text = element_blank(),
# 			axis.ticks = element_blank(),
# 			legend.title = element_text(size = 9),
# 			legend.text = element_text(size = 6),
# 			legend.key.width = grid::unit(0.24, 'cm'),
# 			legend.key.height = grid::unit(0.57, 'cm'),
# 			legend.background = element_rect(fill = alpha('white', 0))
# 		)

# 	### standalone models: future predictions and deltas
# 	####################################################

# 	futs <- deltas <- list()
# 	for (fut in 1:4) {

# 		if (fut == 1) {
# 			var_fut <- var_fut1
# 			title_fut <- 'Biomass-only: SSP245 2041-2070'
# 			title_delta <- 'Change: 2041-2070 - Present'
# 		} else if (fut == 2) {
# 			var_fut <- var_fut2
# 			title_fut <- 'Biomass-only: SSP245 2071-2100'
# 			title_delta <- 'Change: 2071-2100 - Present'
# 		} else if (fut == 3) {
# 			var_fut <- var_fut3
# 			title_fut <- 'Biomass-only: SSP370 2041-2070'
# 			title_delta <- 'Change: 2041-2070 - Present'
# 		} else if (fut == 4) {
# 			var_fut <- var_fut4
# 			title_fut <- 'Biomass-only: SSP370 2071-2100'
# 			title_delta <- 'Change: 2071-2100 - Present'
# 		}

# 		legend_title <- 'Biomass (g)'
# 		futs[[fut]] <- ggplot() +
# 			layer_spatial(uni, aes(fill = .data[[var_fut]]), color = NA) +
# 			layer_spatial(nam, color = 'gray40', fill = NA, linewidth = 0.3) +
# 			layer_spatial(counties_with_ag, pch = 16, color = alpha('gray20', 0.3), size = 0.15) +
# 			layer_spatial(range_core_uni_fut[[fut]], color = 'orange2', fill = NA, linewidth = 0.4) +
# 			layer_spatial(range_anticore_uni_fut[[fut]], color = 'red', fill = NA, linewidth = 0.4) +
# 			layer_spatial(site_vect, pch = 3, size = 1) +
# 			scale_fill_gradientn(
# 				name = legend_title,
# 				colors = c('#edf8e9', '#74c476', '#005a32'),
# 				limits = resp_limits
# 			) +
# 			ggtitle(title_fut) +
# 			coord_sf(xlim = c(extent_coords[1], extent_coords[2]), ylim = c(extent_coords[3], extent_coords[4]), expand = FALSE) +
# 			theme(
# 				plot.title = element_text(size = 10),
# 				axis.text = element_blank(),
# 				axis.ticks = element_blank(),
# 				legend.title = element_text(size = 9),
# 				legend.text = element_text(size = 6),
# 				legend.key.width = grid::unit(0.24, 'cm'),
# 				legend.key.height = grid::unit(0.57, 'cm'),
# 				legend.background = element_rect(fill = alpha('white', 0))
# 			)

# 		legend_title <- 'Change'
# 		deltas[[fut]] <- ggplot() +
# 			layer_spatial(uni, aes(fill = .data[[var_fut]] - .data[[var_sq]]), color = NA) +
# 			layer_spatial(nam, color = 'gray40', fill = NA, linewidth = 0.3) +
# 			layer_spatial(counties_with_ag, pch = 16, color = alpha('gray20', 0.3), size = 0.15) +
# 			layer_spatial(range_core_uni_fut[[fut]], color = 'orange2', fill = NA, linewidth = 0.4) +
# 			layer_spatial(range_anticore_uni_fut[[fut]], color = 'red', fill = NA, linewidth = 0.4) +
# 			layer_spatial(site_vect, pch = 3, size = 1) +
# 			scale_fill_gradientn(
# 				name = legend_title,
# 				colors = c('#8c510a', '#d8b365', '#f6e8c3', '#f5f5f5', '#c7eae5', '#5ab4ac', '#01665e'),
# 				limits = delta_resp_limits
# 			) +
# 			ggtitle(title_delta) +
# 			coord_sf(xlim = c(extent_coords[1], extent_coords[2]), ylim = c(extent_coords[3], extent_coords[4]), expand = FALSE) +
# 			theme(
# 				plot.title = element_text(size = 10),
# 				axis.text = element_blank(),
# 				axis.ticks = element_blank(),
# 				legend.title = element_text(size = 9),
# 				legend.text = element_text(size = 6),
# 				legend.key.width = grid::unit(0.24, 'cm'),
# 				legend.key.height = grid::unit(0.57, 'cm'),
# 				legend.background = element_rect(fill = alpha('white', 0))
# 			)

# 	}

# 	maps_futs <- plot_grid(plotlist = futs, ncol = 1, align = 'hv')
# 	maps_deltas <- plot_grid(plotlist = deltas, ncol = 1, align = 'hv')
# 	maps_futs_deltas <- plot_grid(maps_futs, maps_deltas, ncol = 2, align = 'hv')
# 	maps <- plot_grid(sq_uni, maps_futs_deltas, ncol = 1, align = 'hv', rel_heights = c(1.2, 4))

# 	ggsave(maps, filename = './outputs_loretta/integrated_sdm_pdm/maps_biomass_all_maps_non_integrated.png', width = 6, height = 9, dpi = 600, bg = 'white')

# 	### integrated models: future predictions and deltas
# 	####################################################

# 	futs <- deltas <- list()
# 	for (fut in 1:4) {

# 		if (fut == 1) {
# 			var_fut <- var_fut1
# 			title_fut <- 'Integrated: SSP245 2041-2070'
# 			title_delta <- 'Change: 2041-2070 - Present'
# 		} else if (fut == 2) {
# 			var_fut <- var_fut2
# 			title_fut <- 'Integrated: SSP245 2071-2100'
# 			title_delta <- 'Change: 2071-2100 - Present'
# 		} else if (fut == 3) {
# 			var_fut <- var_fut3
# 			title_fut <- 'Integrated: SSP370 2041-2070'
# 			title_delta <- 'Change: 2041-2070 - Present'
# 		} else if (fut == 4) {
# 			var_fut <- var_fut4
# 			title_fut <- 'Integrated: SSP370 2071-2100'
# 			title_delta <- 'Change: 2071-2100 - Present'
# 		}

# 		legend_title <- 'Biomass (g)'
# 		futs[[fut]] <- ggplot() +
# 			layer_spatial(multi, aes(fill = .data[[var_fut]]), color = NA) +
# 			layer_spatial(nam, color = 'gray40', fill = NA, linewidth = 0.3) +
# 			layer_spatial(counties_with_ag, pch = 16, color = alpha('gray20', 0.3), size = 0.15) +
# 			layer_spatial(range_core_multi_fut[[fut]], color = 'orange2', fill = NA, linewidth = 0.4) +
# 			layer_spatial(range_anticore_multi_fut[[fut]], color = 'red', fill = NA, linewidth = 0.4) +
# 			layer_spatial(site_vect, pch = 3, size = 1) +
# 			scale_fill_gradientn(
# 				name = legend_title,
# 				colors = c('#edf8e9', '#74c476', '#005a32'),
# 				limits = resp_limits
# 			) +
# 			ggtitle(title_fut) +
# 			coord_sf(xlim = c(extent_coords[1], extent_coords[2]), ylim = c(extent_coords[3], extent_coords[4]), expand = FALSE) +
# 			theme(
# 				plot.title = element_text(size = 10),
# 				axis.text = element_blank(),
# 				axis.ticks = element_blank(),
# 				legend.title = element_text(size = 9),
# 				legend.text = element_text(size = 6),
# 				legend.key.width = grid::unit(0.24, 'cm'),
# 				legend.key.height = grid::unit(0.57, 'cm'),
# 				legend.background = element_rect(fill = alpha('white', 0))
# 			)

# 		legend_title <- 'Change'
# 		deltas[[fut]] <- ggplot() +
# 			layer_spatial(multi, aes(fill = .data[[var_fut]] - .data[[var_sq]]), color = NA) +
# 			layer_spatial(nam, color = 'gray40', fill = NA, linewidth = 0.3) +
# 			layer_spatial(counties_with_ag, pch = 16, color = alpha('gray20', 0.3), size = 0.15) +
# 			layer_spatial(range_core_multi_fut[[fut]], color = 'orange2', fill = NA, linewidth = 0.4) +
# 			layer_spatial(range_anticore_multi_fut[[fut]], color = 'red', fill = NA, linewidth = 0.4) +
# 			layer_spatial(site_vect, pch = 3, size = 1) +
# 			scale_fill_gradientn(
# 				name = legend_title,
# 				colors = c('#8c510a', '#d8b365', '#f6e8c3', '#f5f5f5', '#c7eae5', '#5ab4ac', '#01665e'),
# 				limits = delta_resp_limits
# 			) +
# 			ggtitle(title_delta) +
# 			coord_sf(xlim = c(extent_coords[1], extent_coords[2]), ylim = c(extent_coords[3], extent_coords[4]), expand = FALSE) +
# 			theme(
# 				plot.title = element_text(size = 10),
# 				axis.text = element_blank(),
# 				axis.ticks = element_blank(),
# 				legend.title = element_text(size = 9),
# 				legend.text = element_text(size = 6),
# 				legend.key.width = grid::unit(0.24, 'cm'),
# 				legend.key.height = grid::unit(0.57, 'cm'),
# 				legend.background = element_rect(fill = alpha('white', 0))
# 			)

# 	}

# 	maps_futs <- plot_grid(plotlist = futs, ncol = 1, align = 'hv')
# 	maps_deltas <- plot_grid(plotlist = deltas, ncol = 1, align = 'hv')
# 	maps_futs_deltas <- plot_grid(maps_futs, maps_deltas, ncol = 2, align = 'hv')
# 	maps <- plot_grid(sq_multi, maps_futs_deltas, ncol = 1, align = 'hv', rel_heights = c(1.2, 4))

# 	ggsave(maps, filename = './outputs_loretta/integrated_sdm_pdm/maps_biomass_all_maps_integrated.png', width = 6, height = 9, dpi = 600, bg = 'white')

# say('##################################################################################')
# say('### maps of non-biomass traits and change in non-biomass traits for supplement ###')
# say('##################################################################################')

# 	# ### layout
# 	# 	standalone sq (just one map)
# 	# 	fut1	delta
# 	# 	fut2	delta
# 	# 	fut3	delta
# 	# 	fut4	delta

# 	# 	(repeat on next page for integrated)

# 	# best models
# 	models <- read_xlsx('./outputs_loretta/integrated_sdm_pdm/summary_of_ALL_top_models_manual_assessment.xlsx', sheet = 'summary_of_ALL_top_models')
# 	models <- as.data.table(models)
# 	models <- models[models$selected]

# 	# states/provinces
# 	nam <- vect('./data_from_gadm/gadm_4pt1_level_1_north_america_sans_alaska_lambert.gpkg')
# 	nam <- simplifyGeom(nam, tolerance = 1000)

# 	# multivariate models
# 	multi_morph <- vect('./outputs_loretta/integrated_sdm_pdm/models_integrated/non_centered_MVN([occs~hurdlepoisson_offset]_[biomass~hurdleln]_[bla_can_hei])_eta=2/prediction_vector_nam_force_presence.gpkg')
	
# 	multi_phys <- vect('./outputs_loretta/integrated_sdm_pdm/models_integrated/non_centered_MVN([occs~hurdlepoisson_offset]_[biomass~hurdleln]_[cnr_pho_spa_sto_tra])_eta=2/prediction_vector_nam_force_presence.gpkg')

# 	counties_with_ag <- multi_morph[multi_morph$n_andropogon_gerardi > 0]
# 	counties_with_ag <- centroids(counties_with_ag)

# 	# sample sites
# 	sites <- readRDS('./outputs_loretta/integrated_sdm_pdm/models_integrated/data_biomass_sites_bla_can_hei.rds')
# 	site_vect <- vect(sites$site_data_raw, geom = c('LONGITUDE', 'LATITUDE'), crs = getCRS('WGS84'))
# 	site_vect <- project(site_vect, nam)

# 	# extent
# 	extent <- ext(site_vect)
# 	extent <- as.polygons(extent, crs = nam)
# 	extent <- buffer(extent, width = 320 * 1000) # nominal plot extent
# 	extent <- ext(extent)
# 	extent_coords <- as.vector(extent)

# 	# corop to display area
# 	nam <- crop(nam, extent)
# 	multi_morph <- crop(multi_morph, extent)
# 	multi_phys <- crop(multi_phys, extent)

# 	facets <- sort(unique(models$facet))
# 	facets <- facets[facets != 'biomass']

# 	identity_trans <- function() {
# 		trans_new(
# 			'identity',
# 			transform = function(x) x,
# 			inverse   = function(x) x,
# 			domain    = c(-Inf, Inf)
# 		)
# 	}

# 	log10p1_trans <- function() {
# 		trans_new(
# 			'log10p1',
# 			transform = function(x) log10(x + 1),
# 			inverse   = function(x) 10^x - 1,
# 			domain    = c(0, Inf)
# 		)
# 	}

# 	log2p1_trans <- function() {
# 		trans_new(
# 			'log2p1',
# 			transform = function(x) log2(x + 1),
# 			inverse   = function(x) 2^x - 1,
# 			domain    = c(0, Inf)
# 		)
# 	}

# 	asinh_trans <- function() {
# 		trans_new(
# 			'asinh_t',
# 			transform = function(x) asinh(x),
# 			inverse   = function(x) sinh(x),
# 			domain    = c(0, Inf)
# 		)
# 	}

# 	legend_x <- 1.5
# 	legend_x_delta <- 1.4
# 	legend_y <- 0
# 	legend_title_size <- 8

# 	for (f in seq_along(facets)) {

# 		facet <- facets[f]
# 		say(facet)
# 		nice <- get_nice_trait(facet)

# 		type <- nonbiomass_facets[[facet]]$type

# 		if (type == 'morphological') {
# 			multi <- multi_morph
# 		} else {
# 			multi <- multi_phys
# 		}

# 		if (facet == 'blade_width') {
# 			site_mean_var <- 'mBladeWidth'
# 			trans <- asinh_trans()
# 		} else if (facet == 'canopy_diameter') {
# 			site_mean_var <- 'mCanopyDiam'
# 			trans <- identity_trans()
# 		} else if (facet == 'cn_ratio') {
# 			site_mean_var <- 'mCN_ratio'
# 			trans <- identity_trans()
# 		} else if (facet == 'height') {
# 			site_mean_var <- 'mHeight'
# 			trans <- identity_trans()
# 		} else if (facet == 'photosynthetic_rate') {
# 			site_mean_var <- 'mPhotoRate'
# 			trans <- identity_trans()
# 		} else if (facet == 'spad') {
# 			site_mean_var <- 'mSPAD'
# 			trans <- identity_trans()
# 		} else if (facet == 'stomatal_conductance') {
# 			site_mean_var <- 'mStomCond'
# 			trans <- identity_trans()
# 		} else if (facet == 'transpiration_rate') {
# 			site_mean_var <- 'mTranspRate'
# 			trans <- identity_trans()
# 		}

# 		# response variables
# 		var_sq <- paste0(facet, '_mean_sq')
# 		var_fut1 <- paste0(facet, '_mean_ssp245_2041_2070')
# 		var_fut2 <- paste0(facet, '_mean_ssp245_2071_2100')
# 		var_fut3 <- paste0(facet, '_mean_ssp370_2041_2070')
# 		var_fut4 <- paste0(facet, '_mean_ssp370_2071_2100')

# 		i <- which(models$facet == facet)
# 		form <- nonbiomass_facets[[facet]]$filename
# 		resp_distrib <- tolower(nonbiomass_facets[[facet]]$resp_distrib)

# 		uni <- vect(paste0('./outputs_loretta/integrated_sdm_pdm/models_', facet, '/[', facet, '~', resp_distrib, '(', form, ')]/prediction_vector_nam.gpkg'))

# 		uni <- crop(uni, extent)

# 		# delineate current range "core"
# 		range_core_uni_sq <- delineate_range_core(uni, column = var_sq, core_quant = core_quant)
# 		range_core_multi_sq <- delineate_range_core(multi, column = var_sq, core_quant = core_quant)

# 		range_anticore_uni_sq <- delineate_range_core(uni, column = var_sq, core_quant = anticore_quant, rule = '<')
# 		range_anticore_multi_sq <- delineate_range_core(multi, column = var_sq, core_quant = anticore_quant, rule = '<')

# 		range_core_uni_fut <- range_core_multi_fut <- list()
# 		range_anticore_uni_fut <- range_anticore_multi_fut <- list()
# 		for (fut in 1:4) {
		
# 			range_core_uni_fut[[fut]] <- delineate_range_core(uni, column = get(paste0('var_fut', fut)), core_quant = core_quant)
# 			range_core_multi_fut[[fut]] <- delineate_range_core(multi, column = get(paste0('var_fut', fut)), core_quant = core_quant)

# 			range_anticore_uni_fut[[fut]] <- delineate_range_core(uni, column = get(paste0('var_fut', fut)), core_quant = anticore_quant, rule = '<')
# 			range_anticore_multi_fut[[fut]] <- delineate_range_core(multi, column = get(paste0('var_fut', fut)), core_quant = anticore_quant, rule = '<')

# 		}

# 		# across both models, get max value of predictions for scaling
# 		vars <- c(var_sq, var_fut1, var_fut2, var_fut3, var_fut4)
# 		max_val <- -Inf
# 		min_val <- Inf
# 		for (this_vect in c('multi', 'uni')) {
# 			x <- get(this_vect)
# 			for (var in vars) {
# 				x_var <- x[[var]]
# 				x_var <- unlist(x_var)
# 				min_val <- min(min_val, min(x_var, na.rm = TRUE))
# 				max_val <- max(max_val, max(x_var, na.rm = TRUE))
# 			}
# 		}

# 		min_val <- min(min_val, min(unlist(site_vect[[site_mean_var]])))
# 		max_val <- max(max_val, max(unlist(site_vect[[site_mean_var]])))

# 		min_val <- 0.975 * min_val
# 		max_val <- 1.025 * max_val

# 		max_site_mean <- max(unlist(site_vect[[site_mean_var]]))
# 		if (max_val > max_site_mean) {
# 			max_val_plotted <- 1.25 * roundTo(max_site_mean, 0.1, ceiling)
# 		} else {
# 			max_val_plotted <- Inf
# 		}
		
# 		if (!is.infinite(max_val_plotted)) max_val <- max_val_plotted
# 		resp_limits <- c(min_val, max_val)
# 		uni[[var_sq]][uni[[var_sq]] > max_val_plotted] <- max_val_plotted
# 		multi[[var_sq]][multi[[var_sq]] > max_val_plotted] <- max_val_plotted

# 		# across both models, get min/max value of change in predictions for scaling
# 		vars <- c(var_fut1, var_fut2, var_fut3, var_fut4)
# 		min_val <- Inf
# 		max_val <- -Inf
# 		for (this_vect in c('multi', 'uni')) {
# 			x <- get(this_vect)
# 			sq <- x[[var_sq]]
# 			sq <- unlist(sq)
# 			for (var in vars) {
# 				x_var <- x[[var]]
# 				x_var <- unlist(x_var)
# 				delta <- x_var / sq
# 				min_val <- min(min_val, min(delta, na.rm = TRUE))
# 				max_val <- max(max_val, max(delta, na.rm = TRUE))
# 			}
# 		}

# 		# max_abs_val <- max(abs(c(min_val, max_val)))
# 		# delta_resp_limits <- c(-max_abs_val, max_abs_val)
# 		# if (!is.infinite(max_val_plotted)) delta_resp_limits <- c(-max_val_plotted, max_val_plotted)

# 		delta_resp_limits <- c(min_val, max_val)

# 		title <- paste0(nice$short, '-only: Present')
# 		legend_title <- nice$legend_title
# 		if (!is.infinite(max_val_plotted)) legend_title <- paste0(legend_title, '\n(max = ', round(max_val_plotted, 1), ')')
# 		sq_uni <- ggplot() +
# 			layer_spatial(uni, aes(fill = .data[[var_sq]]), color = NA) +
# 			layer_spatial(nam, color = 'gray40', fill = NA, linewidth = 0.3) +
# 			layer_spatial(counties_with_ag, pch = 16, color = alpha('gray20', 0.3), size = 0.15) +
# 			layer_spatial(range_core_uni_sq, color = 'orange2', fill = NA, linewidth = 0.4) +
# 			layer_spatial(range_anticore_uni_sq, color = 'red', fill = NA, linewidth = 0.4) +
# 			layer_spatial(site_vect, aes(fill = .data[[site_mean_var]]), pch = 21, size = 2) +
# 			scale_fill_gradientn(
# 				name = legend_title,
# 				colors = c('#edf8e9', '#74c476', '#005a32'),
# 				trans = trans,
# 				limits = resp_limits
# 			) +
# 			ggtitle(title) +
# 			coord_sf(xlim = c(extent_coords[1], extent_coords[2]), ylim = c(extent_coords[3], extent_coords[4]), expand = FALSE) +
# 			theme(
# 				plot.title = element_text(size = 10),
# 				axis.text = element_blank(),
# 				axis.ticks = element_blank(),
# 				legend.title = element_text(size = legend_title_size),
# 				legend.text = element_text(size = 6),
# 				legend.position = c(legend_x, legend_y),
# 				legend.justification = c(1, 0),
# 				legend.key.width = grid::unit(0.24, 'cm'),
# 				legend.key.height = grid::unit(0.35, 'cm'),
# 				legend.background = element_rect(fill = alpha('white', 0))
# 			)

# 		title <- 'Integrated: Present'
# 		legend_title <- nice$legend_title
# 		if (!is.infinite(max_val_plotted)) legend_title <- paste0(legend_title, '\n(max = ', round(max_val_plotted, 1), ')')
# 		sq_multi <- ggplot() +
# 			layer_spatial(multi, aes(fill = .data[[var_sq]]), color = NA) +
# 			layer_spatial(nam, color = 'gray40', fill = NA, linewidth = 0.3) +
# 			layer_spatial(counties_with_ag, pch = 16, color = alpha('gray20', 0.3), size = 0.15) +
# 			layer_spatial(range_core_multi_sq, color = 'orange2', fill = NA, linewidth = 0.4) +
# 			layer_spatial(range_anticore_multi_sq, color = 'red', fill = NA, linewidth = 0.4) +
# 			layer_spatial(site_vect, aes(fill = .data[[site_mean_var]]), pch = 21, size = 2) +
# 			scale_fill_gradientn(
# 				name = legend_title,
# 				colors = c('#edf8e9', '#74c476', '#005a32'),
# 				trans = trans,
# 				limits = resp_limits
# 			) +
# 			ggtitle(title) +
# 			coord_sf(xlim = c(extent_coords[1], extent_coords[2]), ylim = c(extent_coords[3], extent_coords[4]), expand = FALSE) +
# 			theme(
# 				plot.title = element_text(size = 10),
# 				axis.text = element_blank(),
# 				axis.ticks = element_blank(),
# 				legend.title = element_text(size = legend_title_size),
# 				legend.text = element_text(size = 6),
# 				legend.position = c(legend_x, legend_y),
# 				legend.justification = c(1, 0),
# 				legend.key.width = grid::unit(0.24, 'cm'),
# 				legend.key.height = grid::unit(0.35, 'cm'),
# 				legend.background = element_rect(fill = alpha('white', 0))
# 			)

# 		### standalone models: future predictions and deltas
# 		####################################################

# 		futs <- deltas <- list()
# 		for (fut in 1:4) {

# 			if (fut == 1) {
# 				var_fut <- var_fut1
# 				title_fut <- paste0(nice$short, '-only: SSP245 2041-2070')
# 				title_delta <- 'Change: 2041-2070 - Present'
# 			} else if (fut == 2) {
# 				var_fut <- var_fut2
# 				title_fut <- paste0(nice$short, '-only: SSP245 2071-2100')
# 				title_delta <- 'Change: 2071-2100 - Present'
# 			} else if (fut == 3) {
# 				var_fut <- var_fut3
# 				title_fut <- paste0(nice$short, '-only: SSP370 2041-2070')
# 				title_delta <- 'Change: 2041-2070 - Present'
# 			} else if (fut == 4) {
# 				var_fut <- var_fut4
# 				title_fut <- paste0(nice$short, '-only: SSP370 2071-2100')
# 				title_delta <- 'Change: 2071-2100 - Present'
# 			}

# 			uni[[var_fut]][uni[[var_fut]] > max_val_plotted] <- max_val_plotted

# 			legend_title <- nice$legend_title
# 			if (!is.infinite(max_val_plotted)) legend_title <- paste0(legend_title, '\n(max = ', round(max_val_plotted, 1), ')')
# 			futs[[fut]] <- ggplot() +
# 				layer_spatial(uni, aes(fill = .data[[var_fut]]), color = NA) +
# 				layer_spatial(nam, color = 'gray40', fill = NA, linewidth = 0.3) +
# 				layer_spatial(counties_with_ag, pch = 16, color = alpha('gray20', 0.3), size = 0.15) +
# 				layer_spatial(range_core_uni_fut[[fut]], color = 'orange2', fill = NA, linewidth = 0.4) +
# 				layer_spatial(range_anticore_uni_fut[[fut]], color = 'red', fill = NA, linewidth = 0.4) +
# 				layer_spatial(site_vect, pch = 3, size = 1) +
# 				scale_fill_gradientn(
# 					name = legend_title,
# 					colors = c('#edf8e9', '#74c476', '#005a32'),
# 					trans = trans,
# 					limits = resp_limits
# 				) +
# 				ggtitle(title_fut) +
# 				coord_sf(xlim = c(extent_coords[1], extent_coords[2]), ylim = c(extent_coords[3], extent_coords[4]), expand = FALSE) +
# 				theme(
# 					plot.title = element_text(size = 10),
# 					axis.text = element_blank(),
# 					axis.ticks = element_blank(),
# 					legend.title = element_text(size = legend_title_size),
# 					legend.text = element_text(size = 6),
# 					legend.key.width = grid::unit(0.24, 'cm'),
# 					legend.key.height = grid::unit(0.35, 'cm'),
# 					legend.position = c(legend_x, legend_y),
# 					legend.justification = c(1, 0),
# 					legend.background = element_rect(fill = alpha('white', 0))
# 				)

# 			legend_title <- 'log10-ratio'
# 			if (!is.infinite(max_val_plotted)) legend_title <- paste0(legend_title, '*')
# 			deltas[[fut]] <- ggplot() +
# 				layer_spatial(uni, aes(fill = .data[[var_fut]] / .data[[var_sq]]), color = NA) +
# 				layer_spatial(nam, color = 'gray40', fill = NA, linewidth = 0.3) +
# 				layer_spatial(counties_with_ag, pch = 16, color = alpha('gray20', 0.3), size = 0.15) +
# 				layer_spatial(range_core_uni_fut[[fut]], color = 'orange2', fill = NA, linewidth = 0.4) +
# 				layer_spatial(range_anticore_uni_fut[[fut]], color = 'red', fill = NA, linewidth = 0.4) +
# 				layer_spatial(site_vect, pch = 3, size = 1) +
# 				scale_fill_gradientn(
# 					name = legend_title,
# 					colors = c('#8c510a', '#d8b365', '#f6e8c3', '#f5f5f5', '#c7eae5', '#5ab4ac', '#01665e'),
# 					limits = delta_resp_limits,
# 					trans = 'log10',
# 					values = scales::rescale(c(delta_resp_limits[1], 1, delta_resp_limits[2]))
# 				) +
# 				ggtitle(title_delta) +
# 				coord_sf(xlim = c(extent_coords[1], extent_coords[2]), ylim = c(extent_coords[3], extent_coords[4]), expand = FALSE) +
# 				theme(
# 					plot.title = element_text(size = 10),
# 					axis.text = element_blank(),
# 					axis.ticks = element_blank(),
# 					legend.title = element_text(size = legend_title_size),
# 					legend.text = element_text(size = 6),
# 					legend.key.width = grid::unit(0.24, 'cm'),
# 					legend.key.height = grid::unit(0.35, 'cm'),
# 					legend.position = c(legend_x_delta, legend_y),
# 					legend.justification = c(1, 0),
# 					legend.background = element_rect(fill = alpha('white', 0))
# 				)

# 		}

# 		maps_futs <- plot_grid(plotlist = futs, ncol = 1, align = 'hv')
# 		maps_deltas <- plot_grid(plotlist = deltas, ncol = 1, align = 'hv')
# 		maps_futs_deltas <- plot_grid(maps_futs, maps_deltas, ncol = 2, align = 'hv')
# 		maps <- plot_grid(sq_uni, maps_futs_deltas, ncol = 1, align = 'hv', rel_heights = c(1.2, 4))

# 		ggsave(maps, filename = paste0('./outputs_loretta/integrated_sdm_pdm/maps_', facet, '_all_maps_non_integrated.png'), width = 6, height = 9, dpi = 600, bg = 'white')

# 		### integrated models: future predictions and deltas
# 		####################################################

# 		futs <- deltas <- list()
# 		for (fut in 1:4) {

# 			if (fut == 1) {
# 				var_fut <- var_fut1
# 				title_fut <- 'Integrated: SSP245 2041-2070'
# 				title_delta <- 'Change: 2041-2070 - Present'
# 			} else if (fut == 2) {
# 				var_fut <- var_fut2
# 				title_fut <- 'Integrated: SSP245 2071-2100'
# 				title_delta <- 'Change: 2071-2100 - Present'
# 			} else if (fut == 3) {
# 				var_fut <- var_fut3
# 				title_fut <- 'Integrated: SSP370 2041-2070'
# 				title_delta <- 'Change: 2041-2070 - Present'
# 			} else if (fut == 4) {
# 				var_fut <- var_fut4
# 				title_fut <- 'Integrated: SSP370 2071-2100'
# 				title_delta <- 'Change: 2071-2100 - Present'
# 			}

# 			multi[[var_fut]][multi[[var_fut]] > max_val_plotted] <- max_val_plotted

# 			legend_title <- nice$legend_title
# 			if (!is.infinite(max_val_plotted)) legend_title <- paste0(legend_title, '\n(max = ', round(max_val_plotted, 1), ')')
# 			futs[[fut]] <- ggplot() +
# 				layer_spatial(multi, aes(fill = .data[[var_fut]]), color = NA) +
# 				layer_spatial(nam, color = 'gray40', fill = NA, linewidth = 0.3) +
# 				layer_spatial(counties_with_ag, pch = 16, color = alpha('gray20', 0.3), size = 0.15) +
# 				layer_spatial(range_core_multi_fut[[fut]], color = 'orange2', fill = NA, linewidth = 0.4) +
# 				layer_spatial(range_anticore_multi_fut[[fut]], color = 'red', fill = NA, linewidth = 0.4) +
# 				layer_spatial(site_vect, pch = 3, size = 1) +
# 				scale_fill_gradientn(
# 					name = legend_title,
# 					colors = c('#edf8e9', '#74c476', '#005a32'),
# 					trans = trans,
# 				limits = resp_limits
# 				) +
# 				ggtitle(title_fut) +
# 				coord_sf(xlim = c(extent_coords[1], extent_coords[2]), ylim = c(extent_coords[3], extent_coords[4]), expand = FALSE) +
# 				theme(
# 					plot.title = element_text(size = 10),
# 					axis.text = element_blank(),
# 					axis.ticks = element_blank(),
# 					legend.title = element_text(size = legend_title_size),
# 					legend.text = element_text(size = 6),
# 					legend.key.width = grid::unit(0.24, 'cm'),
# 					legend.key.height = grid::unit(0.35, 'cm'),
# 					legend.position = c(legend_x, legend_y),
# 					legend.justification = c(1, 0),
# 					legend.background = element_rect(fill = alpha('white', 0))
# 				)

# 			legend_title <- 'log10-ratio'
# 			if (!is.infinite(max_val_plotted)) legend_title <- paste0(legend_title, '*')
# 			deltas[[fut]] <- ggplot() +
# 				layer_spatial(multi, aes(fill = .data[[var_fut]] / .data[[var_sq]]), color = NA) +
# 				layer_spatial(nam, color = 'gray40', fill = NA, linewidth = 0.3) +
# 				layer_spatial(counties_with_ag, pch = 16, color = alpha('gray20', 0.3), size = 0.15) +
# 				layer_spatial(range_core_multi_fut[[fut]], color = 'orange2', fill = NA, linewidth = 0.4) +
# 				layer_spatial(range_anticore_multi_fut[[fut]], color = 'red', fill = NA, linewidth = 0.4) +
# 				layer_spatial(site_vect, pch = 3, size = 1) +
# 				scale_fill_gradientn(
# 					name = legend_title,
# 					colors = c('#8c510a', '#d8b365', '#f6e8c3', '#f5f5f5', '#c7eae5', '#5ab4ac', '#01665e'),
# 					limits = delta_resp_limits,
# 					trans = 'log10',
# 					values = scales::rescale(c(delta_resp_limits[1], 1, delta_resp_limits[2]))
# 				) +
# 				ggtitle(title_delta) +
# 				coord_sf(xlim = c(extent_coords[1], extent_coords[2]), ylim = c(extent_coords[3], extent_coords[4]), expand = FALSE) +
# 				theme(
# 					plot.title = element_text(size = 10),
# 					axis.text = element_blank(),
# 					axis.ticks = element_blank(),
# 					legend.title = element_text(size = legend_title_size),
# 					legend.text = element_text(size = 6),
# 					legend.key.width = grid::unit(0.24, 'cm'),
# 					legend.key.height = grid::unit(0.35, 'cm'),
# 					legend.position = c(legend_x_delta, legend_y),
# 					legend.justification = c(1, 0),
# 					legend.background = element_rect(fill = alpha('white', 0))
# 				)

# 		}

# 		maps_futs <- plot_grid(plotlist = futs, ncol = 1, align = 'hv')
# 		maps_deltas <- plot_grid(plotlist = deltas, ncol = 1, align = 'hv')
# 		maps_futs_deltas <- plot_grid(maps_futs, maps_deltas, ncol = 2, align = 'hv')
# 		maps <- plot_grid(sq_multi, maps_futs_deltas, ncol = 1, align = 'hv', rel_heights = c(1.2, 4))

# 		ggsave(maps, filename = paste0('./outputs_loretta/integrated_sdm_pdm/maps_', facet, '_all_maps_integrated.png'), width = 6, height = 9, dpi = 600, bg = 'white')

# 	} # next facet

# say('#####################################################################')
# say('### maps of abundance in Dust Bowl region in 1930s for supplement ###')
# say('#####################################################################')

# # NB in code, "_db" == > Dust Bowl; "_sq" == > status quo (present)

# 	# US states
# 	nam <- vect('./data_from_gadm/gadm_4pt1_level_1_north_america_sans_alaska_lambert.gpkg')
# 	nam <- simplifyGeom(nam, tolerance = 1000)

# 	# Dust Bowl region
# 	dust_bowl <- vect('./data_from_others/dust_bowl_counties_with_most_severe_wind_erosion.gpkg')
# 	dust_bowl <- aggregate(dust_bowl)

# 	# sites
# 	data_traits <- prepare_nonbiomass_data(facet = 'height', formula = ~ 1, n_response_curve_values = n_response_curve_values, calib = calib)
# 	site_vect <- vect(data_traits$site_data_raw, geom = c('LONGITUDE', 'LATITUDE'), crs = getCRS('WGS84'))
# 	site_vect <- project(site_vect, nam)

# 	# response vectors
# 	uni_db <- vect('./outputs_loretta/integrated_sdm_pdm/models_occurrence/[occs~overdispersed_hurdlepoisson(bio1^2_log(bio12)^2_bio15^2)]_[bias~poaceae_offset_plus_exp_neg_1_WITH_zero_correction_covariate]/prediction_vector_conus_1930s.gpkg')
# 	uni_sq <- vect('./outputs_loretta/integrated_sdm_pdm/models_occurrence/[occs~overdispersed_hurdlepoisson(bio1^2_log(bio12)^2_bio15^2)]_[bias~poaceae_offset_plus_exp_neg_1_WITH_zero_correction_covariate]/prediction_vector_nam.gpkg')

# 	multi_db <- vect('./outputs_loretta/integrated_sdm_pdm/models_integrated/non_centered_MVN([occs~hurdlepoisson_offset]_[biomass~hurdleln]_[bla_can_hei])_eta=2/prediction_vector_1930s.gpkg')
# 	multi_sq <- vect('./outputs_loretta/integrated_sdm_pdm/models_integrated/non_centered_MVN([occs~hurdlepoisson_offset]_[biomass~hurdleln]_[bla_can_hei])_eta=2/prediction_vector_nam.gpkg')

# 	uni_db <- uni_db[uni_db$state_province %in% c('Colorado', 'Nebraska', 'Kansas', 'Oklahoma', 'Texas', 'New Mexico')]
# 	uni_sq <- uni_sq[uni_sq$state_province %in% c('Colorado', 'Nebraska', 'Kansas', 'Oklahoma', 'Texas', 'New Mexico')]

# 	multi_db <- multi_db[multi_db$state_province %in% c('Colorado', 'Nebraska', 'Kansas', 'Oklahoma', 'Texas', 'New Mexico')]
# 	multi_sq <- multi_sq[multi_sq$state_province %in% c('Colorado', 'Nebraska', 'Kansas', 'Oklahoma', 'Texas', 'New Mexico')]

# 	uni_db <- project(uni_db, nam)
# 	uni_sq <- project(uni_sq, nam)
# 	multi_db <- project(multi_db, nam)
# 	multi_sq <- project(multi_sq, nam)

# 	extent <- ext(dust_bowl)
# 	extent <- as.polygons(extent, crs = nam)
# 	extent <- buffer(extent, width = 20 * 1000) # nominal plot extent
# 	extent_display <- buffer(extent, width = 30 * 1000) # larger than plot extent
# 	extent <- ext(extent)
# 	extent <- as.vector(extent)

# 	# change
# 	uni_db$delta_abund <- (uni_db$N_ag_mean_1930s - uni_sq$N_ag_mean_sq) / uni_sq$N_ag_mean_sq
# 	uni_db$delta_psi <- uni_db$psi_1930s - uni_sq$psi_sq
	
# 	multi_db$delta_abund <- (multi_db$N_ag_mean_1930s - multi_sq$N_ag_mean_sq) / multi_sq$N_ag_mean_sq
# 	multi_db$delta_psi <- multi_db$psi_1930s - multi_sq$psi_sq

# 	# crop to display area
# 	nam <- crop(nam, extent)
# 	uni_db <- crop(uni_db, extent)
# 	uni_sq <- crop(uni_sq, extent)
# 	multi_db <- crop(multi_db, extent)
# 	multi_sq <- crop(multi_sq, extent)

# 	# response limits for color scales
# 	min_val <- min(uni_db$delta_abund, multi_db$delta_abund)
# 	max_val <- max(uni_db$delta_abund, multi_db$delta_abund)
# 	resp_limits_abund <- c(min_val, max_val)

# 	min_val <- min(uni_db$delta_psi, multi_db$delta_psi)
# 	max_val <- max(uni_db$delta_psi, multi_db$delta_psi)
# 	resp_limits_psi <- c(min_val, max_val)

# 	# state labels at state centroids
# 	centroids <- crds(centroids(nam))
# 	centroids <- as.data.frame(centroids)
# 	centroids$y[nam$NAME_1 == 'Oklahoma'] <- centroids$y[nam$NAME_1 == 'Oklahoma'] - 90000

# 	# maps
# 	plot_title_size <- 14
# 	state_labels_size <- 4
	
# 	map_uni_abund <- ggplot() +
# 		layer_spatial(nam, fill = 'gainsboro', color = NA) +
# 		layer_spatial(
# 			uni_db,
# 			aes(fill = delta_abund),
# 			color = NA
# 		) +
# 		scale_fill_gradient2(
# 			name = 'Change',
# 			low = '#c51b7d',
# 			mid = '#f7f7f7',
# 			high = '#4d9221',
# 			midpoint = 0,
# 			limits = resp_limits_abund,
# 			labels = scales::percent_format(accuracy = 1)
# 		) +
# 		layer_spatial(nam, color = 'gray40', fill = NA) +
# 		geom_text(data = centroids, aes(x = x, y = y, label = nam$NAME_1), size = state_labels_size) +
# 		layer_spatial(dust_bowl, color = 'gray20', fill = NA, linewidth = 1) +
# 		layer_spatial(site_vect, pch = 3, size = 4) +
# 		annotation_scale(location = 'br') +
# 		coord_sf(xlim = c(extent[1], extent[2]), ylim = c(extent[3], extent[4]), expand = FALSE) +
# 		ggtitle('Abundance: Abundance-only Model') +
# 		theme(
# 			plot.title = element_text(size = plot_title_size),
# 			axis.title = element_blank()
# 		)

# 	map_uni_psi <- ggplot() +
# 		layer_spatial(nam, fill = 'gainsboro', color = NA) +
# 		layer_spatial(
# 			uni_db,
# 			aes(fill = delta_psi),
# 			color = NA
# 		) +
# 		scale_fill_gradient2(
# 			name = 'Change',
# 			low = '#c51b7d',
# 			mid = '#f7f7f7',
# 			high = '#4d9221',
# 			midpoint = 0,
# 			limits = resp_limits_psi,
# 			labels = scales::percent_format(accuracy = 1)
# 		) +
# 		layer_spatial(nam, color = 'gray40', fill = NA) +
# 		geom_text(data = centroids, aes(x = x, y = y, label = nam$NAME_1), size = state_labels_size) +
# 		layer_spatial(dust_bowl, color = 'gray20', fill = NA, linewidth = 1) +
# 		layer_spatial(site_vect, pch = 3, size = 4) +
# 		annotation_scale(location = 'br') +
# 		coord_sf(xlim = c(extent[1], extent[2]), ylim = c(extent[3], extent[4]), expand = FALSE) +
# 		ggtitle('Occupancy: Abundance-only Model') +
# 		theme(
# 			plot.title = element_text(size = plot_title_size),
# 			axis.title = element_blank()
# 		)

# 	map_multi_abund <- ggplot() +
# 		layer_spatial(nam, fill = 'gainsboro', color = NA) +
# 		layer_spatial(
# 			multi_db,
# 			aes(fill = delta_abund),
# 			color = NA
# 		) +
# 		scale_fill_gradient2(
# 			name = 'Change',
# 			low = '#c51b7d',
# 			mid = '#f7f7f7',
# 			high = '#4d9221',
# 			midpoint = 0,
# 			limits = resp_limits_abund,
# 			labels = scales::percent_format(accuracy = 1)
# 		) +
# 		layer_spatial(nam, color = 'gray40', fill = NA) +
# 		geom_text(data = centroids, aes(x = x, y = y, label = nam$NAME_1), size = state_labels_size) +
# 		layer_spatial(dust_bowl, color = 'gray20', fill = NA, linewidth = 1) +
# 		layer_spatial(site_vect, pch = 3, size = 4) +
# 		annotation_scale(location = 'br') +
# 		coord_sf(xlim = c(extent[1], extent[2]), ylim = c(extent[3], extent[4]), expand = FALSE) +
# 		ggtitle('Abundance: Integrated Model') +
# 		theme(
# 			plot.title = element_text(size = plot_title_size),
# 			axis.title = element_blank()
# 		)

# 	map_multi_psi <- ggplot() +
# 		layer_spatial(nam, fill = 'gainsboro', color = NA) +
# 		layer_spatial(
# 			multi_db,
# 			aes(fill = delta_psi),
# 			color = NA
# 		) +
# 		scale_fill_gradient2(
# 			name = 'Change',
# 			low = '#c51b7d',
# 			mid = '#f7f7f7',
# 			high = '#4d9221',
# 			midpoint = 0,
# 			limits = resp_limits_psi,
# 			labels = scales::percent_format(accuracy = 1)
# 		) +
# 		layer_spatial(nam, color = 'gray40', fill = NA) +
# 		geom_text(data = centroids, aes(x = x, y = y, label = nam$NAME_1), size = state_labels_size) +
# 		layer_spatial(dust_bowl, color = 'gray20', fill = NA, linewidth = 1) +
# 		layer_spatial(site_vect, pch = 3, size = 4) +
# 		annotation_scale(location = 'br') +
# 		coord_sf(xlim = c(extent[1], extent[2]), ylim = c(extent[3], extent[4]), expand = FALSE) +
# 		ggtitle('Occupancy: Integrated Model') +
# 		theme(
# 			plot.title = element_text(size = plot_title_size),
# 			axis.title = element_blank()
# 		)

# 	map_grobs <- list(map_uni_abund, map_multi_abund, map_uni_psi, map_multi_psi)
# 	names(map_grobs) <- c('delta_univariate', 'delta_multivariate', 'delta_psi_univariate', 'delta_psi_multivariate')
# 	saveRDS(map_grobs, './outputs_loretta/integrated_sdm_pdm/maps_1930s_dust_bowl_abundance_change.rds')

# 	maps <- plot_grid(plotlist = map_grobs, ncol = 2, align = 'hv')
# 	ggsave(maps, filename = './outputs_loretta/integrated_sdm_pdm/maps_1930s_dust_bowl_abundance_change.png', width = 10, height = 11, dpi = 600, bg = 'white')

# say('###################################################################')
# say('### maps of biomass in Dust Bowl region in 1930s for supplement ###')
# say('###################################################################')

# # NB in code, "_db" == > Dust Bowl; "_sq" == > status quo (present)

# 	# US states
# 	nam <- vect('./data_from_gadm/gadm_4pt1_level_1_north_america_sans_alaska_lambert.gpkg')
# 	nam <- simplifyGeom(nam, tolerance = 1000)

# 	# Dust Bowl region
# 	dust_bowl <- vect('./data_from_others/dust_bowl_counties_with_most_severe_wind_erosion.gpkg')
# 	dust_bowl <- aggregate(dust_bowl)

# 	# sites
# 	data_traits <- prepare_nonbiomass_data(facet = 'height', formula = ~ 1, n_response_curve_values = n_response_curve_values, calib = calib)
# 	site_vect <- vect(data_traits$site_data_raw, geom = c('LONGITUDE', 'LATITUDE'), crs = getCRS('WGS84'))
# 	site_vect <- project(site_vect, nam)

# 	# response vectors
# 	uni_db <- vect('./outputs_loretta/integrated_sdm_pdm/models_biomass/[biomass~hurdleln_log(bio12)]/prediction_vector_conus_1930s.gpkg')
# 	uni_sq <- vect('./outputs_loretta/integrated_sdm_pdm/models_biomass/[biomass~hurdleln_log(bio12)]/prediction_vector_nam.gpkg')

# 	multi_db <- vect('./outputs_loretta/integrated_sdm_pdm/models_integrated/non_centered_MVN([occs~hurdlepoisson_offset]_[biomass~hurdleln]_[bla_can_hei])_eta=2/prediction_vector_1930s.gpkg')
# 	multi_sq <- vect('./outputs_loretta/integrated_sdm_pdm/models_integrated/non_centered_MVN([occs~hurdlepoisson_offset]_[biomass~hurdleln]_[bla_can_hei])_eta=2/prediction_vector_nam.gpkg')

# 	uni_db <- uni_db[uni_db$state_province %in% c('Colorado', 'Nebraska', 'Kansas', 'Oklahoma', 'Texas', 'New Mexico')]
# 	uni_sq <- uni_sq[uni_sq$state_province %in% c('Colorado', 'Nebraska', 'Kansas', 'Oklahoma', 'Texas', 'New Mexico')]

# 	multi_db <- multi_db[multi_db$state_province %in% c('Colorado', 'Nebraska', 'Kansas', 'Oklahoma', 'Texas', 'New Mexico')]
# 	multi_sq <- multi_sq[multi_sq$state_province %in% c('Colorado', 'Nebraska', 'Kansas', 'Oklahoma', 'Texas', 'New Mexico')]

# 	uni_db <- project(uni_db, nam)
# 	uni_sq <- project(uni_sq, nam)
# 	multi_db <- project(multi_db, nam)
# 	multi_sq <- project(multi_sq, nam)

# 	extent <- ext(dust_bowl)
# 	extent <- as.polygons(extent, crs = nam)
# 	extent <- buffer(extent, width = 20 * 1000) # nominal plot extent
# 	extent_display <- buffer(extent, width = 30 * 1000) # larger than plot extent
# 	extent <- ext(extent)
# 	extent <- as.vector(extent)

# 	# change and response limits
# 	uni_sq$delta_y <- (uni_db$biomass_mean_1930s - uni_sq$biomass_mean_sq) / uni_sq$biomass_mean_sq
# 	uni_sq$delta_psi <- uni_db$psi_1930s - uni_sq$psi_sq
	
# 	multi_sq$delta_y <- (multi_db$biomass_mean_1930s - multi_sq$biomass_mean_sq) / multi_sq$biomass_mean_sq
# 	multi_sq$delta_psi <- multi_db$psi_1930s - multi_sq$psi_sq

# 	# crop to display area
# 	nam <- crop(nam, extent)
# 	uni_db <- crop(uni_db, extent)
# 	uni_sq <- crop(uni_sq, extent)
# 	multi_db <- crop(multi_db, extent)
# 	multi_sq <- crop(multi_sq, extent)

# 	# response limits for color scales
# 	min_val <- min(uni_sq$delta_y, multi_sq$delta_y)
# 	max_val <- max(uni_sq$delta_y, multi_sq$delta_y)
# 	resp_limits_y <- c(min_val, max_val)

# 	min_val <- min(uni_sq$delta_psi, multi_sq$delta_psi)
# 	max_val <- max(uni_sq$delta_psi, multi_sq$delta_psi)
# 	resp_limits_psi <- c(min_val, max_val)

# 	# state labels at state centroids
# 	centroids <- crds(centroids(nam))
# 	centroids <- as.data.frame(centroids)
# 	centroids$y[nam$NAME_1 == 'Oklahoma'] <- centroids$y[nam$NAME_1 == 'Oklahoma'] - 90000

# 	# maps
# 	plot_title_size <- 14
# 	state_labels_size <- 4
	
# 	map_uni_y <- ggplot() +
# 		layer_spatial(nam, fill = 'gainsboro', color = NA) +
# 		layer_spatial(
# 			uni_sq,
# 			aes(fill = delta_y),
# 			color = NA
# 		) +
# 		scale_fill_gradient2(
# 			name = 'Change',
# 			low = '#c51b7d',
# 			mid = '#f7f7f7',
# 			high = '#4d9221',
# 			midpoint = 0,
# 			limits = resp_limits_y,
# 			labels = scales::percent_format(accuracy = 1)
# 		) +
# 		layer_spatial(nam, color = 'gray40', fill = NA) +
# 		geom_text(data = centroids, aes(x = x, y = y, label = nam$NAME_1), size = state_labels_size) +
# 		layer_spatial(dust_bowl, color = 'gray20', fill = NA, linewidth = 1) +
# 		layer_spatial(site_vect, pch = 3, size = 4) +
# 		annotation_scale(location = 'br') +
# 		coord_sf(xlim = c(extent[1], extent[2]), ylim = c(extent[3], extent[4]), expand = FALSE) +
# 		ggtitle('Biomass: Biomass-only Model') +
# 		theme(
# 			plot.title = element_text(size = plot_title_size),
# 			axis.title = element_blank()
# 		)

# 	map_uni_psi <- ggplot() +
# 		layer_spatial(nam, fill = 'gainsboro', color = NA) +
# 		layer_spatial(
# 			uni_sq,
# 			aes(fill = delta_psi),
# 			color = NA
# 		) +
# 		scale_fill_gradient2(
# 			name = 'Change',
# 			low = '#c51b7d',
# 			mid = '#f7f7f7',
# 			high = '#4d9221',
# 			midpoint = 0,
# 			limits = resp_limits_psi,
# 			labels = scales::percent_format(accuracy = 1)
# 		) +
# 		layer_spatial(nam, color = 'gray40', fill = NA) +
# 		geom_text(data = centroids, aes(x = x, y = y, label = nam$NAME_1), size = state_labels_size) +
# 		layer_spatial(dust_bowl, color = 'gray20', fill = NA, linewidth = 1) +
# 		layer_spatial(site_vect, pch = 3, size = 4) +
# 		annotation_scale(location = 'br') +
# 		coord_sf(xlim = c(extent[1], extent[2]), ylim = c(extent[3], extent[4]), expand = FALSE) +
# 		ggtitle('Occupancy: Biomass-only Model') +
# 		theme(
# 			plot.title = element_text(size = plot_title_size),
# 			axis.title = element_blank()
# 		)

# 	map_multi_y <- ggplot() +
# 		layer_spatial(nam, fill = 'gainsboro', color = NA) +
# 		layer_spatial(
# 			multi_sq,
# 			aes(fill = delta_y),
# 			color = NA
# 		) +
# 		scale_fill_gradient2(
# 			name = 'Change',
# 			low = '#c51b7d',
# 			mid = '#f7f7f7',
# 			high = '#4d9221',
# 			midpoint = 0,
# 			limits = resp_limits_y,
# 			labels = scales::percent_format(accuracy = 1)
# 		) +
# 		layer_spatial(nam, color = 'gray40', fill = NA) +
# 		geom_text(data = centroids, aes(x = x, y = y, label = nam$NAME_1), size = state_labels_size) +
# 		layer_spatial(dust_bowl, color = 'gray20', fill = NA, linewidth = 1) +
# 		layer_spatial(site_vect, pch = 3, size = 4) +
# 		annotation_scale(location = 'br') +
# 		coord_sf(xlim = c(extent[1], extent[2]), ylim = c(extent[3], extent[4]), expand = FALSE) +
# 		ggtitle('Biomass: Integrated Model') +
# 		theme(
# 			plot.title = element_text(size = plot_title_size),
# 			axis.title = element_blank()
# 		)

# 	map_multi_psi <- ggplot() +
# 		layer_spatial(nam, fill = 'gainsboro', color = NA) +
# 		layer_spatial(
# 			multi_sq,
# 			aes(fill = delta_psi),
# 			color = NA
# 		) +
# 		scale_fill_gradient2(
# 			name = 'Change',
# 			low = '#c51b7d',
# 			mid = '#f7f7f7',
# 			high = '#4d9221',
# 			midpoint = 0,
# 			limits = resp_limits_psi,
# 			labels = scales::percent_format(accuracy = 1)
# 		) +
# 		layer_spatial(nam, color = 'gray40', fill = NA) +
# 		geom_text(data = centroids, aes(x = x, y = y, label = nam$NAME_1), size = state_labels_size) +
# 		layer_spatial(dust_bowl, color = 'gray20', fill = NA, linewidth = 1) +
# 		layer_spatial(site_vect, pch = 3, size = 4) +
# 		annotation_scale(location = 'br') +
# 		coord_sf(xlim = c(extent[1], extent[2]), ylim = c(extent[3], extent[4]), expand = FALSE) +
# 		ggtitle('Occupancy: Integrated Model') +
# 		theme(
# 			plot.title = element_text(size = plot_title_size),
# 			axis.title = element_blank()
# 		)

# 	map_grobs <- list(map_uni_y, map_multi_y, map_uni_psi, map_multi_psi)
# 	names(map_grobs) <- c('delta_univariate', 'delta_multivariate', 'delta_psi_univariate', 'delta_psi_multivariate')
# 	saveRDS(map_grobs, file = './outputs_loretta/integrated_sdm_pdm/maps_1930s_dust_bowl_biomass_change.rds')

# 	maps <- plot_grid(plotlist = map_grobs, ncol = 2, align = 'hv')
# 	ggsave(maps, filename = './outputs_loretta/integrated_sdm_pdm/maps_1930s_dust_bowl_biomass_change.png', width = 10, height = 11, dpi = 600, bg = 'white')

# say('#######################################################################')
# say('### maps of non-biomass in Dust Bowl region in 1930s for supplement ###')
# say('#######################################################################')

# 	# NB in code, "_db" == > Dust Bowl; "_sq" == > status quo (present)

# 	# US states
# 	nam <- vect('./data_from_gadm/gadm_4pt1_level_1_north_america_sans_alaska_lambert.gpkg')
# 	nam <- simplifyGeom(nam, tolerance = 1000)

# 	# Dust Bowl region
# 	dust_bowl <- vect('./data_from_others/dust_bowl_counties_with_most_severe_wind_erosion.gpkg')
# 	dust_bowl <- aggregate(dust_bowl)

# 		# plot extent
# 		extent <- ext(dust_bowl)
# 		extent <- as.polygons(extent, crs = nam)
# 		extent <- buffer(extent, width = 20 * 1000) # nominal plot extent
# 		extent <- ext(extent)
# 		extent <- as.vector(extent)

# 	# sites
# 	data_traits <- prepare_nonbiomass_data(facet = 'height', formula = ~ 1, n_response_curve_values = n_response_curve_values, calib = calib)
# 	site_vect <- vect(data_traits$site_data_raw, geom = c('LONGITUDE', 'LATITUDE'), crs = getCRS('WGS84'))
# 	site_vect <- project(site_vect, nam)

# 	# best models
# 	models <- read_xlsx('./outputs_loretta/integrated_sdm_pdm/summary_of_ALL_top_models_manual_assessment.xlsx', sheet = 'summary_of_ALL_top_models')
# 	models <- as.data.table(models)
# 	models <- models[models$selected]
# 	facets <- models$facet
# 	facets <- facets[facets != 'biomass']

# 	# multivariate models
# 	multi_morph_db <- vect('./outputs_loretta/integrated_sdm_pdm/models_integrated/non_centered_MVN([occs~hurdlepoisson_offset]_[biomass~hurdleln]_[bla_can_hei])_eta=2/prediction_vector_1930s.gpkg')
# 	multi_morph_sq <- vect('./outputs_loretta/integrated_sdm_pdm/models_integrated/non_centered_MVN([occs~hurdlepoisson_offset]_[biomass~hurdleln]_[bla_can_hei])_eta=2/prediction_vector_nam.gpkg')
	
# 	multi_phys_db <- vect('./outputs_loretta/integrated_sdm_pdm/models_integrated/non_centered_MVN([occs~hurdlepoisson_offset]_[biomass~hurdleln]_[cnr_pho_spa_sto_tra])_eta=2/prediction_vector_1930s.gpkg')
# 	multi_phys_sq <- vect('./outputs_loretta/integrated_sdm_pdm/models_integrated/non_centered_MVN([occs~hurdlepoisson_offset]_[biomass~hurdleln]_[cnr_pho_spa_sto_tra])_eta=2/prediction_vector_nam.gpkg')

# 	for (f in seq_along(facets)) {

# 		facet <- facets[f]
# 		say(facet)
	
# 		resp_distrib <- tolower(nonbiomass_facets[[facet]]$resp_distrib)
# 		form_filename <- nonbiomass_facets[[facet]]$filename

# 		type <- nonbiomass_facets[[facet]]$type
# 		if (type == 'morphological') {
# 			multi_db <- multi_morph_db
# 			multi_sq <- multi_morph_sq
# 		} else {
# 			multi_db <- multi_phys_db
# 			multi_sq <- multi_phys_sq
# 		}

# 		# response vectors
# 		uni_db <- vect(paste0('./outputs_loretta/integrated_sdm_pdm/models_', facet, '/[', facet, '~', resp_distrib, '(', form_filename, ')]/prediction_vector_conus_1930s.gpkg'))
# 		uni_sq <- vect(paste0('./outputs_loretta/integrated_sdm_pdm/models_', facet, '/[', facet, '~', resp_distrib, '(', form_filename, ')]/prediction_vector_nam.gpkg'))

# 		uni_db <- uni_db[uni_db$state_province %in% c('Colorado', 'Nebraska', 'Kansas', 'Oklahoma', 'Texas', 'New Mexico')]
# 		uni_sq <- uni_sq[uni_sq$state_province %in% c('Colorado', 'Nebraska', 'Kansas', 'Oklahoma', 'Texas', 'New Mexico')]

# 		multi_db <- multi_db[multi_db$state_province %in% c('Colorado', 'Nebraska', 'Kansas', 'Oklahoma', 'Texas', 'New Mexico')]
# 		multi_sq <- multi_sq[multi_sq$state_province %in% c('Colorado', 'Nebraska', 'Kansas', 'Oklahoma', 'Texas', 'New Mexico')]

# 		uni_db <- project(uni_db, nam)
# 		uni_sq <- project(uni_sq, nam)
# 		multi_db <- project(multi_db, nam)
# 		multi_sq <- project(multi_sq, nam)

# 		# state labels at state centroids
# 		centroids <- crds(centroids(nam))
# 		centroids <- as.data.frame(centroids)
# 		centroids$y[nam$NAME_1 == 'Oklahoma'] <- centroids$y[nam$NAME_1 == 'Oklahoma'] - 90000

# 		# change and response limits
# 		uni_sq$delta_y <- (uni_db[[paste0(facet, '_mean_1930s')]] - uni_sq[[paste0(facet, '_mean_sq')]]) / uni_sq[[paste0(facet, '_mean_sq')]]
# 		uni_sq$delta_psi <- uni_db$psi_1930s - uni_sq$psi_sq
		
# 		multi_sq$delta_y <- (multi_db[[paste0(facet, '_mean_1930s')]] - multi_sq[[paste0(facet, '_mean_sq')]]) / multi_sq[[paste0(facet, '_mean_sq')]]
# 		multi_sq$delta_psi <- multi_db$psi_1930s - multi_sq$psi_sq

# 		# crop to display area
# 		nam <- crop(nam, extent)
# 		uni_db <- crop(uni_db, extent)
# 		uni_sq <- crop(uni_sq, extent)
# 		multi_db <- crop(multi_db, extent)
# 		multi_sq <- crop(multi_sq, extent)

# 		# response limits for color scales
# 		min_val <- min(uni_sq$delta_y, multi_sq$delta_y)
# 		max_val <- max(uni_sq$delta_y, multi_sq$delta_y)
# 		resp_limits_y <- c(min_val, max_val)

# 		min_val <- min(uni_sq$delta_psi, multi_sq$delta_psi)
# 		max_val <- max(uni_sq$delta_psi, multi_sq$delta_psi)
# 		resp_limits_psi <- c(min_val, max_val)

# 		# state labels at state centroids
# 		centroids <- crds(centroids(nam))
# 		centroids <- as.data.frame(centroids)
# 		centroids$y[nam$NAME_1 == 'Oklahoma'] <- centroids$y[nam$NAME_1 == 'Oklahoma'] - 90000

# 		# maps
# 		plot_title_size <- 14
# 		state_labels_size <- 4
		
# 		nice <- get_nice_trait(facet)
		
# 		title <- paste0(nice$short, ': ', nice$short, '-only Model')
# 		map_uni_y <- ggplot() +
# 			layer_spatial(nam, fill = 'gainsboro', color = NA) +
# 			layer_spatial(
# 				uni_sq,
# 				aes(fill = delta_y),
# 				color = NA
# 			) +
# 			scale_fill_gradient2(
# 				name = 'Change',
# 				low = '#c51b7d',
# 				mid = '#f7f7f7',
# 				high = '#4d9221',
# 				midpoint = 0,
# 				limits = resp_limits_y,
# 				labels = scales::percent_format(accuracy = 1)
# 			) +
# 			layer_spatial(nam, color = 'gray40', fill = NA) +
# 			geom_text(data = centroids, aes(x = x, y = y, label = nam$NAME_1), size = state_labels_size) +
# 			layer_spatial(dust_bowl, color = 'gray20', fill = NA, linewidth = 1) +
# 			layer_spatial(site_vect, pch = 3, size = 4) +
# 			annotation_scale(location = 'br') +
# 			coord_sf(xlim = c(extent[1], extent[2]), ylim = c(extent[3], extent[4]), expand = FALSE) +
# 			ggtitle(title) +
# 			theme(
# 				plot.title = element_text(size = plot_title_size),
# 				axis.title = element_blank()
# 			)

# 		title <- paste0('Occupancy: ', nice$short, '-only Model')
# 		map_uni_psi <- ggplot() +
# 			layer_spatial(nam, fill = 'gainsboro', color = NA) +
# 			layer_spatial(
# 				uni_sq,
# 				aes(fill = delta_psi),
# 				color = NA
# 			) +
# 			scale_fill_gradient2(
# 				name = 'Change',
# 				low = '#c51b7d',
# 				mid = '#f7f7f7',
# 				high = '#4d9221',
# 				midpoint = 0,
# 				limits = resp_limits_psi,
# 				labels = scales::percent_format(accuracy = 1)
# 			) +
# 			layer_spatial(nam, color = 'gray40', fill = NA) +
# 			geom_text(data = centroids, aes(x = x, y = y, label = nam$NAME_1), size = state_labels_size) +
# 			layer_spatial(dust_bowl, color = 'gray20', fill = NA, linewidth = 1) +
# 			layer_spatial(site_vect, pch = 3, size = 4) +
# 			annotation_scale(location = 'br') +
# 			coord_sf(xlim = c(extent[1], extent[2]), ylim = c(extent[3], extent[4]), expand = FALSE) +
# 			ggtitle(title) +
# 			theme(
# 				plot.title = element_text(size = plot_title_size),
# 				axis.title = element_blank()
# 			)

# 		title <- paste0(nice$short, ': Integrated Model')
# 		map_multi_y <- ggplot() +
# 			layer_spatial(nam, fill = 'gainsboro', color = NA) +
# 			layer_spatial(
# 				multi_sq,
# 				aes(fill = delta_y),
# 				color = NA
# 			) +
# 			scale_fill_gradient2(
# 				name = 'Change',
# 				low = '#c51b7d',
# 				mid = '#f7f7f7',
# 				high = '#4d9221',
# 				midpoint = 0,
# 				limits = resp_limits_y,
# 				labels = scales::percent_format(accuracy = 1)
# 			) +
# 			layer_spatial(nam, color = 'gray40', fill = NA) +
# 			geom_text(data = centroids, aes(x = x, y = y, label = nam$NAME_1), size = state_labels_size) +
# 			layer_spatial(dust_bowl, color = 'gray20', fill = NA, linewidth = 1) +
# 			layer_spatial(site_vect, pch = 3, size = 4) +
# 			annotation_scale(location = 'br') +
# 			coord_sf(xlim = c(extent[1], extent[2]), ylim = c(extent[3], extent[4]), expand = FALSE) +
# 			ggtitle(title) +
# 			theme(
# 				plot.title = element_text(size = plot_title_size),
# 				axis.title = element_blank()
# 			)

# 		title <- paste0('Occupancy: Integrated Model')
# 		map_multi_psi <- ggplot() +
# 			layer_spatial(nam, fill = 'gainsboro', color = NA) +
# 			layer_spatial(
# 				multi_sq,
# 				aes(fill = delta_psi),
# 				color = NA
# 			) +
# 			scale_fill_gradient2(
# 				name = 'Change',
# 				low = '#c51b7d',
# 				mid = '#f7f7f7',
# 				high = '#4d9221',
# 				midpoint = 0,
# 				limits = resp_limits_psi,
# 				labels = scales::percent_format(accuracy = 1)
# 			) +
# 			layer_spatial(nam, color = 'gray40', fill = NA) +
# 			geom_text(data = centroids, aes(x = x, y = y, label = nam$NAME_1), size = state_labels_size) +
# 			layer_spatial(dust_bowl, color = 'gray20', fill = NA, linewidth = 1) +
# 			layer_spatial(site_vect, pch = 3, size = 4) +
# 			annotation_scale(location = 'br') +
# 			coord_sf(xlim = c(extent[1], extent[2]), ylim = c(extent[3], extent[4]), expand = FALSE) +
# 			ggtitle(title) +
# 			theme(
# 				plot.title = element_text(size = plot_title_size),
# 				axis.title = element_blank()
# 			)

# 		map_grobs <- list(map_uni_y, map_multi_y, map_uni_psi, map_multi_psi)
# 		names(map_grobs) <- c('delta_univariate', 'delta_multivariate', 'delta_psi_univariate', 'delta_psi_multivariate')
# 		saveRDS(map_grobs, file = paste0('./outputs_loretta/integrated_sdm_pdm/maps_1930s_dust_bowl_', facet, '_change.rds'))

# 		maps <- plot_grid(plotlist = map_grobs, ncol = 2, align = 'hv')
# 		ggsave(maps, filename = paste0('./outputs_loretta/integrated_sdm_pdm/maps_1930s_dust_bowl_', facet, '_change.png'), width = 10, height = 11, dpi = 600, bg = 'white')

# 	} # next facet

# say('#################################################################')
# say('### maps of change in Dust Bowl region in 1930s for main text ###')
# say('#################################################################')

# 	abund <- readRDS('./outputs_loretta/integrated_sdm_pdm/maps_1930s_dust_bowl_abundance_change.rds')
# 	biomass <- readRDS('./outputs_loretta/integrated_sdm_pdm/maps_1930s_dust_bowl_biomass_change.rds')
# 	morph <- readRDS('./outputs_loretta/integrated_sdm_pdm/maps_1930s_dust_bowl_height_change.rds')
# 	phys <- readRDS('./outputs_loretta/integrated_sdm_pdm/maps_1930s_dust_bowl_photosynthetic_rate_change.rds')

# 	maps <- list(
# 		abund$delta_univariate,
# 		abund$delta_multivariate,
# 		abund$delta_psi_univariate,
# 		abund$delta_psi_multivariate,
# 		biomass$delta_univariate,
# 		biomass$delta_multivariate,
# 		morph$delta_univariate,
# 		morph$delta_multivariate,
# 		phys$delta_univariate,
# 		phys$delta_multivariate
# 	)
	
# 	for (i in seq_along(maps)) {
		
# 		maps[[i]] <- maps[[i]] +
# 			theme(
# 				plot.title = element_text(size = 7),
# 				axis.text = element_text(size = 4),
# 				legend.title = element_text(size = 7),
# 				legend.text = element_text(size = 6),
# 				legend.key.width = grid::unit(0.24, 'cm')
# 			)

# 		# change size of state names
# 		maps[[i]]$layers$geom_text$aes_params$size <- 2

# 		# change size of sample site points
# 		maps[[i]]$layers[[6]]$aes_params$size <- 1.6

# 		# change size of Dust Bowl region outline
# 		maps[[i]]$layers[[5]]$aes_params$linewidth <- 0.4

# 	}

# 	maps[[1]] <- maps[[1]] + ggtitle('a) Abundance:\nAbundance-only Model')
# 	maps[[2]] <- maps[[2]] + ggtitle('b) Abundance:\nIntegrated Model')

# 	maps[[3]] <- maps[[3]] + ggtitle('c) Occupancy:\nAbundance-only Model')
# 	maps[[4]] <- maps[[4]] + ggtitle('d) Occupancy:\nIntegrated Model')

# 	maps[[5]] <- maps[[5]] + ggtitle('e) Biomass:\nBiomass-only Model')
# 	maps[[6]] <- maps[[6]] + ggtitle('f) Biomass:\nIntegrated Model')

# 	maps[[7]] <- maps[[7]] + ggtitle('g) Height:\nHeight-only Model')
# 	maps[[8]] <- maps[[8]] + ggtitle('h) Height:\nIntegrated Model')

# 	maps[[9]] <- maps[[9]] + ggtitle('i) Photosynthetic Rate:\nPhotosynthetic Rate-only Model')
# 	maps[[10]] <- maps[[10]] + ggtitle('j) Photosynthetic Rate:\nIntegrated Model')

# 	map_grid <- plot_grid(plotlist = maps, ncol = 2, align = 'hv')
	
# 	ggsave(map_grid, filename = './outputs_loretta/integrated_sdm_pdm/maps_1930s_dust_bowl_change.png', width = 4.25, height = 9.5, dpi = 600, bg = 'white')

say(date())
say('FINIS!', deco = '+', level = 1)

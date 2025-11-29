# source('C:/Kaji/R/andropogon_integratedEcology/sdm_pdm_for_sem/TEMP.r')

	rm(list = ls())

	drive <- 'C:/Kaji/'

	setwd(paste0(drive, '/Research/Andropogon/Andropogon'))
	source(paste0(drive, '/R/andropogon_integratedEcology/sdm_pdm_for_sem/sdm_pdm_for_sem_00_shared_functions_and_variables.r'))


	# trial <- TRUE # TRUE for testing
	trial <- FALSE # TRUE for testing

	# do cross-validation?
	# crossvalidate <- FALSE
	crossvalidate <- TRUE

	### formula for how aspects of species responds to environment

	# formula_biomass_mu <- ~ 1 + bio12 # response of biomass to environment
	# preds_filename <- 'bio12'

	# formula_biomass_mu <- ~ 1 + bio12 + I(bio12^2) # response of biomass to environment
	# preds_filename <- 'bio12^2'

	# formula_biomass_mu <- ~ 1 + bio1 # response of biomass to environment
	# preds_filename <- 'bio1'

	# formula_biomass_mu <- ~ 1 + bio1 + I(bio1^2) # response of biomass to environment
	# preds_filename <- 'bio1^2'

	# formula_biomass_mu <- ~ 1 + bio1 + bio12 # response of biomass to environment
	# preds_filename <- 'bio1_bio12'

	# formula_biomass_mu <- ~ 1 + bio1 + bio12 + bio1:bio12 # response of biomass to environment
	# preds_filename <- 'bio1_x_bio12'

	# formula_biomass_mu <- ~ 1 + bio1 + bio12 + I(bio1^2) # response of biomass to environment
	# preds_filename <- 'bio1^2_bio12'

	# formula_biomass_mu <- ~ 1 + bio1 + bio12 + I(bio12^2) # response of biomass to environment
	# preds_filename <- 'bio1_bio12^2'

	# formula_biomass_mu <- ~ 1 + bio1 + bio12 + I(bio1^2) + I(bio12^2) # response of biomass to environment
	# preds_filename <- 'bio1^2_bio12^2'

	formula_biomass_mu <- ~ 1 + bio1 + bio12 + I(bio1^2) + bio1:bio12 # response of biomass to environment
	preds_filename <- 'bio1^2_x_bio12'

	# out_dir <- paste0('./outputs_loretta/sdm_pdm_for_sem/models_biomass/', ifelse(trial, 'TRIAL_', ''), '[biomass_gamma~normal_homoscedastic_', preds_filename, ']/')
	# zero_inflated <- FALSE # SPECIFIC TO THIS SCRIPT--SHOULD NOT BE CHANGED
	# formula_pzero <- NULL

	out_dir <- paste0('./outputs_loretta/sdm_pdm_for_sem/models_biomass/', ifelse(trial, 'TRIAL_', ''), '[biomass_zig~normal_homoscedastic_', preds_filename, ']/')
	zero_inflated <- TRUE # SPECIFIC TO THIS SCRIPT--SHOULD NOT BE CHANGED
	formula_pzero <- formula_biomass_mu

	# calib <- TRUE # use just counties with non-NA Poaceae for calibration region
	calib <- FALSE # use all of North America for calibration region

	homoscedastic <- TRUE # SPECIFIC TO THIS SCRIPT--SHOULD NOT BE CHANGED
	formula_biomass_sigma <- NULL

	data_biomass_mu <- prepare_biomass(formula_biomass = formula_biomass_mu, n_response_curve_values = n_response_curve_values, calib = calib)
	data_occs <- prepare_occurrences(formula_occs = ~ 1, formula_occs_bias = ~ 1, n_response_curve_values = n_response_curve_values, psa_quant = psa_quant, calib = calib)

	chains <- readRDS(paste0(out_dir, '/chains.rds'))

	pred_vect_nam <- burn_biomass_into_vector(demesne = 'nam', chains = chains, formula_biomass_mu = formula_biomass_mu, formula_biomass_sigma = formula_biomass_sigma, formula_pzero = formula_pzero)
	pred_vect_1930s <- burn_biomass_into_vector(demesne = '1930s', chains = chains, formula_biomass_mu = formula_biomass_mu, formula_biomass_sigma = formula_biomass_sigma, formula_pzero = formula_pzero)

	### BIOMASS: current map
	########################
	say('BIOMASS: current map', level = 2)

		mean_form <- paste(as.character(formula_biomass_mu), collapse = ' ')
		mean_form <- gsub(mean_form, pattern = 'I\\(', replacement = '')
		mean_form <- gsub(mean_form, pattern = '\\^2\\)', replacement = '²')
		mean_form <- gsub(mean_form, pattern = '*)', replacement = '×')

		if (homoscedastic) {
			subtitle <- paste0('Biomass μ ', mean_form, ' | 1991-2020')
		} else {
			sd_form <- paste(as.character(formula_biomass_sigma), collapse = ' ')
			subtitle <- paste0('Biomass μ ', mean_form, ' & biomass σ ', sd_form, ' | 1991-2020')
		}

		map <- map_biomass(
			out_dir = out_dir,
			filename_append = 'present',
			pred_vect_nam = pred_vect_nam,
			response_var = 'mu_biomass_county_mean_sq',
			response_var_type = 'mu',
			data_occs = data_occs,
			data_biomass = data_biomass_mu,
			title = bquote('Present-day distribution of mean ' * italic('Andropogon gerardi') * ' biomass '),
			subtitle = subtitle,
			legend_title = 'Biomass (g)',
			plot_range_core = TRUE,
			# plot_range_core = FALSE,
			ag_core_quant = 0.95
		)
print(NON)
		if (!is.null(formula_pzero)) {
			
			facet <- 'biomass'
			map <- map_pzero(
				out_dir = out_dir,
				filename_append = 'present',
				pred_vect_nam = pred_vect_nam,
				facet = facet,
				response_var = paste0('pzero_', facet, '_county_sq'),
				data_traits = data_biomass_mu,
				title = bquote('Present-day probability of zero biomass '),
				subtitle = subtitle,
				plot_range_core = TRUE,
				ag_core_quant = 0.95
			)

		}

	### BIOMASS: future maps
	########################
	say('BIOMASS: future maps', level = 2)

		for (fut in futs) {

			say(fut)

			response_var <- paste0('mu_biomass_county_mean_', fut)

			mean_form <- paste(as.character(formula_biomass_mu), collapse = ' ')
			if (homoscedastic) {
				subtitle <- paste0('Biomass μ ', mean_form, '  | SSP ', substr(fut, 4, 6), ' ', substr(fut, 8, 11), '-', substr(fut, 13, 16))
			} else {

				sd_form <- paste(as.character(formula_biomass_sigma), collapse = ' ')
				subtitle <- paste0('Biomass μ ', mean_form, ' & biomass σ ', sd_form, ' | SSP ', substr(fut, 4, 6), ' ', substr(fut, 8, 11), '-', substr(fut, 13, 16))
			}

			map_biomass_fut <- map_biomass(
				out_dir = out_dir,
				filename_append = fut,
				pred_vect_nam = pred_vect_nam,
				response_var = response_var,
				response_var_type = 'mu',
				data_occs = data_occs,
				data_biomass = data_biomass_mu,
				title = bquote('Future distribution of mean ' * italic('Andropogon gerardi') * ' biomass '),
				subtitle = subtitle,
				legend_title = 'Biomass (g)',
				plot_range_core = TRUE,
				ag_core_quant = 0.95
			)

			if (!is.null(formula_pzero)) {

				face <- 'biomass'				
				map_fut_pzero <- map_pzero(
					out_dir = out_dir,
					filename_append = fut,
					pred_vect_nam = pred_vect_nam,
					facet = facet,
					response_var = paste0('pzero_', facet, '_county_', fut),
					data_traits = data_biomass_mu,
					title = bquote('Future probability of zero biomass '),
					subtitle = subtitle,
					plot_range_core = TRUE,
					ag_core_quant = 0.95
				)

			}

		}

	# ### BIOMASS: future change maps
	# ###############################
	# say('BIOMASS: future change maps', level = 2)

	# 	maps_biomass_change <- list()
	# 	for (fut in futs) {

	# 		say(fut)

	# 		response_var <- paste0('mu_biomass_county_mean_', fut)
	# 		title <- bquote('Change in ' * italic('Andropogon gerardi') * ' biomass ')

	# 		mean_form <- paste(as.character(formula_biomass_mu), collapse = ' ')
	# 		if (homoscedastic) {
	# 			subtitle <- paste0('Biomass μ ', mean_form, '  | SSP ', substr(fut, 4, 6), ' ', substr(fut, 8, 11), '-', substr(fut, 13, 16))
	# 		} else {

	# 			sd_form <- paste(as.character(formula_biomass_sigma), collapse = ' ')
	# 			subtitle <- paste0('Biomass μ ', mean_form, ' & biomass σ ', sd_form, ' | SSP ', substr(fut, 4, 6), ' ', substr(fut, 8, 11), '-', substr(fut, 13, 16))
	# 		}

	# 		maps_biomass_change[[length(maps_biomass_change) + 1]] <- map_biomass_change(
	# 			out_dir = out_dir,
	# 			filename_append = fut,
	# 			fut = fut,
	# 			pred_vect_nam = pred_vect_nam,
	# 			response_var = response_var,
	# 			data_occs = data_occs,
	# 			data_biomass = data_biomass_mu,
	# 			title = title,
	# 			subtitle = subtitle,
	# 			plot_range_core = TRUE,
	# 			ag_core_quant = ag_core_quant
	# 		)

	# 		if (!is.null(formula_pzero)) {

	# 			title <- bquote('Change in probability of zero biomass ')

	# 			facet <- 'biomass'				
	# 			map_pzero_change <- map_pzero_change(
	# 				out_dir = out_dir,
	# 				filename_append = fut,
	# 				fut = fut,
	# 				pred_vect_nam = pred_vect_nam,
	# 				facet = facet,
	# 				response_var = paste0('pzero_', facet, '_county_', fut),
	# 				response_var_sq = paste0('pzero_', facet, '_county_sq'),
	# 				data_traits = data_biomass_mu,
	# 				title = title,
	# 				subtitle = subtitle,
	# 				plot_range_core = TRUE,
	# 				ag_core_quant = 0.95
	# 			)

	# 		}

	# 	}

	### OCCURRENCE: 1930s change maps
	#################################
	say('BIOMASS: 1930s change maps', level = 2)

	map_thirties <- map_biomass_traits_change_1930s(facet = 'biomass', data_biomass_traits = data_biomass_mu, pred_vect_nam = pred_vect_nam, pred_vect_1930s = pred_vect_1930s)



#' Burn predictions from an occurrence model into a SpatVector of North America or Dust Bowl map
#'
#' demesne					'nam' (North America) or '1930s' (Dust Bowl area) or '1950s' (post-Dust Bowl)
#' chains					Chains from NIMBLE
#' formula_occs_mu			Formula for occurrences
#' formula_occs_sigma		Formula for occurrences s.d. or `NULL`
#' formula_biomas_mu		Formula for biomass mean
#' formula_biomas_sigma		Formula for biomass sigma
#' formula_pzero			Formula for inflated zero or `NULL`
#' log_precip				Logical--if `TRUE`, log bios 12-14 and 16-19
burn_occs_biomass_into_vector <- function(demesne, chains, formula_occs, formula_occs_sigma, formula_biomass_mu, formula_biomass_sigma, formula_pzero, log_precip) {

	homoscedastic_occs <- is.null(formula_occs_sigma)
	homoscedastic_biomass <- is.null(formula_biomass_sigma)
	zero_inflated <- !is.null(formula_pzero)

	### occurrence data
	###################

	data_occs <- prepare_occurrences(formula_occs = formula_occs, formula_occs_bias = ~ 1, log_precip = log_precip, n_response_curve_values = n_response_curve_values, psa_quant = psa_quant, calib = calib)

	data_biomass <- prepare_occurrences(formula_occs = formula_biomass_mu, formula_occs_bias = ~ 1, log_precip = log_precip, n_response_curve_values = n_response_curve_values, psa_quant = psa_quant, calib = calib)

	if (zero_inflated) data_occs_pzero <- prepare_occurrences(formula_occs = formula_pzero, formula_occs_bias = ~ 1, log_precip = loc_precip, n_response_curve_values = n_response_curve_values, psa_quant = psa_quant, calib = calib)

	pred_vect <- if (demesne == 'nam') {
		pred_vect <- data_occs$ag_vect_sq
	} else if (demesne == '1930s') {
		pred_vect <- vect(paste0('./data_from_adam_and_loretta/andropogon_gerardi_occurrences_with_environment_1931_1940_prism.gpkg'))
	} else if (demesne == '1950s') {
		pred_vect <- vect(paste0('./data_from_adam_and_loretta/andropogon_gerardi_occurrences_with_environment_1952_1961_prism.gpkg'))
	}

	### predict mu
	##############

	x_occs <- if (demesne == 'nam') {
		data_occs$counties_x_occs_sq
	} else if (demesne == '1930s') {
		data_occs$counties_x_occs_thirties
	} else if (demesne == '1950s') {
		data_occs$counties_x_occs_fifties
	}

	x_biomass <- if (demesne == 'nam') {
		data_biomass$counties_x_occs_sq
	} else if (demesne == '1930s') {
		data_biomass$counties_x_occs_thirties
	} else if (demesne == '1950s') {
		data_biomass$counties_x_occs_fifties
	}

	x_pzero <- if (!zero_inflated) {
		NULL
	} else if (demesne == 'nam') {
		data_occs_pzero$counties_x_occs_sq
	} else if (demesne == '1930s') {
		data_occs_pzero$counties_x_occs_thirties
	} else if (demesne == '1950s') {
		data_occs_pzero$counties_x_occs_fifties
	}

	preds <- predict_occs_biomass(chains = chains, x_occs = x_occs, x_biomass = x_biomass, homoscedastic_occs = homoscedastic_occs, homoscedastic_biomass = homoscedastic_biomass, zero_inflated = zero_inflated, type = 'mu', x_pzero = x_pzero)

	preds_occs <- preds$preds_occs
	preds_biomass <- preds$preds_biomass

	preds_occs_mean <- colMeans(preds_occs)
	preds_occs_sd <- apply(preds_occs, 2, sd)

	preds_biomass_mean <- colMeans(preds_biomass)
	preds_biomass_sd <- apply(preds_biomass, 2, sd)

	pred_vect$DUMMY1 <- preds_occs_mean
	pred_vect$DUMMY2 <- preds_occs_sd

	pred_vect$DUMMY3 <- preds_biomass_mean
	pred_vect$DUMMY4 <- preds_biomass_sd

	index <- (ncol(pred_vect) - 3):ncol(pred_vect)
	if (demesne == 'nam') {
		names(pred_vect)[index] <- c('N_ag_county_mean_sq', 'N_ag_county_sd_sq', 'mu_biomass_county_mean_sq', 'mu_biomass_county_sd_sq')
	} else if (demesne == '1930s') {
		names(pred_vect)[index] <- c('N_ag_county_mean_1930s', 'N_ag_county_sd_1930s', 'mu_biomass_county_mean_1930s', 'mu_biomass_county_sd_1930s')
	} else if (demesne == '1950s') {
		names(pred_vect)[index] <- c('N_ag_county_mean_1950s', 'N_ag_county_sd_1950s', 'mu_biomass_county_mean_1930s', 'mu_biomass_county_sd_1930s')
	}

	# ### predict sigma
	# #################

	# if (!is.null(formula_occs_sigma)) {
	
	# 	preds <- predict_biomass(chains = chains, x = x, homoscedastic = homoscedastic, type = 'sigma')
	# 	# preds_mean <- colMeans(preds)
	# 	preds_median <- apply(preds, 2, median)
	# 	# pred_vect$DUMMY1 <- preds_mean
	# 	pred_vect$DUMMY1 <- preds_median

	# 	index <- ncol(pred_vect)
	# 	if (demesne == 'nam') {
	# 		names(pred_vect)[index] <- c('sigma_occs_county_median_sq')
	# 	} else if (demesne == '1930s') {
	# 		names(pred_vect)[index] <- c('sigma_occs_county_median_1930s')
	# 	} else if (demesne == '1950s') {
	# 		names(pred_vect)[index] <- c('sigma_occs_county_median_1950s')
	# 	}

	# }

	### predict probability of zero
	###############################

	if (!is.null(formula_pzero)) {
		
		x <- if (demesne == 'nam') {
			data_occs_pzero$counties_x_occs_sq
		} else if (demesne == '1930s') {
			data_occs_pzero$counties_x_occs_thirties
		} else if (demesne == '1950s') {
			data_occs_pzero$counties_x_occs_fifties
		}

		preds <- predict_pzero(chains = chains, x = x)
		preds_mean <- colMeans(preds)
		pred_vect$DUMMY1 <- preds_mean

		index <- ncol(pred_vect)
		if (demesne == 'nam') {
			names(pred_vect)[index] <- c('pzero_county_sq')
		} else if (demesne == '1930s') {
			names(pred_vect)[index] <- c('pzero_county_1930s')
		} else if (demesne == '1950s') {
			names(pred_vect)[index] <- c('pzero_county_1950s')
		}

	}

	### future occurrences--only for entire North American domain
	##############################################################
	if (demesne == 'nam') {

		### predict mean to each future
		###############################
		for (fut in futs) {
			
			x_occs <- data_occs[paste0('counties_x_occs_', fut)]
			x_occs <- x_occs[[1]]
		
			x_biomass <- data_biomass[paste0('counties_x_occs_', fut)]
			x_biomass <- x_biomass[[1]]
		
			 if (!zero_inflated) {
				x_pzero <- NULL
			} else {
				x_pzero <- data_occs_pzero[paste0('counties_x_occs_', fut)]
				x_pzero <- x_pzero[[1]]
			}
		
			preds <- predict_occs_biomass(chains = chains, x_occs = x_occs, x_biomass = x_biomass, homoscedastic_occs = homoscedastic_occs, homoscedastic_biomass = homoscedastic_biomass, zero_inflated = zero_inflated, type = 'mu', x_pzero = x_pzero)

			preds_occs <- preds$preds_occs
			preds_biomass <- preds$preds_biomass

			preds_occs_mean <- colMeans(preds_occs)
			preds_occs_sd <- apply(preds_occs, 2, sd)

			preds_biomass_mean <- colMeans(preds_biomass)
			preds_biomass_sd <- apply(preds_biomass, 2, sd)

			pred_vect$DUMMY1 <- preds_occs_mean
			pred_vect$DUMMY2 <- preds_occs_sd

			pred_vect$DUMMY3 <- preds_biomass_mean
			pred_vect$DUMMY4 <- preds_biomass_sd

			index <- (ncol(pred_vect) - 3):ncol(pred_vect)
			names(pred_vect)[index] <- c(
				paste0('N_ag_county_mean_', fut),
				paste0('N_ag_county_sd_', fut),
				paste0('mu_biomass_county_mean_', fut),
				paste0('mu_biomass_county_sd_', fut)
			)

		} # next future

		# ### predict sd of site-level abundance to each future--only if heteroscedastic model
		# ##################################################################################
		# if (!is.null(formula_occs_sigma)) {
		
		# 	data_occs_sigma <- prepare_biomass(formula_biomass = formula_occs_sigma, calib = calib)

		# 	# predict to each future
		# 	for (fut in futs) {
				
		# 		x <- data_occs_sigma[[paste0('counties_x_occs_', fut)]]
		# 		preds <- predict_biomass(chains = chains, x = x, homoscedastic = homoscedastic, type = 'sigma')

		# 		# preds_mean <- colMeans(preds)
		# 		preds_median <- apply(preds, 2, median)
		# 		preds_sd <- apply(preds, 2, sd)

		# 		# pred_vect$DUMMY1 <- preds_mean
		# 		pred_vect$DUMMY1 <- preds_median
		# 		# pred_vect$DUMMY2 <- preds_sd

		# 		# index <- (ncol(pred_vect) - 1):ncol(pred_vect)
		# 		index <- ncol(pred_vect)
		# 		# names(pred_vect)[index] <- paste0(c('sigma_occs_county_median_', 'sigma_occs_county_sd_'), fut)
		# 		names(pred_vect)[index] <- paste0('sigma_occs_county_median_', fut)

		# 	} # next future

		# } # if heteroscedastic

		### predict probability of zero AG to future
		###################################################

		if (!is.null(formula_pzero)) {
		
			for (fut in futs) {

				x <- data_occs_pzero[paste0('counties_x_occs_', fut)]
				x <- x[[1]]

				preds <- predict_pzero(chains = chains, x = x)
				preds_mean <- colMeans(preds)
				pred_vect$DUMMY1 <- preds_mean

				index <- ncol(pred_vect)
				names(pred_vect)[index] <- paste0('pzero_county_', fut)

			} # next future

		} # if zero-inflated

	} # if demesne is North America

	pred_vect

}

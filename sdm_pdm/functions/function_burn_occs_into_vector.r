#' Burn predictions from an occurrence model into a SpatVector of North America or Dust Bowl map
#'
#' demesne					'nam' (North America) or '1930s' (Dust Bowl area) or '1950s' (post-Dust Bowl)
#' chains					Chains from NIMBLE
#' formula_occs_mu			Formula for occurrences
#' formula_occs_sigma		Formula for occurrences s.d. or `NULL`
#' formula_occs_pzero		Formula for occurrence probability of inflated zero or `NULL`
burn_occs_into_vector <- function(demesne, chains, formula_occs, formula_occs_sigma, formula_occs_pzero) {

	homoscedastic <- is.null(formula_occs_sigma)
	zero_inflated <- !is.null(formula_occs_pzero)

	### occurrence data
	data_occs <- prepare_occurrences(formula_occs = formula_occs, formula_occs_bias = ~ 1, n_response_curve_values = n_response_curve_values, psa_quant = psa_quant, calib = calib)

	if (zero_inflated) data_occs_pzero <- prepare_occurrences(formula_occs = formula_occs_pzero, formula_occs_bias = ~ 1, n_response_curve_values = n_response_curve_values, psa_quant = psa_quant, calib = calib)

	pred_vect <- if (demesne == 'nam') {
		pred_vect <- data_occs$ag_vect_sq
	} else if (demesne == '1930s') {
		pred_vect <- vect(paste0('./data_from_adam_and_loretta/andropogon_gerardi_occurrences_with_environment_1931_1940_prism.gpkg'))
	} else if (demesne == '1950s') {
		pred_vect <- vect(paste0('./data_from_adam_and_loretta/andropogon_gerardi_occurrences_with_environment_1952_1961_prism.gpkg'))
	}

	### predict mu
	##############

	x <- if (demesne == 'nam') {
		data_occs$counties_x_occs_sq
	} else if (demesne == '1930s') {
		data_occs$counties_x_occs_thirties
	} else if (demesne == '1950s') {
		data_occs$counties_x_occs_fifties
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

	preds <- predict_occs(chains = chains, x = x, homoscedastic = homoscedastic, zero_inflated = zero_inflated, type = 'mu', x_pzero = x_pzero)

	preds_mean <- colMeans(preds)
	preds_sd <- apply(preds, 2, sd)

	pred_vect$DUMMY1 <- preds_mean
	pred_vect$DUMMY2 <- preds_sd

	index <- (ncol(pred_vect) - 1):ncol(pred_vect)
	if (demesne == 'nam') {
		names(pred_vect)[index] <- c('N_ag_county_mean_sq', 'N_ag_county_sd_sq')
	} else if (demesne == '1930s') {
		names(pred_vect)[index] <- c('N_ag_county_mean_1930s', 'N_ag_county_sd_1930s')
	} else if (demesne == '1950s') {
		names(pred_vect)[index] <- c('N_ag_county_mean_1950s', 'N_ag_county_sd_1950s')
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

	if (!is.null(formula_occs_pzero)) {
		
		x <- if (demesne == 'nam') {
			data_occs_pzero$counties_x_occs_sq
		} else if (demesne == '1930s') {
			data_occs_pzero$counties_x_occs_thirties
		} else if (demesne == '1950s') {
			data_occs_pzero$counties_x_occs_fifties
		}

		preds <- predict_occs(chains = chains, x = x, homoscedastic = homoscedastic, zero_inflated = zero_inflated, type = 'pzero')
		preds_mean <- colMeans(preds)
		pred_vect$DUMMY1 <- preds_mean

		index <- ncol(pred_vect)
		if (demesne == 'nam') {
			names(pred_vect)[index] <- c('pzero_occs_county_sq')
		} else if (demesne == '1930s') {
			names(pred_vect)[index] <- c('pzero_occs_county_1930s')
		} else if (demesne == '1950s') {
			names(pred_vect)[index] <- c('pzero_occs_county_1950s')
		}

	}

	### future occurrences--only for entire North American domain
	##############################################################
	if (demesne == 'nam') {

		### predict mean to each future
		###############################
		for (fut in futs) {
			
			x <- data_occs[paste0('counties_x_occs_', fut)]
			x <- x[[1]]
		
			x_pzero <- if (!zero_inflated) {
				NULL
			} else {
				data_occs_pzero[paste0('counties_x_occs_', fut)]
				x_pzero <- x_pzero[[1]]
			}
		
			preds <- predict_occs(chains = chains, x = x, homoscedastic = homoscedastic, zero_inflated = zero_inflated, type = 'mu', x_pzero = x_pzero)

			preds_mean <- colMeans(preds)
			preds_sd <- apply(preds, 2, sd)

			pred_vect$DUMMY1 <- preds_mean
			pred_vect$DUMMY2 <- preds_sd

			index <- (ncol(pred_vect) - 1):ncol(pred_vect)
			names(pred_vect)[index] <- paste0(c('N_ag_county_mean_', 'N_ag_county_sd_'), fut)

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

		### predict probability of zero abundance to future
		###################################################

		if (!is.null(formula_occs_pzero)) {
		
			for (fut in futs) {

				x <- data_occs_pzero[paste0('counties_x_occs_', fut)]
				x <- x[[1]]

				preds <- predict_occs(chains = chains, homoscedastic = homoscedastic, zero_inflated = zero_inflated, x = x, type = 'pzero')
				preds_mean <- colMeans(preds)
				pred_vect$DUMMY1 <- preds_mean

				index <- ncol(pred_vect)
				names(pred_vect)[index] <- paste0('pzero_occs_county_', fut)

			} # next future

		} # if zero-inflated

	} # if demesne is North America

	pred_vect

}

#' Burn predictions from a biomass model into a SpatVector of North America or Dust Bowl map
#'
#' demesne					'nam' (North America) or '1930s' (Dust Bowl area) or '1950s' (post-Dust Bowl)
#' chains					Chains from NIMBLE
#' formula_biomass_mu		Formula for mean of site-level biomass
#' formula_biomass_sigma	Formula for sd of site-level biomass. Ignored if `NULL`.
#' formula_pzero	Formula for probability of zero biomass. Ignored if `NULL`.
burn_biomass_into_vector <- function(demesne, chains, formula_biomass_mu, formula_biomass_sigma, formula_pzero) {

	homoscedastic <- is.null(formula_biomass_sigma)
	zero_inflated <- is.null(formula_pzero)

	### mean site-level biomass	
	data_biomass_mu <- prepare_biomass(formula_biomass = formula_biomass_mu, calib = calib)
	data_occs <- prepare_occurrences(formula_occs = ~ 1, formula_occs_bias = ~ 1, n_response_curve_values = n_response_curve_values, psa_quant = psa_quant, calib = calib)

	pred_vect <- if (demesne == 'nam') {
		pred_vect <- data_occs$ag_vect_sq
	} else if (demesne == '1930s') {
		pred_vect <- vect(paste0('./data_from_adam_and_loretta/andropogon_gerardi_occurrences_with_environment_1931_1940_prism.gpkg'))
	} else if (demesne == '1950s') {
		pred_vect <- vect(paste0('./data_from_adam_and_loretta/andropogon_gerardi_occurrences_with_environment_1952_1961_prism.gpkg'))
	}

	# predict to present
	x <- if (demesne == 'nam') {
		data_biomass_mu$counties_x_biomass_sq
	} else if (demesne == '1930s') {
		data_biomass_mu$counties_x_biomass_thirties
	} else if (demesne == '1950s') {
		data_biomass_mu$counties_x_biomass_fifties
	}

	### predict mean biomass
	########################

	preds <- predict_biomass(chains = chains, x = x, homoscedastic = homoscedastic, type = 'mu')

	preds_mean <- colMeans(preds)
	preds_sd <- apply(preds, 2, sd)

	pred_vect$DUMMY1 <- preds_mean
	pred_vect$DUMMY2 <- preds_sd

	index <- (ncol(pred_vect) - 1):ncol(pred_vect)
	if (demesne == 'nam') {
		names(pred_vect)[index] <- c('mu_biomass_county_mean_sq', 'mu_biomass_county_sd_sq')
	} else if (demesne == '1930s') {
		names(pred_vect)[index] <- c('mu_biomass_county_mean_1930s', 'mu_biomass_county_sd_1930s')
	} else if (demesne == '1950s') {
		names(pred_vect)[index] <- c('mu_biomass_county_mean_1950s', 'mu_biomass_county_sd_1950s')
	}

	### predict sigma
	#################

	if (!is.null(formula_biomass_sigma)) {
	
		preds <- predict_biomass(chains = chains, x = x, homoscedastic = homoscedastic, type = 'sigma')
		# preds_mean <- colMeans(preds)
		preds_median <- apply(preds, 2, median)
		# pred_vect$DUMMY1 <- preds_mean
		pred_vect$DUMMY1 <- preds_median

		index <- ncol(pred_vect)
		if (demesne == 'nam') {
			names(pred_vect)[index] <- c('sigma_biomass_county_median_sq')
		} else if (demesne == '1930s') {
			names(pred_vect)[index] <- c('sigma_biomass_county_median_1930s')
		} else if (demesne == '1950s') {
			names(pred_vect)[index] <- c('sigma_biomass_county_median_1950s')
		}

	}

	### predict probability of zero biomass
	#######################################

	if (!is.null(formula_pzero)) {
	
		preds <- predict_biomass(chains = chains, x = x, homoscedastic = homoscedastic, type = 'pzero')
		preds_mean <- colMeans(preds)
		pred_vect$DUMMY1 <- preds_mean

		index <- ncol(pred_vect)
		if (demesne == 'nam') {
			names(pred_vect)[index] <- c('pzero_biomass_county_sq')
		} else if (demesne == '1930s') {
			names(pred_vect)[index] <- c('pzero_biomass_county_1930s')
		} else if (demesne == '1950s') {
			names(pred_vect)[index] <- c('pzero_biomass_county_1950s')
		}

	}

	### future mean biomass--only for entire North American domain
	##############################################################
	if (demesne == 'nam') {

		### predict mean to each future
		###############################
		for (fut in futs) {
			
			x <- data_biomass_mu[paste0('counties_x_biomass_', fut)]
			x <- x[[1]]
		
			preds <- predict_biomass(chains = chains, x = x, homoscedastic = homoscedastic, type = 'mu')

			preds_mean <- colMeans(preds)
			preds_sd <- apply(preds, 2, sd)

			pred_vect$DUMMY1 <- preds_mean
			pred_vect$DUMMY2 <- preds_sd

			index <- (ncol(pred_vect) - 1):ncol(pred_vect)
			names(pred_vect)[index] <- paste0(c('mu_biomass_county_mean_', 'mu_biomass_county_sd_'), fut)

		} # next future

		### predict sd of site-level biomass to each future--only if heteroscedastic model
		##################################################################################
		if (!is.null(formula_biomass_sigma)) {
		
			data_biomass_sigma <- prepare_biomass(formula_biomass = formula_biomass_sigma, calib = calib)

			# predict to each future
			for (fut in futs) {
				
				x <- data_biomass_sigma[[paste0('counties_x_biomass_', fut)]]
				preds <- predict_biomass(chains = chains, x = x, homoscedastic = homoscedastic, type = 'sigma')

				# preds_mean <- colMeans(preds)
				preds_median <- apply(preds, 2, median)
				preds_sd <- apply(preds, 2, sd)

				# pred_vect$DUMMY1 <- preds_mean
				pred_vect$DUMMY1 <- preds_median
				# pred_vect$DUMMY2 <- preds_sd

				# index <- (ncol(pred_vect) - 1):ncol(pred_vect)
				index <- ncol(pred_vect)
				# names(pred_vect)[index] <- paste0(c('sigma_biomass_county_median_', 'sigma_biomass_county_sd_'), fut)
				names(pred_vect)[index] <- paste0('sigma_biomass_county_median_', fut)

			} # next future

		} # if heteroscedastic

		### predict probability of zero biomass to future
		#######################################

		if (!is.null(formula_pzero)) {
		
			for (fut in futs) {

				x <- data_biomass_mu[paste0('counties_x_biomass_', fut)]
				x <- x[[1]]

				preds <- predict_biomass(chains = chains, homoscedastic = homoscedastic, x = x, type = 'pzero')
				preds_mean <- colMeans(preds)
				pred_vect$DUMMY1 <- preds_mean

				index <- ncol(pred_vect)
				names(pred_vect)[index] <- paste0('pzero_biomass_county_', fut)

			} # next future

		} # if zero-inflated

	} # if demesne is North America

	pred_vect

}

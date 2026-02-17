#' Burn predictions from a biomass model into a SpatVector of North America or Dust Bowl map
#'
#' demesne					'nam' (North America) or '1930s' (Dust Bowl area) or '1950s' (post-Dust Bowl)
#' chains					Chains from NIMBLE
#' formula_biomass			Formula for mean of site-level biomass
#' formula_psi				Formula for probability of presence. Ignored if `NULL`.
#' resp_distrib 			Response distribution: 'gamma' or 'lognormal' or 'ZIG' or 'ZILN'.
#' transform 				Either 'exponential' or 'identity'.
#' log_precip				If `TRUE`, use log of BIOs 12-14 and 16-19.
#'
#' Returns a SpatVector.
burn_biomass_into_vector <- function(demesne, chains, formula_biomass, formula_psi, resp_distrib, transform, log_precip) {

	zero_inflated <- !is.null(formula_psi)

	### mean site-level biomass	
	data_biomass <- prepare_biomass_data(formula_biomass = formula_biomass, log_precip = log_precip, calib = calib)
	data_occs <- prepare_occurrence_data(formula_occs = ~ 1, formula_occs_bias = ~ 1, log_precip = log_precip, n_response_curve_values = n_response_curve_values, psa_quant = psa_quant, calib = calib)

	pred_vect <- if (demesne == 'nam') {
		pred_vect <- data_occs$ag_vect_sq
	} else if (demesne == '1930s') {
		pred_vect <- vect(paste0('./data_from_adam_and_loretta/andropogon_gerardi_occurrences_with_environment_1931_1940_prism.gpkg'))
	} else if (demesne == '1950s') {
		pred_vect <- vect(paste0('./data_from_adam_and_loretta/andropogon_gerardi_occurrences_with_environment_1952_1961_prism.gpkg'))
	}

	x <- if (demesne == 'nam') {
		data_biomass$counties_x_biomass_sq
	} else if (demesne == '1930s') {
		data_biomass$counties_x_biomass_thirties
	} else if (demesne == '1950s') {
		data_biomass$counties_x_biomass_fifties
	}

	### predict mean biomass
	########################

	say('   burning present...')
	preds <- predict_biomass(chains = chains, x = x, resp_distrib = resp_distrib, transform = transform)

	pred_vect$DUMMY1 <- colMeans(preds)
	pred_vect$DUMMY2 <- apply(preds, 2, median)
	pred_vect$DUMMY3 <- apply(preds, 2, inner_quant)

	index <- (ncol(pred_vect) - 2):ncol(pred_vect)
	if (demesne == 'nam') {
		names(pred_vect)[index] <- c('mu_biomass_county_mean_sq', 'mu_biomass_county_median_sq', 'mu_biomass_county_inner_quant_sq')
	} else if (demesne == '1930s') {
		names(pred_vect)[index] <- c('mu_biomass_county_mean_1930s', 'mu_biomass_county_median_1930s', 'mu_biomass_county_inner_quant_1930s')
	} else if (demesne == '1950s') {
		names(pred_vect)[index] <- c('mu_biomass_county_mean_1950s', 'mu_biomass_county_emdian_1950s', 'mu_biomass_county_inner_quant_1950s')
	}

	### predict probability of zero biomass
	#######################################

	if (!is.null(formula_psi)) {
	
		preds <- predict_psi(chains = chains, x = x)
		preds_mean <- colMeans(preds)
		pred_vect$DUMMY1 <- preds_mean

		index <- ncol(pred_vect)
		if (demesne == 'nam') {
			names(pred_vect)[index] <- c('psi_county_sq')
		} else if (demesne == '1930s') {
			names(pred_vect)[index] <- c('psi_county_1930s')
		} else if (demesne == '1950s') {
			names(pred_vect)[index] <- c('psi_county_1950s')
		}

	}

	### future mean biomass--only for entire North American domain
	##############################################################
	if (demesne == 'nam') {

		### predict mean to each future
		###############################
		for (fut in futs) {

			say('   burning ', fut, '...')

			x <- data_biomass[paste0('counties_x_biomass_', fut)]
			x <- x[[1]]
		
			preds <- predict_biomass(chains = chains, x = x, resp_distrib = resp_distrib, transform = transform)

			pred_vect$DUMMY1 <- colMeans(preds)
			pred_vect$DUMMY2 <- apply(preds, 2, median)
			pred_vect$DUMMY3 <- apply(preds, 2, inner_quant)

			index <- (ncol(pred_vect) - 2):ncol(pred_vect)
			names(pred_vect)[index] <- paste0(c('mu_biomass_county_mean_', 'mu_biomass_county_median_', 'mu_biomass_county_inner_quant_'), fut)

		} # next future

		### predict probability of zero biomass to future
		#################################################

		if (!is.null(formula_psi)) {
		
			for (fut in futs) {

				x <- data_biomass[paste0('counties_x_biomass_', fut)]
				x <- x[[1]]

				preds <- predict_psi(chains = chains, x = x)
				preds_mean <- colMeans(preds)
				pred_vect$DUMMY1 <- preds_mean

				index <- ncol(pred_vect)
				names(pred_vect)[index] <- paste0('psi_county_', fut)

			} # next future

		} # if zero-inflated

	} # if demesne is North America

	pred_vect

}

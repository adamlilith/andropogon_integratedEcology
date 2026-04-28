#' Burn predictions from a non-biomass facet model into a SpatVector of North America or Dust Bowl map
#'
#' facet					Name of facet.
#' demesne					'nam' (North America) or '1930s' (Dust Bowl area) or '1950s' (post-Dust Bowl)
#' chains					Chains from NIMBLE
#' formula_facet			Formula for mean of site-level mean value of facet
#' formula_psi				Formula for probability of presence. Ignored if `NULL`.
#' resp_distrib 			Response distribution: 'hGamma' or 'hurdleLN'
#' transform 				Either 'exponential' or 'identity' (depending on value of resp_distrib).
#'
#' Returns a SpatVector.
burn_nonbiomass_into_vector <- function(facet, demesne, chains, formula_facet, formula_psi, resp_distrib, transform) {

	zero_inflated <- !is.null(formula_psi)

	### mean site-level biomass	
	data_facet <- prepare_nonbiomass_data(facet = facet, formula_facet = formula_facet, calib = calib)
	data_occs <- prepare_occurrence_data(formula_occs = ~ 1, formula_bias = ~ 1, n_response_curve_values = n_response_curve_values, psa_quant = psa_quant, calib = calib)

	pred_vect <- if (demesne == 'nam') {
		pred_vect <- data_occs$ag_vect_sq
	} else if (demesne == '1930s') {
		pred_vect <- vect(paste0('./data_from_adam_and_loretta/andropogon_gerardi_occurrences_with_environment_1931_1940_prism.gpkg'))
	}
	pred_vect <- simplifyGeom(pred_vect, tolerance = 1000)

	x <- if (demesne == 'nam') {
		data_facet$counties_x_sq
	} else if (demesne == '1930s') {
		data_facet$counties_x_thirties
	}

	### predict mean biomass
	########################

	say('   burning present...')
	preds <- predict_nonbiomass_single_trait(chains = chains, x = x, resp_distrib = resp_distrib, transform = transform)

	pred_vect$DUMMY1 <- colMeans(preds)
	pred_vect$DUMMY2 <- apply(preds, 2, median)
	pred_vect$DUMMY3 <- apply(preds, 2, inner_quant)

	index <- (ncol(pred_vect) - 2):ncol(pred_vect)
	if (demesne == 'nam') {
		names(pred_vect)[index] <- c(paste0(facet, '_mean_sq'), paste0(facet, '_median_sq'), paste0(facet, '_inner_quant_sq'))
	} else if (demesne == '1930s') {
		names(pred_vect)[index] <- c(paste0(facet, '_mean_1930s'), paste0(facet, '_median_1930s'), paste0(facet, '_inner_quant_1930s'))
	}

	### predict probability of zero biomass
	#######################################

	if (!is.null(formula_psi)) {
	
		preds <- predict_psi(chains = chains, x = x)
		preds_mean <- colMeans(preds)
		pred_vect$DUMMY1 <- preds_mean

		index <- ncol(pred_vect)
		if (demesne == 'nam') {
			names(pred_vect)[index] <- c('psi_sq')
		} else if (demesne == '1930s') {
			names(pred_vect)[index] <- c('psi_1930s')
		}

	}

	### future mean biomass--only for entire North American domain
	##############################################################
	if (demesne == 'nam') {

		### predict mean to each future
		###############################
		for (fut in futs) {

			say('   burning ', fut, '...')

			x <- data_facet[paste0('counties_x_', fut)]
			x <- x[[1]]
		
			preds <- predict_nonbiomass_single_trait(chains = chains, x = x, resp_distrib = resp_distrib, transform = transform)

			pred_vect$DUMMY1 <- colMeans(preds)
			pred_vect$DUMMY2 <- apply(preds, 2, median)
			pred_vect$DUMMY3 <- apply(preds, 2, inner_quant)

			index <- (ncol(pred_vect) - 2):ncol(pred_vect)
			names(pred_vect)[index] <- paste0(c(paste0(facet, '_mean_'), paste0(facet, '_median_'), paste0(facet, '_inner_quant_')), fut)

		} # next future

		### predict probability of zero biomass to future
		#################################################

		if (!is.null(formula_psi)) {
		
			for (fut in futs) {

				x <- data_facet[paste0('counties_x_', fut)]
				x <- x[[1]]

				preds <- predict_psi(chains = chains, x = x)
				preds_mean <- colMeans(preds)
				pred_vect$DUMMY1 <- preds_mean

				index <- ncol(pred_vect)
				names(pred_vect)[index] <- paste0('psi_', fut)

			} # next future

		} # if zero-inflated

	} # if demesne is North America

	pred_vect

}

#' Burn predictions from an occurrence model into a SpatVector of North America or Dust Bowl map
#'
#' source('C:/Kaji/R/andropogon_integratedEcology/sdm_pdm/functions/function_burn_occs_biomass_into_vector.r')
#'
#' demesne				'nam' (North America) or '1930s' (Dust Bowl area) or '1950s' (post-Dust Bowl)
#' chains				Chains from NIMBLE
#' formula_occs_mu		Formula for occurrences
#' formula_biomas_mu	Formula for biomass mean
#' formula_psi			Formula for inflated zero or `NULL`
#' resp_distrib 		Named vector of response distributions. For occurrence, this can be 'Poisson' or 'ZIP'. For biomass this can be 'gamma', 'ZIG' (zero-inflated gamma), 'lognormal', or 'ZILN' (zero-inflated lognormal)
#' transform			Named vector of transformations to translate MVN to mean occurrence intensity or biomass: 'identity', 'softplus' or 'exponential'.
#' log_precip			Logical--if `TRUE`, log bios 12-14 and 16-19
burn_occs_biomass_into_vector <- function(
	demesne,
	chains,
	formula_occs,
	formula_biomass,
	formula_psi,
	resp_distrib,
	transform,
	log_precip
) {

	zero_inflated <- !is.null(formula_psi)

	### occurrence data
	###################

	data_occs <- prepare_occurrence_data(formula_occs = formula_occs, formula_occs_bias = ~ 1, log_precip = log_precip, n_response_curve_values = n_response_curve_values, psa_quant = psa_quant, calib = calib)

	data_biomass <- prepare_occurrence_data(formula_occs = formula_biomass, formula_occs_bias = ~ 1, log_precip = log_precip, n_response_curve_values = n_response_curve_values, psa_quant = psa_quant, calib = calib)

	if (zero_inflated) data_occs_psi <- prepare_occurrence_data(formula_occs = formula_psi, formula_occs_bias = ~ 1, log_precip = log_precip, n_response_curve_values = n_response_curve_values, psa_quant = psa_quant, calib = calib)

	pred_vect <- if (demesne == 'nam') {
		pred_vect <- data_occs$ag_vect_sq
	} else if (demesne == '1930s') {
		pred_vect <- vect(paste0('./data_from_adam_and_loretta/andropogon_gerardi_occurrences_with_environment_1931_1940_prism.gpkg'))
	} else if (demesne == '1950s') {
		pred_vect <- vect(paste0('./data_from_adam_and_loretta/andropogon_gerardi_occurrences_with_environment_1952_1961_prism.gpkg'))
	}

	### predict present
	###################
	say('   predicting present occurrence and biomass...')

	x_occs <- if (demesne == 'nam') {
		data_occs$counties_x_occs_sq
	} else if (demesne == '1930s') {
		data_occs$counties_x_occs_thirties
	}

	x_biomass <- if (demesne == 'nam') {
		data_biomass$counties_x_occs_sq
	} else if (demesne == '1930s') {
		data_biomass$counties_x_occs_thirties
	}

	x_psi <- if (!zero_inflated) {
		NULL
	} else if (demesne == 'nam') {
		data_occs_psi$counties_x_occs_sq
	} else if (demesne == '1930s') {
		data_occs_psi$counties_x_occs_thirties
	}

	preds <- predict_occs_biomass(chains = chains, resp_distrib = resp_distrib, transform = transform, x_occs = x_occs, x_biomass = x_biomass, x_psi = x_psi)

	preds_occs <- preds$preds_occs
	preds_biomass <- preds$preds_biomass

	pred_vect$DUMMY1 <- colMeans(preds_occs)
	pred_vect$DUMMY2 <- apply(preds_occs, 2, median)
	pred_vect$DUMMY3 <- apply(preds_occs, 2, inner_quant)

	pred_vect$DUMMY4 <- colMeans(preds_biomass)
	pred_vect$DUMMY5 <- apply(preds_biomass, 2, median)
	pred_vect$DUMMY6 <- apply(preds_biomass, 2, inner_quant)

	index <- (ncol(pred_vect) - 5):ncol(pred_vect)
	if (demesne == 'nam') {
		names(pred_vect)[index] <- c('N_ag_county_mean_sq', 'N_ag_county_median_sq', 'N_ag_county_inner_quant_sq', 'mu_biomass_county_mean_sq', 'mu_biomass_county_median_sq', 'mu_biomass_county_inner_quant_sq')
	} else if (demesne == '1930s') {
		names(pred_vect)[index] <- c('N_ag_county_mean_1930s', 'N_ag_county_median_1930s', 'N_ag_county_inner_quant_1930s', 'mu_biomass_county_mean_1930s', 'mu_biomass_county_median_1930s', 'mu_biomass_county_inner_quant_1930s')
	}

	### predict probability of zero
	###############################

	if (!is.null(formula_psi)) {

		say('   predicting present probability of presence...')
		
		x <- if (demesne == 'nam') {
			data_occs_psi$counties_x_occs_sq
		} else if (demesne == '1930s') {
			data_occs_psi$counties_x_occs_thirties
		}

		preds <- predict_psi(chains = chains, x = x)
		preds_mean <- colMeans(preds)
		pred_vect$DUMMY1 <- preds_mean

		index <- ncol(pred_vect)
		if (demesne == 'nam') {
			names(pred_vect)[index] <- c('psi_county_sq')
		} else if (demesne == '1930s') {
			names(pred_vect)[index] <- c('psi_county_1930s')
		}

	}

	### future occurrences--only for entire North American domain
	##############################################################
	if (demesne == 'nam') {

		### predict mean to each future
		###############################
		for (fut in futs) {

			say('   predicting ', fut, ' occurrence and biomass...')
			
			x_occs <- data_occs[paste0('counties_x_occs_', fut)]
			x_occs <- x_occs[[1]]
		
			x_biomass <- data_biomass[paste0('counties_x_occs_', fut)]
			x_biomass <- x_biomass[[1]]
		
			 if (!zero_inflated) {
				x_psi <- NULL
			} else {
				x_psi <- data_occs_psi[paste0('counties_x_occs_', fut)]
				x_psi <- x_psi[[1]]
			}
			
			preds <- predict_occs_biomass(chains = chains, resp_distrib = resp_distrib, transform = transform, x_occs = x_occs, x_biomass = x_biomass, x_psi = x_psi)

			preds_occs <- preds$preds_occs
			preds_biomass <- preds$preds_biomass

			pred_vect$DUMMY1 <- colMeans(preds_occs)
			pred_vect$DUMMY2 <- apply(preds_occs, 2, median)
			pred_vect$DUMMY3 <- apply(preds_occs, 2, inner_quant)

			pred_vect$DUMMY4 <- colMeans(preds_biomass)
			pred_vect$DUMMY5 <- apply(preds_biomass, 2, median)
			pred_vect$DUMMY6 <- apply(preds_biomass, 2, inner_quant)

			index <- (ncol(pred_vect) - 5):ncol(pred_vect)
			names(pred_vect)[index] <- paste0(c('N_ag_county_mean_', 'N_ag_county_median_', 'N_ag_county_inner_quant_', 'mu_biomass_county_mean_', 'mu_biomass_county_median_', 'mu_biomass_county_inner_quant_'), fut)

		} 

		### predict probability of presence to future
		#############################################

		if (zero_inflated) {
		
			for (fut in futs) {

				say('   predicting ', fut, ' probability of presence...')

				x <- data_occs_psi[paste0('counties_x_occs_', fut)]
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

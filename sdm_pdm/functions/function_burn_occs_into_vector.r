#' Burn predictions from an occurrence model into a SpatVector of North America or Dust Bowl map
#'
#' demesne					'nam' (North America) or '1930s' (Dust Bowl area) or '1950s' (post-Dust Bowl)
#' chains					Chains from NIMBLE
#' formula_occs_mu			Formula for occurrences
#' formula_psi				Formula for occurrence probability of inflated zero or `NULL`
#'
#' Returns a SpatVector.
burn_occs_into_vector <- function(demesne, chains, formula_occs, formula_psi = NULL) {

	zero_inflated <- !is.null(formula_psi)

	### occurrence data
	data_occs <- prepare_occurrence_data(formula_occs = formula_occs, formula_occs_bias = ~ 1, log_precip = log_precip, n_response_curve_values = n_response_curve_values, psa_quant = psa_quant, calib = calib)

	if (zero_inflated) data_occs_psi <- prepare_occurrence_data(formula_occs = formula_psi, formula_occs_bias = ~ 1, log_precip = log_precip, n_response_curve_values = n_response_curve_values, psa_quant = psa_quant, calib = calib)

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

	x_psi <- if (!zero_inflated) {
		NULL
	} else if (demesne == 'nam') {
		data_occs_psi$counties_x_occs_sq
	} else if (demesne == '1930s') {
		data_occs_psi$counties_x_occs_thirties
	} else if (demesne == '1950s') {
		data_occs_psi$counties_x_occs_fifties
	}

	say('   burning present lambda...')
	preds <- predict_occs(chains = chains, x = x, zero_inflated = zero_inflated, x_psi = x_psi)

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

	### predict probability of zero
	###############################

	if (!is.null(formula_psi)) {
		
		x <- if (demesne == 'nam') {
			data_occs_psi$counties_x_occs_sq
		} else if (demesne == '1930s') {
			data_occs_psi$counties_x_occs_thirties
		} else if (demesne == '1950s') {
			data_occs_psi$counties_x_occs_fifties
		}

		say('   burning present psi...')
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

	### future occurrences--only for entire North American domain
	##############################################################
	if (demesne == 'nam') {

		### predict mean to each future
		###############################
		for (fut in futs) {
			
			x <- data_occs[paste0('counties_x_occs_', fut)]
			x <- x[[1]]
		
			 if (!zero_inflated) {
				x_psi <- NULL
			} else {
				x_psi <- data_occs_psi[paste0('counties_x_occs_', fut)]
				x_psi <- x_psi[[1]]
			}
		
			say('   burning ', fut, ' lambda...')
			preds <- predict_occs(chains = chains, x = x, x_psi = x_psi, zero_inflated = zero_inflated)

			preds_mean <- colMeans(preds)
			preds_sd <- apply(preds, 2, sd)

			pred_vect$DUMMY1 <- preds_mean
			pred_vect$DUMMY2 <- preds_sd

			index <- (ncol(pred_vect) - 1):ncol(pred_vect)
			names(pred_vect)[index] <- paste0(c('N_ag_county_mean_', 'N_ag_county_sd_'), fut)

		} # next future

		### predict probability of zero abundance to future
		###################################################

		if (!is.null(formula_psi)) {
		
			for (fut in futs) {

				x <- data_occs_psi[paste0('counties_x_occs_', fut)]
				x <- x[[1]]

				say('   burning ', fut, ' psi...')
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

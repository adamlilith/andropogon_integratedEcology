#' Burn predictions from a biomass model into a SpatVector of North America or Dust Bowl map
#'
#' chains					Chains from NIMBLE
#' formula_biomass_mu		Formula for mean of site-level biomass
#' formula_biomass_sigma	Formula for sd of site-level biomass. Ignored if `NULL`.
#' formula_pzero	Formula for probability of zero biomass. Ignored if `NULL`.
burn_biomass_into_vector <- function(chains, formula_biomass_mu, formula_pzero) {

	zero_inflated <- !is.null(formula_pzero)

	### mean site-level biomass	
	data_biomass_mu <- prepare_biomass(formula_biomass = formula_biomass_mu, calib = calib)
	data_occs <- prepare_occurrences(formula_occs = ~ 1, formula_occs_bias = ~ 1, n_response_curve_values = n_response_curve_values, psa_quant = psa_quant, calib = calib)

	pred_vect <- data_occs$ag_vect_sq

	# predict to present
	x <- data_biomass_mu$counties_x_biomass_sq

	### predict mean biomass
	########################

	preds <- predict_biomass(chains = chains, x = x, zero_inflated = zero_inflated, type = 'mu')

	preds_mean <- colMeans(preds)
	preds_sd <- apply(preds, 2, sd)

	pred_vect$DUMMY1 <- preds_mean
	pred_vect$DUMMY2 <- preds_sd

	index <- (ncol(pred_vect) - 1):ncol(pred_vect)
	names(pred_vect)[index] <- c('mu_biomass_county_mean_sq', 'mu_biomass_county_sd_sq')

	### predict probability of zero biomass
	#######################################

	if (!is.null(formula_pzero)) {
	
		preds <- predict_biomass(chains = chains, x = x, zero_inflated = zero_inflated, type = 'pzero')
		preds_mean <- colMeans(preds)
		pred_vect$DUMMY1 <- preds_mean

		index <- ncol(pred_vect)
		names(pred_vect)[index] <- c('pzero_biomass_county_sq')

	}

	pred_vect

}

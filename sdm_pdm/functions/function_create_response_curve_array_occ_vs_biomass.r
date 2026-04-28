#' This function creates a response curve array for the relationship between the number of occurrences and biomass.
#'
#' @param formula The formula for the response of abundance to biomass. May or may not include an intercept. If not explicilitly included, the intercept will be added.
#' @param data_biomass Output for `prepare_biomass_data()` function.
#' @param n_response_curve_values A numeric value specifying the number of values to use for the response curve. Default is 200.
#' @return A list with an unscaled version of biomass and a model matrix with scaled biomass.
create_response_curve_array_occ_vs_biomass <- function(formula, data_biomass, n_response_curve_values = 200) {

	# create matrix for predicting response of SDM lambda to biomass
	site_biomass_mean <- data_biomass$x_centers # mean biomass for scaling
	site_biomass_sd <- data_biomass$x_scales # sd of biomass for scaling

	biomass_max <- omnibus::roundTo(1.1 * max(data_biomass$y_biomass), 10, ceiling) # max value of biomass for response curve

	biomass_seq <- seq(0, biomass_max, length.out = n_response_curve_values)
	response_curves_x <- scale(biomass_seq, center = site_biomass_mean, scale = site_biomass_sd)
	response_curves_x <- data.frame(biomass = response_curves_x)
	response_curves_x <- model.matrix(formula, response_curves_x)

	list(
		unscaled_biomass = biomass_seq,	# unscaled sequence of biomasses
		response_curves_x_occs_vs_biomass = response_curves_x	# scaled model matrix for response curve
	)

}
	

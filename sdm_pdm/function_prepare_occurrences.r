#' This function loads and prepares the occurrence data for distribution modeling. See the "return()` line for details on the output.
#' 
#' @param formula_occs A formula object specifying the model to be fit for the response of AG abundance to the environment. Include an intercept and the RHS only.
#' @param formula_occs_bias A formula object specifying the model to be fit for the sampling bias in the number of AG observed occurrences. Include an intercept and the RHS only.
#' @param n_response_curve_values A numeric value specifying the number of values to use for the response curve. Default is 200.
#' @param psa_quant A numeric value between 0 and 1 specifying the quantile of Poaceae records to use for pseudo-absences. Default is 0.99.
#' @param calib If `TRUE`, then the training data is subset only to counties with non-`NA` for Poaceae.
#'
prepare_occurrences <- function(formula_occs, formula_occs_bias, n_response_curve_values = 200, psa_quant = 0.99, calib = TRUE) {

	# formulae
	terms_occs <- terms(formula_occs)
	terms_occs <- attr(terms_occs, 'term.labels')

	covariates <- terms_occs
	covariates <- covariates[!grepl(covariates, pattern = '\\^2')]
	covariates <- covariates[!grepl(covariates, pattern = '\\:')]
	covariates <- covariates[!grepl(covariates, pattern = '\\*')]

	n_covariates <- length(covariates)

	terms_occs_bias <- terms(formula_occs_bias)
	terms_occs_bias <- attr(terms_occs_bias, 'term.labels')

	covariates_bias <- terms_occs_bias
	covariates_bias <- covariates_bias[!grepl(covariates_bias, pattern = '\\^2')]
	covariates_bias <- covariates_bias[!grepl(covariates_bias, pattern = '\\:')]
	covariates_bias <- covariates_bias[!grepl(covariates_bias, pattern = '\\*')]

	n_covariates_bias <- length(covariates_bias)

	# load county-level environmental data
	ag_vect_sq <- vect('./outputs_loretta/integrated_sdm_pdm/andropogon_gerardi_occurrences_with_environment_1961_2020_for_integration.gpkg')
	ag_vect_sq <- ag_vect_sq[ , c('focal_region', 'geofold', 'country', 'state_province', 'county', 'area_km2', 'elevation_m', 'n_andropogon_gerardi', 'n_poaceae', covariates)]

	# remove occurrences in largest parts of Ontario and Manitoba on basis that they are just too big to indicate species-environment relationships
	redacts <- c(
		which(ag_vect_sq$state_province == 'Ontario' & ag_vect_sq$county == 'Kenora'),
		which(ag_vect_sq$state_province == 'Ontario' & ag_vect_sq$county == 'Cochrane'),
		which(ag_vect_sq$state_province == 'Manitoba' & ag_vect_sq$county == 'Division No. 19')
	)

	for (redact in redacts) {
		ag_vect_sq$n_andropogon_gerardi[redact] <- 0
	}

	# define focal region for defining range of environment for response curves, and region for calibration of model
	if (calib) {
		ag_vect_sq$calib_region <- !is.na(ag_vect_sq$n_poaceae)
	} else {
		ag_vect_sq$calib_region <- TRUE
	}

	# for counties with NA Poaceae and AG, assign a maximal number of Poaceae and 0 AG
	n_pseudoabs <- quantile(ag_vect_sq$n_poaceae, psa_quant, na.rm = TRUE)
	ag_vect_sq$n_poaceae[is.na(ag_vect_sq$n_poaceae)] <- n_pseudoabs
	ag_vect_sq$n_andropogon_gerardi[is.na(ag_vect_sq$n_andropogon_gerardi)] <- 0

	ag_sq <- as.data.frame(ag_vect_sq)
	ag_sq$area_km2_log10 <- log10(ag_sq$area_km2)
	ag_sq$n_poaceae_log10p1 <- log10(ag_sq$n_poaceae + 1)

	bias_frame <- data.frame(area_km2_log10 = ag_sq$area_km2_log10, n_poaceae_log10p1 = ag_sq$n_poaceae_log10p1)
	bias_frame <- scale(bias_frame)
	bias_centers <- attributes(bias_frame)$`scaled:center`
	bias_scales <- attributes(bias_frame)$`scaled:scale`
	bias_frame <- as.data.frame(bias_frame)
	w_occs_bias <- model.matrix(formula_occs_bias, bias_frame)

	terms_occs_bias <- terms(formula_occs_bias)
	terms_occs_bias <- attr(terms_occs_bias, 'term.labels')
	n_terms_occs_bias <- ncol(w_occs_bias)

	# scale covariates... only use calibration region for calculating centers and scales
	covars <- ag_sq[ag_vect_sq$calib_region, c(covariates)]
	covars_scaled <- scale(covars)

	x_centers <- attributes(covars_scaled)$`scaled:center`
	x_scales <- attributes(covars_scaled)$`scaled:scale`

	ag_sq <- ag_sq[ , covariates, drop = FALSE]
	ag_sq <- scale(ag_sq, center = x_centers[covariates], scale = x_scales[covariates])
	ag_sq <- as.data.frame(ag_sq)
	counties_x_sq <- model.matrix(formula_occs, ag_sq)

	n_terms <- ncol(counties_x_sq)

	for (fut in futs) {
	
		this_fut <- vect(paste0('./data_from_adam_and_loretta/andropogon_gerardi_occurrences_with_environment_ensemble_8GCMs_', fut, '_climatena.gpkg'))

		this_fut <- this_fut[ , c('country', 'state_province', 'county', covariates)]

		assign(paste0('counties_', fut), this_fut)

		for (covariate in covariates) {
			
			x <- this_fut[[covariate]]
			x <- unlist(x)
			x <- scale(x, center = x_centers[covariate], scale = x_scales[covariate])
			x <- as.numeric(x)
			this_fut[ , covariate] <- x

		}

		mm <- as.data.frame(this_fut)[ , covariates, drop = FALSE] 
		mm <- model.matrix(formula_occs, mm)

		assign(paste0('counties_x_', fut), mm)
	
	}

	# 20th century climate
	thirties <- vect(paste0('./data_from_adam_and_loretta/andropogon_gerardi_occurrences_with_environment_1931_1940_prism.gpkg'))
	fifties <- vect(paste0('./data_from_adam_and_loretta/andropogon_gerardi_occurrences_with_environment_1952_1961_prism.gpkg'))
	
	for (clim in c('thirties', 'fifties')) {

		x <- get(clim)
		x <- x[ , c('country', 'state_province', 'county', covariates)]

		for (covariate in covariates) {
			
			xx <- x[[covariate]]
			xx <- unlist(xx)
			xx <- scale(xx, center = x_centers[covariate], scale = x_scales[covariate])
			xx <- as.numeric(xx)
			x[ , covariate] <- xx

		}

		mm <- as.data.frame(x)[ , covariates, drop = FALSE] 
		mm <- model.matrix(formula_occs, mm)

		assign(paste0('counties_x_occs_', clim), mm)
	
	}
	n_counties_20th_cent <- nrow(thirties)

	### response curve array: occurrence vs environment
	resp_arrays <- create_response_curve_array(
		formula = formula_occs,
		centers = x_centers,
		scales = x_scales,
		ag_vect_sq = ag_vect_sq[ag_vect_sq$focal_region],
		vects = list(
			counties_ssp245_2041_2070[ag_vect_sq$focal_region],
			counties_ssp245_2071_2100[ag_vect_sq$focal_region],
			counties_ssp370_2041_2070[ag_vect_sq$focal_region],
			counties_ssp370_2071_2100[ag_vect_sq$focal_region],
			thirties, fifties
		),
		n_response_curve_values = n_response_curve_values
	)

	# response curve array: thinning vs bias covariates
	ag_vect_sq_bias <- ag_vect_sq
	ag_vect_sq_bias$area_km2_log10 <- log10(ag_vect_sq_bias$area_km2)
	ag_vect_sq_bias$n_poaceae_log10p1 <- log10(1 + ag_vect_sq_bias$n_poaceae)
	resp_arrays_bias <- create_response_curve_array(
		formula = formula_occs_bias,
		centers = bias_centers,
		scales = bias_scales,
		ag_vect_sq = ag_vect_sq_bias[ag_vect_sq_bias$focal_region],
		vects = NULL,
		n_response_curve_values = n_response_curve_values
	)

	y_n_ag <- ag_vect_sq$n_andropogon_gerardi[ag_vect_sq$calib_region]

	list(
		terms_occs = terms_occs,				# terms_occs in formula_occs
		covariates_occs = covariates,			# covariates in formula_occs (can appear in >1 term)
		n_terms_occs = n_terms,					# number of terms_occs in formula_occs
		n_covariates_occs = n_covariates,		# number of variables in formula_occs (which can appear in >1 term)
	
		terms_occs_bias = terms_occs_bias,
		covariates_occs_bias = covariates_bias,
		n_terms_occs_bias = n_terms_occs_bias,
		n_covariates_occs_bias = n_covariates_bias,		# number of variables in formula_occs (which can appear in >1 term)
		w_occs_bias = w_occs_bias,

		y_n_ag = y_n_ag,						# number of AG records in each county
		n_counties_occs_calib = sum(ag_vect_sq$calib_region),	# number of counties in calibration region
		n_counties = nrow(ag_vect_sq),			# number of counties in dataset

		x_centers_occs = x_centers,				# vector of means for each covariate
		x_scales_occs = x_scales,				# vector of standard deviations for each covariate

		w_centers_occs_bias = bias_centers,				# vector of means for each covariate
		w_scales_occs_bias = bias_scales,				# vector of standard deviations for each covariate

		ag_vect_sq = ag_vect_sq,				# SpatVector of counties for plotting

		counties_x_occs_sq = counties_x_sq,		# model matrix of county-level environmental data for current conditions

		counties_occs_ssp245_2041_2070 = counties_ssp245_2041_2070, # SpatVector of county-level environmental data for future conditions (SSP245, 2041-2070)
		counties_occs_ssp245_2071_2100 = counties_ssp245_2071_2100, # SpatVector of county-level environmental data for future conditions (SSP245, 2071-2100)
		counties_occs_ssp370_2041_2070 = counties_ssp370_2041_2070, # SpatVector of county-level environmental data for future conditions (SSP370, 2041-2070)
		counties_occs_ssp370_2071_2100 = counties_ssp370_2071_2100, # SpatVector of county-level environmental data for future conditions (SSP370, 2071-2100)

		counties_x_occs_ssp245_2041_2070 = counties_x_ssp245_2041_2070, # model matrix of county-level environmental data for future conditions (SSP245, 2041-2070)
		counties_x_occs_ssp245_2071_2100 = counties_x_ssp245_2071_2100, # model matrix of county-level environmental data for future conditions (SSP245, 2071-2100)
		counties_x_occs_ssp370_2041_2070 = counties_x_ssp370_2041_2070, # model matrix of county-level environmental data for future conditions (SSP370, 2041-2070)
		counties_x_occs_ssp370_2071_2100 = counties_x_ssp370_2071_2100, # model matrix of county-level environmental data for future conditions (SSP370, 2071-2100),

		n_counties_20th_cent = n_counties_20th_cent, # SpatVector of county-level environmental data for 20th century
		counties_thirties = thirties, # SpatVector of county-level environmental data for 20th century
		counties_fifties = fifties, # SpatVector of county-level environmental data for 20th century

		counties_x_occs_thirties = counties_x_occs_thirties, # model matrix of county-level environmental data for 20th century
		counties_x_occs_fifties = counties_x_occs_fifties, # model matrix of county-level environmental data for 20th century

		resp_curves_x_occs = resp_arrays$response_curves_x_scaled, # response curve array
		resp_curves_x_occs_unscaled = resp_arrays$resp_curves_x_unscaled, # response curve array (unscaled)

		resp_curves_w_occs = resp_arrays_bias$response_curves_x_scaled, # response curve array
		resp_curves_w_occs_bias_unscaled = resp_arrays_bias$resp_curves_x_unscaled # response curve array (unscaled)

	)

}

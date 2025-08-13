#' This function loads and prepares the occurrence data for distribution modeling. See the "return()` line for details on the output.
#' 
#' @param formula A formula object specifying the model to be fit. Include an intercept and the RHS only.
#' @param n_response_curve_values A numeric value specifying the number of values to use for the response curve. Default is 200.
#' @param psa_quant A numeric value between 0 and 1 specifying the quantile of Poaceae records to use for pseudo-absences. Default is 0.99.
#' @param calib If `TRUE`, then the training data is subset only to counties with non-`NA` for Poaceae.
#' @param site_vect A points SpatVector with locations of sites where phenotypes were measured.
#'
prepare_occurrences_vis_a_vis_phenotyped_sites <- function(formula, n_response_curve_values = 200, psa_quant = 0.99, calib = TRUE, site_vect = NULL) {

	# formula
	terms <- terms(formula)
	terms <- attr(terms, 'term.labels')

	covariates <- terms
	covariates <- covariates[!grepl(covariates, pattern = '\\^2')]
	covariates <- covariates[!grepl(covariates, pattern = '\\:')]
	covariates <- covariates[!grepl(covariates, pattern = '\\*')]

	n_covariates <- length(covariates)

	# load county-level environmental data
	ag_vect_sq <- vect('./outputs_loretta/integrated_sdm_pdm/andropogon_gerardi_occurrences_with_environment_1961_2020_for_integration.gpkg')
	ag_vect_sq <- ag_vect_sq[ , c('focal_region', 'geofold', 'country', 'state_province', 'county', 'area_km2', 'elevation_m', 'n_andropogon_gerardi', 'n_poaceae', covariates)]

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

	# scale covariates... only use calibration region for calculating centers and scales
	covars <- ag_sq[ag_vect_sq$calib_region, c('area_km2_log10', 'n_poaceae_log10p1', covariates)]
	covars_scaled <- scale(covars)

	x_centers <- attributes(covars_scaled)$`scaled:center`
	x_scales <- attributes(covars_scaled)$`scaled:scale`

	area_km2_log10 <- covars_scaled[ , 'area_km2_log10']
	n_poaceae_log10p1 <- covars_scaled[ , 'n_poaceae_log10p1']

	ag_sq <- ag_sq[ , covariates, drop = FALSE]
	ag_sq <- scale(ag_sq, center = x_centers[covariates], scale = x_scales[covariates])
	ag_sq <- as.data.frame(ag_sq)
	counties_x_sq <- model.matrix(formula, ag_sq)

	n_terms <- ncol(counties_x_sq)

	for (fut in futs) {
	
		this_fut <- vect(paste0('./data_from_adam_and_loretta/andropogon_gerardi_occurrences_with_environment_ensemble_8GCMs_', fut, '.gpkg'))

		this_fut <- this_fut[ , c('country', 'stateProvince', 'county', 'area_km2', covariates)]
		names(this_fut)[names(this_fut) == 'stateProvince'] <- 'state_province'

		assign(paste0('counties_', fut), this_fut)

		for (covariate in covariates) {
			
			x <- this_fut[[covariate]]
			x <- unlist(x)
			x <- scale(x, center = x_centers[covariate], scale = x_scales[covariate])
			x <- as.numeric(x)
			this_fut[ , covariate] <- x

		}

		mm <- as.data.frame(this_fut)[ , covariates, drop = FALSE] 
		mm <- model.matrix(formula, mm)

		assign(paste0('counties_x_', fut), mm)
	
	}

	### response curve array
	resp_arrays <- create_response_curve_array(
		formula = formula,
		centers = x_centers,
		scales = x_scales,
		ag_vect_sq = ag_vect_sq[ag_vect_sq$focal_region],
		vects = list(
			counties_ssp245_2041_2070[ag_vect_sq$focal_region],
			counties_ssp245_2071_2100[ag_vect_sq$focal_region],
			counties_ssp370_2041_2070[ag_vect_sq$focal_region],
			counties_ssp370_2071_2100[ag_vect_sq$focal_region]
		),
		n_response_curve_values = n_response_curve_values
	)

	y_n_ag <- ag_vect_sq$n_andropogon_gerardi[ag_vect_sq$calib_region]

	### separate values at phenotyped sites from other values
	ag_vect_sq_TEMP <- ag_vect_sq
	ag_vect_sq_TEMP$id <- 1:nrow(ag_vect_sq_TEMP)
	pheno_counties <- extract(ag_vect_sq_TEMP, site_vect)
	ag_vect_sq$pheno_site <- NA_character_
	ag_vect_sq$pheno_site[pheno_counties$id] <- site_vect$site_id

	list(
		terms_occs = terms,						# terms in formula
		covariates_occs = covariates,			# covariates in formula (can appear in >1 term)
		n_terms_occs = n_terms,					# number of terms in formula
		n_covariates_occs = n_covariates,		# number of variables in formula (which can appear in >1 term)
	
		y_n_ag = y_n_ag,						# number of AG records in each county
		y_n_ag_pheno_sites = y_n_ag_pheno_sites,	# observed number of AG records in each county with a phenotype site
		area_km2_log10 = area_km2_log10,	# area of each county (log10)
		n_poaceae_log10p1 = n_poaceae_log10p1,	# number of Poaceae records in each county (log10 + 1)
		n_counties_occs_calib = sum(ag_vect_sq$calib_region),	# number of counties in calibration region
		n_counties = nrow(ag_vect_sq),			# number of counties in dataset

		x_centers_occs = x_centers,				# vector of means for each covariate
		x_scales_occs = x_scales,				# vector of standard deviations for each covariate

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

		resp_curves_x_occs = resp_arrays$response_curves_x_scaled, # response curve array
		resp_curves_x_occs_unscaled = resp_arrays$resp_curves_x_unscaled # response curve array (unscaled)

	)

}

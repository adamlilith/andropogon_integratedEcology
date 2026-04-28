#' This function loads and prepares the biomass data for distribution modeling. See the "return()` line for details on the output.
#' 
#' source('C:/Kaji/R/andropogon_integratedEcology/sdm_pdm/functions/function_prepare_biomass_data.r')
#'
#' formula_biomass A formula_biomass object specifying the model to be fit. Include an intercept and the RHS only.
#' n_response_curve_values A numeric value specifying the number of values to use for the response curve. Default is 200.
#' calib If `TRUE`, then the training data is subset only to counties with non-`NA` for Poaceae.
#'
prepare_biomass_data <- function(formula_biomass, n_response_curve_values = 200, calib = TRUE) {

	# formula
	terms <- terms(formula_biomass)
	terms <- attr(terms, 'term.labels')

	covariates <- terms
	covariates <- covariates[!grepl(covariates, pattern = '\\^2')]
	covariates <- covariates[!grepl(covariates, pattern = '\\:')]
	covariates <- covariates[!grepl(covariates, pattern = '\\*')]

	n_covariates <- length(covariates)

	# load site data
	site_data_raw <- readRDS('./data_from_loretta/sdm_pdm_00_merged_site_data_with_climate/sites.rds')
	biomass_data_raw <- readRDS('./data_from_loretta/sdm_pdm_00_merged_site_data_with_climate/biomass.rds')

	site_data_raw <- calculate_logged_vars(site_data_raw)
	biomass_data_raw <- calculate_logged_vars(biomass_data_raw)

	# ensure sites in each array appear in the same order
	site_data_raw <- site_data_raw[order(site_id)]
	biomass_data_raw <- biomass_data_raw[order(SITE)]

	stopifnot(all(site_data_raw$site_id == unique(biomass_data_raw$SITE)))

	# add soil chemistry variables from field samples
	soil_data_raw <- openxlsx::read.xlsx('./data_from_loretta/!plant_sitelevel_data_11NOV2024 - Google Sheets [aggregated by Erica].xlsx', sheet = 'original_soil_field_data')

	soil_data_raw_site_plant <- tolower(soil_data_raw$X1)

	soil_data_raw_site_plant_spaces <- unlist(gregexpr(soil_data_raw_site_plant, pattern = ' '))
	soil_data_raw_site <- substr(soil_data_raw_site_plant, 1, soil_data_raw_site_plant_spaces - 1)
	last_char <- substr(soil_data_raw_site, nchar(soil_data_raw_site), nchar(soil_data_raw_site))
	append_num <- rep('', length(soil_data_raw_site))
	append_num[last_char %notin% c('1', '2', '3')] <- '1'
	soil_data_raw_site <- paste0(soil_data_raw_site, append_num)

	site_data_raw_site <- gsub(tolower(site_data_raw$site_id), pattern = '_', replacement = '')
	
	matches <- match(site_data_raw_site, soil_data_raw_site)
	site_data_raw$site_ph <- soil_data_raw$pH[matches]
	site_data_raw$site_ph_sikora <- soil_data_raw$`Sikora_pH`[matches]
	# site_data_raw$site_nitrogen <- soil_data_raw$`Total.N.%`[matches]
	site_data_raw$nitrogen <- soil_data_raw$`Total.N.%`[matches]

	biomass_data_raw_site <- gsub(tolower(biomass_data_raw$SITE), pattern = '_', replacement = '')
	matches <- match(biomass_data_raw_site, soil_data_raw_site)
	biomass_data_raw$site_ph <- soil_data_raw$pH[matches]
	biomass_data_raw$site_ph_sikora <- soil_data_raw$`Sikora_pH`[matches]
	# biomass_data_raw$site_nitrogen <- soil_data_raw$`Total.N.%`[matches]
	biomass_data_raw$nitrogen <- soil_data_raw$`Total.N.%`[matches]

	# add soil texture variables from field samples

	site_data_raw$sand <- site_data_raw$SAND
	site_data_raw$silt <- site_data_raw$SILT
	site_data_raw$clay <- site_data_raw$CLAY

	biomass_data_raw[ , c('sand', 'silt', 'clay') := NA_real_]
	for (i in 1:nrow(site_data_raw)) {

		biomass_data_raw$sand[biomass_data_raw$SITE == site_data_raw$site_id[i]] <- site_data_raw$sand[i]
		biomass_data_raw$silt[biomass_data_raw$SITE == site_data_raw$site_id[i]] <- site_data_raw$silt[i]
		biomass_data_raw$clay[biomass_data_raw$SITE == site_data_raw$site_id[i]] <- site_data_raw$clay[i]

	}

	n_pheno_sites <- length(unique(site_data_raw$site_id))
	
	# create spatial versions of site and phenotype data for plotting
	site_data <- merge(site_data_raw[ , c('site_id', 'LONGITUDE', 'LATITUDE')], biomass_data_raw, by.x = 'site_id', by.y = 'SITE')
	
	# means and sd
	site_data_agg <- site_data[ , .(
		biomass_mean = mean(Biomass, na.rm = TRUE),
		biomass_median = median(Biomass, na.rm = TRUE),
		biomass_sd = sd(Biomass, na.rm = TRUE)
	), by = site_id]
	site_data_agg <- merge(site_data_agg, site_data_raw[ , c('site_id', 'LONGITUDE', 'LATITUDE')], by.x = 'site_id', by.y = 'site_id')
	site_vect <- vect(site_data_agg, geom = c('LONGITUDE', 'LATITUDE'), crs = getCRS('WGS84'))
	
	# center and scale at site level
	x <- site_data_raw
	x <- x[ , ..covariates]

	centers_scales <- fread('./outputs_loretta/integrated_sdm_pdm/centers_and_scales_for_covariates.csv')
	x_centers <- centers_scales$center_sites[match(covariates, centers_scales$variable)]
	x_scales <- centers_scales$scale_sites[match(covariates, centers_scales$variable)]
	names(x_centers) <- covariates
	names(x_scales) <- covariates

	x <- scale(x, center = x_centers, scale = x_scales)
	x <- as.data.frame(x)
	x_by_site <- model.matrix(formula_biomass, x)

	n_terms <- ncol(x_by_site)

	site_biomass_mean <- biomass_data_raw[ , .(val = mean(Biomass)), by = SITE][['val']]
	site_biomass_sd <- biomass_data_raw[ , .(val = sd(Biomass)), by = SITE][['val']]

	across_sites_biomass_mean <- mean(site_biomass_mean)
	across_sites_biomass_sd <- sd(site_biomass_mean)

	# make vector of which sampled site matches each row in the biomass and morphology/physiology data
	n <- nrow(biomass_data_raw)

	site_index <- rep(NA, n)
	for (i in 1:n) site_index[i] <- which(site_data_raw$site_id == biomass_data_raw$SITE[i])

	### county-level environmental data
	ag_vect_sq <- vect('./outputs_loretta/integrated_sdm_pdm/andropogon_gerardi_occurrences_with_environment_1961_2020_for_integration.gpkg')
	ag_vect_sq <- simplifyGeom(ag_vect_sq, tolerance = 1000)
	# names(ag_vect_sq)[names(ag_vect_sq) == 'nitrogen'] <- 'site_nitrogen'
	ag_vect_sq <- calculate_logged_vars(ag_vect_sq)
	ag_sq <- as.data.frame(ag_vect_sq)
	ag_sq <- ag_sq[ , covariates, drop = FALSE]

	ag_sq <- scale(ag_sq, center = x_centers[covariates], scale = x_scales[covariates])
	ag_sq <- as.data.frame(ag_sq)
	counties_x_sq <- model.matrix(formula_biomass, ag_sq)

	n_counties <- nrow(ag_sq)

	site_vect <- project(site_vect, ag_vect_sq)

	# define focal region for defining range of environment for response curves, and region for calibration of model
	if (calib) {
		ag_vect_sq$calib_region <- !is.na(ag_vect_sq$n_poaceae)
	} else {
		ag_vect_sq$calib_region <- TRUE
	}
	index_biomass_county_in_calib_county  <- which(ag_vect_sq$calib_region)

	### future environmental data
	for (fut in futs) {
	
		this_fut <- vect(paste0('./data_from_adam_and_loretta/andropogon_gerardi_occurrences_with_environment_ensemble_8GCMs_', fut, '.gpkg'))
		this_fut <- simplifyGeom(this_fut, tolerance = 1000)
		this_fut <- calculate_logged_vars(this_fut)
		# names(this_fut)[names(this_fut) == 'nitrogen'] <- 'site_nitrogen'

		if (n_covariates > 0) this_fut <- this_fut[ , covariates]
		assign(paste0('counties_', fut), this_fut)

		this_fut <- as.data.frame(this_fut)[ , covariates, drop = FALSE]

		this_fut <- scale(this_fut, center = x_centers[covariates], scale = x_scales[covariates])
		this_fut <- as.data.frame(this_fut)
		this_fut <- model.matrix(formula_biomass, this_fut)

		assign(paste0('counties_x_', fut), this_fut)
	
	}

	# 20th century climate
	thirties <- vect(paste0('./data_from_adam_and_loretta/andropogon_gerardi_occurrences_with_environment_1931_1940_prism.gpkg'))
	thirties <- calculate_logged_vars(thirties)

	# names(thirties)[names(thirties) == 'nitrogen'] <- 'site_nitrogen'

	# for (clim in c('thirties', 'fifties')) {
	for (clim in c('thirties')) {

		x <- get(clim)
		x <- x[ , c('country', 'state_province', 'county', covariates)]

		if (length(covariates) > 0) {
			
			for (covariate in covariates) {
				
				xx <- x[[covariate]]
				xx <- unlist(xx)

				xx <- scale(xx, center = x_centers[covariate], scale = x_scales[covariate])
				xx <- as.numeric(xx)
				x[ , covariate] <- xx

			}

			mm <- as.data.frame(x)[ , covariates, drop = FALSE] 
			mm <- model.matrix(formula_biomass, mm)

		} else {
		
			mm <- NA

		}

		assign(paste0('counties_x_', clim), mm)
		
		
	}
	n_counties_1930s <- nrow(thirties)
	
	### response array
	resp_arrays <- create_response_curve_array(
		formula = formula_biomass,
		centers = x_centers,
		scales = x_scales,
		site_data_raw = site_data_raw,
		ag_vect_sq = ag_vect_sq[ag_vect_sq$n_andropogon_gerardi > 0],
		vects = list(
			counties_ssp245_2041_2070[ag_vect_sq$n_andropogon_gerardi > 0],
			counties_ssp245_2071_2100[ag_vect_sq$n_andropogon_gerardi > 0],
			counties_ssp370_2041_2070[ag_vect_sq$n_andropogon_gerardi > 0],
			counties_ssp370_2071_2100[ag_vect_sq$n_andropogon_gerardi > 0],
			thirties#, fifties
		),
		n_response_curve_values = n_response_curve_values
	)

	list(

		facet = 'biomass',

		site_data_raw = site_data_raw,			# raw site data
		raw_data_biomass = biomass_data_raw,	# raw data for biomass

		terms = terms,					# terms in formula_biomass
		covariates = covariates,		# covariates in formula_biomass (which can appear in >1 term)
		n_terms = n_terms,				# number of terms in formula_biomass
		n_covariates = n_covariates,	# number of variables in formula_biomass (which can appear in >1 term)
	
		y_biomass = biomass_data_raw$Biomass,		# response variable (biomass)
		n_biomass = length(biomass_data_raw$Biomass),	# number of observations in response variable
		site_index_biomass = site_index,		# index of which values in the response variable correspond to which sampled site
		n_pheno_sites = n_pheno_sites,			# number of sites at which phenotypes were measured

		site_mean = site_biomass_mean,	# mean site-level biomass
		site_sd = site_biomass_sd,		# sd of site-level biomass

		across_sites_mean = across_sites_biomass_mean,	# mean biomass across sites
		across_sites_sd = across_sites_biomass_mean,	# s.d. across sites

		x_by_site = x_by_site,			# MM with covariates (scaled)
		x_centers = x_centers,			# vector of means for each covariate
		x_scales = x_scales,			# vector of standard deviations for each covariate

		site_vect = site_vect,			# SpatVector of site data for plotting

		index_biomass_county_in_calib_county = index_biomass_county_in_calib_county, # index of counties in SDM calibration region

		n_counties = n_counties,				# number of counties in dataset
		# counties_x_sq_calib = counties_x_sq_calib, # MM of county-level environmental data for present (SDM calibration region)
		counties_x_sq = counties_x_sq,		# MM of county-level environmental data for current conditions
		counties_x_ssp245_2041_2070 = counties_x_ssp245_2041_2070, # MM of county-level environmental data for future conditions (SSP245, 2041-2070)
		counties_x_ssp245_2071_2100 = counties_x_ssp245_2071_2100, # MM of county-level environmental data for future conditions (SSP245, 2071-2100)
		counties_x_ssp370_2041_2070 = counties_x_ssp370_2041_2070, # MM of county-level environmental data for future conditions (SSP370, 2041-2070)
		counties_x_ssp370_2071_2100 = counties_x_ssp370_2071_2100, # MM of county-level environmental data for future conditions (SSP370, 2071-2100)

		n_counties_1930s = n_counties_1930s, # SpatVector of county-level environmental data for 20th century
		counties_thirties = thirties, # SpatVector of county-level environmental data for 20th century
		# counties_fifties = fifties, # SpatVector of county-level environmental data for 20th century

		counties_x_thirties = counties_x_thirties, # model matrix of county-level environmental data for 20th century
		# counties_x_fifties = counties_x_fifties, # model matrix of county-level environmental data for 20th century

		resp_curves_x = resp_arrays$response_curves_x_scaled,			# model matrices for response curves
		resp_curves_x_unscaled = resp_arrays$resp_curves_x_unscaled	# matrices for response curves

	)

}

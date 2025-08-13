#' This function loads and prepares non-biomass trait (morphological and physiological) data for distribution modeling. See the "return()` line for details on the output.
#' 
#' @param trait One of:
#' * delta13c (will fail bc some values are NA!!!)
#' * n_concentration
#' * cn_ratio
#' * height
#' * blade_width
#' * leaf_thickness
#' * spad
#' * canopy_diameter
#' * water_potential
#' * photosynthetic_rate
#' * stomatal_conductance
#' * internal_co2
#' * transpiration_rate
#' @param formula_traits A formula_traits object specifying the model to be fit. Include an intercept and the RHS only.
#' @param n_response_curve_values A numeric value specifying the number of values to use for the response curve. Default is 200.
#' @param calib If `TRUE`, then the training data is subset only to counties with non-`NA` for Poaceae.
#'
prepare_nonbiomass_traits <- function(trait, formula_traits, n_response_curve_values = 200, calib = TRUE) {

	if (trait == 'delta13c') stop('Cannot use `delta13c` because some values are NA.')

	# get traits
	trait_columns <- c('Delta13C', 'N_conc', 'CN_ratio', 'Height', 'BladeWidth', 'LeafThick', 'SPAD', 'CanopyDiam', 'WatPot', 'PhotoRate', 'StomCond', 'IntCO2', 'TranspRate')
	nice_trait <- c('delta13c', 'n_concentration', 'cn_ratio', 'height', 'blade_width', 'leaf_thickness', 'spad', 'canopy_diameter', 'water_potential', 'photosynthetic_rate', 'stomatal_conductance', 'internal_co2', 'transpiration_rate')

	traits_columns <- trait_columns[match(trait, nice_trait)]

	# formula
	terms <- terms(formula_traits)
	terms <- attr(terms, 'term.labels')

	covariates <- terms
	covariates <- covariates[!grepl(covariates, pattern = '\\^2')]
	covariates <- covariates[!grepl(covariates, pattern = '\\:')]
	covariates <- covariates[!grepl(covariates, pattern = '\\*')]

	n_covariates <- length(covariates)

	# load site data
	site_data_raw <- readRDS('./data_from_loretta/sdm_pdm_00_merged_site_data_with_climate/sites.rds')
	biomass_data_raw <- readRDS('./data_from_loretta/sdm_pdm_00_merged_site_data_with_climate/biomass.rds')
	data_raw <- readRDS('./data_from_loretta/sdm_pdm_00_merged_site_data_with_climate/morpho_phys.rds')

	# ensure sites in each array appear in the same order
	site_data_raw <- site_data_raw[order(site_id)]
	biomass_data_raw <- biomass_data_raw[order(SITE)]
	data_raw <- data_raw[order(SITE)]

	# get just columns of interest--some traits have NAs for some plants
	columns <- c('SITE', 'PLANT', traits_columns, paste0('bio', 1:19), 'aridity')
	data_raw <- data_raw[ , ..columns]
	data_raw <- data_raw[complete.cases(data_raw)]

	stopifnot(all(site_data_raw$site_id == unique(biomass_data_raw$SITE)))
	stopifnot(all(site_data_raw$site_id == unique(data_raw$SITE)))

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
	site_data_raw$site_soil_pH <- soil_data_raw$pH[matches]
	site_data_raw$site_soil_pH_sikora <- soil_data_raw$`Sikora_pH`[matches]
	site_data_raw$site_soil_N <- soil_data_raw$`Total.N.%`[matches]

	traits_data_raw_site <- gsub(tolower(data_raw$SITE), pattern = '_', replacement = '')
	matches <- match(traits_data_raw_site, soil_data_raw_site)
	data_raw$site_soil_pH <- soil_data_raw$pH[matches]
	data_raw$site_soil_pH_sikora <- soil_data_raw$`Sikora_pH`[matches]
	data_raw$site_soil_N <- soil_data_raw$`Total.N.%`[matches]

	# add soil texture variables from field samples
	site_data_raw$sand <- site_data_raw$SAND / 100
	site_data_raw$silt <- site_data_raw$SILT / 100
	site_data_raw$clay <- site_data_raw$CLAY / 100

	data_raw[ , c('sand', 'silt', 'clay') := NA_real_]
	for (i in 1:nrow(site_data_raw)) {

		data_raw$sand[data_raw$SITE == site_data_raw$site_id[i]] <- site_data_raw$sand[i]
		data_raw$silt[data_raw$SITE == site_data_raw$site_id[i]] <- site_data_raw$silt[i]
		data_raw$clay[data_raw$SITE == site_data_raw$site_id[i]] <- site_data_raw$clay[i]

	}

	n_pheno_sites <- length(unique(site_data_raw$site_id))
	
	### create spatial versions of site and phenotype data for plotting
	site_data <- merge(site_data_raw[ , c('site_id', 'LONGITUDE', 'LATITUDE')], data_raw, by.x = 'site_id', by.y = 'SITE')
	
	# means and sd
	site_data_agg <- data.table(site = unique(site_data$site_id))
	for (i in seq_along(traits_columns)) {
	
		trait_column <- traits_columns[i]

		site_data_agg[ , c('this_mean', 'this_sd') := NA_real_]
		for (site in site_data_agg$site) {
		
			x <- data_raw[SITE == site, ..trait_column]
			x <- x[[1]]
			mu <- mean(x)
			sigma <- sd(x)

			site_data_agg$this_mean[site_data_agg$site == site] <- mu
			site_data_agg$this_sd[site_data_agg$site == site] <- sigma


		}

		names(site_data_agg)[names(site_data_agg) == 'this_mean'] <- paste0(trait[i], '_mean')
		names(site_data_agg)[names(site_data_agg) == 'this_sd'] <- paste0(trait[i], '_sd')

	}
	site_data_agg <- merge(site_data_agg, site_data_raw[ , c('site_id', 'LONGITUDE', 'LATITUDE')], by.x = 'site', by.y = 'site_id')
	site_vect <- terra::vect(site_data_agg, geom = c('LONGITUDE', 'LATITUDE'), crs = enmSdmX::getCRS('WGS84'))
	
	### means and means of sds for model initialization
	cols <- paste0(trait, '_mean')
	site_means <- colMeans(site_data_agg[ , ..cols])
	cols <- paste0(trait, '_sd')
	site_sds <- apply(site_data_agg[ , ..cols], 2, sd)

	### center and scale at site level
	x <- site_data_raw
	names(x)[names(x) == 'site_soil_pH'] <- 'ph'
	names(x)[names(x) == 'SAND'] <- 'sand'
	names(x)[names(x) == 'SILT'] <- 'silt'
	names(x)[names(x) == 'CLAY'] <- 'clay'
	x <- x[ , ..covariates]
	x <- scale(x)
	x_centers <- attr(x, 'scaled:center')
	x_scales <- attr(x, 'scaled:scale')
	x <- as.data.frame(x)
	x_by_site <- model.matrix(formula_traits, x)

	n_terms <- ncol(x_by_site)

	# make vector of which sampled site matches each row in the morphology/physiology data
	n <- nrow(data_raw)

	site_index <- rep(NA, n)
	for (i in 1:n) site_index[i] <- which(site_data_raw$site_id == data_raw$SITE[i])

	### county-level environmental data
	ag_vect_sq <- vect('./outputs_loretta/integrated_sdm_pdm/andropogon_gerardi_occurrences_with_environment_1961_2020_for_integration.gpkg')
	ag_sq <- as.data.frame(ag_vect_sq)
	ag_sq <- ag_sq[ , covariates, drop = FALSE]
	ag_sq <- scale(ag_sq, center = x_centers[covariates], scale = x_scales[covariates])
	ag_sq <- as.data.frame(ag_sq)
	counties_x_sq <- model.matrix(formula_traits, ag_sq)

	n_counties <- nrow(ag_sq)

	site_vect <- project(site_vect, ag_vect_sq)

	# define focal region for defining range of environment for response curves, and region for calibration of model
	if (calib) {
		ag_vect_sq$calib_region <- !is.na(ag_vect_sq$n_poaceae)
	} else {
		ag_vect_sq$calib_region <- TRUE
	}
	index_traits_county_in_calib_county  <- which(ag_vect_sq$calib_region)

	### future environmental data
	for (fut in futs) {
	
		this_fut <- vect(paste0('./data_from_adam_and_loretta/andropogon_gerardi_occurrences_with_environment_ensemble_8GCMs_', fut, '_climatena.gpkg'))

		if (n_covariates > 0) this_fut <- this_fut[ , covariates]
		assign(paste0('counties_', fut), this_fut)

		this_fut <- as.data.frame(this_fut)[ , covariates, drop = FALSE]
		this_fut <- scale(this_fut, center = x_centers[covariates], scale = x_scales[covariates])
		this_fut <- as.data.frame(this_fut)
		this_fut <- model.matrix(formula_traits, this_fut)

		assign(paste0('counties_x_', fut), this_fut)
	
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
		mm <- model.matrix(formula_traits, mm)

		assign(paste0('counties_x_', clim), mm)
	
	}
	n_counties_20th_cent <- nrow(thirties)

	### response array
	resp_arrays <- create_response_curve_array(
		formula = formula_traits,
		centers = x_centers,
		scales = x_scales,
		ag_vect_sq = ag_vect_sq[ag_vect_sq$n_andropogon_gerardi > 0],
		vects = list(
			counties_ssp245_2041_2070[ag_vect_sq$n_andropogon_gerardi > 0],
			counties_ssp245_2071_2100[ag_vect_sq$n_andropogon_gerardi > 0],
			counties_ssp370_2041_2070[ag_vect_sq$n_andropogon_gerardi > 0],
			counties_ssp370_2071_2100[ag_vect_sq$n_andropogon_gerardi > 0],
			thirties, fifties
		),
		n_response_curve_values = n_response_curve_values
	)

	### responses
	y <- data_raw[[traits_columns]]

	list(

		trait = trait,

		site_data_raw = site_data_raw,			# raw site data
		raw_data_traits = data_raw,		# raw data for traits

		terms_traits = terms,					# terms in formula_traits
		covariates_traits = covariates,			# covariates in formula_traits (which can appear in >1 term)
		n_terms_traits = n_terms,				# number of terms in formula_traits
		n_covariates_traits = n_covariates,		# number of variables in formula_traits (which can appear in >1 term)
	
		y_traits = y,							# response variable (traits)
		n_trait_values = length(y),				# number of observations in response variable
		site_index_traits = site_index,			# index of which values in the response variable correspond to which sampled site
		plant = data_raw$PLANT,
		n_pheno_sites = n_pheno_sites,			# number of sites at which phenotypes were measured

		site_means = site_means,				# mean site-level trait values
		site_sds = site_sds,					# mean of site-level traits' s.d.

		x_by_site_traits = x_by_site,			# MM with covariates (scaled)
		x_centers_traits = x_centers,			# vector of means for each covariate
		x_scales_traits = x_scales,				# vector of standard deviations for each covariate

		site_vect_traits = site_vect,			# SpatVector of site data for plotting

		index_traits_county_in_calib_county = index_traits_county_in_calib_county, # index of counties in SDM calibration region

		n_counties = n_counties,				# number of counties in dataset
		# counties_x_traits_sq_calib = counties_x_sq_calib, # MM of county-level environmental data for present (SDM calibration region)
		counties_x_traits_sq = counties_x_sq,		# MM of county-level environmental data for current conditions
		counties_x_traits_ssp245_2041_2070 = counties_x_ssp245_2041_2070, # MM of county-level environmental data for future conditions (SSP245, 2041-2070)
		counties_x_traits_ssp245_2071_2100 = counties_x_ssp245_2071_2100, # MM of county-level environmental data for future conditions (SSP245, 2071-2100)
		counties_x_traits_ssp370_2041_2070 = counties_x_ssp370_2041_2070, # MM of county-level environmental data for future conditions (SSP370, 2041-2070)
		counties_x_traits_ssp370_2071_2100 = counties_x_ssp370_2071_2100, # MM of county-level environmental data for future conditions (SSP370, 2071-2100)

		n_counties_20th_cent = n_counties_20th_cent, # SpatVector of county-level environmental data for 20th century
		counties_thirties = thirties, # SpatVector of county-level environmental data for 20th century
		counties_fifties = fifties, # SpatVector of county-level environmental data for 20th century

		counties_x_thirties = counties_x_thirties, # model matrix of county-level environmental data for 20th century
		counties_x_fifties = counties_x_fifties, # model matrix of county-level environmental data for 20th century

		resp_curves_x_traits = resp_arrays$response_curves_x_scaled,			# model matrices for response curves
		resp_curve_x_traits_unscaled = resp_arrays$resp_curves_x_unscaled	# matrices for response curves

	)

}

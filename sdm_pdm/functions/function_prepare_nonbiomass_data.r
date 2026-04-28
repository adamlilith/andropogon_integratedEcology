#' This function loads and prepares non-biomass facet (morphological and physiological) data for distribution modeling. See the "return()` line for details on the output.
#' 
#' source('C:/Kaji/R/andropogon_integratedEcology/sdm_pdm/functions/function_prepare_nonbiomass_data.r')
#' 
#' facet One of:
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
#' formula_facet			A formula_facet object specifying the model to be fit. Include an intercept and the RHS only.
#' n_response_curve_values	A numeric value specifying the number of values to use for the response curve. Default is 200.
#' calib					If `TRUE`, then the training data is subset only to counties with non-`NA` for Poaceae.
#'
prepare_nonbiomass_data <- function(facet, formula_facet, n_response_curve_values = 200, calib = FALSE) {

	if (facet == 'delta13c') stop('Cannot use `delta13c` because some values are NA.')

	# get traits
	trait_columns <- c('Delta13C', 'N_conc', 'CN_ratio', 'Height', 'BladeWidth', 'LeafThick', 'SPAD', 'CanopyDiam', 'WatPot', 'PhotoRate', 'StomCond', 'IntCO2', 'TranspRate')
	nice_trait <- c('delta13c', 'n_concentration', 'cn_ratio', 'height', 'blade_width', 'leaf_thickness', 'spad', 'canopy_diameter', 'water_potential', 'photosynthetic_rate', 'stomatal_conductance', 'internal_co2', 'transpiration_rate')

	traits_columns <- trait_columns[match(facet, nice_trait)]
	facet_raw <- get_raw_trait_name_from_rfriendly(facet)

	# formula
	terms <- terms(formula_facet)
	terms <- attr(terms, 'term.labels')

	covariates <- terms
	covariates <- covariates[!grepl(covariates, pattern = '\\^2')]
	covariates <- covariates[!grepl(covariates, pattern = '\\:')]
	covariates <- covariates[!grepl(covariates, pattern = '\\*')]

	n_covariates <- length(covariates)

	# load site data
	site_data_raw <- readRDS('./data_from_loretta/sdm_pdm_00_merged_site_data_with_climate/sites.rds')
	raw_data_facet <- readRDS('./data_from_loretta/sdm_pdm_00_merged_site_data_with_climate/morpho_phys.rds')

	site_data_raw <- calculate_logged_vars(site_data_raw)
	raw_data_facet <- calculate_logged_vars(raw_data_facet)

	raw_data_facet[ , insolation_2000_growing_season_kWh_per_m2 := site_data_raw$insolation_2000_growing_season_kWh_per_m2[match(raw_data_facet$SITE, site_data_raw$site_id)]]
	raw_data_facet[ , insolation_2023_growing_season_kWh_per_m2 := site_data_raw$insolation_2023_growing_season_kWh_per_m2[match(raw_data_facet$SITE, site_data_raw$site_id)]]

	# ensure sites in each array appear in the same order
	site_data_raw <- site_data_raw[order(site_id)]
	raw_data_facet <- raw_data_facet[order(SITE)]

	# get just columns of interest--some traits have NAs for some plants
	columns <- c('SITE', 'PLANT', traits_columns, paste0('bio', 1:19), paste0('bio', c(12:14, 16:19), '_log10p1'), 'aridity', 'insolation_2000_growing_season_kWh_per_m2', 'insolation_2023_growing_season_kWh_per_m2', 'geofold')
	raw_data_facet <- raw_data_facet[ , ..columns]
	raw_data_facet <- raw_data_facet[complete.cases(raw_data_facet)]

	stopifnot(all(site_data_raw$site_id == unique(raw_data_facet$SITE)))
	stopifnot(all(site_data_raw$site_id == unique(raw_data_facet$SITE)))

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

	traits_data_raw_site <- gsub(tolower(raw_data_facet$SITE), pattern = '_', replacement = '')
	matches <- match(traits_data_raw_site, soil_data_raw_site)
	raw_data_facet$site_ph <- soil_data_raw$pH[matches]
	raw_data_facet$site_ph_sikora <- soil_data_raw$`Sikora_pH`[matches]
	# raw_data_facet$site_nitrogen <- soil_data_raw$`Total.N.%`[matches]
	raw_data_facet$nitrogen <- soil_data_raw$`Total.N.%`[matches]

	# add soil texture variables from field samples
	site_data_raw$sand <- site_data_raw$SAND
	site_data_raw$silt <- site_data_raw$SILT
	site_data_raw$clay <- site_data_raw$CLAY

	raw_data_facet[ , c('sand', 'silt', 'clay') := NA_real_]
	for (i in 1:nrow(site_data_raw)) {

		raw_data_facet$sand[raw_data_facet$SITE == site_data_raw$site_id[i]] <- site_data_raw$sand[i]
		raw_data_facet$silt[raw_data_facet$SITE == site_data_raw$site_id[i]] <- site_data_raw$silt[i]
		raw_data_facet$clay[raw_data_facet$SITE == site_data_raw$site_id[i]] <- site_data_raw$clay[i]

	}

	n_pheno_sites <- length(unique(site_data_raw$site_id))
	
	### create spatial versions of site and phenotype data for plotting
	site_data <- merge(site_data_raw[ , c('site_id', 'LONGITUDE', 'LATITUDE')], raw_data_facet, by.x = 'site_id', by.y = 'SITE')
	
	# means and sd
	site_data_agg <- data.table(site = unique(site_data$site_id))
	for (i in seq_along(traits_columns)) {
	
		trait_column <- traits_columns[i]

		site_data_agg[ , c('this_mean', 'this_sd') := NA_real_]
		for (site in site_data_agg$site) {
		
			x <- raw_data_facet[SITE == site, ..trait_column]
			x <- x[[1]]
			mu <- mean(x)
			sigma <- sd(x)

			site_data_agg$this_mean[site_data_agg$site == site] <- mu
			site_data_agg$this_sd[site_data_agg$site == site] <- sigma


		}

		names(site_data_agg)[names(site_data_agg) == 'this_mean'] <- paste0(facet[i], '_mean')
		names(site_data_agg)[names(site_data_agg) == 'this_sd'] <- paste0(facet[i], '_sd')

	}
	site_data_agg <- merge(site_data_agg, site_data_raw[ , c('site_id', 'LONGITUDE', 'LATITUDE')], by.x = 'site', by.y = 'site_id')
	site_vect <- terra::vect(site_data_agg, geom = c('LONGITUDE', 'LATITUDE'), crs = enmSdmX::getCRS('WGS84'))
	
	### means and sds
	cols <- paste0(facet, '_mean')
	site_means <- site_data_agg[[paste0(facet, '_mean')]]
	site_sds <- site_data_agg[[paste0(facet, '_sd')]]

	### center and scale at site level
	x <- site_data_raw
	names(x)[names(x) == 'site_ph'] <- 'ph'
	names(x)[names(x) == 'SAND'] <- 'sand'
	names(x)[names(x) == 'SILT'] <- 'silt'
	names(x)[names(x) == 'CLAY'] <- 'clay'
	x <- calculate_logged_vars(x)
	x <- x[ , ..covariates]

	centers_scales <- fread('./outputs_loretta/integrated_sdm_pdm/centers_and_scales_for_covariates.csv')
	x_centers <- centers_scales$center_sites[match(covariates, centers_scales$variable)]
	x_scales <- centers_scales$scale_sites[match(covariates, centers_scales$variable)]
	names(x_centers) <- covariates
	names(x_scales) <- covariates

	x <- scale(x, center = x_centers, scale = x_scales)
	x <- as.data.frame(x)
	x_by_site <- model.matrix(formula_facet, x)

	n_terms <- ncol(x_by_site)

	# make vector of which sampled site matches each row in the morphology/physiology data
	n <- nrow(raw_data_facet)

	site_index <- rep(NA, n)
	for (i in 1:n) site_index[i] <- which(site_data_raw$site_id == raw_data_facet$SITE[i])

	### county-level environmental data
	ag_vect_sq <- vect('./outputs_loretta/integrated_sdm_pdm/andropogon_gerardi_occurrences_with_environment_1961_2020_for_integration.gpkg')
	ag_vect_sq <- simplifyGeom(ag_vect_sq, tolerance = 1000)
	# names(ag_vect_sq)[names(ag_vect_sq) == 'nitrogen'] <- 'site_nitrogen'
	ag_vect_sq <- calculate_logged_vars(ag_vect_sq)

	ag_sq <- as.data.frame(ag_vect_sq)
	ag_sq <- ag_sq[ , covariates, drop = FALSE]

	ag_sq <- scale(ag_sq, center = x_centers[covariates], scale = x_scales[covariates])
	ag_sq <- as.data.frame(ag_sq)
	counties_x_sq <- model.matrix(formula_facet, ag_sq)

	n_counties <- nrow(ag_sq)
	site_vect <- project(site_vect, ag_vect_sq)

	# define focal region for defining range of environment for response curves, and region for calibration of model
	if (calib) {
		ag_vect_sq$calib_region <- !is.na(ag_vect_sq$n_poaceae)
	} else {
		ag_vect_sq$calib_region <- TRUE
	}
	index_facet_county_in_calib_county  <- which(ag_vect_sq$calib_region)

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
		this_fut <- model.matrix(formula_facet, this_fut)

		assign(paste0('counties_x_', fut), this_fut)
	
	}

	# 1930s climate
	thirties <- vect(paste0('./data_from_adam_and_loretta/andropogon_gerardi_occurrences_with_environment_1931_1940_prism.gpkg'))	
	thirties <- simplifyGeom(thirties, tolerance = 1000)
	thirties <- calculate_logged_vars(thirties)
	# names(thirties)[names(thirties) == 'nitrogen'] <- 'site_nitrogen'

	thirties <- thirties[ , c('country', 'state_province', 'county', covariates)]

	for (covariate in covariates) {
		
		x <- thirties[[covariate]]
		x <- unlist(x)

		x <- scale(x, center = x_centers[covariate], scale = x_scales[covariate])
		x <- as.numeric(x)
		thirties[ , covariate] <- x

	}

	mm <- as.data.frame(thirties)[ , covariates, drop = FALSE] 
	mm <- model.matrix(formula_facet, mm)

	counties_x_thirties <- mm
	n_counties_1930s <- nrow(thirties)

	### response array
	resp_arrays <- create_response_curve_array(
		formula = formula_facet,
		centers = x_centers,
		scales = x_scales,
		site_data_raw = site_data_raw,
		ag_vect_sq = ag_vect_sq[ag_vect_sq$n_andropogon_gerardi > 0],
		vects = list(
			counties_ssp245_2041_2070[ag_vect_sq$n_andropogon_gerardi > 0],
			counties_ssp245_2071_2100[ag_vect_sq$n_andropogon_gerardi > 0],
			counties_ssp370_2041_2070[ag_vect_sq$n_andropogon_gerardi > 0],
			counties_ssp370_2071_2100[ag_vect_sq$n_andropogon_gerardi > 0]
		),
		n_response_curve_values = n_response_curve_values
	)

	### responses
	y <- raw_data_facet[[traits_columns]]

	list(

		facet = facet,
		facet_raw = facet_raw,					# Name of facet in original data sheets

		site_data_raw = site_data_raw,			# raw site data
		raw_data_facet = raw_data_facet,		# raw data for traits

		terms = terms,							# terms in formula_facet
		covariates = covariates,				# covariates in formula_facet (which can appear in >1 term)
		n_terms = n_terms,						# number of terms in formula_facet
		n_covariates = n_covariates,			# number of variables in formula_facet (which can appear in >1 term)
	
		y_facet = y,							# response variable (traits)
		n_plants = length(y),					# number of observations (plants) in response variable
		site_index_facet = site_index,			# index of which values in the response variable correspond to which sampled site
		plant = raw_data_facet$PLANT,
		n_plants_per_site_facet = length(unique(raw_data_facet$PLANT)), # number of phenotyped plants per site (ASSUMING SAME ACROSS SITES!)
		n_pheno_sites = n_pheno_sites,			# number of sites at which phenotypes were measured

		site_means = site_means,				# mean site-level facet values
		site_sds = site_sds,					# mean of site-level traits' s.d.

		x_by_site = x_by_site,			# MM with covariates (scaled)
		x_centers = x_centers,			# vector of means for each covariate
		x_scales = x_scales,				# vector of standard deviations for each covariate

		site_vect = site_vect,			# SpatVector of site data for plotting

		index_facet_county_in_calib_county = index_facet_county_in_calib_county, # index of counties in SDM calibration region

		n_counties = n_counties,				# number of counties in dataset
		# counties_x_facet_sq_calib = counties_x_sq_calib, # MM of county-level environmental data for present (SDM calibration region)
		counties_x_sq = counties_x_sq,		# MM of county-level environmental data for current conditions
		counties_x_ssp245_2041_2070 = counties_x_ssp245_2041_2070, # MM of county-level environmental data for future conditions (SSP245, 2041-2070)
		counties_x_ssp245_2071_2100 = counties_x_ssp245_2071_2100, # MM of county-level environmental data for future conditions (SSP245, 2071-2100)
		counties_x_ssp370_2041_2070 = counties_x_ssp370_2041_2070, # MM of county-level environmental data for future conditions (SSP370, 2041-2070)
		counties_x_ssp370_2071_2100 = counties_x_ssp370_2071_2100, # MM of county-level environmental data for future conditions (SSP370, 2071-2100)

		n_counties_1930s = n_counties_1930s, # SpatVector of county-level environmental data for 20th century
		counties_thirties = thirties, # SpatVector of county-level environmental data for 20th century
		counties_x_thirties = counties_x_thirties, # model matrix of county-level environmental data for 20th century

		resp_curves_x = resp_arrays$response_curves_x_scaled,			# model matrices for response curves
		resp_curves_x_unscaled = resp_arrays$resp_curves_x_unscaled	# matrices for response curves

	)

}

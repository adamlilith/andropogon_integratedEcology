### MODELING ANDROPOGON GERARDI DISTRIBUTION, PHENOTYPE, PHYSIOLOGY, GENOTYPE, and ASSOCIATED MICROBIAL COMMUNITIES
### Erica Newman | Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2023-12
###
### This script constructs an integrated SDM-genetics model using occurrence data from Smith et al. (2017 GEB) for the SDM and results from an ADMIXTURE analysis by Jack Sytsma and Loretta Johnson.
###
### source('C:/Ecology/R/andropogon_integratedEcology/sdm_admixture/admixture_01_sdm_species_level_priors.r')
### source('C:/Subarashi/R/andropogon_integratedEcology/sdm_admixture/admixture_01_sdm_species_level_priors.r')
###
### CONTENTS ###
### setup ###
### integrated SDM-ADMIXTURE model ###
### model diagnostics for SDM component ###
### model diagnostics for ADMIXTURE component ###
### response curves ###
### instill predictions into spatial vectors ###
### map of CURRENT species distribution ###
### maps of CURRENT species and ADMIXTURE distributions ###
### maps of FUTURE species and ADMIXTURE distributions ###
### parameter estimates for SDM ###
### parameter estimates for ADMIXTURE ###
### response curves ###

#############
### setup ###
#############

	rm(list = ls())

	drive <- 'C:/Ecology/'
	# drive <- 'C:/Subarashi/'

	setwd(paste0(drive, '/Research/Andropogon/Andropogon'))

	library(bayesplot) # graphing
	library(coda) # Bayesian diagnostics
	library(cowplot) # combining ggplots
	library(enmSdmX) # SDMing
	library(ggplot2) # plotting
	library(ggnewscale) # change scales in ggplots
	library(ggspatial) # plotting spatial things
	library(scatterpie) # pie charts
	library(nimble) # Bayes
	library(nimbleHMC) # Hamiltonian Monte Carlo samplers
	library(omnibus)
	library(scales) # for plotting transparency
	library(terra) # spatial objects

	### user-defined values
	#######################

		# should we use a subsample of the data (for development--to make things go faster)
		abbreviated <- FALSE
		# abbreviated <- TRUE

		ok <- function() {
			say('', pre = 2)
			x <- readline(paste0('`abbreviated` is ', abbreviated, '. Is this OK? (y/n) '))
			if (x != 'y') stop('All your base are belong to us.')
		}
		ok()

		### output folder (base... parts added depending on decisions below)
		if (abbreviated) {

			out_dir <- paste0('./outputs_loretta/admixture_[TEST]')

		} else {
			
			out_dir <- paste0('./outputs_loretta/admixture_[sdm_nmixture]_[admixture_priors_from_species]')

		}

		psa_quant <- 0.95 # to define pseudoabsences, use this quantile of Poaceae occurrences across counties with Poaceae occurrences 

		K_start <- 7 # total number of ADMIXTURE populations (before removing rare genotypes)

		### include the THREE rare genotypes mainly found at just on site each?
		# include_rare_geno <- TRUE
		include_rare_geno <- FALSE

		if (!include_rare_geno) {
			out_dir <- paste0(out_dir, '_[remove_rare_geno]')
		} else {
			out_dir <- paste0(out_dir, '_[include_rare_geno]')
		}

		### include soil predictors?
		# include_soil <- TRUE
		include_soil <- FALSE

		if (include_soil) {
			out_dir <- paste0(out_dir, '_[climate2_soil2]/')
		} else {
			out_dir <- paste0(out_dir, '_[climate2]/')
		}

		### MCMC settings
		if (!abbreviated) {

			niter <- 130000
			nburnin <- 10000
			thin <- 120
			nchains <- 4

		} else {

			# for testing
			niter <- 220
			nburnin <- 20
			thin <- 1
			nchains <- 2

			say('ABBREVIATED!!!', level = 2)

		}

	dirCreate(out_dir)
	dirCreate(out_dir, 'diagnostics')

	sink(paste0(out_dir, '/runtime_log.txt'), split = TRUE)
	say()

say('######################################')
say('### integrated SDM-ADMIXTURE model ###')
say('######################################')
say(date(), post = 1)

	say('MCMC settings:', level = 2)
	say('abbreviated ......... ', abbreviated)
	say('niter ............... ', niter)
	say('nburnin ............. ', nburnin)
	say('thin ................ ', thin)
	say('nchains ............. ', nchains)
	say('psa_quant ........... ', psa_quant)
	say('include_soil ........ ', include_soil)
	say('include_rare_geno ... ', include_soil, post = 1)

	### model predictors and terms
	say('Predictors and model formula:', level = 2)

	climate_predictor_names <- c('aridity', 'bio7', 'bio12', 'bio15')
	soil_predictor_names <- c('ph', 'sand', 'silt')

	predictor_names <- if (include_soil) {
		c(climate_predictor_names, soil_predictor_names)
	} else {
		climate_predictor_names
	}

	say('We are using predictors: ', paste(predictor_names, collapse = ' '), post = 1)

	form <- c(1, predictor_names, paste0('I(', predictor_names, '^2)'))

	# combos <- combn(predictor_names, 2, simplify = FALSE)

	# for (i in seq_along(combos)) form <- c(form, paste(combos[[i]], collapse = ':'))
	form <- paste(form, collapse = ' + ')
	form <- paste0('~ ', form)

	# create model formula
	# form <- paste0(' ~ 1 + ', paste(c(linear_terms, terms), collapse = ' + '))
	say('Model formula:')
	say(form, post = 1, breaks = 80)
	form <- as.formula(form)

	### collate occurrence data for SDM
	###################################

	### load AG data
	ag_vect_sq <- vect('./data_from_adam_and_loretta/andropogon_gerardi_occurrences_with_environment_1961_2020.gpkg')

	### spatial vectors with future values of predictors
	ag_vect_ssp245_2041_2070 <- vect('./data_from_adam_and_loretta/andropogon_gerardi_occurrences_with_environment_ensemble_8GCMs_ssp245_2041_2070.gpkg')
	ag_vect_ssp245_2071_2100 <- vect('./data_from_adam_and_loretta/andropogon_gerardi_occurrences_with_environment_ensemble_8GCMs_ssp245_2071_2100.gpkg')
	ag_vect_ssp370_2041_2070 <- vect('./data_from_adam_and_loretta/andropogon_gerardi_occurrences_with_environment_ensemble_8GCMs_ssp370_2041_2070.gpkg')
	ag_vect_ssp370_2071_2100 <- vect('./data_from_adam_and_loretta/andropogon_gerardi_occurrences_with_environment_ensemble_8GCMs_ssp370_2071_2100.gpkg')

	### soil predictors
	if (include_soil) {

		ag_vect_soil <- vect('./data_from_adam_and_loretta/andropogon_gerardi_occurrences_with_environment_soil.gpkg')
	
		soil_columns <- as.data.frame(ag_vect_soil)
		soil_columns <- soil_columns[ , soil_predictor_names]
		
		ag_vect_sq <- cbind(ag_vect_sq, soil_columns)
		ag_vect_ssp245_2041_2070 <- cbind(ag_vect_ssp245_2041_2070, soil_columns)
		ag_vect_ssp245_2071_2100 <- cbind(ag_vect_ssp245_2071_2100, soil_columns)
		ag_vect_ssp370_2041_2070 <- cbind(ag_vect_ssp370_2041_2070, soil_columns)
		ag_vect_ssp370_2071_2100 <- cbind(ag_vect_ssp370_2071_2100, soil_columns)

	}

	fields <- c('country', 'stateProvince', 'county', 'area_km2', 'any_ag_quality1to3', 'num_poaceae_records', 'elevation_m', predictor_names)
	ag_vect_sq <- ag_vect_sq[ , fields]

	### pseudoabsences
	# for counties with NA Poaceae and AG, assign a maximal number of Poaceae and 0 AG
	n_pseudoabs <- quantile(ag_vect_sq$num_poaceae_records, psa_quant, na.rm = TRUE)
	ag_vect_sq$num_poaceae_records[is.na(ag_vect_sq$num_poaceae_records)] <- n_pseudoabs
	ag_vect_sq$any_ag_quality1to3[is.na(ag_vect_sq$any_ag_quality1to3)] <- 0

	if (abbreviated) {
	
		sampled <- sample(nrow(ag_vect_sq), 500)
		ag_vect_sq <- ag_vect_sq[sampled]
		
		ag_vect_ssp245_2041_2070 <- ag_vect_ssp245_2041_2070[sampled]
		ag_vect_ssp245_2071_2100 <- ag_vect_ssp245_2071_2100[sampled]

		ag_vect_ssp370_2041_2070 <- ag_vect_ssp370_2041_2070[sampled]
		ag_vect_ssp370_2071_2100 <- ag_vect_ssp370_2071_2100[sampled]

	}

	ag_sq <- as.data.frame(ag_vect_sq)

	### county area... used as an offset to make underlying model fit an IPP
	area_km2 <- ag_sq$area_km2
	log_area_km2 <- log(area_km2)
	log_area_km2_scaled <- scale(log_area_km2)
	log_area_center <- attributes(log_area_km2_scaled)$`scaled:center`
	log_area_scale <- attributes(log_area_km2_scaled)$`scaled:scale`
	log_area_km2_scaled <- log_area_km2_scaled[ , 1]

	### number of Poaceae records... used to model sampling bias
	log_num_poaceae_records <- log1p(ag_sq$num_poaceae_records) # log(occs_x_sq + 1)
	log_num_poaceae_records_scaled <- scale(log_num_poaceae_records)
	log_num_poaceae_records_scaled <- log_num_poaceae_records_scaled[ , 1]

	### subset, scale, and manipulate predictors into model frame
	# status quo: counties with training data
	occs_x_raw <- as.data.frame(ag_sq[ , c('area_km2', predictor_names)])
	log_area_km2 <- log(occs_x_raw$area_km2)
	log_area_km2 <- data.frame(log_area_km2 = log_area_km2)
	occs_x_raw$area_km2 <- NULL
	occs_x_raw <- cbind(log_area_km2, occs_x_raw)
	occs_x_raw_scaled <- scale(occs_x_raw)
	occs_x_centers <- attr(occs_x_raw_scaled, 'scaled:center')
	occs_x_scales <- attr(occs_x_raw_scaled, 'scaled:scale')
	occs_x_raw_scaled <- as.data.frame(occs_x_raw_scaled)
	occs_x_sq <- model.matrix(form, as.data.frame(occs_x_raw_scaled))

	# ssp245_2041_2070
	x_fut <- ag_vect_ssp245_2041_2070
	x_fut <- as.data.frame(x_fut[ , c('area_km2', predictor_names)])
	log_area_km2 <- log(x_fut$area_km2)
	log_area_km2 <- data.frame(log_area_km2 = log_area_km2)
	x_fut$area_km2 <- NULL
	x_fut <- cbind(log_area_km2, x_fut)
	x_fut <- scale(x_fut, center = occs_x_centers, scale = occs_x_scales)
	x_fut <- as.data.frame(x_fut)
	x_ssp245_2041_2070 <- model.matrix(form, as.data.frame(x_fut))

	# ssp245_2071_2100
	x_fut <- ag_vect_ssp245_2071_2100
	x_fut <- as.data.frame(x_fut[ , c('area_km2', predictor_names)])
	log_area_km2 <- log(x_fut$area_km2)
	log_area_km2 <- data.frame(log_area_km2 = log_area_km2)
	x_fut$area_km2 <- NULL
	x_fut <- cbind(log_area_km2, x_fut)
	x_fut <- scale(x_fut, center = occs_x_centers, scale = occs_x_scales)
	x_fut <- as.data.frame(x_fut)
	x_ssp245_2071_2100 <- model.matrix(form, as.data.frame(x_fut))

	# ssp370_2041_2070
	x_fut <- ag_vect_ssp370_2041_2070
	x_fut <- as.data.frame(x_fut[ , c('area_km2', predictor_names)])
	log_area_km2 <- log(x_fut$area_km2)
	log_area_km2 <- data.frame(log_area_km2 = log_area_km2)
	x_fut$area_km2 <- NULL
	x_fut <- cbind(log_area_km2, x_fut)
	x_fut <- scale(x_fut, center = occs_x_centers, scale = occs_x_scales)
	x_fut <- as.data.frame(x_fut)
	x_ssp370_2041_2070 <- model.matrix(form, as.data.frame(x_fut))

	# ssp245_2071_2100
	x_fut <- ag_vect_ssp370_2071_2100
	x_fut <- as.data.frame(x_fut[ , c('area_km2', predictor_names)])
	log_area_km2 <- log(x_fut$area_km2)
	log_area_km2 <- data.frame(log_area_km2 = log_area_km2)
	x_fut$area_km2 <- NULL
	x_fut <- cbind(log_area_km2, x_fut)
	x_fut <- scale(x_fut, center = occs_x_centers, scale = occs_x_scales)
	x_fut <- as.data.frame(x_fut)
	x_ssp370_2071_2100 <- model.matrix(form, as.data.frame(x_fut))

	### collate data for ADMIXTURE component
	########################################

	admix <- read.csv('./data_from_loretta/!admixture_data_collated/admixture_with_env.csv')

	# CRITICAL: to ensure we keep plants of the same population together, we sort by population before doing anything else
	admix <- admix[order(admix$population_proper), ]

	# ### remove rare genotypes
	# # This chunk removes the genetic populations that are rare by finding genetic populations where the highest-ranked proportional representation across sites is >0.9 more than the second-ranked proportional representation across sites (e.g., if the site where that genetic population is most dominant has a proportion of that genetic population that is 0.9 more than the site where the genetic population is second-most abundant, then this genetic population is primarily found in one site). The problem with this approach is that for the sites where a rare genotype dominates (basically, is 100% of the genotypes there, minus a small rounding factor), the remainder of the genetic populations have proportions that are equal and extremely small. Rescaling them causes the resultant proportions to be the same across genetic populations at those sites.
	# if (!include_rare_geno) {
	
		# site_level_geno <- aggregate(admix, by = list(admix$population_proper), FUN = mean)
		# site_level_geno$population_proper <- NULL
		# names(site_level_geno)[1] <- 'population_proper'
	
		# removes <- integer()
		# for (k in 1:K_start) {
			
		# 	props <- site_level_geno[ , paste0('k', k)]
		# 	ranks <- rank(1 - props)

		# 	remove <- props[ranks == 1] - props[ranks == 2] > 0.9
		# 	if (remove) removes <- c(removes, k)

		# }

		# # index of each genetic population
		# K <- {1:K_start}[-removes]

		# # remove rare genotypes
		# removes <- paste0('k', removes)
		# admix <- admix[ , names(admix) %notin% removes]

		# # NB rescaling genotype proportions is done below

	# } else {
	# 	# index of each genetic population
	# 	K <- 1:K_start
	# }

	### remove sites dominated by rare genotypes
	# This chunk 1) identifies sites that are dominated by a single genotype; and 2) genotypes that are found almost exclusively at one site. It then removes *sites* that are dominated by these rare genoytpes.
	if (!include_rare_geno) {

		site_level_geno <- aggregate(admix, by = list(admix$population_proper), FUN = mean)
		site_level_geno$population_proper <- NULL
		names(site_level_geno)[1] <- 'population_proper'

		# identify sites dominated by one genotype
		maxs <- apply(site_level_geno[ , paste0('k', 1:K_start)], 1, max)
		dominated_sites <- which(maxs > 0.99)

		# identify rare genotypes (present in one site)
		rares <- integer()
		for (k in 1:K_start) {
			
			props <- site_level_geno[ , paste0('k', k)]
			ranks <- rank(1 - props)

			is_rare <- props[ranks == 1] - props[ranks == 2] > 0.9
			if (is_rare) rares <- c(rares, k)

		}

		# identify sites dominated by a rare genotype
		remove_sites <- integer()
		for (i in dominated_sites) {
			most_abund_geno <- which.max(site_level_geno[i, paste0('k', 1:K_start)])
			if (most_abund_geno %in% rares) remove_sites <- c(remove_sites, i)
		}
		
		population_proper <- unique(admix$population_proper)
		keep_population_proper <- population_proper[-remove_sites]

		admix <- admix[admix$population_proper %in% keep_population_proper, ]
	
		# index of each genetic population
		K <- 1:K_start

		# NB rescaling genotype proportions is done below

	} else {
		# index of each genetic population
		K <- 1:K_start
	}

	nK <- length(K)

	# rename predictor columns so they match that of occurrence data
	# The ADMIXTURE data is paired with soil variables taken on-site, while the occurrence data is paired with county-wide soil data from SoilGrids 2.0. We assume these represent the same variables in the model (i.e., some coefficients, etc. will be shared).

	# rename predictors so they match columns in the the spatial vectors
	names_admix <- names(admix)
	if (any(names_admix == 'sand_field')) names_admix[names_admix == 'sand_field'] <- 'sand'
	if (any(names_admix == 'silt_field')) names_admix[names_admix == 'silt_field'] <- 'silt'
	if (any(names_admix == 'clay_field')) names_admix[names_admix == 'clay_field'] <- 'clay'
	if (any(names_admix == 'ph_field')) names_admix[names_admix == 'ph_field'] <- 'ph'
	if (any(names_admix == 'nitrogen_field_perc')) names_admix[names_admix == 'nitrogen_field_perc'] <- 'nitrogen'
	if (any(names_admix == 'soc_field_perc')) names_admix[names_admix == 'soc_field_perc'] <- 'soc'
	if (any(names_admix == 'soc_field_perc')) names_admix[names_admix == 'soc_field_perc'] <- 'soc'
	names(admix) <- names_admix

	admix_x_sq <- admix[ , predictor_names]

	# center and scale
	admix_x_centers <- colMeans(admix_x_sq)
	admix_x_scales <- apply(admix_x_sq, 2, sd)

	admix_x_sq <- sweep(admix_x_sq, 2, admix_x_centers, '-')
	admix_x_sq <- sweep(admix_x_sq, 2, admix_x_scales, '/')
	admix_x_sq <- model.matrix(form, admix_x_sq)

	# calculate offsets between centers/scales of occurrence data and ADMIXTURE data for rescaling responses
	delta_occ_admix_centers <- occs_x_centers[predictor_names] - admix_x_centers
	delta_occ_admix_scales <- occs_x_scales[predictor_names] / admix_x_scales

	### response matrix for ADMIXTURE
	admix_y <- admix[ , paste0('k', K)]
	admix_y <- as.matrix(admix_y)
	row_sums <- rowSums(admix_y)
	admix_y <- sweep(admix_y, 1, row_sums, '/')

	n_admix_samples <- nrow(admix_y)

	### study design for ADMIXTURE
	admix_design <- admix[ , 'population_proper', drop = FALSE]
	admix_design$population_proper <- as.factor(admix_design$population_proper)

	n_admix_sites <- length(unique(admix_design$population_proper))

	admix_pop <- model.matrix(~ -1 + population_proper, admix_design)

	### response array
	##################

	# This is an array for plotting the response curves of the SDM and ADMIXTURE proportions
	# columns: variables, all but one held at mean value across counties with AG
	# rows: focal variable increments, others held constant
	# depths: identity of focal variable changes

	# number of rows in array--represents number of points along which we'll be able to make the response curve
	n_predictors <- length(predictor_names)
	n_array_rows <- 200

	env_array <- array(
		NA,
		dim = c(n_array_rows, n_predictors, n_predictors),
		dimnames = list(
			1:n_array_rows,
			predictor_names,
			predictor_names
		)
	)

	# calculate means of each variable across counties with AG presences
	ag_pres <- ag_vect_sq[ag_vect_sq$any_ag_quality1to3 > 0]
	ag_pres <- as.data.frame(ag_pres)[ , predictor_names]
	means <- colMeans(ag_pres)

	for (i in seq_along(predictor_names)) {
		for (j in seq_along(predictor_names)) {
			env_array[ , i, j] <- means[i]
		}
	}

	mins <- apply(ag_pres, 2, min)
	maxs <- apply(ag_pres, 2, max)

	for (i in seq_along(predictor_names)) {
		env_array[ , i, i] <- seq(mins[i], maxs[i], length.out = n_array_rows)
	}

	# scale
	env_array_scaled <- env_array
	for (i in seq_along(predictor_names)) {
		env_array_scaled[ , , i] <- sweep(env_array_scaled[ , , i], 2, occs_x_centers[predictor_names], '-')
		env_array_scaled[ , , i] <- sweep(env_array_scaled[ , , i], 2, occs_x_scales[predictor_names], '/')
	}

	# make into model matrix
	mm_array <- list()
	for (i in seq_along(predictor_names)) {
		mm_array[[i]] <- model.matrix(form, as.data.frame(env_array_scaled[ , , i, drop = TRUE]))
	}
	nc <- ncol(mm_array[[1]])
	response_curve_x <- array(
		as.numeric(unlist(mm_array)),
		dim = c(n_array_rows, nc, n_predictors),
		dimnames = list(1:n_array_rows, colnames(mm_array[[1]]), predictor_names)
	)

	n_response_curve_rows <- nrow(response_curve_x)

	### model inputs
	################

	say('Inputs:', level = 2)
	data <- list(
		occs_y = ag_sq$any_ag_quality1to3, # number of AG records in each county
		admix_y = admix_y # for each plant, proportion of ancestry assignable to each lineage
	)

	n_counties <- nrow(ag_sq)
	n_counties_future <- nrow(x_ssp370_2071_2100)

	n_beta_terms <- ncol(occs_x_sq)

	if (any(names(occs_x_sq) != names(admix_x_sq))) stop('Mismatch between SDM and ADMIXTURE column names in `x`.')

	constants <- list(

		n_counties = n_counties,
		n_counties_future = n_counties_future,

		# admixture data
		admix_pop = admix_pop,
		n_admix_sites = n_admix_sites,
		n_admix_samples = n_admix_samples,
		nK = nK,

		# predictors
		n_beta_terms = n_beta_terms,

		# model matrices
		occs_x_sq = occs_x_sq,
		admix_x_sq = admix_x_sq,
		
		# response curves
		n_predictors = n_predictors,
		n_response_curve_rows = n_response_curve_rows,
		response_curve_x = response_curve_x,

		# sampling bias
		log_area_km2_scaled = log_area_km2_scaled,
		log_num_poaceae_records_scaled = log_num_poaceae_records_scaled,

		# future model matrices
		x_ssp245_2041_2070 = x_ssp245_2041_2070,
		x_ssp245_2071_2100 = x_ssp245_2071_2100,

		x_ssp370_2041_2070 = x_ssp370_2041_2070,
		x_ssp370_2071_2100 = x_ssp370_2071_2100

	)

	# lambda_fut_inits <- rep(mean(ag_sq$any_ag_quality1to3), n_counties_future)
	# inits <- list()
	# for (i in seq_len(nchains)) {

	# 	# N, lambdas, and betas for SDM
	# 	N_inits <- 10 * ag_sq$any_ag_quality1to3
	# 	lambda_sq_inits <- 1 + ag_sq$any_ag_quality1to3
	# 	beta_occs_inits <- runif(n_beta_terms, -0.5, 0.5)

	# 	# betas for ADMIXTURE
	# 	n <- nK * n_beta_terms
	# 	beta_admix_inits <- runif(n, -0.5, 0.5)
	# 	beta_admix_inits <- matrix(beta_admix_inits, nrow = nK, ncol = n_beta_terms)

	# 	# alphas for ADMIXTURE
	# 	n <- nK * n_admix_sites
	# 	alpha_admix_pop_inits <- runif(n, -0.5, 0.5)
	# 	alpha_admix_pop_inits <- matrix(alpha_admix_pop_inits, nrow = nK, ncol = n_admix_sites)

	# 	# ADMIXTURE post posterior sampler for present day
	# 	admix_y_sq_pps_inits <- runifMatrix(nrow = n_counties, ncol = nK, stand = 'rows')

	# 	# ADMIXTURE post posterior sampler for future
	# 	admix_y_future_pps_inits <- runifMatrix(nrow = n_counties_future, ncol = nK, stand = 'rows')

	# 	# SDM response curves
	# 	response_curves_sdm_inits <- runifMatrix(nrow = n_response_curve_rows, ncol = n_predictors, min = 0, max = 20)
	# 	response_curves_sdm_inits <- round(response_curves_sdm_inits)

	# 	# ADMIXTURE response curves
	# 	nr <- n_response_curve_rows
	# 	nc <- nK
	# 	nd <- n_predictors
	# 	n <- nr * nc * nd
		
	# 	r <- round(runif(n, 0, 10))
	# 	response_curves_admix_inits <- array(r, dim = c(nr, nc, nd))

	# 	inits[[i]] <- list(
			
	# 		beta_occs = beta_occs_inits,
			
	# 		beta_admix = beta_admix_inits,
	# 		alpha0_occs_sampling = runif(1, -0.5, 0.5),
	# 		alpha_occs_area = runif(1, -0.5, 0.5),
	# 		alpha_occs_poaceae = runif(1, -0.5, 0.5),

	# 		alpha_admix = alpha_admix_pop_inits,
	# 		admix_y_sq_pps = admix_y_sq_pps_inits,
			
	# 		N = N_inits,

	# 		lambda_sq = lambda_sq_inits,

	# 		lambda_ssp245_2041_2070 = lambda_fut_inits,
	# 		lambda_ssp245_2071_2100 = lambda_fut_inits,
	# 		lambda_ssp370_2041_2070 = lambda_fut_inits,
	# 		lambda_ssp370_2071_2100 = lambda_fut_inits,

	# 		admix_y_ssp245_2041_2070_pps = admix_y_future_pps_inits,
	# 		admix_y_ssp245_2071_2100_pps = admix_y_future_pps_inits,

	# 		admix_y_ssp370_2041_2070_pps = admix_y_future_pps_inits,
	# 		admix_y_ssp370_2071_2100_pps = admix_y_future_pps_inits,

	# 		response_curves_sdm = response_curves_sdm_inits,
	# 		response_curves_admix = response_curves_admix_inits

	# 	)

	# }

	# N, lambdas, and betas for SDM
	N_inits <- 10 * ag_sq$any_ag_quality1to3
	lambda_sq_inits <- 1 + ag_sq$any_ag_quality1to3
	lambda_fut_inits <- 1 + ag_sq$any_ag_quality1to3
	beta_occs_inits <- rep(0, n_beta_terms)

	# betas for ADMIXTURE
	n <- nK * n_beta_terms
	beta_admix_inits <- rep(0, n)
	beta_admix_inits <- matrix(beta_admix_inits, nrow = nK, ncol = n_beta_terms)

	# alphas for ADMIXTURE
	n <- nK * n_admix_sites
	alpha_admix_pop_inits <- rep(0, n)
	alpha_admix_pop_inits <- matrix(alpha_admix_pop_inits, nrow = nK, ncol = n_admix_sites)

	# ADMIXTURE post posterior sampler for present day
	admix_y_sq_pps_inits <- runifMatrix(nrow = n_counties, ncol = nK, stand = 'rows')

	# ADMIXTURE post posterior sampler for future
	admix_y_future_pps_inits <- runifMatrix(nrow = n_counties_future, ncol = nK, stand = 'rows')

	# SDM response curves
	response_curves_sdm_inits <- runifMatrix(nrow = n_response_curve_rows, ncol = n_predictors, min = 0, max = 20)
	response_curves_sdm_inits <- round(response_curves_sdm_inits)

	# ADMIXTURE response curves
	nr <- n_response_curve_rows
	nc <- nK
	nd <- n_predictors
	n <- nr * nc * nd
	
	r <- round(runif(n, 0, 10))
	response_curves_admix_inits <- array(r, dim = c(nr, nc, nd))

	inits <- list(
		
		beta_occs = rep(0, n_beta_terms),
		
		beta_admix = beta_admix_inits,
		alpha0_occs_sampling = 0,
		alpha_occs_area = 0,
		alpha_occs_poaceae = 0,

		alpha_admix = alpha_admix_pop_inits,
		admix_y_sq_pps = admix_y_sq_pps_inits,
		
		N = N_inits,

		lambda_sq = lambda_sq_inits,

		lambda_ssp245_2041_2070 = lambda_fut_inits,
		lambda_ssp245_2071_2100 = lambda_fut_inits,
		lambda_ssp370_2041_2070 = lambda_fut_inits,
		lambda_ssp370_2071_2100 = lambda_fut_inits,

		admix_y_ssp245_2041_2070_pps = admix_y_future_pps_inits,
		admix_y_ssp245_2071_2100_pps = admix_y_future_pps_inits,

		admix_y_ssp370_2041_2070_pps = admix_y_future_pps_inits,
		admix_y_ssp370_2071_2100_pps = admix_y_future_pps_inits,

		response_curves_sdm = response_curves_sdm_inits,
		response_curves_admix = response_curves_admix_inits

	)

	say('Data:')
	print(str(data))

	say('Constants:', pre = 1)
	print(str(constants))

	say('Initializations:', pre = 1)
	print(str(inits))

	### define model
	say('nimbleCode():', level = 2)
	say('For the SDM, we assume an N-mixture model (latent, real abundance ~ Poisson, and observations ~ binomial draws from latent abundance.', post = 2)
	say('For the ADMIXTURE data, we assume a multinomial model (proportion of ancestry from a given population ~ dmulti(environment).', post = 2)
	say('Priors for the SDM and ADMIXTURE component are drawn from the same distributions.', post = 2)
	
	model_code <- nimbleCode({
	  
		# priors for SDM sampling bias
		alpha0_occs_sampling ~ dnorm(0, sd = 10)
		alpha_occs_area ~ dnorm(0, sd = 10)
		alpha_occs_poaceae ~ dnorm(0, sd = 10)

		# priors for SDM relationship to environment
		beta_occs[1] ~ dnorm(0, sd = 10) # intercept... not regularized
		for (j in 2:n_beta_terms) {
			# beta_occs[j] ~ ddexp(0, rate = 1) # regularization toward 0
			beta_occs[j] ~ dnorm(0, sd = 10) # broad prior
		}

		# priors for ADMIXTURE relationship to environment
		# betas_admix are a matrix with rows as populations and columns as betas
		for (k in 1:nK) {
		
			beta_admix[k, 1] ~ dnorm(0, sd = 10) # intercept
			for (j in 2:n_beta_terms) {
				# beta_admix[k, j] ~ dnorm(beta_occs[j], sd = 10) # broad prior
				beta_admix[k, j] ~ ddexp(beta_occs[j], rate = 1) # regularized prior bc small number of sites
			}
		
		}

		# likelihood of SDM
		for (i in 1:n_counties) {    # this specifies estimates be made for each county
			
			### actual abundance (latent--unobserved)
			N[i] ~ dpois(lambda_sq[i])

			# relationship between latent abundance and environment
			# rescale betas so they are the same between SDM and ADMIXTURE
			log(lambda_sq[i]) <- inprod(beta_occs[1:n_beta_terms], occs_x_sq[i, 1:n_beta_terms])

			### observed number of AG
			occs_y[i] ~ dbin(prob = occs_p[i], size = N[i])

			# sampling bias
			logit(occs_p[i]) <- alpha0_occs_sampling + alpha_occs_area * log_area_km2_scaled[i] + alpha_occs_poaceae * log_num_poaceae_records_scaled[i]

		}

		# prior for ADMIXTURE site-level effects
		for (k in 1:nK) {
			for (i in 1:n_admix_sites) {
				alpha_admix[k, i] ~ dnorm(0, sd = 10)
			}
		}

		# likelihood of ADMIXTURE data
		# admix_y rows are individuals, columns are ancestry proportions
		# admix_pop rows are individuals, columns are sites
		for (i in 1:n_admix_samples) {
		
			admix_y[i, 1:nK] ~ ddirch(admix_psi[i, 1:nK])
			for (k in 1:nK) {

				log(admix_psi[i, k]) <-
					inprod(beta_admix[k, 1:n_beta_terms], admix_x_sq[i, 1:n_beta_terms]) + # environment
					inprod(alpha_admix[k, 1:n_admix_sites], admix_pop[i, 1:n_admix_sites]) # site offset
				

			}
		
		}

		# posterior predictive sampler (PPS) for predicting ADMIXTURE proportions to present-day
		for (i in 1:n_counties) {    # this specifies estimates be made for each county
			
			# predict ADMIXTURE to present
			for (k in 1:nK) {
				log(admix_y_sq[i, k]) <-
					inprod(beta_admix[k, 1:n_beta_terms], occs_x_sq[i, 1:n_beta_terms])
			}
			admix_y_sq_pps[i, 1:nK] ~ ddirch(admix_y_sq[i, 1:nK])

		}

		# posterior predictive samplers (PPS) for future SDM and ADMIXTURE predictions
		for (i in 1:n_counties_future) {

			# SDM: AG occurrence ~ environment
			log(lambda_ssp245_2041_2070[i]) <- inprod(beta_occs[1:n_beta_terms], x_ssp245_2041_2070[i, 1:n_beta_terms])
			log(lambda_ssp245_2071_2100[i]) <- inprod(beta_occs[1:n_beta_terms], x_ssp245_2071_2100[i, 1:n_beta_terms])
			
			log(lambda_ssp370_2041_2070[i]) <- inprod(beta_occs[1:n_beta_terms], x_ssp370_2041_2070[i, 1:n_beta_terms])
			log(lambda_ssp370_2071_2100[i]) <- inprod(beta_occs[1:n_beta_terms], x_ssp370_2071_2100[i, 1:n_beta_terms])

			# ADMIXTURE: ancestry proportion ~ environment
			for (k in 1:nK) {
				
				log(admix_y_ssp245_2041_2070[i, k]) <-
					inprod(beta_admix[k, 1:n_beta_terms], x_ssp245_2041_2070[i, 1:n_beta_terms])
				log(admix_y_ssp245_2071_2100[i, k]) <-
					inprod(beta_admix[k, 1:n_beta_terms], x_ssp245_2071_2100[i, 1:n_beta_terms])

				log(admix_y_ssp370_2041_2070[i, k]) <-
					inprod(beta_admix[k, 1:n_beta_terms], x_ssp370_2041_2070[i, 1:n_beta_terms])
				log(admix_y_ssp370_2071_2100[i, k]) <-
					inprod(beta_admix[k, 1:n_beta_terms], x_ssp370_2071_2100[i, 1:n_beta_terms])


			}

			admix_y_ssp245_2041_2070_pps[i, 1:nK] ~ ddirch(admix_y_ssp245_2041_2070[i, 1:nK])
			admix_y_ssp245_2071_2100_pps[i, 1:nK] ~ ddirch(admix_y_ssp245_2071_2100[i, 1:nK])

			admix_y_ssp370_2041_2070_pps[i, 1:nK] ~ ddirch(admix_y_ssp370_2041_2070[i, 1:nK])
			admix_y_ssp370_2071_2100_pps[i, 1:nK] ~ ddirch(admix_y_ssp370_2071_2100[i, 1:nK])

		}

		# SDM: posterior predictive sampler for response curves
		for (j in 1:n_predictors) {
			for (i in 1:n_response_curve_rows) {
				
				log(response_curves_sdm_lambda[i, j]) <-
					inprod(beta_occs[1:n_beta_terms], response_curve_x[i, 1:n_beta_terms, j])

				response_curves_sdm[i, j] ~ dpois(response_curves_sdm_lambda[i, j])

			}
		}

		# ADMIXTURE: posterior predictive sampler for response curves
		# response_curves_admix_psi and admix_response_curves
		# 	rows: 1:number of values of x we want to estimate response at {1, 2, 3, ...}
		#	columns: ancestral populations {1, 2, 3, ..., nK}
		#	depths: predictors {1, 2, 3, ..., number of predictors} (note: number of *predictors*, not terms in models)
		for (i in 1:n_response_curve_rows) {
			for (j in 1:n_predictors) {
				for (k in 1:nK) {

					log(response_curves_admix_psi[i, k, j]) <-
						inprod(beta_admix[k, 1:n_beta_terms], response_curve_x[i, 1:n_beta_terms, j])
		
				}

				response_curves_admix[i, 1:nK, j] ~ ddirch(response_curves_admix_psi[i, 1:nK, j])

			}
		}
	  
	})

	print(model_code)

	say('nimbleModel():', level = 2)
	model <- nimbleModel(
		code = model_code, # our model
		constants = constants, # constants
		data = data, # data
		inits = inits, # initialization values
		check = TRUE, # any errors?
		calculate = FALSE,
		buildDerivs = TRUE # need for Hamiltonian Monte Carlo
	)

	say('$initializeInfo() and $calculate():', level = 2)
	model$initializeInfo()
	model$calculate()

	say('configureMCMC():', level = 2)
	monitors <- c(
		'beta_occs', 'alpha0_occs_sampling', 'alpha_occs_area', 'alpha_occs_poaceae', 'lambda_sq',
		'lambda_ssp245_2041_2070', 'lambda_ssp245_2071_2100', 'lambda_ssp370_2041_2070', 'lambda_ssp370_2071_2100',
		'response_curves_sdm',
		'beta_admix', 'alpha_admix', 'admix_y_sq_pps',
		'admix_y_ssp245_2041_2070_pps', 'admix_y_ssp245_2071_2100_pps', 'admix_y_ssp370_2041_2070_pps', 'admix_y_ssp370_2071_2100_pps',
		'response_curves_admix'
	)
	
	conf <- configureMCMC(
		model,
		monitors = monitors,
		print = TRUE,
		enableWAIC = FALSE
	)

	# add no U-turn sampler (Hamiltonian Monte Carlo)
	target <- monitors
	conf$addSampler(target = monitors, type = 'NUTS')

	### compile/build/run model/save MCMC
	build <- buildMCMC(conf)

	compiled <- compileNimble(model, build, showCompilerOutput = FALSE)

	chains <- runMCMC(
		compiled$build,
		niter = niter,
		nburnin = nburnin,
		thin = thin,
		nchains = nchains,
		inits = inits,
		progressBar = TRUE,
		samplesAsCodaMCMC = TRUE,
		summary = TRUE,
		WAIC = FALSE,
		perChainWAIC = FALSE
	)

	saveRDS(chains, paste0(out_dir, '/chains.rds'))

	say('session info', level = 2)
	print(sessionInfo())

	say(date(), pre = 1)
	sink()

say('###########################################')
say('### model diagnostics for SDM component ###')
say('###########################################')

	# chains <- readRDS(paste0(out_dir, '/chains.rds'))
	mcmc <- chains$samples

	cols <- c(paste0('beta_occs[', 1:n_beta_terms, ']'), 'alpha0_occs_sampling', 'alpha_occs_poaceae', 'alpha_occs_area')
	for (i in 1:nchains) mcmc[[i]] <- mcmc[[i]][ , cols]

	# graphing trace plots for all betas
	pars <- paste0('beta_occs[', 1:n_beta_terms, ']')
	file <- paste0(out_dir, '/diagnostics/admixture_occs_beta_trace.png')
	ggsave(mcmc_trace(mcmc, pars = pars), file = file, width = 10, height = 8, dpi = 450, bg = 'white')

	# graphing trace plots for all "extra" betas
	pars <- c('alpha0_occs_sampling', 'alpha_occs_poaceae', 'alpha_occs_area')
	file <- paste0(out_dir, '/diagnostics/admixture_occs_alpha_trace.png')
	ggsave(mcmc_trace(mcmc, pars = pars), file = file, width = 10, height = 4, dpi = 450, bg = 'white')

	# graphing density plots for all betas
	pars <- paste0('beta_occs[', 1:n_beta_terms, ']')
	file <- paste0(out_dir, '/diagnostics/admixture_occs_beta_density.png')
	ggsave(mcmc_dens_overlay(mcmc, pars = pars), file = file, width = 12, height = 8, dpi = 450, bg = 'white')

	# graphing trace plots for all alphas
	pars <- c('alpha0_occs_sampling', 'alpha_occs_poaceae', 'alpha_occs_area')
	file <- paste0(out_dir, '/diagnostics/admixture_occs_alpha_density.png')
	ggsave(mcmc_dens_overlay(mcmc, pars = pars), file = file, width = 10, height = 4, dpi = 450, bg = 'white')

	# Gelman-Rubin statistics for convergence
	rhats <- gelman.diag(mcmc, autoburnin = FALSE, multivariate = TRUE)

	sink(paste0(out_dir, '/diagnostics/convergence_sdm_gelman_rubin.txt'), split = TRUE)
	say('GELMAN-RUBIN STATISTICS')
	say(date(), post = 2)
	print(rhats)
	sink()

	# Geweke statistics for convergence
	gew <- geweke.diag(mcmc)
	for (i in 1:nchains) gew[[i]] <- gew[[i]]$z
	gew <- do.call(cbind, gew)
	write.csv(gew, paste0(out_dir, '/diagnostics/convergence_sdm_geweke.csv'), row.names = TRUE)

	# effective sample size
	ess <- effectiveSize(mcmc)
	ess <- cbind(ess)
	write.csv(ess, paste0(out_dir, '/diagnostics/effective_sample_sizes_sdm.csv'), row.names = TRUE)

say('#################################################')
say('### model diagnostics for ADMIXTURE component ###')
say('#################################################')

	# chains <- readRDS(paste0(out_dir, '/chains.rds'))
	mcmc <- chains$samples

	expand <- expand.grid(k = K, betas = 1:n_beta_terms)
	expand <- apply(expand, 1, paste, collapse = ', ')
	cols_1 <- paste0('beta_admix[', expand, ']')

	expand <- expand.grid(k = K, n_admix_sites = 1:n_admix_sites)
	expand <- apply(expand, 1, paste, collapse = ', ')
	cols_2 <- paste0('alpha_admix[', expand, ']')

	cols <- c(cols_1, cols_2)

	for (i in 1:nchains) {
		mcmc[[i]] <- mcmc[[i]][ , cols]
	}

	# graphing trace and density plots for all betas
	for (i in 1:n_beta_terms) {

		pars <- paste0('beta_admix[', 1:nK, ', ', i, ']')
		file <- paste0(out_dir, '/diagnostics/admixture_admix_beta_trace_beta_', i, '.png')
		ggsave(mcmc_trace(mcmc, pars = pars), file = file, width = 12, height = 8, dpi = 450, bg = 'white')

		file <- paste0(out_dir, '/diagnostics/admixture_admix_beta_density_beta_', i, '.png')
		ggsave(mcmc_dens_overlay(mcmc, pars = pars), file = file, width = 12, height = 8, dpi = 450, bg = 'white')

	}

	# graphing trace plots for all alphas
	for (i in 1:n_admix_sites) {

		pars <- paste0('alpha_admix[', 1:nK, ', ', i, ']')

		file <- paste0(out_dir, paste0('/diagnostics/admixture_admix_alpha_trace_alpha_', i, '.png'))
		ggsave(mcmc_trace(mcmc, pars = pars), file = file, width = 12, height = 8, dpi = 450, bg = 'white')

		file <- paste0(out_dir, paste0('/diagnostics/admixture_admix_alpha_density_alpha_', i, '.png'))
		ggsave(mcmc_dens_overlay(mcmc, pars = pars), file = file, width = 12, height = 8, dpi = 450, bg = 'white')

	}

	# Gelman-Rubin statistics for convergence
	rhats <- gelman.diag(mcmc, autoburnin = FALSE, multivariate = TRUE)

	sink(paste0(out_dir, '/diagnostics/convergence_admix_gelman_rubin.txt'), split = TRUE)
	say('GELMAN-RUBIN STATISTICS')
	say(date(), post = 2)
	print(rhats)
	sink()

	# Geweke statistics for convergence
	gew <- geweke.diag(mcmc)
	for (i in 1:nchains) gew[[i]] <- gew[[i]]$z
	gew <- do.call(cbind, gew)
	write.csv(gew, paste0(out_dir, '/diagnostics/convergence_admix_geweke.csv'), row.names = TRUE)

	# effective sample size
	ess <- effectiveSize(mcmc)
	ess <- cbind(ess)
	write.csv(ess, paste0(out_dir, '/diagnostics/effective_sample_sizes_admix.csv'), row.names = TRUE)

say('#######################', pre = 1)
say('### response curves ###')
say('#######################')

	# chains <- readRDS(paste0(out_dir, '/chains.rds'))
	mcmc <- chains$samples

	# 1st is species
	colors <- c('black', 'red', 'green', 'cyan', 'orange', 'yellow', 'violet', 'floralwhite')

	# make a plot for how AG SDM and each ancestral population responds to each predictor
	responses <- list()
	for (pred in 1:n_predictors) {

		# unscaled predictor value
		x <- env_array[ , pred, pred]

		# create data frame with SDM response
		pars <- paste0('response_curves_sdm[', 1:n_response_curve_rows, ', ', pred, ']')
		sdm_response_mean <- chains$summary$all.chains[pars, 'Mean']
		sdm_response_lower <- chains$summary$all.chains[pars, '95%CI_low']
		sdm_response_upper <- chains$summary$all.chains[pars, '95%CI_upp']

		# data frames to hold predictions in long format
		df_mean <- data.frame(
			x = x,
			response = sdm_response_mean,
			type = 'Species-level'
		)

		df_lower <- data.frame(
			x = x,
			response = sdm_response_lower,
			type = 'Species-level'
		)

		df_upper <- data.frame(
			x = x,
			response = sdm_response_upper,
			type = 'Species-level'
		)

		# add to data frames ADMIXTURE responses
		admix_response_mean <- admix_response_lower <- admix_response_upper <- list()
		for (k in 1:nK) {

			# response of this ancestral population to this predictor
			pars <- paste0('response_curves_admix[', 1:n_response_curve_rows, ', ', k, ', ', pred, ']')
			admix_response_mean <- chains$summary$all.chains[pars, 'Mean']
			admix_response_lower <- chains$summary$all.chains[pars, '95%CI_low']
			admix_response_upper <- chains$summary$all.chains[pars, '95%CI_upp']

			# data frames to hold predictions in long format
			df_mean <- rbind(
				df_mean,
				data.frame(
					x = x,
					response = admix_response_mean,
					type = paste0('K', k)
				)
			)

			df_lower <- rbind(
				df_lower,
				data.frame(
					x = x,
					response = admix_response_lower,
					type = paste0('K', k)
				)
			)

			df_upper <- rbind(
				df_upper,
				data.frame(
					x = x,
					response = admix_response_upper,
					type = paste0('K', k)
				)
			)

		}


		responses[[pred]] <- ggplot()
		types <- unique(df_mean$type)

		# plot CI first
		for (i in seq_along(types)) {

			type <- types[i]

			this_df_upper <- df_upper[df_mean$type == type, ]
			this_df_lower <- df_lower[df_mean$type == type, ]

			max_score <- max(this_df_upper$response, na.rm = TRUE)

			this_df_upper$response <- this_df_upper$response / max_score
			this_df_lower$response <- this_df_lower$response / max_score

			this_df_ci <- data.frame(
				x = c(x, rev(x)),
				y = c(this_df_upper$response, rev(this_df_lower$response)),
				type = type
			)

			responses[[pred]] <- responses[[pred]] +
				geom_polygon(
					data = this_df_ci,
					mapping = aes(x = x, y = y),
					color = NA,
					fill = alpha(colors[i], 0.1)
				)

		}

		# plot trendlines second
		for (i in seq_along(types)) {

			type <- types[i]

			this_df_mean <- df_mean[df_mean$type == type, ]
			this_df_upper <- df_upper[df_mean$type == type, ]

			max_score <- max(this_df_upper$response, na.rm = TRUE)

			this_df_mean$response <- this_df_mean$response / max_score

			responses[[pred]] <- responses[[pred]] +
				geom_line(
					data = this_df_mean,
					mapping = aes(x = x, y = response),
					color = colors[i]
				)

		}

		if (predictor_names[pred] == 'aridity') {
			nice_title <- 'Aridity (temperature / precipitation)'
			nice_axis <- 'Aridity (°C / mm)'
		} else if (predictor_names[pred] == 'bio7') {
			nice_title <- 'Temperature annual range (BIO 07)'
			nice_axis <- 'Temperature annual range (°C)'
		} else if (predictor_names[pred] == 'bio12') {
			nice_title <- 'Total annual precipitation (BIO 12)'
			nice_axis <- 'Total annual precipitation (mm)'
		} else if (predictor_names[pred] == 'bio15') {
			nice_title <- 'Precipitation seasonality (BIO 15)'
			nice_axis <- 'Precipitation seasonality'
		} else if (predictor_names[pred] == 'bio15') {
			nice_title <- 'Precipitation seasonality (BIO 15)'
			nice_axis <- 'Precipitation seasonality'
		} else if (predictor_names[pred] == 'ph') {
			nice_title <- 'Soil pH'
			nice_axis <- 'pH'
		} else if (predictor_names[pred] == 'sand') {
			nice_title <- 'Soil proportion sand'
			nice_axis <- 'Proportion sand'
		} else if (predictor_names[pred] == 'silt') {
			nice_title <- 'Soil proportion silt'
			nice_axis <- 'Proportion silt'
		} else if (predictor_names[pred] == 'soc') {
			nice_title <- 'Soil organic matter'
			nice_axis <- 'Soil organic matter'
		}

		responses[[pred]] <- responses[[pred]] +
			ylim(0, 1) +
			scale_color_discrete(name = 'Entity') +
			scale_fill_discrete(name = 'Entity') +
			xlab(nice_axis) +
			ylab('Scaled ancestry proportion or\nSDM-predicted abundance') +
			ggtitle(nice_title) +
			theme(
				plot.title = element_text(size = 10)
			)

	} # next predictor

	if (length(predictor_names) <= 4) {
		nrow <- 1
		width <- 14
		height <- 4
	} else {
		nrow <- 2
		width <- 16
		height <- 10	
	}

	responses <- plot_grid(plotlist = responses, nrow = nrow)

	ggsave(plot = responses, filename = paste0(out_dir, '/response_curves.png'), width = width, height = height, dpi = 600)

say('################################################')
say('### instill predictions into spatial vectors ###')
say('################################################')

	# chains <- readRDS(paste0(out_dir, '/chains.rds'))

	# subset chain summary to just the lambdas
	summary <- chains$summary$all.chains

	### present-day SDM predictions
	###############################

	# just counties with data
	which_lambda <- grepl(rownames(summary), pattern = 'lambda_sq')
	lambda <- summary[which_lambda, ]

	lambda_mean <- lambda[ , 'Mean']
	lambda_0.05ci <- lambda[ , '95%CI_low']
	lambda_0.95ci <- lambda[ , '95%CI_upp']
	lambda_90ci <- lambda_0.95ci - lambda_0.05ci

	ag_vect_sq$lambda_mean <- lambda_mean
	ag_vect_sq$lambda_90ci <- lambda_90ci

	### present-day ADMIXTURE predictions
	#####################################

	### make columns with predicted ancestry proportions
	for (k in 1:nK) {

		# just counties with data
		pars <- paste0('admix_y_sq_pps[', 1:n_counties, ', ', k, ']')
		matches <- which(rownames(summary) %in% pars)
		x <- summary[matches, ]

		ag_vect_sq$MEANDUMMY <- x[ , 'Mean']
		lower <- x[ , '95%CI_low']
		upper <- x[ , '95%CI_upp']
		ag_vect_sq$CIDUMMY <- upper - lower

		names(ag_vect_sq)[(ncol(ag_vect_sq) - 1):ncol(ag_vect_sq)] <-
			c(paste0('k', k, '_admix_mean'), paste0('k', k, '_admix_90ci'))

	}

	writeVector(ag_vect_sq, paste0(out_dir, '/sdm_admixture_1961_2020.gpkg'), overwrite = TRUE)

	### future SDM and ADMIXTURE predictions
	########################################

	futs <- c(
		'ssp245_2041_2070',
		'ssp245_2071_2100',
		'ssp370_2041_2070',
		'ssp370_2071_2100'
	)

	# add predictions to vectors
	ag_vect_futs <- list()
	for (fut in futs) {

		ag_vect_fut <- get(paste0('ag_vect_', fut))

		### SDM

		# just counties with data
		pars <- paste0('lambda_', fut, '[', 1:n_counties, ']')
		matches <- which(rownames(summary) %in% pars)
		x <- summary[matches, ]

		means <- x[ , 'Mean']
		ag_vect_fut$MEANDUMMY <- means

		lower <- x[ , '95%CI_low']
		upper <- x[ , '95%CI_upp']
		ci <- upper - lower

		ag_vect_fut$CIDUMMY <- ci

		names(ag_vect_fut)[(ncol(ag_vect_fut) - 1):ncol(ag_vect_fut)] <-
			c('lambda_mean', 'lambda_90ci')

		### ADMIXTURE
		for (k in 1:nK) {

			# just counties with data
			pars <- paste0('admix_y_', fut, '_pps[', 1:n_counties, ', ', k, ']')
			matches <- which(rownames(summary) %in% pars)
			x <- summary[matches, ]

			ag_vect_fut$MEANDUMMY <- x[ , 'Mean']
			lower <- x[ , '95%CI_low']
			upper <- x[ , '95%CI_upp']
			ag_vect_fut$CIDUMMY <- upper - lower

			names(ag_vect_fut)[(ncol(ag_vect_fut) - 1):ncol(ag_vect_fut)] <-
				c(paste0('k', k, '_admix_mean'), paste0('k', k, '_admix_90ci'))

		}

		writeVector(ag_vect_fut, paste0(out_dir, '/sdm_admixture_', fut, '.gpkg'), overwrite = TRUE)
		ag_vect_futs[[length(ag_vect_futs) + 1]] <- ag_vect_fut
		names(ag_vect_futs)[length(ag_vect_futs)] <- fut

	}

say('###########################################', pre = 1)
say('### map of CURRENT species distribution ###')
say('###########################################')

	# North America
	nam <- vect(paste0(drive, '/Research Data/GADM/Version 4.1/High Res North America Level 1 sans Great Lakes SpatVector WGS84.gpkg'))
	nam <- project(nam, ag_vect_sq)

	### SDM map
	###########

	# prediction values for color scaling
	lambda_mean <- ag_vect_sq$lambda_mean 
	lambda_90ci <- ag_vect_sq$lambda_90ci

	lambda_sq_quants <- quantile(lambda_mean, c(0.25, 0.5, 0.75, 0.95))
	quant_labels <- c('[0 - 0.25)', '[0.25 - 0.5)', '[0.5 - 0.75)', '[0.75 - 0.95)', '[0.95-1]')
	
	lambda_quant_col <- NA
	lambda_quant_col[lambda_mean < lambda_sq_quants[1]] <- quant_labels[1]
	lambda_quant_col[lambda_mean >= lambda_sq_quants[1] & lambda_mean < lambda_sq_quants[2]] <- quant_labels[2]
	lambda_quant_col[lambda_mean >= lambda_sq_quants[2] & lambda_mean < lambda_sq_quants[3]] <- quant_labels[3]
	lambda_quant_col[lambda_mean >= lambda_sq_quants[3] & lambda_mean < lambda_sq_quants[4]] <- quant_labels[4]
	lambda_quant_col[lambda_mean >= lambda_sq_quants[4]] <- quant_labels[5]

	# set plot extent
	ag_vect_sq_pres <- ag_vect_sq[ag_vect_sq$any_ag_quality1to3 > 0]
	extent <- ext(ag_vect_sq_pres)
	extent <- as.vector(extent)
	x_range <- (extent[2] - extent[1])
	y_range <- (extent[4] - extent[3])
	extent[1] <- extent[1] + 0.01 * x_range
	extent[3] <- extent[3] + 0.125 * y_range
	extent[4] <- extent[4] - 0.2 * y_range

	# points for counties with at least one AG record
	cents_with_ag <- ag_vect_sq[ag_vect_sq$any_ag_quality1to3 > 0]
	cents_with_ag <- centroids(cents_with_ag)

	# color scale for SDM predictions
	fill_scale <- c('[0 - 0.25)' = 'gray85', '[0.25 - 0.5)' = alpha('forestgreen', 0.25), '[0.5 - 0.75)' = alpha('forestgreen', 0.45), '[0.75 - 0.95)' = alpha('forestgreen', 0.7), '[0.95-1]' = 'forestgreen')

	# location of plots sampled for ADMIXTURE
	admix_vect <- vect(admix, geom = c('longitude', 'latitude'), crs = getCRS('WGS84'))
	admix_vect <- project(admix_vect, ag_vect_sq)

	sdm_sq_map <- ggplot() +
		layer_spatial(ag_vect_sq, aes(fill = lambda_quant_col), color = NA) +
		layer_spatial(nam, color = 'gray30', fill = NA, linewidth = 0.1) +
		scale_fill_manual(
			name = 'Quantile\n of λ',
			values = fill_scale
		) +
		layer_spatial(cents_with_ag, pch = 16, color = alpha('gray20', 0.5), size = 0.35) +
		layer_spatial(admix_vect, fill = NA, color = 'orange', pch = 1, size = 6) +
		guides(fill = guide_legend(reverse = TRUE)) +
		xlim(extent[1], extent[2]) + ylim(extent[3], extent[4]) +
		ggtitle(expression('Present-day distribution of ' * italic('Andropogon gerardi')), subtitle = '1961-2020 | N-mixture model') +
		theme(
			plot.title = element_text(size = 16),
			plot.subtitle = element_text(size = 14)
		)		

	ggsave(plot = sdm_sq_map, filename = paste0(out_dir, '/map_sdm_1961_2020.png'), width = 12, height = 10, dpi = 600)

say('###########################################################', pre = 1)
say('### maps of CURRENT species and ADMIXTURE distributions ###')
say('###########################################################')

	# plot ADMIXTURE predictions only in counties that have predicted AG abundance >= max TSS threshold OR are within this many km of such counties
	focal_region_buffer_m <- 200000

	# use same extent as above
	summary <- chains$summary$all.chains

	admix_agg <- aggregate(admix, by = list(admix$population_proper), 'mean')
	admix_agg$population_proper <- NULL
	names(admix_agg)[1] <- 'site'

	### defining plot region:
	# 1. select counties with recorded AG
	# 2. add buffer
	# 3. add buffer to ADMIXTURE sample sites
	# 4. merge buffers
	# 5. select whole counties within this merged buffer
	# 6. remove counties above a given elevation (montane counties)
	# 7. select set of counties that form the largest contiguous block of counties

	# buffers around ADMIXTURE samples
	admix_vect <- vect(admix_agg, geom = c('longitude', 'latitude'), crs = getCRS('WGS84'))
	admix_vect <- project(admix_vect, ag_vect_sq)
	admix_vect_buffer <- buffer(admix_vect, focal_region_buffer_m)
	admix_vect_buffer <- aggregate(admix_vect_buffer)

	# buffers around AG records
	ag_vect_sq$id <- 1:nrow(ag_vect_sq)
	ag_vect_sq$area_m2_ALL_AG <- expanse(ag_vect_sq)
	ag_vect_focal <- ag_vect_sq[ag_vect_sq$any_ag_quality1to3 > 0]
	ag_vect_buffer <- buffer(ag_vect_focal, focal_region_buffer_m)
	ag_vect_buffer <- aggregate(ag_vect_buffer)
	
	# union buffers and select
	buffer <- union(ag_vect_buffer, admix_vect_buffer)
	buffer <- aggregate(buffer)
	ag_vect_focal <- crop(ag_vect_sq, buffer)
	ag_vect_focal$area_m2_FOCAL <- expanse(ag_vect_focal)
	ag_vect_focal <- ag_vect_focal[ag_vect_focal$area_m2_ALL_AG == ag_vect_focal$area_m2_FOCAL]

	ag_vect_sq$id <- NULL
	ag_vect_sq$area_m2_ALL_AG <- NULL

	# mask out counties with elevations > some value
	max_elevation_m <- quantile(ag_vect_focal$elevation_m[ag_vect_focal$any_ag_quality1to3 > 0], 0.99)
	ag_vect_focal <- ag_vect_focal[ag_vect_focal$elevation_m <= max_elevation_m]

	ag_vect_focal_agg <- aggregate(ag_vect_focal)
	ag_vect_focal_agg <- disagg(ag_vect_focal_agg)
	agg_areas <- expanse(ag_vect_focal_agg)
	ag_vect_focal_agg <- ag_vect_focal_agg[which.max(agg_areas)]

	ag_vect_focal <- crop(ag_vect_focal, ext(ag_vect_focal_agg))

	fill_color <- c('red', 'green', 'cyan', 'orange', 'yellow', 'violet', 'floralwhite')
	admix_sq_maps <- list()
	for (k in 1:nK) {

		# color each county by (scaled) predicted ancestry proportions
		fill <- unlist(as.data.frame(ag_vect_focal[ , paste0('k', k, '_admix_mean')]))
		ag_vect_focal$fill <- fill

		admix_sq_maps[[k]] <- ggplot() +
			layer_spatial(nam, color = NA, fill = 'gray80') +
			layer_spatial(ag_vect_focal, aes(fill = fill), color = NA) +
			layer_spatial(nam, color = 'gray30', fill = NA, linewidth = 0.1) +
			guides(fill = guide_legend(reverse = TRUE)) +
			scale_fill_gradientn(
				name = paste0('K', k),
				# limits = c(0, 0.8),
				# colors = c('gray20', fill_color[k])
				colors = c('gray30', fill_color[k]),
				label_number(accuracy = 0.01)
			) +
			xlim(extent[1], extent[2]) + ylim(extent[3], extent[4]) +
			ggtitle(
				paste0('K', k),
				subtitle = '1961-2020'
			) +
			theme(
				plot.title = element_text(size = 16),
				plot.subtitle = element_text(size = 14),
				legend.background = element_rect(color = 'black', fill = 'white'),
				legend.position = c(0.92, 0.15)
			)

		# add pie charts for observed ancestry proportions
		pies <- admix_agg[ , c('site', 'longitude', 'latitude', paste0('k', 1:nK))]
		pies <- vect(pies, geom = c('longitude', 'latitude'), crs = getCRS('WGS84'))
		pies <- project(pies, ag_vect_sq)
		pie_coords <- crds(pies)
		pies <- cbind(pies, pie_coords)
		pies <- as.data.frame(pies)

		# pie_colors <- rep(alpha('white', 0), K)
		pie_colors <- rep('gray20', nK)
		pie_colors[k] <- fill_color[k]
		names(pie_colors) <- paste0('k', K)

		admix_sq_maps[[k]] <- admix_sq_maps[[k]] +
			new_scale_fill() +
			geom_scatterpie(
				data = pies,
				mapping = aes(
					x = x, y = y,
					group = site
				),
				cols = paste0('k', K),
				pie_scale = 1
			) +
			scale_fill_manual(
				guide = 'none',
				values = pie_colors,
				label_number(accuracy = 0.01)
			)

	} # next k

	sdm_sq_map_leg_adjust <- sdm_sq_map
	sdm_sq_map_leg_adjust <- sdm_sq_map_leg_adjust +
		ggtitle('SDM') +
		theme(
			legend.background = element_rect(color = 'black', fill = 'white'),
			legend.position = c(0.92, 0.15)
		)

	combo_map_list <- c(list(sdm_sq_map_leg_adjust), admix_sq_maps)
	combined_sq_maps <- plot_grid(plotlist = combo_map_list, nrow = 2, ncol = 4)

	ggsave(plot = combined_sq_maps, filename = paste0(out_dir, '/map_sdm_admixture_1961_2020.png'), width = 19.2, height = 10.8, dpi = 600, bg = 'white')

say('##########################################################')
say('### maps of FUTURE species and ADMIXTURE distributions ###')
say('##########################################################')

	# use same extent as above
	summary <- chains$summary$all.chains

	fill_color <- c('red', 'green', 'cyan', 'orange', 'yellow', 'violet', 'floralwhite')

	futs <- c(
		'ssp245_2041_2070',
		'ssp245_2071_2100',
		'ssp370_2041_2070',
		'ssp370_2071_2100'
	)

	focal_geog <- as.data.frame(ag_vect_focal)[ , c('country', 'stateProvince', 'county')]
	focal_geog <- apply(focal_geog, 1, paste, collapse = ' ')
	for (fut in futs) {

		say(fut)

		ag_vect_fut <- ag_vect_futs[[fut]]

		# crop to area we want to map
		fut_geog <- as.data.frame(ag_vect_fut)[ , c('country', 'stateProvince', 'county')]
		fut_geog <- apply(fut_geog, 1, paste, collapse = ' ')
		keeps <- fut_geog %in% focal_geog
		ag_vect_fut <- ag_vect_fut[keeps]

		ssp <- substr(fut, 4, 6)
		start_year <- substr(fut, 8, 11)
		end_year <- substr(fut, 13, 16)

		subtitle <- paste0(start_year, '-', end_year, ' | SSP ', ssp)

		### SDM map
		###########

		# prediction values for color scaling
		lambda_mean <- ag_vect_fut$lambda_mean
		lambda_90ci <- ag_vect_fut$lambda_90ci

		# using quantile breaks and labels from creation of SQ SDM map above
		lambda_quant_col <- NA
		lambda_quant_col[lambda_mean < lambda_sq_quants[1]] <- quant_labels[1]
		lambda_quant_col[lambda_mean >= lambda_sq_quants[1] & lambda_mean < lambda_sq_quants[2]] <- quant_labels[2]
		lambda_quant_col[lambda_mean >= lambda_sq_quants[2] & lambda_mean < lambda_sq_quants[3]] <- quant_labels[3]
		lambda_quant_col[lambda_mean >= lambda_sq_quants[3] & lambda_mean < lambda_sq_quants[4]] <- quant_labels[4]
		lambda_quant_col[lambda_mean >= lambda_sq_quants[4]] <- quant_labels[5]
		ag_vect_fut$lambda_quant_col <- lambda_quant_col

		# location of plots sampled for ADMIXTURE
		admix_vect <- vect(admix, geom = c('longitude', 'latitude'), crs = getCRS('WGS84'))
		admix_vect <- project(admix_vect, ag_vect_fut)

		sdm_fut_map <- ggplot() +
			layer_spatial(ag_vect_fut, aes(fill = lambda_quant_col), color = NA) +
			layer_spatial(nam, color = 'gray30', fill = NA, linewidth = 0.1) +
			scale_fill_manual(
				name = 'Quantile\n of λ',
				values = fill_scale
			) +
			layer_spatial(cents_with_ag, pch = 16, color = alpha('gray20', 0.5), size = 0.35) +
			layer_spatial(admix_vect, fill = NA, color = 'orange', pch = 1, size = 6) +
			guides(fill = guide_legend(reverse = TRUE)) +
			xlim(extent[1], extent[2]) + ylim(extent[3], extent[4]) +
			ggtitle('SDM', subtitle = subtitle) +
			theme(
				plot.title = element_text(size = 16),
				plot.subtitle = element_text(size = 14),
				legend.background = element_rect(color = 'black', fill = 'white'),
				legend.position = c(0.92, 0.15)
			)

		### ADMIXTURE MAPS

		admix_fut_maps <- list()
		for (k in 1:nK) {

			# color each county by (scaled) predicted ancestry proportions
			fill <- unlist(as.data.frame(ag_vect_fut[ , paste0('k', k, '_admix_mean')]))
			ag_vect_fut$fill <- fill

			admix_fut_maps[[k]] <- ggplot() +
				layer_spatial(nam, color = NA, fill = 'gray80') +
				layer_spatial(ag_vect_fut, aes(fill = fill), color = NA) +
				layer_spatial(nam, color = 'gray30', fill = NA, linewidth = 0.1) +
				guides(fill = guide_legend(reverse = TRUE)) +
				scale_fill_gradientn(
					name = paste0('K', k),
					# limits = c(0, 0.8),
					# colors = c('gray20', fill_color[k])
					colors = c('gray30', fill_color[k])
				) +
				xlim(extent[1], extent[2]) + ylim(extent[3], extent[4]) +
				ggtitle(
					paste0('K', k),
					subtitle = subtitle
				) +
				theme(
					plot.title = element_text(size = 16),
					plot.subtitle = element_text(size = 14),
					legend.background = element_rect(color = 'black', fill = 'white'),
					legend.position = c(0.92, 0.15)
				)

			# add pie charts for observed ancestry proportions
			pies <- admix_agg[ , c('site', 'longitude', 'latitude', paste0('k', 1:nK))]
			pies <- vect(pies, geom = c('longitude', 'latitude'), crs = getCRS('WGS84'))
			pies <- project(pies, ag_vect_fut)
			pie_coords <- crds(pies)
			pies <- cbind(pies, pie_coords)
			pies <- as.data.frame(pies)

			# pie_colors <- rep(alpha('white', 0), K)
			pie_colors <- rep('gray20', nK)
			pie_colors[k] <- fill_color[k]
			names(pie_colors) <- paste0('k', 1:nK)

			admix_fut_maps[[k]] <- admix_fut_maps[[k]] +
				new_scale_fill() +
				geom_scatterpie(
					data = pies,
					mapping = aes(
						x = x, y = y,
						group = site
					),
					cols = paste0('k', 1:nK),
					pie_scale = 1
				) +
				scale_fill_manual(
					guide = 'none',
					values = pie_colors,
					label_number(accuracy = 0.01)
				)

		} # next k

		combo_map_list <- c(list(sdm_fut_map), admix_fut_maps)
		combined_fut_maps <- plot_grid(plotlist = combo_map_list, nrow = 2, ncol = 4)

		ggsave(plot = combined_fut_maps, filename = paste0(out_dir, '/map_sdm_admixture_', fut, '.png'), width = 19.2, height = 10.8, dpi = 600, bg = 'white')

	} # next future

say('###################################')
say('### parameter estimates for SDM ###')
say('###################################')

	# chains <- readRDS(paste0(out_dir, '/chains.rds'))

	# subset chain summary to just the lambdas
	summary <- chains$summary$all.chains

	model_term <- attr(terms(form), 'term.labels')
	model_term_nice <- model_term
	model_term_nice <- gsub(model_term_nice, pattern = 'I\\(', replacement = '')
	model_term_nice <- gsub(model_term_nice, pattern = '\\)', replacement = '')
	model_term_nice <- c('Intercept', model_term_nice)

	pars <- paste0('beta_occs[', 1:(1 + length(model_term)), ']')

	coeff_chains <- chains$samples
	for (i in 1:nchains) {
		coeff_chains[[i]] <- coeff_chains[[i]][ , pars]
		colnames(coeff_chains[[i]]) <- model_term_nice
	}

	caterpillars <- mcmc_intervals(coeff_chains, pars = model_term_nice)

	# plot posterior of coefficient estimates
	caterpillars <- caterpillars +
		xlab('Estimated value') +
		theme(plot.background = element_rect(fill = 'white'))

	ggsave(caterpillars, filename = paste0(out_dir, '/admixture_occs_beta_coefficients.png'), width = 10, height = 12, dpi = 300)

say('#########################################')
say('### parameter estimates for ADMIXTURE ###')
say('#########################################')

	# chains <- readRDS(paste0(out_dir, '/chains.rds'))

	model_term <- attr(terms(form), 'term.labels')
	model_terms_nice <- model_term
	model_terms_nice <- gsub(model_terms_nice, pattern = 'I\\(', replacement = '')
	model_terms_nice <- gsub(model_terms_nice, pattern = '\\)', replacement = '')
	model_terms_nice <- c('Intercept', model_terms_nice)

	coeffs <- list()
	for (i in seq_along(model_terms_nice)) {

		model_term_nice <- model_terms_nice[i]

		pars <- paste0('beta_admix[', 1:nK, ', ', i, ']')

		coeff_chains <- chains$samples
		for (j in 1:nchains) {
			coeff_chains[[j]] <- coeff_chains[[j]][ , pars]
			colnames(coeff_chains[[j]]) <- paste0('K', 1:nK)
		}

		coeffs[[i]] <- mcmc_intervals(
				coeff_chains,
				pars = paste0('K', 1:nK)
			) + 
			geom_vline(xintercept = 0) +
			xlab('Estimated value') +
			ggtitle(model_term_nice) +
			theme(
				text = element_text(family = 'TT Arial'),
				title = element_text(size = 14),
				axis.text = element_text(size = 10),
        		axis.title = element_text(size = 11),
				plot.background = element_rect(fill = 'white')
			)

	} # next ADMIXTURE model term

	caterpillars <- plot_grid(plotlist = coeffs, nrow = 3)

	ggsave(caterpillars, filename = paste0(out_dir, '/admixture_admix_beta_coefficients.png'), width = 19.2, height = 10.8, dpi = 600, bg = 'white')

say('#######################', pre = 1)
say('### response curves ###')
say('#######################')

	# chains <- readRDS(paste0(out_dir, '/chains.rds'))
	mcmc <- chains$samples

	# 1st is species
	colors <- c('black', 'red', 'green', 'cyan', 'orange', 'yellow', 'violet', 'floralwhite')

	# make a plot for how AG SDM and each ancestral population responds to each predictor
	responses <- list()
	for (pred in 1:n_predictors) {

		# unscaled predictor value
		x <- env_array[ , pred, pred]

		# create data frame with SDM response
		pars <- paste0('response_curves_sdm[', 1:n_response_curve_rows, ', ', pred, ']')
		sdm_response_mean <- chains$summary$all.chains[pars, 'Mean']
		sdm_response_lower <- chains$summary$all.chains[pars, '95%CI_low']
		sdm_response_upper <- chains$summary$all.chains[pars, '95%CI_upp']

		# data frames to hold predictions in long format
		df_mean <- data.frame(
			x = x,
			response = sdm_response_mean / max(sdm_response_mean),
			type = 'Species-level'
		)

		df_lower <- data.frame(
			x = x,
			response = sdm_response_lower,
			type = 'Species-level'
		)

		df_upper <- data.frame(
			x = x,
			response = sdm_response_upper,
			type = 'Species-level'
		)

		# add to data frames ADMIXTURE responses
		admix_response_mean <- admix_response_lower <- admix_response_upper <- list()
		for (k in 1:nK) {

			# response of this ancestral population to this predictor
			pars <- paste0('response_curves_admix[', 1:n_response_curve_rows, ', ', k, ', ', pred, ']')
			admix_response_mean <- chains$summary$all.chains[pars, 'Mean']
			admix_response_lower <- chains$summary$all.chains[pars, '95%CI_low']
			admix_response_upper <- chains$summary$all.chains[pars, '95%CI_upp']

			# data frames to hold predictions in long format
			df_mean <- rbind(
				df_mean,
				data.frame(
					x = x,
					response = admix_response_mean,
					type = paste0('K', k)
				)
			)

			df_lower <- rbind(
				df_lower,
				data.frame(
					x = x,
					response = admix_response_lower,
					type = paste0('K', k)
				)
			)

			df_upper <- rbind(
				df_upper,
				data.frame(
					x = x,
					response = admix_response_upper,
					type = paste0('K', k)
				)
			)

		}


		responses[[pred]] <- ggplot()
		types <- unique(df_mean$type)

		# plot CI first
		for (i in seq_along(types)) {

			type <- types[i]

			this_df_upper <- df_upper[df_mean$type == type, ]
			this_df_lower <- df_lower[df_mean$type == type, ]

			max_score <- max(this_df_upper$response)

			this_df_upper$response <- this_df_upper$response / max_score
			this_df_lower$response <- this_df_lower$response / max_score

			this_df_ci <- data.frame(
				x = c(x, rev(x)),
				y = c(this_df_upper$response, rev(this_df_lower$response)),
				type = type
			)

			responses[[pred]] <- responses[[pred]] +
				geom_polygon(
					data = this_df_ci,
					mapping = aes(x = x, y = y),
					color = NA,
					fill = alpha(colors[i], 0.1)
				)

		}

		# plot trendlines second
		for (i in seq_along(types)) {

			type <- types[i]

			this_df_mean <- df_mean[df_mean$type == type, ]
			this_df_upper <- df_upper[df_mean$type == type, ]
			this_df_lower <- df_lower[df_mean$type == type, ]

			max_score <- max(this_df_mean$response, this_df_upper$response)

			this_df_mean$response <- this_df_mean$response / max_score

			responses[[pred]] <- responses[[pred]] +
				geom_line(
					data = this_df_mean,
					mapping = aes(x = x, y = response),
					color = colors[i]
				)

		}

		if (predictor_names[pred] == 'aridity') {
			nice_title <- 'Aridity (temperature / precipitation)'
			nice_axis <- 'Aridity (°C / mm)'
		} else if (predictor_names[pred] == 'bio7') {
			nice_title <- 'Temperature annual range (BIO 07)'
			nice_axis <- 'Temperature annual range (°C)'
		} else if (predictor_names[pred] == 'bio12') {
			nice_title <- 'Total annual precipitation (BIO 12)'
			nice_axis <- 'Total annual precipitation (mm)'
		} else if (predictor_names[pred] == 'bio15') {
			nice_title <- 'Precipitation seasonality (BIO 15)'
			nice_axis <- 'Precipitation seasonality'
		} else if (predictor_names[pred] == 'bio15') {
			nice_title <- 'Precipitation seasonality (BIO 15)'
			nice_axis <- 'Precipitation seasonality'
		} else if (predictor_names[pred] == 'ph') {
			nice_title <- 'Soil pH'
			nice_axis <- 'pH'
		} else if (predictor_names[pred] == 'sand') {
			nice_title <- 'Soil proportion sand'
			nice_axis <- 'Proportion sand'
		} else if (predictor_names[pred] == 'silt') {
			nice_title <- 'Soil proportion silt'
			nice_axis <- 'Proportion silt'
		} else if (predictor_names[pred] == 'soc') {
			nice_title <- 'Soil organic matter'
			nice_axis <- 'Soil organic matter'
		}

		responses[[pred]] <- responses[[pred]] +
			ylim(0, 1) +
			scale_color_discrete(name = 'Entity') +
			scale_fill_discrete(name = 'Entity') +
			xlab(nice_axis) +
			ylab('Scaled ancestry proportion or\nSDM-predicted abundance') +
			ggtitle(nice_title) +
			theme(
				plot.title = element_text(size = 10)
			)

	} # next predictor

	if (length(predictor_names) <= 4) {
		nrow <- 1
		width <- 12
		height <- 4
	} else {
		nrow <- 2
		width <- 12
		height <- 10	
	}

	responses <- plot_grid(plotlist = responses, nrow = nrow)

	ggsave(plot = responses, filename = paste0(out_dir, '/response_curves.png'), width = width, height = height, dpi = 600)


say('FINIS!', deco = '+', level = 1)

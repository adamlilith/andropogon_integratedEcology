### MODELING ANDROPOGON GERARDI DISTRIBUTION, PHENOTYPE, PHYSIOLOGY, GENOTYPE, and ASSOCIATED MICROBIAL COMMUNITIES
### Erica Newman | Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2023-12
###
### This script constructs an integrated species distribution model for Andropogon gerardi, where occurrence is a function of environmental variables and biomass, which is in turn also a function of environmental variables.
###
### source('C:/Ecology/R/andropogon_integratedEcology/sdm_pdm/sdm_pdm_04_integrated_model_occurrence_fx_of_biomass.r')
### source('C:/Subarashi/R/andropogon_integratedEcology/sdm_pdm/sdm_pdm_04_integrated_model_occurrence_fx_of_biomass.r')
###
### CONTENTS ###
### setup ###

#############
### setup ###
#############

	rm(list = ls())

	drive <- 'C:/Ecology/'
	# drive <- 'C:/Subarashi/'

	setwd(paste0(drive, '/Research/Andropogon/Andropogon'))

	library(data.table) # fast data tables
	library(enmSdmX) # GIS & SDMing
	library(nimble) # Bayesian modeling
	library(nimbleHMC) # Bayesian modeling
	library(omnibus) # utilities
	library(predicts) # GIS & SDMing
	library(readxl) # Excel
	library(terra) # spatial objects

	### user-defined values
	#######################

		predictor_names <- c('bio1', 'bio12', 'bio15') # should include predictors from SDM *and* phenotype model!!!

		occ_form <- ~ 1 + bio1 + bio12 + bio15 + I(bio1^2) + I(bio12^2) + I(bio15^2)
		pheno_form_mu <- ~ 1 + bio12
		pheno_form_sigma <- ~ 1 + bio12 + I(bio12^2)

		psa_quant <- 0.99 # to define pseudoabsences, use this quantile of Poaceae occurrences across counties with Poaceae occurrences 
		# psa_quant <- 0 # to define pseudoabsences, use this quantile of Poaceae occurrences across counties with Poaceae occurrences 

		# ### MCMC settings
		# niter <- 280000
		# nburnin <- 40000
		# thin <- 200
		# nchains <- 4
		# trial <- FALSE

		# ### MCMC settings FOR TESTING
		niter <- 220
		nburnin <- 20
		thin <- 2
		nchains <- 2
		trial <- TRUE

		out_dir <- paste0('./outputs_loretta/sdm_pdm_occurence_conditional_on_biomass', ifelse(trial, '_TRIAL', NULL), '/')
		dirCreate(out_dir)

		# number of values in response curve array used to depict responses of SDM and PDM
		n_response_curve_values <- 200
		
#################################################
### model occurrence conditional on phenotype ###
#################################################

	sink(paste0(out_dir, '/runtime_log.txt'), split = TRUE)
	say('MODELING OCCURRENCE CONDITIONAL ON PHENOTYPE - SIMPLE MODEL WITH JUST ONE PHENOTYPIC VARIABLE')
	say('Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | ', date(), post = 1)

	### data collation
	##################
	say('data collation', post = 1)

	say('This model is for the spatial distribution of AG conditional on a single phenotypic variable. Occurrences are at the county level, so county area is used as an offset of occurrence, as is the number of Poaceae collections in each county.', breaks = 60, post = 1)

	say('MCMC settings:', level = 2)
	say('niter ...... .', niter)
	say('nburnin ..... ', nburnin)
	say('thin ........ ', thin)
	say('nchains ..... ', nchains)
	say('psa_quant ... ', psa_quant, post = 1)

	##################################
	### model predictors and terms ###
	##################################

	say('Predictors and model formula:', level = 2)
	say('Predictors: ', paste(predictor_names, collapse = ', '))
	say('SDM formula: ', as.character(occ_form))
	say('Phenotype formula for mean biomass: ', as.character(pheno_form_mu))
	say('Phenotype formula for SD of biomass: ', as.character(pheno_form_sigma))

	terms <- terms(occ_form)
	terms <- attr(terms, 'term.labels')
	linear_terms <- terms
	linear_terms <- linear_terms[!grepl(linear_terms, pattern = '\\^2')]
	linear_terms <- linear_terms[!grepl(linear_terms, pattern = '\\:')]
	linear_terms <- linear_terms[!grepl(linear_terms, pattern = '\\*')]
	occ_predictors <- linear_terms
	
	terms <- terms(pheno_form_mu)
	terms <- attr(terms, 'term.labels')
	linear_terms <- terms
	linear_terms <- linear_terms[!grepl(linear_terms, pattern = '\\^2')]
	linear_terms <- linear_terms[!grepl(linear_terms, pattern = '\\:')]
	linear_terms <- linear_terms[!grepl(linear_terms, pattern = '\\*')]
	biomass_predictors_mu <- linear_terms
	
	terms <- terms(pheno_form_sigma)
	terms <- attr(terms, 'term.labels')
	linear_terms <- terms
	linear_terms <- linear_terms[!grepl(linear_terms, pattern = '\\^2')]
	linear_terms <- linear_terms[!grepl(linear_terms, pattern = '\\:')]
	linear_terms <- linear_terms[!grepl(linear_terms, pattern = '\\*')]
	biomass_predictors_sigma <- linear_terms
	
	occ_and_pheno_predictors <- sort(unique(c(occ_predictors, biomass_predictors_mu, biomass_predictors_sigma)))
	n_occ_predictors <- length(occ_predictors)
	n_biomass_predictors_mu <- length(biomass_predictors_mu)
	n_biomass_predictors_sigma <- length(biomass_predictors_sigma)

	### collate data for PDM
	########################
	say('load and collate site and biomass data for phenotype distribution model', level = 3)

	site_data_raw <- readRDS('./data_from_loretta/sdm_pdm_00_merged_site_data_with_climate/sites.rds')
	biomass_data_raw <- readRDS('./data_from_loretta/sdm_pdm_00_merged_site_data_with_climate/biomass.rds')

	# ensure sites in each array appear in the same order
	site_data_raw <- site_data_raw[order(site_id)]
	biomass_data_raw <- biomass_data_raw[order(SITE)]

	stopifnot(all(site_data_raw$site_id == unique(biomass_data_raw$SITE)))

	n_pheno_sites <- length(unique(site_data_raw$site_id))
	
	# center and scale at site level
	x <- site_data_raw[ , ..biomass_predictors_mu]
	x <- scale(x)
	biomass_x_centers <- attr(x, 'scaled:center')
	biomass_x_scales <- attr(x, 'scaled:scale')
	x <- as.data.frame(x)
	biomass_by_site_x_mu <- model.matrix(pheno_form_mu, x)
	biomass_by_site_x_sigma <- model.matrix(pheno_form_sigma, x)

	# make vector of which sampled site matches each row in the biomass and morphology/physiology data
	n_biomass <- nrow(biomass_data_raw)

	biomass_site_index <- rep(NA, n_biomass)
	for (i in 1:n_biomass) {
		biomass_site_index[i] <- which(site_data_raw$site_id == biomass_data_raw$SITE[i])
	}
	
	### load occurrence data
	########################
	say('load occurrences and present-day + future environmental data', level = 3)
	ag_vect_sq <- vect('./data_from_adam_and_loretta/andropogon_gerardi_occurrences_with_environment_1961_2020.gpkg')
	
	### future climates
	ag_vect_ssp245_2041_2070 <- vect('./data_from_adam_and_loretta/andropogon_gerardi_occurrences_with_environment_ensemble_8GCMs_ssp245_2041_2070.gpkg')
	ag_vect_ssp245_2071_2100 <- vect('./data_from_adam_and_loretta/andropogon_gerardi_occurrences_with_environment_ensemble_8GCMs_ssp245_2071_2100.gpkg')
	ag_vect_ssp370_2041_2070 <- vect('./data_from_adam_and_loretta/andropogon_gerardi_occurrences_with_environment_ensemble_8GCMs_ssp370_2041_2070.gpkg')
	ag_vect_ssp370_2071_2100 <- vect('./data_from_adam_and_loretta/andropogon_gerardi_occurrences_with_environment_ensemble_8GCMs_ssp370_2071_2100.gpkg')

	fields <- c('area_km2', 'any_ag_quality_1_to_3', 'num_poaceae_records', occ_and_pheno_predictors)
	ag_vect_sq <- ag_vect_sq[ , fields]

	# for counties with NA Poaceae and AG, assign a maximal number of Poaceae and 0 AG
	n_pseudoabs <- quantile(ag_vect_sq$num_poaceae_records, psa_quant, na.rm = TRUE)
	ag_vect_sq$num_poaceae_records[is.na(ag_vect_sq$num_poaceae_records)] <- n_pseudoabs
	ag_vect_sq$any_ag_quality_1_to_3[is.na(ag_vect_sq$any_ag_quality_1_to_3)] <- 0

	### collate data for SDM (and PDM)
	##################################
	say('collate SDM data', level = 3)

	ag_sq <- as.data.frame(ag_vect_sq)

	### county area... used as an offset to account for sampling bias
	area_km2 <- ag_sq$area_km2
	log_area_km2 <- log(area_km2)
	log_area_km2_scaled <- scale(log_area_km2)
	log_area_center <- attributes(log_area_km2_scaled)$`scaled:center`
	log_area_scale <- attributes(log_area_km2_scaled)$`scaled:scale`
	log_area_km2_scaled <- log_area_km2_scaled[ , 1]

	### number of Poaceae records... used as an offset to account for sampling bias
	log_num_poaceae_records <- log1p(ag_sq$num_poaceae_records) # log(occ_x_counties_sq + 1)
	log_num_poaceae_records_scaled <- scale(log_num_poaceae_records)
	log_num_poaceae_records_scaled <- log_num_poaceae_records_scaled[ , 1]

	### subset, scale, and manipulate predictors into model frame for SDM
	# status quo: counties with training data
	log_area_km2 <- log(ag_sq$area_km2)
	log_area_km2 <- data.frame(log_area_km2 = log_area_km2)

	x_raw <- ag_sq[ , occ_predictors, drop = FALSE]
	x_raw <- as.data.frame(x_raw)
	
	x_raw_scaled <- scale(x_raw)
	occ_x_centers <- attr(x_raw_scaled, 'scaled:center')
	occ_x_scales <- attr(x_raw_scaled, 'scaled:scale')
	x_raw_scaled <- as.data.frame(x_raw_scaled)
	occ_x_counties_sq <- model.matrix(occ_form, x_raw_scaled)
	
	### subset, scale, and manipulate predictors into model frame for phenotype
	x_raw <- ag_sq[ , biomass_predictors_mu, drop = FALSE]
	x_raw_scaled <- scale(x_raw, scale = biomass_x_scales, center = biomass_x_centers)
	x_raw_scaled <- as.data.frame(x_raw_scaled)
	biomass_x_counties_mu_sq <- model.matrix(pheno_form_mu, x_raw_scaled)
	biomass_x_counties_sigma_sq <- model.matrix(pheno_form_sigma, x_raw_scaled)

	# futures
	for (ssp in c('ssp245_2041_2070', 'ssp245_2071_2100', 'ssp370_2041_2070', 'ssp370_2071_2100')) {

		# future: counties with future data
		x_fut <- get(paste0('ag_vect_', ssp))
		x_fut <- as.data.frame(x_fut)

		# SDM
		this_x_fut <- x_fut[ , occ_predictors, drop = FALSE]
		x_fut_scaled <- scale(this_x_fut, center = occ_x_centers, scale = occ_x_scales)
		x_fut_scaled <- as.data.frame(x_fut_scaled)
		mm <- model.matrix(occ_form, x_fut_scaled)
		assign(paste0('occ_x_counties_', ssp), mm)

		# PDM: mean of biomass
		this_x_fut <- x_fut[ , biomass_predictors_mu, drop = FALSE]
		x_fut_scaled <- scale(this_x_fut, center = biomass_x_centers, scale = biomass_x_scales)
		x_fut_scaled <- as.data.frame(x_fut_scaled)
		mm <- model.matrix(pheno_form_mu, x_fut_scaled)
		assign(paste0('pheno_x_counties_mu_', ssp), mm)

		mm <- model.matrix(pheno_form_sigma, x_fut_scaled)
		assign(paste0('pheno_x_counties_sigma_', ssp), mm)

	}

	### response array for SDM
	##########################
	say('response array for occurrence', level = 3)

	# This is an array for plotting the response curves of the SDM
	# columns: variables, all but one held at mean value across counties with AG
	# rows: focal variable increments, others held constant
	# depths: identity of focal variable changes

	# generalization
	form <- occ_form
	x_centers <- occ_x_centers
	x_scales <- occ_x_scales
	linear_terms <- occ_predictors

	n_linear_terms <- length(linear_terms)

	env_array <- array(
		NA,
		dim = c(n_response_curve_values, n_linear_terms, n_linear_terms),
		dimnames = list(
			1:n_response_curve_values,
			linear_terms,
			linear_terms
		)
	)

	# calculate means of each variable across counties with AG presences
	ag_pres <- ag_vect_sq[ag_vect_sq$any_ag_quality_1_to_3 > 0]
	ag_pres <- as.data.frame(ag_pres)[ , linear_terms, drop = FALSE]
	means <- colMeans(ag_pres)

	for (i in seq_along(linear_terms)) {
		for (j in seq_along(linear_terms)) {
			env_array[ , i, j] <- means[i]
		}
	}

	mins <- apply(ag_pres[ , linear_terms, drop = FALSE], 2, min)
	maxs <- apply(ag_pres[ , linear_terms, drop = FALSE], 2, max)

	for (i in seq_along(linear_terms)) {
		env_array[ , i, i] <- seq(mins[i], maxs[i], length.out = n_response_curve_values)
	}

	# scale
	env_array_scaled <- env_array
	for (i in seq_along(linear_terms)) {
		env_array_scaled[ , , i] <- sweep(env_array_scaled[ , , i], 2, x_centers[linear_terms], '-')
		env_array_scaled[ , , i] <- sweep(env_array_scaled[ , , i], 2, x_scales[linear_terms], '/')
	}

	# make into model matrix
	mm_array <- list()
	for (i in seq_along(linear_terms)) {
		mm_array[[i]] <- model.matrix(form, as.data.frame(env_array_scaled[ , , i, drop = TRUE]))
	}

	nc <- ncol(mm_array[[1]])
	response_curve_x <- array(
		as.numeric(unlist(mm_array)),
		dim = c(n_response_curve_values, nc, n_linear_terms),
		dimnames = list(1:n_response_curve_values, colnames(mm_array[[1]]), linear_terms)
	)

	# remember
	occ_x_response_curve <- response_curve_x

	### response array for PDM: 1 phenotypic response
	#################################################
	say('response array for phenotype: 1 phenotypic response', level = 3)

	# This is an array for plotting the response curves of the PDM
	# For a single phenotypic response, the output will be 2D matrix.
	# columns: variables, all but one held at mean value across *sampled sites*
	# rows: focal variable increments, others held constant
	# depths: identity of focal variable changes

	# generalization
	form_mu <- pheno_form_mu
	form_sigma <- pheno_form_sigma
	x_centers <- biomass_x_centers
	x_scales <- biomass_x_scales
	linear_terms_mu <- biomass_predictors_mu
	linear_terms_sigma <- biomass_predictors_sigma

	n_linear_terms_mu <- length(linear_terms_mu)
	n_linear_terms_sigma <- length(linear_terms_sigma)

	env_array <- matrix(
		NA,
		nrow = n_response_curve_values,
		ncol = n_linear_terms,
		dimnames = list(
			1:n_response_curve_values,
			linear_terms
		)
	)

	# calculate means of each variable across sampled sites
	means <- colMeans(site_data_raw[ , ..biomass_predictors_mu])

	for (i in seq_along(linear_terms)) {
		env_array[ , i] <- means[i]
	}

	mins <- apply(ag_pres[ , linear_terms, drop = FALSE], 2, min)
	maxs <- apply(ag_pres[ , linear_terms, drop = FALSE], 2, max)

	for (i in seq_along(linear_terms)) {
		env_array[ , i] <- seq(mins[i], maxs[i], length.out = n_response_curve_values)
	}

	# scale
	env_array_scaled <- env_array
	for (i in seq_along(linear_terms)) {
		env_array_scaled <- sweep(env_array_scaled, 2, x_centers[linear_terms], '-')
		env_array_scaled <- sweep(env_array_scaled, 2, x_scales[linear_terms], '/')
	}

	# make into model matrix
	response_curve_x_mu <- model.matrix(form_mu, as.data.frame(env_array_scaled))
	response_curve_x_sigma <- model.matrix(form_sigma, as.data.frame(env_array_scaled))

	pheno_x_response_curve_mu <- response_curve_x_mu
	pheno_x_response_curve_sigma <- response_curve_x_sigma

	### inputs for nimble
	say('Inputs:', level = 2)
	data <- list(
		y_occ = ag_sq$any_ag_quality_1_to_3, # number of AG records in each county
		y_biomass = biomass_data_raw$Biomass
	)

	n_counties <- nrow(ag_sq)

	n_occ_terms <- ncol(occ_x_counties_sq)
	n_pheno_terms_mu <- ncol(biomass_x_counties_mu_sq)
	n_pheno_terms_sigma <- ncol(biomass_x_counties_sigma_sq)

	constants <- list(
		
		n_pheno_sites = n_pheno_sites, # number of phenotype sample sites
		n_biomass = n_biomass, # number of biomass records

		n_counties = n_counties,

		n_occ_terms = n_occ_terms, # number of terms in SDM, including intercept and non-linear terms
		n_pheno_terms_mu = n_pheno_terms_mu, # number of terms in PDM of mean biomass, including intercept and non-linear terms
		n_pheno_terms_sigma = n_pheno_terms_sigma, # number of terms in PDM of SD of biomass, including intercept and non-linear terms

		biomass_by_site_x_mu = biomass_by_site_x_mu, # model matrix for mean site-level biomass at sampled sites
		biomass_by_site_x_sigma = biomass_by_site_x_sigma, # model matrix for site-level SD of biomass at sampled sites
		biomass_site_index = biomass_site_index, # index of sampled site for each row in biomass data

		n_occ_predictors = n_occ_predictors, # number of linear terms in SDM, excluding intercept
		# n_biomass_predictors_mu = n_biomass_predictors_mu, # number of linear terms in PDM, excluding intercept
		n_response_curve_values = n_response_curve_values,
		occ_x_response_curve = occ_x_response_curve,
		pheno_x_response_curve_mu = pheno_x_response_curve_mu,
		pheno_x_response_curve_sigma = pheno_x_response_curve_sigma,

		log_area_km2_scaled = log_area_km2_scaled,
		log_num_poaceae_records_scaled = log_num_poaceae_records_scaled,

		occ_x_counties_sq = occ_x_counties_sq,
		biomass_x_counties_mu_sq = biomass_x_counties_mu_sq,
		biomass_x_counties_sigma_sq = biomass_x_counties_sigma_sq,

		occ_x_counties_ssp245_2041_2070 = occ_x_counties_ssp245_2041_2070,
		occ_x_counties_ssp245_2071_2100 = occ_x_counties_ssp245_2071_2100,

		occ_x_counties_ssp370_2041_2070 = occ_x_counties_ssp370_2041_2070,
		occ_x_counties_ssp370_2071_2100 = occ_x_counties_ssp370_2071_2100,

		pheno_x_counties_mu_ssp245_2071_2100 = pheno_x_counties_mu_ssp245_2071_2100,
		pheno_x_counties_mu_ssp245_2041_2070 = pheno_x_counties_mu_ssp245_2041_2070,
		pheno_x_counties_mu_ssp370_2041_2070 = pheno_x_counties_mu_ssp370_2041_2070,
		pheno_x_counties_mu_ssp370_2071_2100 = pheno_x_counties_mu_ssp370_2071_2100,

		pheno_x_counties_sigma_ssp245_2071_2100 = pheno_x_counties_sigma_ssp245_2071_2100,
		pheno_x_counties_sigma_ssp245_2041_2070 = pheno_x_counties_sigma_ssp245_2041_2070,
		pheno_x_counties_sigma_ssp370_2041_2070 = pheno_x_counties_sigma_ssp370_2041_2070,
		pheno_x_counties_sigma_ssp370_2071_2100 = pheno_x_counties_sigma_ssp370_2071_2100

	)

	lambda_fut_inits <- rep(mean(ag_sq$any_ag_quality_1_to_3), n_counties)
	N_inits <- 10 * ag_sq$any_ag_quality_1_to_3 + 1
	lambda_sq_inits <- 1 + ag_sq$any_ag_quality_1_to_3

	occ_beta_inits <- rep(0, n_occ_terms)
	occ_beta_inits[grepl(colnames(occ_x_counties_sq), pattern = '\\^2')] <- -1
	occ_beta_inits[grepl(colnames(occ_x_counties_sq), pattern = '\\:')] <- 0

	pheno_beta_mu_inits <- rep(0, n_pheno_terms_mu)
	pheno_beta_mu_inits[grepl(colnames(biomass_x_counties_mu_sq), pattern = '\\^2')] <- -1
	pheno_beta_mu_inits[grepl(colnames(biomass_x_counties_mu_sq), pattern = '\\:')] <- 0

	pheno_beta_sigma_inits <- rep(0, n_pheno_terms_sigma)
	pheno_beta_sigma_inits[grepl(colnames(biomass_x_counties_sigma_sq), pattern = '\\^2')] <- 1
	pheno_beta_sigma_inits[grepl(colnames(biomass_x_counties_sigma_sq), pattern = '\\:')] <- 0

	pheno_biomass_sq_inits <- pheno_biomass_fut_inits <- rep(mean(biomass_data_raw$Biomass), n_counties)

	inits <- list(

		y_biomass = rep(mean(biomass_data_raw$Biomass), n_biomass), # biomass

		occ_beta = occ_beta_inits, # SDM coefficients
		pheno_beta_mu = pheno_beta_mu_inits, # PDM coefficients
		pheno_beta_sigma = pheno_beta_sigma_inits, # PDM coefficients
		occ_by_pheno_beta = 0, # coefficient for relationship between SDM and PDM
		
		occ_alpha0_sampling = 0, # SDM sampling bias coefficient
		occ_alpha_area = 1, # SDM sampling bias coefficient
		occ_alpha_poaceae = 1, # SDM sampling bias coefficient
		
		occ_N = N_inits, # numler of AG in a county

		occ_lambda_sq = lambda_sq_inits, # latent abundance of AG in counties

		pheno_biomass_site_mu_hat = rep(mean(biomass_data_raw$Biomass), n_pheno_sites), # latent mean biomass at a site
		pheno_biomass_site_sigma_hat = rep(sd(biomass_data_raw$Biomass), n_pheno_sites), # latent sd of biomass at a site

		occ_lambda_ssp245_2041_2070 = lambda_fut_inits, # future AG abundance in counties
		occ_lambda_ssp245_2071_2100 = lambda_fut_inits,
		occ_lambda_ssp370_2041_2070 = lambda_fut_inits,
		occ_lambda_ssp370_2071_2100 = lambda_fut_inits,

		biomass_county_mu_ssp245_2041_2070 = pheno_biomass_fut_inits, # future biomass in counties
		biomass_county_mu_ssp245_2071_2100 = pheno_biomass_fut_inits,
		biomass_county_mu_ssp370_2041_2070 = pheno_biomass_fut_inits,
		biomass_county_mu_ssp370_2071_2100 = pheno_biomass_fut_inits

	)

	say('Data:')
	print(str(data))

	say('Constants:', pre = 1)
	print(str(constants))

	say('Initializations:', pre = 1)
	print(str(inits))

	### define model
	say('nimbleCode():', level = 2)
	say('For the occurrence component, we assume an N-mixture model (latent, real abundance ~ Poisson, and observations ~ binomial draws from latent abundance.', post = 1)
	say('For the biomass component, We assume a log-linear model model (log(biomass) ~ environment * coefficients.', post = 2)
	model_code <- nimbleCode({
	  
		# phenotype: priors for relationship of mean and variance to environment
		pheno_beta_mu[1] ~ dnorm(0, sd = 10) # intercept... not regularized
		for (i in 2:n_pheno_terms_mu) {
			# pheno_beta_mu[i] ~ ddexp(0, rate = 1) # regularization toward 0
			pheno_beta_mu[i] ~ dnorm(0, sd = 10) # broad prior
		}

		pheno_beta_sigma[1] ~ dnorm(0, sd = 10) # intercept... not regularized
		for (i in 2:n_pheno_terms_sigma) {
			# pheno_beta_sigma[i] ~ ddexp(0, rate = 1) # regularization toward 0
			pheno_beta_sigma[i] ~ dnorm(0, sd = 10) # broad prior
		}

		# occurrences: priors for relationship to phenotype
		occ_by_pheno_beta ~ dnorm(0, sd = 10)

		# occurrences: priors for relationship to environment
		occ_beta[1] ~ dnorm(0, sd = 10) # intercept... not regularized
		for (j in 2:n_occ_terms) {
			# occ_beta[j] ~ ddexp(0, rate = 1) # regularization toward 0
			occ_beta[j] ~ dnorm(0, sd = 10) # broad prior
		}

		# occurrences: priors for sampling bias
		occ_alpha0_sampling ~ dnorm(0, sd = 10)
		occ_alpha_area ~ dnorm(0, sd = 10)
		occ_alpha_poaceae ~ dnorm(0, sd = 10)

		# likelihood of biomass
		# mean and sd of biomass at a site are latent and functions of environment
		# plant biomasses are samples from the distribution defined by the mean and sd of the given site
		for (i in 1:n_pheno_sites) {
			
			log(pheno_biomass_site_mu_hat[i]) <-  inprod(pheno_beta_mu[1:n_pheno_terms_mu], biomass_by_site_x_mu[i, 1:n_pheno_terms_mu])

			log(pheno_biomass_site_sigma_hat[i]) <-  inprod(pheno_beta_sigma[1:n_pheno_terms_sigma], biomass_by_site_x_sigma[i, 1:n_pheno_terms_sigma])
			
		}

		for (i in 1:n_biomass) {

			log(y_biomass[i]) ~ dnorm(pheno_biomass_site_mu_hat[biomass_site_index[i]], sd = pheno_biomass_site_sigma_hat[biomass_site_index[i]])
	
		}
			
		# likelihood of occurrences
		for (i in 1:n_counties) {
			
			### actual abundance (latent--unobserved)
			occ_N[i] ~ dpois(occ_lambda_sq[i])

			# relationship between latent abundance and environment and phenotype
			log(occ_lambda_sq[i]) <-
				inprod(occ_beta[1:n_occ_terms], occ_x_counties_sq[i, 1:n_occ_terms]) + # environmental dependency
				occ_by_pheno_beta * biomass_county_mu_sq[i] # phenotypic dependency

			### observed number of AG
			y_occ[i] ~ dbin(prob = p[i], size = occ_N[i])

			# sampling bias
			logit(p[i]) <- occ_alpha0_sampling + occ_alpha_area * log_area_km2_scaled[i] + occ_alpha_poaceae * log_num_poaceae_records_scaled[i]

		}

		# posterior samplers for predictions to counties in status quo and future
		for (i in 1:n_counties) {

			# biomass mean and sigma: status quo
			log(biomass_county_mu_sq[i]) <- inprod(pheno_beta_mu[1:n_pheno_terms_mu], biomass_x_counties_mu_sq[i, 1:n_pheno_terms_mu])
			log(biomass_county_sigma_sq[i]) <- inprod(pheno_beta_sigma[1:n_pheno_terms_sigma], biomass_x_counties_sigma_sq[i, 1:n_pheno_terms_sigma])

			# biomass mean: future
			log(biomass_county_mu_ssp245_2041_2070[i]) <- inprod(pheno_beta_mu[1:n_pheno_terms_mu], pheno_x_counties_mu_ssp245_2041_2070[i, 1:n_pheno_terms_mu])
			log(biomass_county_mu_ssp245_2071_2100[i]) <- inprod(pheno_beta_mu[1:n_pheno_terms_mu], pheno_x_counties_mu_ssp245_2071_2100[i, 1:n_pheno_terms_mu])
			
			log(biomass_county_mu_ssp370_2041_2070[i]) <- inprod(pheno_beta_mu[1:n_pheno_terms_mu], pheno_x_counties_mu_ssp370_2041_2070[i, 1:n_pheno_terms_mu])
			log(biomass_county_mu_ssp370_2071_2100[i]) <- inprod(pheno_beta_mu[1:n_pheno_terms_mu], pheno_x_counties_mu_ssp370_2071_2100[i, 1:n_pheno_terms_mu])

			# biomass sigma: future
			log(biomass_county_sigma_ssp245_2041_2070[i]) <- inprod(pheno_beta_sigma[1:n_pheno_terms_sigma], pheno_x_counties_sigma_ssp245_2041_2070[i, 1:n_pheno_terms_sigma])
			log(biomass_county_sigma_ssp245_2071_2100[i]) <- inprod(pheno_beta_sigma[1:n_pheno_terms_sigma], pheno_x_counties_sigma_ssp245_2071_2100[i, 1:n_pheno_terms_sigma])
			
			log(biomass_county_sigma_ssp370_2041_2070[i]) <- inprod(pheno_beta_sigma[1:n_pheno_terms_sigma], pheno_x_counties_sigma_ssp370_2041_2070[i, 1:n_pheno_terms_sigma])
			log(biomass_county_sigma_ssp370_2071_2100[i]) <- inprod(pheno_beta_sigma[1:n_pheno_terms_sigma], pheno_x_counties_sigma_ssp370_2071_2100[i, 1:n_pheno_terms_sigma])

			# occurrence
			log(occ_lambda_ssp245_2041_2070[i]) <- inprod(occ_beta[1:n_occ_terms], occ_x_counties_ssp245_2041_2070[i, 1:n_occ_terms])
			log(occ_lambda_ssp245_2071_2100[i]) <- inprod(occ_beta[1:n_occ_terms], occ_x_counties_ssp245_2071_2100[i, 1:n_occ_terms])
			
			log(occ_lambda_ssp370_2041_2070[i]) <- inprod(occ_beta[1:n_occ_terms], occ_x_counties_ssp370_2041_2070[i, 1:n_occ_terms])
			log(occ_lambda_ssp370_2071_2100[i]) <- inprod(occ_beta[1:n_occ_terms], occ_x_counties_ssp370_2071_2100[i, 1:n_occ_terms])

		}

		# posterior predictive sampler for response curves for occurrence
		for (j in 1:n_occ_predictors) {
			for (i in 1:n_response_curve_values) {
				
				occ_response_curves_lambda[i, j]) <-
					inprod(occ_beta[1:n_occ_terms], occ_x_response_curve[i, 1:n_occ_terms, j]

				# exp_response_curves_sdm_lambda[i, j] <- exp(occ_response_curves_lambda[i, j])

			}
		}

		# posterior predictive sampler for response curves for biomass: assuming one predictor
		for (i in 1:n_response_curve_values) {
				
			pheno_response_curves_biomass_mu[i] <-
				exp(inprod(pheno_beta_mu[1:n_pheno_terms_mu], pheno_x_response_curve_mu[i, 1:n_pheno_terms_mu]))
			
			pheno_response_curves_biomass_sigma[i] <-
				exp(inprod(pheno_beta_sigma[1:n_pheno_terms_sigma], pheno_x_response_curve_sigma[i, 1:n_pheno_terms_sigma]))
			
		}

	})

	print(model_code)

	say('nimbleModel():', level = 2)
	model_species <- nimbleModel(
		code = model_code, # our model
		constants = constants, # constants
		data = data, # data
		inits = inits, # initialization values
		check = TRUE, # any errors?
		calculate = FALSE
		# buildDerivs = TRUE # need for Hamiltonian Monte Carlo
	)

	say('$initializeInfo() and $calculate():', level = 2)
	model_species$initializeInfo()
	model_species$calculate()

	say('configureMCMC():', level = 2)
	monitors <- c(
		'occ_beta', 'occ_alpha0_sampling', 'occ_alpha_area', 'occ_alpha_poaceae',
		'pheno_beta_mu', 'pheno_beta_sigma',
		'occ_by_pheno_beta',
		'occ_response_curves_lambda',
		'pheno_response_curves_biomass_mu', 'pheno_response_curves_biomass_sigma',
		'occ_lambda_sq',
		'pheno_biomass_site_mu_hat', 'pheno_biomass_site_sigma_hat',
		'biomass_county_mu_sq', 'biomass_county_sigma_sq',
		'occ_lambda_ssp245_2041_2070', 'occ_lambda_ssp245_2071_2100', 'occ_lambda_ssp370_2041_2070', 'occ_lambda_ssp370_2071_2100',
		'biomass_county_mu_ssp245_2041_2070', 'biomass_county_mu_ssp245_2071_2100', 'biomass_county_mu_ssp370_2041_2070', 'biomass_county_mu_ssp370_2071_2100',
		'biomass_county_sigma_ssp245_2041_2070', 'biomass_county_sigma_ssp245_2071_2100', 'biomass_county_sigma_ssp370_2041_2070', 'biomass_county_sigma_ssp370_2071_2100'
	)
	
	conf <- configureMCMC(
		model_species,
		monitors = monitors,
		print = TRUE,
		enableWAIC = FALSE
	)

	# # add no U-turn sampler (Hamiltonian Monte Carlo)
	# conf$addSampler(target = c('occ_beta', 'occ_alpha_area', 'occ_alpha_poaceae'), type = 'NUTS')

	### compile/build/run model/save MCMC
	build <- buildMCMC(conf)

	compiled <- compileNimble(model_species, build, showCompilerOutput = FALSE)

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

	saveRDS(chains, paste0(out_dir, '/sdm_nmixture_chains.rds'))

	say('session info', level = 2)
	print(sessionInfo())

	say(date(), pre = 1)
	sink()
print(NON)
say('#########################')
say('### model diagnostics ###')
say('#########################')

	chains <- readRDS(paste0(out_dir, '/sdm_nmixture_chains.rds'))
	mcmc <- chains$samples

	for (i in 1:nchains) {
		cols <- c(paste0('occ_beta[', 1:n_occ_terms, ']'), 'occ_alpha0_sampling', 'occ_alpha_poaceae', 'occ_alpha_area')
		mcmc[[i]] <- mcmc[[i]][ , cols]
	}

	# graphing trace plots for all betas
	pars <- paste0('occ_beta[', 1:n_occ_terms, ']')
	file <- paste0(out_dir, '/beta_trace.png')
	ggsave(mcmc_trace(mcmc, pars = pars), file = file, width = 10, height = 8, dpi = 450, bg = 'white')

	# graphing trace plots for all "extra" betas
	pars <- c('occ_alpha0_sampling', 'occ_alpha_poaceae', 'occ_alpha_area')
	file <- paste0(out_dir, '/alpha_trace.png')
	ggsave(mcmc_trace(mcmc, pars = pars), file = file, width = 10, height = 4, dpi = 450, bg = 'white')

	# graphing density plots for all betas
	pars <- paste0('occ_beta[', 1:n_occ_terms, ']')
	file <- paste0(out_dir, '/beta_density.png')
	ggsave(mcmc_dens_overlay(mcmc, pars = pars), file = file, width = 12, height = 8, dpi = 450, bg = 'white')

	# graphing trace plots for all "extra" betas
	pars <- c('occ_alpha0_sampling', 'occ_alpha_poaceae', 'occ_alpha_area')
	file <- paste0(out_dir, '/alpha_density.png')
	ggsave(mcmc_dens_overlay(mcmc, pars = pars), file = file, width = 10, height = 4, dpi = 450, bg = 'white')

	# Gelman-Rubin statistic
	rhats <- gelman.diag(mcmc, autoburnin = FALSE, multivariate = TRUE)

	sink(paste0(out_dir, '/convergence.txt'), split = TRUE)
	say('GELMAN-RUBIN STATISTICS')
	say(date(), post = 2)
	print(rhats)
	sink()

say('###################################', pre = 1)
say('### map of current distribution ###')
say('###################################')

	chains <- readRDS(paste0(out_dir, '/sdm_nmixture_chains.rds'))

	# subset chain summary to just the lambdas
	summary <- chains$summary$all.chains

	# just counties with data
	which_lambda <- grepl(rownames(summary), pattern = 'occ_lambda_sq')
	lambda <- summary[which_lambda, ]

	ag_vect_sq$lambda_mean <- lambda[ , 'Mean']
	ag_vect_sq$lambda_0.05ci <- lambda[ , '95%CI_low']
	ag_vect_sq$lambda_0.95ci <- lambda[ , '95%CI_upp']
	ag_vect_sq$lambda_ci <- ag_vect_sq$lambda_0.95ci - ag_vect_sq$lambda_0.05ci

	lambda_sq_quants <- quantile(ag_vect_sq$lambda_mean, c(0.25, 0.5, 0.75, 0.90, 0.95))
	quant_labels <- c('[0 - 0.25)', '[0.25 - 0.5)', '[0.5 - 0.75)', '[0.75 - 0.90)', '[0.90 - 0.95)', '[0.95 - 1]')
	
	ag_vect_sq$quant_col <- NA
	ag_vect_sq$quant_col[ag_vect_sq$lambda_mean < lambda_sq_quants[1]] <- quant_labels[1]
	ag_vect_sq$quant_col[ag_vect_sq$lambda_mean >= lambda_sq_quants[1] & ag_vect_sq$lambda_mean < lambda_sq_quants[2]] <- quant_labels[2]
	ag_vect_sq$quant_col[ag_vect_sq$lambda_mean >= lambda_sq_quants[2] & ag_vect_sq$lambda_mean < lambda_sq_quants[3]] <- quant_labels[3]
	ag_vect_sq$quant_col[ag_vect_sq$lambda_mean >= lambda_sq_quants[3] & ag_vect_sq$lambda_mean < lambda_sq_quants[4]] <- quant_labels[4]
	ag_vect_sq$quant_col[ag_vect_sq$lambda_mean >= lambda_sq_quants[4] & ag_vect_sq$lambda_mean < lambda_sq_quants[5]] <- quant_labels[5]
	ag_vect_sq$quant_col[ag_vect_sq$lambda_mean >= lambda_sq_quants[5]] <- quant_labels[6]

	nam <- vect(paste0(drive, '/Research Data/GADM/Version 4.1/High Res North America Level 1 sans Great Lakes SpatVector WGS84.gpkg'))
	nam <- project(nam, ag_vect_sq)

	ag_vect_sq_pres <- ag_vect_sq[ag_vect_sq$any_ag_quality_1_to_3 > 0]
	extent <- ext(ag_vect_sq_pres)
	extent <- as.vector(extent)
	x_range <- (extent[2] - extent[1])
	y_range <- (extent[4] - extent[3])
	extent[1] <- extent[1] + 0.15 * x_range
	extent[3] <- extent[3] + 0.125 * y_range
	extent[4] <- extent[4] - 0.2 * y_range

	cents_with_ag <- ag_vect_sq[ag_vect_sq$any_ag_quality_1_to_3 > 0]
	cents_with_ag <- centroids(cents_with_ag)

	fill_scale <- c(
		'[0 - 0.25)' = 'gray85',
		'[0.25 - 0.5)' = alpha('forestgreen', 0.20),
		'[0.5 - 0.75)' = alpha('forestgreen', 0.40),
		'[0.75 - 0.90)' = alpha('forestgreen', 0.65),
		'[0.90 - 0.95)' = alpha('forestgreen', 0.75),
		'[0.95 - 1]' = 'forestgreen'
	)

	map <- ggplot() +
		layer_spatial(ag_vect_sq, aes(fill = quant_col), color = NA) +
		layer_spatial(nam, color = 'gray30', fill = NA, linewidth = 0.1) +
		scale_fill_manual(
			name = 'Quantile\n of λ',
			values = fill_scale
		) +
		layer_spatial(cents_with_ag, pch = 16, color = alpha('gray20', 0.5), size = 0.35) +
		guides(fill = guide_legend(reverse = TRUE)) +
		xlim(extent[1], extent[2]) + ylim(extent[3], extent[4]) +
		ggtitle(expression('Present-day distribution of ' * italic('Andropogon gerardi')), subtitle = '1961-2020 | occ_N-mixture model') +
		theme(
			plot.title = element_text(size = 16),
			plot.subtitle = element_text(size = 14)
		)

	ggsave(plot = map, filename = paste0(out_dir, '/sdm_nmixture_lambda_status_quo.png'), width = 12, height = 10, dpi = 600)

	writeVector(ag_vect_sq, paste0(out_dir, '/sdm_nmixture_1961_2020_climate_focus.gpkg'), overwrite = TRUE)

say('###################################')
say('### maps of future distribution ###')
say('###################################')

	futs <- c(
		'ssp245_2041_2070',
		'ssp245_2071_2100',
		'ssp370_2041_2070',
		'ssp370_2071_2100'
	)

	nam <- vect(paste0(drive, '/Research Data/GADM/Version 4.1/High Res North America Level 1 sans Great Lakes SpatVector WGS84.gpkg'))
	nam <- project(nam, ag_vect_sq)

	chains <- readRDS(paste0(out_dir, '/sdm_nmixture_chains.rds'))

	# subset chain summary to just the lambdas
	summary <- chains$summary$all.chains

	for (fut in futs) {

		which_lambda <- grepl(rownames(summary), pattern = paste0('lambda_', fut))
		lambda <- summary[which_lambda, ]

		this_ag_vect <- get(paste0('ag_vect_', fut))

		# get just counties with data
		this_ag_vect$lambda_mean <- lambda[ , 'Mean']
		this_ag_vect$lambda_0.05ci <- lambda[ , '95%CI_low']
		this_ag_vect$lambda_0.95ci <- lambda[ , '95%CI_upp']
		this_ag_vect$lambda_ci <- this_ag_vect$lambda_0.95ci - this_ag_vect$lambda_0.05ci

		quant_labels <- c('[0 - 0.25)', '[0.25 - 0.5)', '[0.5 - 0.75)', '[0.75 - 0.90)', '[0.90 - 0.95)', '[0.95 - 1]') # should be same as above
		# NB lambda_sq_quants is from above!!!
		this_ag_vect$quant_col <- NA
		this_ag_vect$quant_col[this_ag_vect$lambda_mean < lambda_sq_quants[1]] <- quant_labels[1]
		this_ag_vect$quant_col[this_ag_vect$lambda_mean >= lambda_sq_quants[1] & this_ag_vect$lambda_mean < lambda_sq_quants[2]] <- quant_labels[2]
		this_ag_vect$quant_col[this_ag_vect$lambda_mean >= lambda_sq_quants[2] & this_ag_vect$lambda_mean < lambda_sq_quants[3]] <- quant_labels[3]
		this_ag_vect$quant_col[this_ag_vect$lambda_mean >= lambda_sq_quants[3] & this_ag_vect$lambda_mean < lambda_sq_quants[4]] <- quant_labels[4]
		this_ag_vect$quant_col[this_ag_vect$lambda_mean >= lambda_sq_quants[4] & this_ag_vect$lambda_mean < lambda_sq_quants[5]] <- quant_labels[5]
		this_ag_vect$quant_col[this_ag_vect$lambda_mean >= lambda_sq_quants[5]] <- quant_labels[6]

		pretty_title <- paste0('SSP ', substr(fut, 4, 6), ' ', substr(fut, 8, 11), '-', substr(fut, 13, 16), ' | occ_N-mixture model')

		fill_scale <- c(
			'[0 - 0.25)' = 'gray85',
			'[0.25 - 0.5)' = alpha('forestgreen', 0.20),
			'[0.5 - 0.75)' = alpha('forestgreen', 0.40),
			'[0.75 - 0.90)' = alpha('forestgreen', 0.65),
			'[0.90 - 0.95)' = alpha('forestgreen', 0.75),
			'[0.95 - 1]' = 'forestgreen'
		)

		
		map <- ggplot() +
			layer_spatial(this_ag_vect, aes(fill = quant_col), color = NA) +
			layer_spatial(nam, color = 'gray50', fill = NA, linewidth = 0.1) +
			scale_fill_manual(
				name = 'Quantile\n of λ',
				values = fill_scale
			) +
			# layer_spatial(cents_with_ag, pch = 16, col = 'black', size = 0.2) +
			guides(fill = guide_legend(reverse = TRUE)) +
			xlim(extent[1], extent[2]) + ylim(extent[3], extent[4]) +
			ggtitle(expression('Future equilibrial distribution of ' * italic('Andropogon gerardi')), subtitle = pretty_title) +
			theme(
				plot.title = element_text(size = 16),
				plot.subtitle = element_text(size = 14)
			)
				

		ggsave(plot = map, filename = paste0(out_dir, '/sdm_nmixture_lambda_', fut, '.png'), width = 12, height = 10, dpi = 600)

		writeVector(this_ag_vect, paste0(out_dir, '/sdm_nmixture_', fut, '.gpkg'), overwrite = TRUE)

	} # next future

say('###########################')
say('### parameter estimates ###')
say('###########################')

	chains <- readRDS(paste0(out_dir, '/sdm_nmixture_chains.rds'))

	# subset chain summary to just the lambdas
	summary <- chains$summary$all.chains

	model_term <- attr(terms(form), 'term.labels')
	model_term_nice <- model_term
	model_term_nice <- gsub(model_term_nice, pattern = 'I\\(', replacement = '')
	model_term_nice <- gsub(model_term_nice, pattern = '\\)', replacement = '')
	model_term_nice <- c('Intercept', model_term_nice)

	pars <- paste0('occ_beta[', 1:(1 + length(model_term)), ']')

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

	ggsave(caterpillars, filename = paste0(out_dir, '/sdm_nmixture_coefficients.png'), width = 10, height = 12, dpi = 300)

say('#######################')
say('### response curves ###')
say('#######################')

	chains <- readRDS(paste0(out_dir, '/sdm_nmixture_chains.rds'))

	# make a plot for how AG SDM and each ancestral population responds to each predictor
	responses <- list()
	for (i in 1:n_occ_predictors) {

		pred <- predictor_names[i]

		# unscaled predictor value
		x <- env_array[ , pred, pred]
		
		# create data frame with SDM response
		pars <- paste0('exp_response_curves_sdm_lambda[', 1:n_response_curve_rows, ', ', i, ']')
		sdm_response_mean <- chains$summary$all.chains[pars, 'Mean']
		sdm_response_lower <- chains$summary$all.chains[pars, '95%CI_low']
		sdm_response_upper <- chains$summary$all.chains[pars, '95%CI_upp']

		# data frames to hold predictions in long format
		df_mean <- data.frame(
			x = x,
			response = log(sdm_response_mean)
		)

		df_lower <- data.frame(
			x = x,
			response = log(sdm_response_lower)
		)

		df_upper <- data.frame(
			x = x,
			response = log(sdm_response_upper)
		)

		this_df_ci <- data.frame(
			x = c(x, rev(x)),
			y_occ = c(df_upper$response, rev(df_lower$response))
		)

		lambda_max <- max(-Inf, df_upper$response[!is.infinite(df_upper$response)])

		if (pred == 'aridity') {
			nice_title <- 'Aridity ((temperature + 10) / (precipitation / 1000))'
			nice_axis <- 'Aridity'
		} else if (pred == 'bio1') {
			nice_title <- 'Mean annual temperature (BIO01)'
			nice_axis <- 'Mean annual temperature (°C)'
		} else if (pred == 'bio5') {
			nice_title <- 'Temperature of the hottest month (BIO05)'
			nice_axis <- 'Temperature of the hottest month (°C)'
		} else if (pred == 'bio6') {
			nice_title <- 'Temperature of the coldest month (BIO06)'
			nice_axis <- 'Temperature of the coldest month (°C)'
		} else if (pred == 'bio7') {
			nice_title <- 'Temperature annual range (BIO07)'
			nice_axis <- 'Temperature annual range (°C)'
		} else if (pred == 'bio12') {
			nice_title <- 'Total annual precipitation (BIO12)'
			nice_axis <- 'Total annual precipitation (mm)'
		} else if (pred == 'bio15') {
			nice_title <- 'Precipitation seasonality (BIO15)'
			nice_axis <- 'Precipitation seasonality'
		} else if (pred == 'bio18') {
			nice_title <- 'Precipitation of warmest quarter (BIO18)'
			nice_axis <- 'Precipitation of warmest quarter (mm)'
		} else if (pred == 'ph') {
			nice_title <- 'Soil pH'
			nice_axis <- 'pH'
		} else if (pred == 'sand') {
			nice_title <- 'Soil proportion sand'
			nice_axis <- 'Proportion sand'
		} else if (pred == 'silt') {
			nice_title <- 'Soil proportion silt'
			nice_axis <- 'Proportion silt'
		} else if (pred == 'soc') {
			nice_title <- 'Soil organic matter'
			nice_axis <- 'Soil organic matter'
		}

		response <- ggplot() +
			geom_polygon(
				data = this_df_ci,
				mapping = aes(x = x, y_occ = y_occ),
				color = NA,
				fill = alpha('blue', 0.1)
			) +
			geom_line(
				data = df_mean,
				mapping = aes(x = x, y_occ = response)
			) +
			xlab(nice_axis) +
			ylab('Lambda') +
			ggtitle(nice_title) +
			theme(
				plot.title = element_text(size = 10)
			)

		responses[[pred]] <- response

	} # next predictor

	for (i in 1:n_occ_predictors) {
		responses[[i]] <- responses[[i]] + ylim(0, lambda_max)
	}

	if (length(predictor_names) == 1) {
		nrow <- 1
		width <- 10
		height <- 8
	} else if (length(predictor_names) <= 4) {
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

say(date())
say('FINIS!', deco = '+', level = 1)






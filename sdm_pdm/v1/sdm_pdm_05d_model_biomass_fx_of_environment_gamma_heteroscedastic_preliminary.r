### MODELING ANDROPOGON GERARDI DISTRIBUTION, PHENOTYPE, PHYSIOLOGY, GENOTYPE, and ASSOCIATED MICROBIAL COMMUNITIES
### Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2023-12
###
### This script constructs a distribution model of site-level mean biomass for Andropogon gerardi. It assumes that the mean site-level response is a function of climate (and perhaps also soil) variables, and that the standard deviation of biomass within a site is constant across sites.
###
### source('C:/Kaji/R/andropogon_integratedEcology/sdm_pdm/sdm_pdm_05d_model_biomass_fx_of_environment_gamma_heteroscedastic_preliminary.r')
### 
### CONTENTS ###
### setup ###
### model ###

#############
### setup ###
#############

	rm(list = ls())

	drive <- 'C:/Kaji/'

	setwd(paste0(drive, '/Research/Andropogon/Andropogon'))
	source(paste0(drive, '/R/andropogon_integratedEcology/sdm_pdm/sdm_pdm_00_shared_functions_and_variables.r'))

###########################
### user-defined values ###
###########################

	trial <- FALSE # TRUE for testing

	# # gamma and lognormal
	# formula_biomass <- ~ 1 + bio12
	# out_dir <- paste0('./outputs_loretta/integrated_sdm_pdm/biomass_alone_lognormal_homoscedastic_bio12', ifelse(trial, '_TRIAL', ''), '/')

	# # gamma and lognormal
	# formula_biomass <- ~ 1 + bio12 + I(bio12^2)
	# out_dir <- paste0('./outputs_loretta/integrated_sdm_pdm/biomass_alone_lognormal_homoscedastic_bio12^2', ifelse(trial, '_TRIAL', ''), '/')

	# # lognormal
	# formula_biomass <- ~ 1 + bio1 + bio12 + I(bio1^2)
	# out_dir <- paste0('./outputs_loretta/integrated_sdm_pdm/biomass_alone_lognormal_homoscedastic_bio1^2_bio12', ifelse(trial, '_TRIAL', ''), '/')

	# # gamma & lognormal
	# formula_biomass <- ~ 1 + bio1 + aridity + I(bio1^2)
	# out_dir <- paste0('./outputs_loretta/integrated_sdm_pdm/biomass_alone_lognormal_homoscedastic_bio1^2_aridity', ifelse(trial, '_TRIAL', ''), '/')

	# # gamma
	# formula_biomass <- ~ 1 + bio12 + aridity
	# out_dir <- paste0('./outputs_loretta/integrated_sdm_pdm/biomass_alone_lognormal_homoscedastic_bio12_aridity', ifelse(trial, '_TRIAL', ''), '/')

	dirCreate(out_dir)

	### MCMC settings
	niter <- 240000
	nburnin <- 40000
	thin <- 200
	nchains <- 4
	waic <- TRUE
	k_folds <- 5
	trial <- FALSE

	# # ### MCMC settings FOR TESTING
	# niter <- 11000
	# nburnin <- 1000
	# thin <- 10
	# nchains <- 2
	# waic <- FALSE
	# k_folds <- 2
	# trial <- TRUE

	# # ### MCMC settings FOR TESTING
	# niter <- 2200
	# nburnin <- 200
	# thin <- 2
	# nchains <- 2
	# # waic <- FALSE
	# waic <- TRUE
	# k_folds <- 2
	# trial <- TRUE

	# number of values in response curve array used to depict responses of SDM and PDM
	n_response_curve_values <- 200
	
	futs <- c(
		'ssp245_2041_2070',
		'ssp245_2071_2100',
		'ssp370_2041_2070',
		'ssp370_2071_2100'
	)

#############
### model ###
#############

	sink(paste0(out_dir, '/runtime_log.txt'), split = TRUE)
	say('MODELING BIOMASS ALONE')
	say('Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | ', date(), post = 1)

	### data collation
	##################
	say('data collation', post = 1)

	say('This script constructs a distribution model of site-level mean biomass for Andropogon gerardi. It assumes that the mean site-level response is a function of climate (and perhaps also soil) variables, and that the standard deviation of biomass within a site is also a function of climate/soil. We use a gamma distribution to model biomass at the site level as a linear function of climate, then use those coefficients to predict biomass at the county level.', breaks = 60, post = 1)

	say('MCMC settings:', level = 2)
	say('niter ...... .', niter)
	say('nburnin ..... ', nburnin)
	say('thin ........ ', thin)
	say('nchains ..... ', nchains, post = 1)

	##################################
	### model predictors and terms ###
	##################################

	say('Predictors and model formula:', level = 2)
	say('Phenotype formula for mean biomass: ', as.character(formula_biomass))

	terms <- terms(formula_biomass)
	terms <- attr(terms, 'term.labels')
	linear_terms <- terms
	linear_terms <- linear_terms[!grepl(linear_terms, pattern = '\\^2')]
	linear_terms <- linear_terms[!grepl(linear_terms, pattern = '\\:')]
	linear_terms <- linear_terms[!grepl(linear_terms, pattern = '\\*')]
	predictors_biomass <- linear_terms
	
	n_predictors_biomass <- length(predictors_biomass)

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
	
	# create spatial versions of site and phenotype data for plotting
	site_data_with_biomass <- merge(site_data_raw[ , c('site_id', 'LONGITUDE', 'LATITUDE')], biomass_data_raw, by.x = 'site_id', by.y = 'SITE')
	site_data_with_biomass_vect <- vect(site_data_with_biomass, geom = c('LONGITUDE', 'LATITUDE'), crs = getCRS('WGS84'))
	
	# center and scale at site level
	x <- site_data_raw[ , ..predictors_biomass]
	x <- scale(x)
	biomass_x_centers <- attr(x, 'scaled:center')
	biomass_x_scales <- attr(x, 'scaled:scale')
	x <- as.data.frame(x)
	x_by_site_biomass <- model.matrix(formula_biomass, x)

	# make vector of which sampled site matches each row in the biomass and morphology/physiology data
	n_biomass <- nrow(biomass_data_raw)

	site_index_biomass <- rep(NA, n_biomass)
	for (i in 1:n_biomass) {
		site_index_biomass[i] <- which(site_data_raw$site_id == biomass_data_raw$SITE[i])
	}
	
	### load occurrence data
	########################
	# We need this to make spatial predictions.
	say('load occurrences and present-day + future environmental data', level = 3)
	ag_vect_sq <- vect('./data_from_adam_and_loretta/andropogon_gerardi_occurrences_with_environment_1961_2020_climatena.gpkg')
	
	### future climates
	ag_vect_ssp245_2041_2070 <- vect('./data_from_adam_and_loretta/andropogon_gerardi_occurrences_with_environment_ensemble_8GCMs_ssp245_2041_2070_climatena.gpkg')
	ag_vect_ssp245_2071_2100 <- vect('./data_from_adam_and_loretta/andropogon_gerardi_occurrences_with_environment_ensemble_8GCMs_ssp245_2071_2100_climatena.gpkg')
	ag_vect_ssp370_2041_2070 <- vect('./data_from_adam_and_loretta/andropogon_gerardi_occurrences_with_environment_ensemble_8GCMs_ssp370_2041_2070_climatena.gpkg')
	ag_vect_ssp370_2071_2100 <- vect('./data_from_adam_and_loretta/andropogon_gerardi_occurrences_with_environment_ensemble_8GCMs_ssp370_2071_2100_climatena.gpkg')

	# fields <- c('area_km2', 'n_andropogon_gerardi', 'n_poaceae', predictors_biomass)
	# ag_vect_sq <- ag_vect_sq[ , fields]

	site_data_with_biomass_vect <- project(site_data_with_biomass_vect, ag_vect_sq)

	### collate county-level environmental data for spatial predictions
	###################################################################
	say('collate county-level environmental data', level = 3)

	ag_sq <- as.data.frame(ag_vect_sq)
	x_raw <- ag_sq[ , predictors_biomass, drop = FALSE]

	# record counties with infinite environmental values
	bads <- numeric() # indices of counties with infinite values for a predictor
	for (i in 1:ncol(x_raw)) {
		x <- x_raw[ , i, drop = TRUE]
		if (any(is.infinite(x))) {
			bads <- c(bads, which(is.infinite(x)))
		}
	}

	x_raw_scaled <- scale(x_raw, scale = biomass_x_scales, center = biomass_x_centers)
	x_raw_scaled <- as.data.frame(x_raw_scaled)
	pheno_x_counties_sq <- model.matrix(formula_biomass, x_raw_scaled)

	# futures
	for (ssp in c('ssp245_2041_2070', 'ssp245_2071_2100', 'ssp370_2041_2070', 'ssp370_2071_2100')) {

		# future: counties with future data
		x_fut <- get(paste0('ag_vect_', ssp))
		x_fut <- as.data.frame(x_fut)
		x_fut <- x_fut[ , predictors_biomass, drop = FALSE]

		# record counties with infinite environmental values
		for (i in 1:ncol(x_fut)) {
			x <- x_fut[ , i, drop = TRUE]
			if (any(is.infinite(x))) {
				bads <- c(bads, which(is.infinite(x)))
			}
		}

		# scale
		this_x_fut <- x_fut[ , predictors_biomass, drop = FALSE]
		x_fut_scaled <- scale(this_x_fut, center = biomass_x_centers, scale = biomass_x_scales)
		x_fut_scaled <- as.data.frame(x_fut_scaled)

		mm <- model.matrix(formula_biomass, x_fut_scaled)
		assign(paste0('pheno_x_counties_', ssp), mm)

	}

	### remove counties with infinite values for a predictor
	bads <- sort(unique(bads))
	if (length(bads) > 0) {
	
		ag_vect_sq <- ag_vect_sq[-bads]
		ag_vect_ssp245_2041_2070 <- ag_vect_ssp245_2041_2070[-bads]
		ag_vect_ssp245_2071_2100 <- ag_vect_ssp245_2071_2100[-bads]
		ag_vect_ssp370_2041_2070 <- ag_vect_ssp370_2041_2070[-bads]
		ag_vect_ssp370_2071_2100 <- ag_vect_ssp370_2071_2100[-bads]
	
		pheno_x_counties_sq <- pheno_x_counties_sq[-bads, , drop = FALSE]
		pheno_x_counties_ssp245_2041_2070 <- pheno_x_counties_ssp245_2041_2070[-bads, , drop = FALSE]
		pheno_x_counties_ssp245_2071_2100 <- pheno_x_counties_ssp245_2071_2100[-bads, , drop = FALSE]
		pheno_x_counties_ssp370_2041_2070 <- pheno_x_counties_ssp370_2041_2070[-bads, , drop = FALSE]
		pheno_x_counties_ssp370_2071_2100 <- pheno_x_counties_ssp370_2071_2100[-bads, , drop = FALSE]

	}

	### response array for biomass: response of mean
	################################################
	say('response array for biomass: response of mu', level = 3)

	response_curve_x_biomass <- create_response_curve_array(
		form = formula_biomass,
		centers = biomass_x_centers,
		scales = biomass_x_scales,
		ag_pres = ag_vect_sq
	)

#####################
### inputs for nimble
#####################

	say('Inputs:', level = 2)
	data <- list(
		y_biomass = biomass_data_raw$Biomass
	)

	n_counties <- nrow(ag_vect_sq)
	n_terms_biomass <- ncol(pheno_x_counties_sq)

	constants <- list(
		
		n_pheno_sites = n_pheno_sites, # number of phenotype sample sites
		n_biomass = n_biomass, # number of biomass records

		n_counties = n_counties,

		n_terms_biomass = n_terms_biomass, # number of *terms* in linear predictor of mean biomass, including intercept and non-linear terms

		n_predictors_biomass = n_predictors_biomass, # number of *predictors* in linear predictor of mean biomass

		x_by_site_biomass = x_by_site_biomass, # model matrix for mean site-level biomass at sampled sites
		site_index_biomass = site_index_biomass, # index of sampled site for each row in biomass data

		n_response_curve_values = n_response_curve_values,
		response_curve_x_biomass = response_curve_x_biomass,

		pheno_x_counties_sq = pheno_x_counties_sq,

		pheno_x_counties_ssp245_2071_2100 = pheno_x_counties_ssp245_2071_2100,
		pheno_x_counties_ssp245_2041_2070 = pheno_x_counties_ssp245_2041_2070,
		pheno_x_counties_ssp370_2041_2070 = pheno_x_counties_ssp370_2041_2070,
		pheno_x_counties_ssp370_2071_2100 = pheno_x_counties_ssp370_2071_2100

	)

	# use frequentist coefficients as initial values for MCMC: mu
	y_raw <- data$y_biomass
	x <- x_by_site_biomass
	y <- rep(NA_real_, n_pheno_sites)
	for (i in seq_along(unique(site_index_biomass))) {
		y_site <- y_raw[site_index_biomass == i]
		y[i] <- mean(y_site)
	}
	
	frequentist_model <- glm(
		y ~ x - 1,
		family = Gamma(link = 'log')
	)
	inits_beta_biomass <- coefficients(frequentist_model)

	pheno_biomass_sq_inits <- pheno_biomass_fut_inits <- rep(mean(biomass_data_raw$Biomass), n_counties)

	# Calculate standard deviation of biomass across rows grouped by site_id
	biomass_sd_by_site <- biomass_data_raw[, .(biomass_sd = sd(Biomass, na.rm = TRUE)), by = SITE]
	init_biomass_sigma <- mean(biomass_sd_by_site$biomass_sd)

	inits <- list(

		beta_biomass_mu_rate = inits_beta_biomass, # PDM coefficients
		beta_biomass_mu_shape = inits_beta_biomass, # PDM coefficients
		
		log_site_mu_hat_biomass = rep(mean(biomass_data_raw$Biomass), n_pheno_sites), # latent mean biomass at a site

		biomass_county_ssp245_2041_2070 = pheno_biomass_fut_inits, # future biomass in counties
		biomass_county_ssp245_2071_2100 = pheno_biomass_fut_inits,
		biomass_county_ssp370_2041_2070 = pheno_biomass_fut_inits,
		biomass_county_ssp370_2071_2100 = pheno_biomass_fut_inits

	)

	say('Data:')
	print(str(data))

	say('Constants:', pre = 1)
	print(str(constants))

	say('Initializations:', pre = 1)
	print(str(inits))

	### define model
	say('nimbleCode():', level = 2)
	say('For the biomass component, we assume a gamma distribution of values at the site level from which the individual-level biomasses are drawn.', post = 2)
	
	if (n_predictors_biomass > 1) {
	
		### model where formulae for mean biomass has >1 predictor: full model WITH posterior predictive nodes
		######################################################################################################
		model_code <- nimbleCode({
		
			# phenotype: priors for relationship of mean and variance to environment
			beta_biomass_mu_rate[1] ~ dnorm(0, sd = 10) # intercept... not regularized
			beta_biomass_mu_shape[1] ~ dnorm(0, sd = 10) # intercept... not regularized
			for (i in 2:n_terms_biomass) {
				# beta_biomass_mu_rate[i] ~ ddexp(0, rate = 1) # regularization toward 0
				beta_biomass_mu_rate[i] ~ dnorm(0, sd = 10) # broad prior
				beta_biomass_mu_shape[i] ~ dnorm(0, sd = 10) # broad prior
			}

			# parameters of biomass distribution are latent and functions of environment
			# individual plant biomasses are samples from the site-level distribution defined by the site-level distribution (next chunk after this one)
			for (i in 1:n_pheno_sites) {

				log_site_mu_hat_biomass[i] <- inprod(beta_biomass_mu_rate[1:n_terms_biomass], x_by_site_biomass[i, 1:n_terms_biomass])
				log_shape_biomass[i] <- inprod(beta_biomass_mu_shape[1:n_terms_biomass], x_by_site_biomass[i, 1:n_terms_biomass])
				
				site_mu_hat_biomass[i] <- exp(log_site_mu_hat_biomass[i])
				shape_biomass[i] <- exp(log_shape_biomass[i])
				
				rate_biomass[i] <- log_shape_biomass[i] / site_mu_hat_biomass[i]

				# # # # simulate mean site biomass for DHARMa residuals
				# # # site_mu_hat_biomass_sim[i] ~dgamma(shape_biomass[i], rate_biomass[i])
				

			}

			# likelihood of biomass of individual plants
			for (i in 1:n_biomass) {

				# likelihood
				y_biomass[i] ~ dgamma(shape = shape_biomass[site_index_biomass[i]], rate = rate_biomass[site_index_biomass[i]])

				# simulated values for unconditional DHARMa residuals
				y_biomass_sim[i] ~ dgamma(shape = shape_biomass[site_index_biomass[i]], rate = rate_biomass[site_index_biomass[i]])
		
			}
				
			# posterior samplers for predictions to counties in status quo and future
			for (i in 1:n_counties) {

				# biomass: status quo
				log_biomass_county_mu_sq[i] <- inprod(beta_biomass_mu_rate[1:n_terms_biomass], pheno_x_counties_sq[i, 1:n_terms_biomass])
				mu_biomass_county_sq[i] <- exp(log_biomass_county_mu_sq[i])

				# biomass: future
				log_biomass_county_ssp245_2041_2070[i] <- inprod(beta_biomass_mu_rate[1:n_terms_biomass], pheno_x_counties_ssp245_2041_2070[i, 1:n_terms_biomass])
				biomass_county_ssp245_2041_2070[i] <- exp(log_biomass_county_ssp245_2041_2070[i])
				
				log_biomass_county_ssp245_2071_2100[i] <- inprod(beta_biomass_mu_rate[1:n_terms_biomass], pheno_x_counties_ssp245_2071_2100[i, 1:n_terms_biomass])
				biomass_county_ssp245_2071_2100[i] <- exp(log_biomass_county_ssp245_2071_2100[i])
				
				log_biomass_county_ssp370_2041_2070[i] <- inprod(beta_biomass_mu_rate[1:n_terms_biomass], pheno_x_counties_ssp370_2041_2070[i, 1:n_terms_biomass])
				biomass_county_ssp370_2041_2070[i] <- exp(log_biomass_county_ssp370_2041_2070[i])

				log_biomass_county_ssp370_2071_2100[i] <- inprod(beta_biomass_mu_rate[1:n_terms_biomass], pheno_x_counties_ssp370_2071_2100[i, 1:n_terms_biomass])
				biomass_county_ssp370_2071_2100[i] <- exp(log_biomass_county_ssp370_2071_2100[i])

			}

			# posterior predictive sampler for response curves for biomass: site-level mean
			for (j in 1:n_predictors_biomass) {

				for (i in 1:n_response_curve_values) {
						
					log_pheno_response_curves_biomass[i, j] <-
						inprod(beta_biomass_mu_rate[1:n_terms_biomass], response_curve_x_biomass[i, 1:n_terms_biomass, j])

					pheno_response_curves_biomass[i, j] <- exp(log_pheno_response_curves_biomass[i, j])
					
				}

			}

		})

		### model where formulae for mean biomass has >1 predictor: full model WITHOUT posterior predictive nodes
		#########################################################################################################
		model_code_simple <- nimbleCode({
		
			# phenotype: priors for relationship of mean and variance to environment
			beta_biomass_mu_rate[1] ~ dnorm(0, sd = 10) # intercept... not regularized
			beta_biomass_mu_shape[1] ~ dnorm(0, sd = 10) # intercept... not regularized
			for (i in 2:n_terms_biomass) {
				# beta_biomass_mu_rate[i] ~ ddexp(0, rate = 1) # regularization toward 0
				beta_biomass_mu_rate[i] ~ dnorm(0, sd = 10) # broad prior
				beta_biomass_mu_shape[i] ~ dnorm(0, sd = 10) # broad prior
			}

			# parameters of biomass distribution are latent and functions of environment
			# individual plant biomasses are samples from the site-level distribution defined by the site-level distribution (next chunk after this one)
			for (i in 1:n_pheno_sites) {

				log_site_mu_hat_biomass[i] <- inprod(beta_biomass_mu_rate[1:n_terms_biomass], x_by_site_biomass[i, 1:n_terms_biomass])
				log_shape_biomass[i] <- inprod(beta_biomass_mu_shape[1:n_terms_biomass], x_by_site_biomass[i, 1:n_terms_biomass])
				
				site_mu_hat_biomass[i] <- exp(log_site_mu_hat_biomass[i])
				shape_biomass[i] <- exp(log_shape_biomass[i])

				rate_biomass[i] <- shape_biomass[i] / site_mu_hat_biomass[i]

			}

			# likelihood of biomass of individual plants
			for (i in 1:n_biomass) {

				# likelihood
				y_biomass[i] ~ dgamma(shape = shape_biomass[site_index_biomass[i]], rate = rate_biomass[site_index_biomass[i]])
		
			}
				
		})

	} else if (n_predictors_biomass == 1) {
	
		### model where formulae for mean biomass has just 1 predictor: full model WITH posterior predictive nodes
		##########################################################################################################
		model_code <- nimbleCode({
		
			# phenotype: priors for relationship of mean and variance to environment
			beta_biomass_mu_rate[1] ~ dnorm(0, sd = 10) # intercept... not regularized
			beta_biomass_mu_shape[1] ~ dnorm(0, sd = 10) # intercept... not regularized
			for (i in 2:n_terms_biomass) {
				# beta_biomass_mu_rate[i] ~ ddexp(0, rate = 1) # regularization toward 0
				beta_biomass_mu_rate[i] ~ dnorm(0, sd = 10) # broad prior
				beta_biomass_mu_shape[i] ~ dnorm(0, sd = 10) # broad prior
			}

			# parameters of biomass distribution are latent and functions of environment
			# individual plant biomasses are samples from the site-level distribution defined by the site-level distribution (next chunk after this one)
			for (i in 1:n_pheno_sites) {

				log_site_mu_hat_biomass[i] <- inprod(beta_biomass_mu_rate[1:n_terms_biomass], x_by_site_biomass[i, 1:n_terms_biomass])
				log_shape_biomass[i] <- inprod(beta_biomass_mu_shape[1:n_terms_biomass], x_by_site_biomass[i, 1:n_terms_biomass])
				
				site_mu_hat_biomass[i] <- exp(log_site_mu_hat_biomass[i])
				shape_biomass[i] <- exp(log_shape_biomass[i])
				
				rate_biomass[i] <- shape_biomass[i] / site_mu_hat_biomass[i]

				# # # # simulate mean site biomass for DHARMa residuals
				# # # site_mu_hat_biomass_sim[i] ~dgamma(shape_biomass[i], rate_biomass[i])
				

			}

			# likelihood of biomass of individual plants
			for (i in 1:n_biomass) {

				# likelihood
				y_biomass[i] ~ dgamma(shape = shape_biomass[site_index_biomass[i]], rate = rate_biomass[site_index_biomass[i]])

				# simulated values for unconditional DHARMa residuals
				y_biomass_sim[i] ~ dgamma(shape = shape_biomass[site_index_biomass[i]], rate = rate_biomass[site_index_biomass[i]])
		
			}
				
			# posterior samplers for predictions to counties in status quo and future
			for (i in 1:n_counties) {

				# biomass: status quo
				log_biomass_county_mu_sq[i] <- inprod(beta_biomass_mu_rate[1:n_terms_biomass], pheno_x_counties_sq[i, 1:n_terms_biomass])
				mu_biomass_county_sq[i] <- exp(log_biomass_county_mu_sq[i])

				# biomass: future
				log_biomass_county_ssp245_2041_2070[i] <- inprod(beta_biomass_mu_rate[1:n_terms_biomass], pheno_x_counties_ssp245_2041_2070[i, 1:n_terms_biomass])
				biomass_county_ssp245_2041_2070[i] <- exp(log_biomass_county_ssp245_2041_2070[i])
				
				log_biomass_county_ssp245_2071_2100[i] <- inprod(beta_biomass_mu_rate[1:n_terms_biomass], pheno_x_counties_ssp245_2071_2100[i, 1:n_terms_biomass])
				biomass_county_ssp245_2071_2100[i] <- exp(log_biomass_county_ssp245_2071_2100[i])
				
				log_biomass_county_ssp370_2041_2070[i] <- inprod(beta_biomass_mu_rate[1:n_terms_biomass], pheno_x_counties_ssp370_2041_2070[i, 1:n_terms_biomass])
				biomass_county_ssp370_2041_2070[i] <- exp(log_biomass_county_ssp370_2041_2070[i])

				log_biomass_county_ssp370_2071_2100[i] <- inprod(beta_biomass_mu_rate[1:n_terms_biomass], pheno_x_counties_ssp370_2071_2100[i, 1:n_terms_biomass])
				biomass_county_ssp370_2071_2100[i] <- exp(log_biomass_county_ssp370_2071_2100[i])

			}

			# posterior predictive sampler for response curves for biomass: site-level mean
			for (i in 1:n_response_curve_values) {
					
				log_pheno_response_curves_biomass[i] <-
					inprod(beta_biomass_mu_rate[1:n_terms_biomass], response_curve_x_biomass[i, 1:n_terms_biomass])
				
				pheno_response_curves_biomass[i] <- exp(log_pheno_response_curves_biomass[i])
				
			}

		})

		### model where formulae for mean biomass has just 1 predictor: full model WITHOUT posterior predictive nodes
		##############################################################################################################
		model_code_simple <- nimbleCode({
		
			# phenotype: priors for relationship of mean and variance to environment
			beta_biomass_mu_rate[1] ~ dnorm(0, sd = 10) # intercept... not regularized
			beta_biomass_mu_shape[1] ~ dnorm(0, sd = 10) # intercept... not regularized
			for (i in 2:n_terms_biomass) {
				# beta_biomass_mu_rate[i] ~ ddexp(0, rate = 1) # regularization toward 0
				beta_biomass_mu_rate[i] ~ dnorm(0, sd = 10) # broad prior
				beta_biomass_mu_shape[i] ~ dnorm(0, sd = 10) # broad prior
			}

			# parameters of biomass distribution are latent and functions of environment
			# individual plant biomasses are samples from the site-level distribution defined by the site-level distribution (next chunk after this one)
			for (i in 1:n_pheno_sites) {

				log_site_mu_hat_biomass[i] <- inprod(beta_biomass_mu_rate[1:n_terms_biomass], x_by_site_biomass[i, 1:n_terms_biomass])
				log_shape_biomass[i] <- inprod(beta_biomass_mu_shape[1:n_terms_biomass], x_by_site_biomass[i, 1:n_terms_biomass])

				site_mu_hat_biomass[i] <- exp(log_site_mu_hat_biomass[i])
				shape_biomass[i] <- exp(log_shape_biomass[i])

				rate_biomass[i] <- shape_biomass[i] / site_mu_hat_biomass[i]


			}

			# likelihood of biomass of individual plants
			for (i in 1:n_biomass) {

				# likelihood
				y_biomass[i] ~ dgamma(shape = shape_biomass[site_index_biomass[i]], rate = rate_biomass[site_index_biomass[i]])

			}

		})

	}

	print(model_code)

	say('nimbleModel():', level = 2)
	model <- nimbleModel(
		code = model_code, # our model
		constants = constants, # constants
		data = data, # data
		inits = inits, # initialization values
		check = TRUE, # any errors?
		calculate = FALSE
		# buildDerivs = TRUE # need for Hamiltonian Monte Carlo
	)

	model_simple <- nimbleModel(
		code = model_code_simple, # our model
		constants = constants, # constants
		data = data, # data
		inits = inits, # initialization values
		check = TRUE, # any errors?
		calculate = FALSE
		# buildDerivs = TRUE # need for Hamiltonian Monte Carlo
	)

	say('$initializeInfo() and $calculate():', level = 2)
	model$initializeInfo()
	calc <- model$calculate()
	say('model$calculate(): ', calc)

	say('configureMCMC():', level = 2)
	
	monitors_simple <- c('beta_biomass_mu_rate', 'beta_biomass_mu_shape')
	monitors_coeffs <- c('beta_biomass_mu_rate', 'beta_biomass_mu_shape', 'site_mu_hat_biomass')
	monitors <- c(
		monitors_coeffs,
		'y_biomass_sim', # 'site_mu_hat_biomass_sim', # for residuals analysis
		'pheno_response_curves_biomass',
		'mu_biomass_county_sq',
		'biomass_county_ssp245_2041_2070', 'biomass_county_ssp245_2071_2100', 'biomass_county_ssp370_2041_2070', 'biomass_county_ssp370_2071_2100'
	)
	
	conf <- configureMCMC(
		model,
		monitors = monitors,
		print = TRUE,
		enableWAIC = waic
	)

	conf_simple <- configureMCMC(
		model_simple,
		monitors = monitors_simple,
		print = FALSE,
		enableWAIC = FALSE
	)

	# # add no U-turn sampler (Hamiltonian Monte Carlo)
	# conf$addSampler(target = c('occ_beta', 'occ_alpha_area', 'occ_alpha_poaceae'), type = 'NUTS')

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
		WAIC = waic,
		perChainWAIC = FALSE
	)

	saveRDS(chains, paste0(out_dir, '/chains.rds'))

	say('session info', level = 2)
	print(sessionInfo())

	say(date(), pre = 1)
	sink()

say('#########################')
say('### model diagnostics ###')
say('#########################')

	# chains <- readRDS(paste0(out_dir, '/chains.rds'))
	mcmc <- chains$samples
	cols <- c(
		paste0('beta_biomass_mu_rate[', 1:n_terms_biomass, ']'),
		paste0('beta_biomass_mu_shape[', 1:n_terms_biomass, ']')
	)
	for (i in 1:nchains) {
		mcmc[[i]] <- mcmc[[i]][ , cols]
	}

	ggs_mcmc <- ggs(mcmc)

	### trace and density plots
	###########################

	# graphing trace and density plots for all betas
	pars <- paste0('beta_biomass_mu_rate')
	file <- paste0(out_dir, '/beta_rate_biomass_mu_trace.png')
	ggsave(ggs_traceplot(ggs_mcmc, family = pars), file = file, width = 10, height = 8, dpi = 450, bg = 'white')

	file <- paste0(out_dir, '/beta_rate_biomass_mu_density.png')
	ggsave(ggs_density(ggs_mcmc, family = pars, rug = TRUE, hpd = TRUE), file = file, width = 6, height = 8, dpi = 450, bg = 'white')

	# graphing trace and density plots for all betas
	pars <- paste0('beta_biomass_mu_shape')
	file <- paste0(out_dir, '/beta_shape_biomass_mu_trace.png')
	ggsave(ggs_traceplot(ggs_mcmc, family = pars), file = file, width = 10, height = 8, dpi = 450, bg = 'white')

	file <- paste0(out_dir, '/beta_shape_biomass_mu_density.png')
	ggsave(ggs_density(ggs_mcmc, family = pars, rug = TRUE, hpd = TRUE), file = file, width = 6, height = 8, dpi = 450, bg = 'white')

	# correlations between samples of beta
	pars <- c('beta_biomass_mu_shape')
	file <- paste0(out_dir, '/beta_shape_biomass_correlations.png')
	ggsave(ggs_crosscorrelation(ggs_mcmc, family = pars), file = file, width = 8, height = 8, dpi = 450, bg = 'white')

	pars <- c('beta_biomass_mu_rate')
	file <- paste0(out_dir, '/beta_rate_biomass_correlations.png')
	ggsave(ggs_crosscorrelation(ggs_mcmc, family = pars), file = file, width = 8, height = 8, dpi = 450, bg = 'white')

	### DHARMa residuals
	####################

	### residuals by plant
	# NB this uses the mean predicted value of a site as a plant-level prediction
	sims_by_plant <- hammer_subset(chains, param = 'y_biomass_sim', j = TRUE)
	sims_by_plant <- hammer_rbind(sims_by_plant)
	sims_by_plant <- t(sims_by_plant)

	estimates <- hammer_subset(chains, param = 'site_mu_hat_biomass', j = TRUE)
	estimates <- hammer_rbind(estimates)

	site_counts <- biomass_data_raw[ , .N, by = SITE]
	fit_by_site <- apply(estimates, 2, median)
	fit_replicated <- numeric()
	for (i in seq_along(fit_by_site)) {
		fit_replicated <- c(fit_replicated, rep(fit_by_site[i], site_counts$N[i]))
	}

	dharma <- createDHARMa(simulatedResponse = sims_by_plant, observedResponse = data$y_biomass, fittedPredictedResponse = fit_replicated, integerResponse = FALSE)

	file <- paste0(out_dir, '/dharma_by_plant_y_biomass.png')
	png(file, width = 1200, height = 800)
		plot(dharma)
	dev.off()
	
	file <- paste0(out_dir, '/dharma_by_plant_y_biomass_residuals.png')
	png(file, width = 1200, height = 800)
		hist(dharma$scaledResiduals, main = 'DHARMa residuals for biomass by plant', xlab = 'Scaled residuals', breaks = 30)
	dev.off()

	# # # ### residuals by site
	# # # # NB this uses the mean predicted and observed values sites
	# # # sims_by_plant <- hammer_subset(chains, param = 'y_biomass_sim', j = TRUE)
	# # # sims_by_plant <- hammer_rbind(sims_by_plant)
	# # # sims_by_site <- rep(NA, n_pheno_sites)
	# # # for (i in seq_along(unique(site_index_biomass))) {
	# # # 	sims_by_site[i] <- mean(sims_by_plant[ , site_index_biomass == i])
	# # # }

	# # # sims_by_site <- t(sims_by_site)
	# # # rownames(sims_by_site) <- site_data_raw$site_id

	# # # estimates <- hammer_subset(chains, param = 'site_mu_hat_biomass', j = TRUE)
	# # # estimates <- hammer_rbind(estimates)
	# # # site_counts <- biomass_data_raw[ , .N, by = SITE]
	# # # fit_by_site <- apply(estimates, 2, median)

	# # # # Calculate mean biomass grouped by site
	# # # mean_biomass_by_site <- biomass_data_raw[ , .(mean_biomass = mean(Biomass, na.rm = TRUE)), by = SITE]

	# # # dharma <- createDHARMa(simulatedResponse = sims_by_site, observedResponse = mean_biomass_by_site$mean_biomass, fittedPredictedResponse = fit_by_site, integerResponse = FALSE)

	# # # file <- paste0(out_dir, '/dharma_by_site_y_biomass.png')
	# # # png(file, width = 1200, height = 800)
	# # # 	plot(dharma)
	# # # dev.off()
	
	# # # file <- paste0(out_dir, '/dharma_by_site_y_biomass_residuals.png')
	# # # png(file, width = 1200, height = 800)
	# # # 	hist(dharma$scaledResiduals, main = 'DHARMa residuals for biomass by site', xlab = 'Scaled residuals', breaks = 20)
	# # # dev.off()

	### residuals vs covariates
	terms <- terms(formula_biomass)
	terms <- attr(terms, 'term.labels')
	linear_terms <- terms
	linear_terms <- linear_terms[!grepl(linear_terms, pattern = '\\^2')]
	linear_terms <- linear_terms[!grepl(linear_terms, pattern = '\\:')]
	linear_terms <- linear_terms[!grepl(linear_terms, pattern = '\\*')]
	n_predictors <- length(linear_terms)

	resids <- dharma$scaledResiduals
	resids_vs_covariates <- list()
	site_counts <- biomass_data_raw[ , .N, by = SITE]
	for (i in seq_along(linear_terms)) {

		this_x <- data.frame(
			x = biomass_data_raw[[linear_terms[i]]],
			y = resids
		)

		resids_vs_covariates[[i]] <- ggplot(this_x, aes(x = x, y = y)) +
			geom_point() +
			geom_smooth(method = loess, formula = y ~ x, se = FALSE) +
			xlab(linear_terms[i]) +
			ylab('Scaled residuals') +
			ggtitle('Biomass-only model')

	}

	if (n_predictors == 1) {
		nrow <- 1
		width <- 10
		height <- 8
	} else if (n_predictors <= 4) {
		nrow <- 1
		width <- 14
		height <- 4
	} else {
		nrow <- 2
		width <- 16
		height <- 10	
	}

	resids_vs_covariates <- plot_grid(plotlist = resids_vs_covariates, nrow = nrow)
	ggsave(plot = resids_vs_covariates, filename = paste0(out_dir, '/dharma_residuals_vs_covariates.png'), width = width, height = height, dpi = 600)

	### Gelman-Rubin convergence statistic
	######################################
	mcmc_coeffs <- chains$samples
	cols <- c(
		paste0('beta_biomass_mu_rate[', 1:n_terms_biomass, ']'),
		paste0('beta_biomass_mu_shape[', 1:n_terms_biomass, ']')
	)
	for (i in 1:nchains) {
		mcmc_coeffs[[i]] <- mcmc_coeffs[[i]][ , cols]
	}

	rhats <- tryCatch(
		gelman.diag(mcmc_coeffs, autoburnin = FALSE, multivariate = TRUE),
		error = function(cond) FALSE
	)
	
	sink(paste0(out_dir, '/convergence.txt'), split = TRUE)
	say('GELMAN-RUBIN STATISTICS')
	say(date(), post = 2)
	print(rhats)
	sink()

	### WAIC
	########
	sink(paste0(out_dir, '/waic.txt'), split = TRUE)
	say('WAIC')
	say(date(), post = 2)
	print(chains$WAIC)
	sink()

say('########################')
say('### cross validation ###')
say('########################')

	cv <- runCrossValidate(
		MCMCconfiguration = conf_simple,
		k = k_folds,
		foldFunction = 'random',
		lossFunction = 'MSE',
		MCMCcontrol = list(niter = niter, nburnin = nburnin),
		returnSamples = FALSE,
		nCores = 1,
		nBootReps = 200,
		silent = FALSE
	)

	sink(paste0(out_dir, '/cross_validation.txt'), split = TRUE)
	say('CROSS VALIDATION', post = 2)
	print(cv)
	sink()

say('###########################')
say('### parameter estimates ###')
say('###########################')

	# chains <- readRDS(paste0(out_dir, '/chains.rds'))
	mcmc <- chains$samples
	cols <- c(
		paste0('beta_biomass_mu_rate[', 1:n_terms_biomass, ']'),
		paste0('beta_biomass_mu_shape[', 1:n_terms_biomass, ']'),
		paste0('site_mu_hat_biomass[', 1:n_pheno_sites, ']')
	)
	for (i in 1:nchains) {
		mcmc[[i]] <- mcmc[[i]][ , cols]
	}

	ggs_mcmc <- ggs(mcmc)

	caterpillars <- list()
	for (i in seq_along(monitors_coeffs)) {

		monitors_coeff <- monitors_coeffs[i]
		caterpillars[[i]] <- ggs_caterpillar(ggs_mcmc, family = monitors_coeff) +
			xlab('Estimated value') +
			ggtitle('Biomass-only model') +
			theme(
				axis.title.y = element_blank()
			)

	
	}

	combo <- plot_grid(plotlist = caterpillars, ncol = 2, align = 'v', axis = 'l')
	ggsave(combo, filename = paste0(out_dir, '/coefficients_biomass.png'), width = 12, height = 8, dpi = 300)

say('#######################')
say('### response curves ###')
say('#######################')

	# chains <- readRDS(paste0(out_dir, '/chains.rds'))

	n_predictors <- n_predictors_biomass
	these_predictors_biomass <- predictors_biomass
	y_lab <- bquote('Biomass ' * ' (g)')

	# make a plot for how biomass mu and sigma respond to each predictor
	responses <- list()
	max_val <- -Inf
	for (i in 1:n_predictors) {

		pred <- these_predictors_biomass[i]

		if (pred == 'aridity') {
			nice_title <- 'Aridity ((temp. + 10) / (precip. / 1000))'
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

		nice_title <- bquote('Biomass ' * ' versus ' * .(nice_title))

		# unscale predictor value
		if (n_predictors_biomass == 1) {
			x <- response_curve_x_biomass[ , pred]
		} else {
			x <- response_curve_x_biomass[ , pred, pred]
		}
		x <- x * biomass_x_scales[[pred]] + biomass_x_centers[[pred]]

		# create data frame with SDM response
		if (n_predictors_biomass == 1) {
			pars <- paste0('pheno_response_curves_biomass[', 1:n_response_curve_values, ']')
		} else {
			pars <- paste0('pheno_response_curves_biomass[', 1:n_response_curve_values, ', ', i, ']')
		}
		response_mean <- chains$summary$all.chains[pars, 'Mean']
		response_lower <- chains$summary$all.chains[pars, '95%CI_low']
		response_upper <- chains$summary$all.chains[pars, '95%CI_upp']

		# data frames to hold predictions in long format
		df_mean <- data.frame(
			x = x,
			response = response_mean
		)

		df_lower <- data.frame(
			x = x,
			response = response_lower
		)

		df_upper <- data.frame(
			x = x,
			response = response_upper
		)

		this_df_ci <- data.frame(
			x = c(x, rev(x)),
			response = c(df_upper$response, rev(df_lower$response))
		)

		response <- ggplot() +
			geom_polygon(
				data = this_df_ci,
				mapping = aes(x = x, y = response),
				color = NA,
				fill = alpha('blue', 0.1)
			) +
			geom_line(
				data = df_mean,
				mapping = aes(x = x, y = response)
			) +
			xlab(nice_axis) +
			ylab(y_lab) +
			ggtitle(nice_title) +
			theme(
				legend.position = 'none',
				plot.title = element_text(size = 10)
			)

			site_data_with_biomass_TEMP <- site_data_with_biomass
			site_data_with_biomass_TEMP$x <- site_data_with_biomass_TEMP[[pred]]
			site_data_with_biomass_TEMP$y <- site_data_with_biomass_TEMP$Biomass

			max_val <- max(
				max_val,
				quantile(df_upper$response[!is.infinite(df_upper$response)], 0.2),
				site_data_with_biomass_TEMP$y
			)

			response <- response + geom_point(
				data = site_data_with_biomass_TEMP,
				mapping = aes(x = x, y = y, color = site_id)
			)

		responses[[pred]] <- response

	} # next predictor

	for (i in 1:n_predictors) {
		responses[[i]] <- responses[[i]] + coord_cartesian(ylim = c(0, max_val)) 
	}

	if (n_predictors == 1) {
		nrow <- 1
		width <- 10
		height <- 8
	} else if (n_predictors <= 4) {
		nrow <- 1
		width <- 14
		height <- 4
	} else {
		nrow <- 2
		width <- 16
		height <- 10	
	}

	responses <- plot_grid(plotlist = responses, nrow = nrow)
	ggsave(plot = responses, filename = paste0(out_dir, '/response_curves_biomass.png'), width = width, height = height, dpi = 600)

say('########################', pre = 1)
say('### map of residuals ###')
say('########################')

	# Make map of standardized residuals (predicted - observed) / (observed_mean)

	# # user-defined
	# chains <- readRDS(paste0(out_dir, '/chains.rds'))

	# calculate residuals
	resp <- hammer_subset(chains, param = 'site_mu_hat_biomass', j = TRUE)
	resp <- hammer_extract_summary(resp)

	mean_biomass_by_site <- biomass_data_raw[ , .(mean_biomass = mean(Biomass, na.rm = TRUE)), by = SITE]
	stand_resid <- (mean_biomass_by_site$mean_biomass - resp$all.chains[ , 'Mean']) / mean_biomass_by_site$mean_biomass

	site_data_with_biomass_vect_resid <- site_data_with_biomass_vect
	site_data_with_biomass_vect_resid$stand_resid <- stand_resid

	# administrative boundaries
	nam <- vect('./data_from_gadm/gadm_4pt1_level_1_north_america_sans_alaska_lambert.gpkg')

	# extent
	ag_vect_sq_pres <- ag_vect_sq[ag_vect_sq$n_andropogon_gerardi > 0]
	extent <- ext(ag_vect_sq_pres)
	extent <- as.vector(extent)
	x_range <- (extent[2] - extent[1])
	y_range <- (extent[4] - extent[3])
	extent[1] <- extent[1] + 0.15 * x_range
	extent[3] <- extent[3] + 0.125 * y_range
	extent[4] <- extent[4] - 0.2 * y_range

	legend_title <- 'Biomass\nmean-standardized\nresidual'

	map <- ggplot() +
		layer_spatial(nam, color = 'gray30', fill = 'white', linewidth = 0.3) +
		layer_spatial(site_data_with_biomass_vect_resid, aes(fill = stand_resid), pch = 21, size = 6) +
		scale_fill_gradient2(
			name = legend_title,
			low = 'red',
			mid = 'beige',
			high = 'blue',
			midpoint = 0
		) +
		xlim(extent[1], extent[2]) + ylim(extent[3], extent[4]) +
		ggtitle(
			bquote('Residuals for ' * italic('Andropogon gerardi') * ' biomass '),
			subtitle = '1991-2020 | biomass-only model') +
		theme(
			plot.title = element_text(size = 16),
			plot.subtitle = element_text(size = 14)
		)

	ggsave(plot = map, filename = paste0(out_dir, '/biomass_residuals_by_site_map.png'), width = 12, height = 10, dpi = 600)

say('##############################', pre = 1)
say('### map of current biomass ###')
say('##############################')

	# # user-defined
	# chains <- readRDS(paste0(out_dir, '/chains.rds'))

	# subset chain summary to just the lambdas
	summary <- chains$summary$all.chains

	nam <- vect('./data_from_gadm/gadm_4pt1_level_1_north_america_sans_alaska_lambert.gpkg')

	# get range of values for plotting
	max_val <- max(site_data_with_biomass_vect$Biomass)
	ag_vect_sq_pres <- ag_vect_sq[ag_vect_sq$n_andropogon_gerardi > 0]
	for (fut in futs) {

		this_ag_vect <- get(paste0('ag_vect_', fut))
		this_ag_vect$num <- 1:nrow(this_ag_vect)
		matches <- this_ag_vect$stateProvince %in% ag_vect_sq_pres$state_province & this_ag_vect$county %in% ag_vect_sq_pres$county
		this_ag_vect <- this_ag_vect[matches]	

		which_resp <- grepl(rownames(summary), pattern = paste0('biomass_county_', fut))
		resp <- summary[which_resp, 'Mean']
		resp <- resp[this_ag_vect$num] # reorder to match this_ag_vect

		# get just counties with data
		max_val <- max(c(max_val, resp + 0.001))

	}

	resp_limits <- c(1, max_val)

	# just counties with data
	which_resp <- grepl(rownames(summary), pattern = paste0('mu_biomass_county_sq'))
	resp <- summary[which_resp, ]

	ag_vect_sq$resp_mean <- resp[ , 'Mean']
	ag_vect_sq$resp_0.05ci <- resp[ , '95%CI_low']
	ag_vect_sq$resp_0.95ci <- resp[ , '95%CI_upp']
	ag_vect_sq$resp_ci <- ag_vect_sq$resp_0.95ci - ag_vect_sq$resp_0.05ci

	# extent
	ag_vect_sq_pres <- ag_vect_sq[ag_vect_sq$n_andropogon_gerardi > 0]
	extent <- ext(ag_vect_sq_pres)
	extent <- as.vector(extent)
	x_range <- (extent[2] - extent[1])
	y_range <- (extent[4] - extent[3])
	extent[1] <- extent[1] + 0.15 * x_range
	extent[3] <- extent[3] + 0.125 * y_range
	extent[4] <- extent[4] - 0.2 * y_range

	cents_with_ag <- ag_vect_sq[ag_vect_sq$n_andropogon_gerardi > 0]
	cents_with_ag <- centroids(cents_with_ag)

	legend_title <- 'Biomass (g)'

	map <- ggplot() +
		layer_spatial(ag_vect_sq, aes(fill = resp_mean), color = NA) +
		layer_spatial(nam, color = 'gray30', fill = NA, linewidth = 0.3) +
		layer_spatial(cents_with_ag, pch = 16, color = alpha('gray20', 0.5), size = 0.35) +
		layer_spatial(site_data_with_biomass_vect, aes(fill = Biomass), pch = 21, size = 4) +
		scale_fill_continuous(
			name = legend_title,
			type = 'viridis',
			limits = resp_limits,
			trans = 'log10'
		) +
		xlim(extent[1], extent[2]) + ylim(extent[3], extent[4]) +
		ggtitle(
			bquote('Present-day distribution of ' * italic('Andropogon gerardi') * ' biomass '),
			subtitle = '1991-2020 | biomass-only model') +
		theme(
			plot.title = element_text(size = 16),
			plot.subtitle = element_text(size = 14)
		)

	names(ag_vect_sq)[names(ag_vect_sq) == 'resp_mean'] <- paste0('biomass_mean')
	names(ag_vect_sq)[names(ag_vect_sq) == 'resp_0.05ci'] <- paste0('biomass_0.05ci')
	names(ag_vect_sq)[names(ag_vect_sq) == 'resp_0.95ci'] <- paste0('biomass_0.95ci')
	names(ag_vect_sq)[names(ag_vect_sq) == 'resp_ci'] <- paste0('biomass_ci')

	ggsave(plot = map, filename = paste0(out_dir, '/biomass_status_quo.png'), width = 12, height = 10, dpi = 600)

	writeVector(ag_vect_sq, paste0(out_dir, '/biomass_status_quo.gpkg'), overwrite = TRUE)

say('##############################')
say('### maps of future biomass ###')
say('##############################')

	# user-defined

	# chains <- readRDS(paste0(out_dir, '/chains.rds'))

	nam <- vect('./data_from_gadm/gadm_4pt1_level_1_north_america_sans_alaska_lambert.gpkg')

	# get range of values for plotting
	max_val <- max(site_data_with_biomass_vect$Biomass)
	ag_vect_sq_pres <- ag_vect_sq[ag_vect_sq$n_andropogon_gerardi > 0]
	for (fut in futs) {

		this_ag_vect <- get(paste0('ag_vect_', fut))
		this_ag_vect$num <- 1:nrow(this_ag_vect)
		matches <- this_ag_vect$stateProvince %in% ag_vect_sq_pres$state_province & this_ag_vect$county %in% ag_vect_sq_pres$county
		this_ag_vect <- this_ag_vect[matches]	

		which_resp <- grepl(rownames(summary), pattern = paste0('biomass_county_', fut))
		resp <- summary[which_resp, 'Mean']
		resp <- resp[this_ag_vect$num] # reorder to match this_ag_vect

		# get just counties with data
		max_val <- max(c(max_val, resp + 0.001))

	}

	resp_limits <- c(1, max_val)

	for (fut in futs) {

		say(fut)

		which_resp <- grepl(rownames(summary), pattern = paste0('biomass_county_', fut))
		resp <- summary[which_resp, ]

		this_ag_vect <- get(paste0('ag_vect_', fut))

		# get just counties with data
		this_ag_vect$resp_mean <- resp[ , 'Mean']
		this_ag_vect$resp_0.05ci <- resp[ , '95%CI_low']
		this_ag_vect$resp_0.95ci <- resp[ , '95%CI_upp']
		this_ag_vect$resp_ci <- this_ag_vect$resp_0.95ci - this_ag_vect$resp_0.05ci

		pretty_title <- paste0('SSP ', substr(fut, 4, 6), ' ', substr(fut, 8, 11), '-', substr(fut, 13, 16), ' | biomass-only model')

		# extent
		ag_vect_sq_pres <- ag_vect_sq[ag_vect_sq$n_andropogon_gerardi > 0]
		extent <- ext(ag_vect_sq_pres)
		extent <- as.vector(extent)
		x_range <- (extent[2] - extent[1])
		y_range <- (extent[4] - extent[3])
		extent[1] <- extent[1] + 0.15 * x_range
		extent[3] <- extent[3] + 0.125 * y_range
		extent[4] <- extent[4] - 0.2 * y_range

		legend_title <- 'Biomass (g)'
		# resp_limits <- range(c(this_ag_vect$resp_mean, site_data_with_biomass_vect$Biomass))
	
		# Calculate standard deviation of biomass across rows grouped by site_id

		map <- ggplot() +
			layer_spatial(this_ag_vect, aes(fill = resp_mean), color = NA) +
			layer_spatial(nam, color = 'gray30', fill = NA, linewidth = 0.3) +
			layer_spatial(cents_with_ag, pch = 16, color = alpha('gray20', 0.5), size = 0.35) +
			layer_spatial(site_data_with_biomass_vect, aes(fill = Biomass), pch = 21, size = 4) +
			scale_fill_continuous(
				name = legend_title,
				# low = 'yellow',
				# high = 'darkorange4',
				type = 'viridis',
				limits = resp_limits,
				trans = 'log10'
			) +
			xlim(extent[1], extent[2]) + ylim(extent[3], extent[4]) +
				ggtitle(
					bquote('Future equilibrial distribution of ' * italic('Andropogon gerardi') * ' biomass '),
					subtitle = pretty_title
				) +
			theme(
				plot.title = element_text(size = 16),
				plot.subtitle = element_text(size = 14)
			)

		names(ag_vect_sq)[names(ag_vect_sq) == 'resp_mean'] <- paste0('biomass_mean')
		names(ag_vect_sq)[names(ag_vect_sq) == 'resp_0.05ci'] <- paste0('biomass_0.05ci')
		names(ag_vect_sq)[names(ag_vect_sq) == 'resp_0.95ci'] <- paste0('biomass_0.95ci')
		names(ag_vect_sq)[names(ag_vect_sq) == 'resp_ci'] <- paste0('biomass_ci')

		ggsave(plot = map, filename = paste0(out_dir, '/biomass_', fut, '.png'), width = 12, height = 10, dpi = 600)

		writeVector(this_ag_vect, paste0(out_dir, '/biomass_', fut, '.gpkg'), overwrite = TRUE)

	} # next future

say('########################################')
say('### maps of change in future biomass ###')
say('########################################')

	# user-defined

	# chains <- readRDS(paste0(out_dir, '/chains.rds'))

	nam <- vect('./data_from_gadm/gadm_4pt1_level_1_north_america_sans_alaska_lambert.gpkg')

	ag_vect_deltas <- ag_vect_sq

	# calculate deltas
	summary <- chains$summary$all.chains
	which_resp <- grepl(rownames(summary), pattern = paste0('mu_biomass_county_sq'))
	resp <- summary[which_resp, 'Mean']
	ag_vect_deltas$present_biomass <- resp

	for (fut in futs) {

		this_ag_vect <- get(paste0('ag_vect_', fut))
		which_resp <- grepl(rownames(summary), pattern = paste0('biomass_county_', fut))
		resp <- summary[which_resp, 'Mean']

		ag_vect_deltas$DUMMY <- resp
		names(ag_vect_deltas)[names(ag_vect_deltas) == 'DUMMY'] <- paste0('biomass_', fut)

		ratio <- resp / ag_vect_deltas$present_biomass
		ag_vect_deltas$DUMMY <- ratio
		names(ag_vect_deltas)[names(ag_vect_deltas) == 'DUMMY'] <- paste0('biomass_delta_ratio_', fut)

	}

	# calculate min/max for plotting
	min_val <- Inf
	max_val <- -Inf

	ag_vect_deltas_focal <- ag_vect_deltas[ag_vect_deltas$n_andropogon_gerardi > 0]
	for (fut in futs) {

		vals <- as.data.frame(ag_vect_deltas_focal)[ , paste0('biomass_delta_ratio_', fut), drop = TRUE]
		min_val <- min(c(min_val, vals))
		max_val <- max(c(max_val, vals))

	}

	resp_limits <- c(min_val, max_val)

	### map
	for (fut in futs) {

		say(fut)

		vals <- as.data.frame(ag_vect_deltas)[ , paste0('biomass_delta_ratio_', fut), drop = TRUE]
		ag_vect_deltas$ratio <- vals

		pretty_title <- paste0('SSP ', substr(fut, 4, 6), ' ', substr(fut, 8, 11), '-', substr(fut, 13, 16), ' | biomass-only model')

		# extent
		ag_vect_sq_pres <- ag_vect_sq[ag_vect_sq$n_andropogon_gerardi > 0]
		extent <- ext(ag_vect_sq_pres)
		extent <- as.vector(extent)
		x_range <- (extent[2] - extent[1])
		y_range <- (extent[4] - extent[3])
		extent[1] <- extent[1] + 0.15 * x_range
		extent[3] <- extent[3] + 0.125 * y_range
		extent[4] <- extent[4] - 0.2 * y_range

		legend_title <- 'Biomass change\n(future / present)'

		map <- ggplot() +
			layer_spatial(ag_vect_deltas, aes(fill = ratio), color = NA) +
			layer_spatial(nam, color = 'gray30', fill = NA, linewidth = 0.3) +
			layer_spatial(cents_with_ag, pch = 16, color = alpha('gray20', 0.5), size = 0.35) +
			layer_spatial(site_data_with_biomass_vect, pch = 3, size = 4) +
			scale_fill_gradient2(
				name = legend_title,
				low = '#7b3294',
				# mid = '#ffffbf',
				mid = 'beige',
				high = '#008837',
				midpoint = 1,
				limits = resp_limits
			) +
			xlim(extent[1], extent[2]) + ylim(extent[3], extent[4]) +
				ggtitle(
					bquote('Change in ' * italic('Andropogon gerardi') * ' biomass '),
					subtitle = pretty_title
				) +
			theme(
				plot.title = element_text(size = 16),
				plot.subtitle = element_text(size = 14)
			)

		ggsave(plot = map, filename = paste0(out_dir, '/biomass_change_ratio_', fut, '.png'), width = 12, height = 10, dpi = 600)

	} # next future

say(date())
say('FINIS!', deco = '+', level = 1)

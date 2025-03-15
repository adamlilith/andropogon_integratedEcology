### MODELING ANDROPOGON GERARDI DISTRIBUTION, PHENOTYPE, PHYSIOLOGY, GENOTYPE, and ASSOCIATED MICROBIAL COMMUNITIES
### Erica Newman | Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2023-12
###
### This script constructs an integrated species distribution model for Andropogon gerardi, where occurrence is a function of environmental variables and biomass, which is in turn also a function of environmental variables.
###
### source('C:/Ecology/R/andropogon_integratedEcology/sdm_pdm/sdm_pdm_03b_model_biomass_alone.r')
### source('C:/Subarashi/R/andropogon_integratedEcology/sdm_pdm/sdm_pdm_03b_model_biomass_alone.r')
###
### CONTENTS ###
### setup ###
### user-defined values ###
### model biomass ###

#############
### setup ###
#############

	rm(list = ls())

	drive <- 'C:/Ecology/'
	# drive <- 'C:/Subarashi/'

	setwd(paste0(drive, '/Research/Andropogon/Andropogon'))
	source(paste0(drive, '/R/andropogon_integratedEcology/sdm_pdm/sdm_pdm_00_shared_functions_and_variables.r'))

###########################
### user-defined values ###
###########################


	# predictor_names <- c('bio1', 'bio12', 'bio15')
	# pheno_form_mu <- ~ 1 + bio1 + bio12 + bio15 + I(bio1^2) + I(bio12^2) + I(bio15^2)
	# pheno_form_sigma <- ~ 1 + bio1 + bio12 + bio15 + I(bio1^2) + I(bio12^2) + I(bio15^2)

	predictor_names <- c('bio1', 'bio12')
	pheno_form_mu <- ~ 1 + bio1 + bio12 + I(bio1^2) + I(bio12^2)
	pheno_form_sigma <- ~ 1 + bio1 + bio12 + I(bio1^2) + I(bio12^2)

	# ### MCMC settings
	# niter <- 240000
	# nburnin <- 40000
	# thin <- 200
	# nchains <- 4
	# waic <- TRUE
	# trial <- FALSE

	# ### MCMC settings FOR TESTING
	niter <- 110000
	nburnin <- 10000
	thin <- 100
	nchains <- 2
	waic <- FALSE
	trial <- TRUE

	# # ### MCMC settings FOR TESTING
	# niter <- 220
	# nburnin <- 20
	# thin <- 2
	# nchains <- 2
	# waic <- FALSE
	# trial <- TRUE

	out_dir <- paste0('./outputs_loretta/sdm_pdm_biomass_alone', ifelse(trial, '_TRIAL', NULL), '/')
	dirCreate(out_dir)

	# number of values in response curve array used to depict responses of SDM and PDM
	n_response_curve_values <- 200
	
#####################
### model biomass ###
#####################

	sink(paste0(out_dir, '/runtime_log.txt'), split = TRUE)
	say('MODELING BIOMASS ALONE')
	say('Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | ', date(), post = 1)

	### data collation
	##################
	say('data collation', post = 1)

	say('This model is for the spatial distribution of biomass alone (without other traits and unconnected to the occurrence). We model the log of biomass at the site level as a linear function of climate, then use those coefficients to predict biomass at the county level.', breaks = 60, post = 1)

	say('MCMC settings:', level = 2)
	say('niter ...... .', niter)
	say('nburnin ..... ', nburnin)
	say('thin ........ ', thin)
	say('nchains ..... ', nchains, post = 1)

	##################################
	### model predictors and terms ###
	##################################

	say('Predictors and model formula:', level = 2)
	say('Predictors: ', paste(predictor_names, collapse = ', '))
	say('Phenotype formula for mean biomass: ', as.character(pheno_form_mu))
	say('Phenotype formula for SD of biomass: ', as.character(pheno_form_sigma))

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
	
	# create spatial versions of site and phenotype data for plotting
	site_data_with_biomass <- merge(site_data_raw[ , c('site_id', 'LONGITUDE', 'LATITUDE')], biomass_data_raw, by.x = 'site_id', by.y = 'SITE')
	site_data_with_biomass_vect <- vect(site_data_with_biomass, geom = c('LONGITUDE', 'LATITUDE'), crs = getCRS('WGS84'))

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
	# We need this to make spatial predictions.
	say('load occurrences and present-day + future environmental data', level = 3)
	ag_vect_sq <- vect('./data_from_adam_and_loretta/andropogon_gerardi_occurrences_with_environment_1961_2020.gpkg')
	
	### future climates
	ag_vect_ssp245_2041_2070 <- vect('./data_from_adam_and_loretta/andropogon_gerardi_occurrences_with_environment_ensemble_8GCMs_ssp245_2041_2070.gpkg')
	ag_vect_ssp245_2071_2100 <- vect('./data_from_adam_and_loretta/andropogon_gerardi_occurrences_with_environment_ensemble_8GCMs_ssp245_2071_2100.gpkg')
	ag_vect_ssp370_2041_2070 <- vect('./data_from_adam_and_loretta/andropogon_gerardi_occurrences_with_environment_ensemble_8GCMs_ssp370_2041_2070.gpkg')
	ag_vect_ssp370_2071_2100 <- vect('./data_from_adam_and_loretta/andropogon_gerardi_occurrences_with_environment_ensemble_8GCMs_ssp370_2071_2100.gpkg')

	fields <- c('area_km2', 'any_ag_quality_1_to_3', 'num_poaceae_records', predictor_names)
	ag_vect_sq <- ag_vect_sq[ , fields]

	site_data_with_biomass_vect <- project(site_data_with_biomass_vect, ag_vect_sq)

	### collate county-level environmental data for spatial predictions
	###################################################################
	say('collate county-level environmental data', level = 3)

	ag_sq <- as.data.frame(ag_vect_sq)
	x_raw <- ag_sq[ , biomass_predictors_mu, drop = FALSE]
	x_raw_scaled <- scale(x_raw, scale = biomass_x_scales, center = biomass_x_centers)
	x_raw_scaled <- as.data.frame(x_raw_scaled)
	pheno_x_counties_mu_sq <- model.matrix(pheno_form_mu, x_raw_scaled)
	pheno_x_counties_sigma_sq <- model.matrix(pheno_form_sigma, x_raw_scaled)

	# futures
	for (ssp in c('ssp245_2041_2070', 'ssp245_2071_2100', 'ssp370_2041_2070', 'ssp370_2071_2100')) {

		# future: counties with future data
		x_fut <- get(paste0('ag_vect_', ssp))

		x_fut <- as.data.frame(x_fut)

		# PDM: mean of biomass
		this_x_fut <- x_fut[ , biomass_predictors_mu, drop = FALSE]
		x_fut_scaled <- scale(this_x_fut, center = biomass_x_centers, scale = biomass_x_scales)
		x_fut_scaled <- as.data.frame(x_fut_scaled)
		mm <- model.matrix(pheno_form_mu, x_fut_scaled)
		assign(paste0('pheno_x_counties_mu_', ssp), mm)

		# PDM: SD of biomass
		this_x_fut <- x_fut[ , biomass_predictors_sigma, drop = FALSE]
		x_fut_scaled <- scale(this_x_fut, center = biomass_x_centers, scale = biomass_x_scales)
		x_fut_scaled <- as.data.frame(x_fut_scaled)
		mm <- model.matrix(pheno_form_sigma, x_fut_scaled)
		assign(paste0('pheno_x_counties_sigma_', ssp), mm)

	}

	### response array for biomass: response of mean
	################################################
	say('response array for biomass: response of mu', level = 3)

	biomass_x_response_curve_mu <- create_response_curve_array(
		form = pheno_form_mu,
		centers = biomass_x_centers,
		scales = biomass_x_scales,
		ag_pres = ag_pres
	)

	biomass_x_response_curve_sigma <- create_response_curve_array(
		form = pheno_form_sigma,
		centers = biomass_x_centers,
		scales = biomass_x_scales,
		ag_pres = ag_pres
	)

#####################
### inputs for nimble
#####################

	say('Inputs:', level = 2)
	data <- list(
		y_biomass = biomass_data_raw$Biomass
	)

	n_counties <- nrow(ag_sq)

	n_biomass_terms_mu <- ncol(pheno_x_counties_mu_sq)
	n_biomass_terms_sigma <- ncol(pheno_x_counties_sigma_sq)

	constants <- list(
		
		n_pheno_sites = n_pheno_sites, # number of phenotype sample sites
		n_biomass = n_biomass, # number of biomass records

		n_counties = n_counties,

		n_biomass_terms_mu = n_biomass_terms_mu, # number of *terms* in linear predictor of mean biomass, including intercept and non-linear terms
		n_biomass_terms_sigma = n_biomass_terms_sigma, # number of *terms* in linear predictor of biomass, including intercept and non-linear terms

		n_biomass_predictors_mu = n_biomass_predictors_mu, # number of *predictors* in linear predictor of mean biomass
		n_biomass_predictors_sigma = n_biomass_predictors_sigma, # number of *predictors* in linear predictor of SD of biomass

		biomass_by_site_x_mu = biomass_by_site_x_mu, # model matrix for mean site-level biomass at sampled sites
		biomass_by_site_x_sigma = biomass_by_site_x_sigma, # model matrix for site-level SD of biomass at sampled sites
		biomass_site_index = biomass_site_index, # index of sampled site for each row in biomass data

		n_response_curve_values = n_response_curve_values,
		biomass_x_response_curve_mu = biomass_x_response_curve_mu,
		biomass_x_response_curve_sigma = biomass_x_response_curve_sigma,

		pheno_x_counties_mu_sq = pheno_x_counties_mu_sq,
		pheno_x_counties_sigma_sq = pheno_x_counties_sigma_sq,

		pheno_x_counties_mu_ssp245_2071_2100 = pheno_x_counties_mu_ssp245_2071_2100,
		pheno_x_counties_mu_ssp245_2041_2070 = pheno_x_counties_mu_ssp245_2041_2070,
		pheno_x_counties_mu_ssp370_2041_2070 = pheno_x_counties_mu_ssp370_2041_2070,
		pheno_x_counties_mu_ssp370_2071_2100 = pheno_x_counties_mu_ssp370_2071_2100,

		pheno_x_counties_sigma_ssp245_2071_2100 = pheno_x_counties_sigma_ssp245_2071_2100,
		pheno_x_counties_sigma_ssp245_2041_2070 = pheno_x_counties_sigma_ssp245_2041_2070,
		pheno_x_counties_sigma_ssp370_2041_2070 = pheno_x_counties_sigma_ssp370_2041_2070,
		pheno_x_counties_sigma_ssp370_2071_2100 = pheno_x_counties_sigma_ssp370_2071_2100

	)

	biomass_beta_mu_inits <- rep(0, n_biomass_terms_mu)
	biomass_beta_mu_inits[grepl(colnames(pheno_x_counties_mu_sq), pattern = '\\^2')] <- -1
	biomass_beta_mu_inits[grepl(colnames(pheno_x_counties_mu_sq), pattern = '\\:')] <- 0

	biomass_beta_sigma_inits <- rep(0, n_biomass_terms_sigma)
	biomass_beta_sigma_inits[grepl(colnames(pheno_x_counties_sigma_sq), pattern = '\\^2')] <- -2
	biomass_beta_sigma_inits[grepl(colnames(pheno_x_counties_sigma_sq), pattern = '\\:')] <- 0

	pheno_biomass_sq_inits <- pheno_biomass_fut_inits <- rep(mean(biomass_data_raw$Biomass), n_counties)

	inits <- list(

		# y_biomass = rep(mean(biomass_data_raw$Biomass), n_biomass), # biomass

		biomass_beta_mu = biomass_beta_mu_inits, # PDM coefficients
		biomass_beta_sigma = biomass_beta_sigma_inits, # PDM coefficients
		
		biomass_site_level_mu_hat = rep(mean(biomass_data_raw$Biomass), n_pheno_sites), # latent mean biomass at a site
		biomass_site_level_sigma_hat = rep(sd(biomass_data_raw$Biomass), n_pheno_sites), # latent sd of biomass at a site

		biomass_county_mu_ssp245_2041_2070 = pheno_biomass_fut_inits, # future biomass in counties
		biomass_county_mu_ssp245_2071_2100 = pheno_biomass_fut_inits,
		biomass_county_mu_ssp370_2041_2070 = pheno_biomass_fut_inits,
		biomass_county_mu_ssp370_2071_2100 = pheno_biomass_fut_inits,

		biomass_county_sigma_ssp245_2041_2070 = pheno_biomass_fut_inits, # future biomass in counties
		biomass_county_sigma_ssp245_2071_2100 = pheno_biomass_fut_inits,
		biomass_county_sigma_ssp370_2041_2070 = pheno_biomass_fut_inits,
		biomass_county_sigma_ssp370_2071_2100 = pheno_biomass_fut_inits

	)

	say('Data:')
	print(str(data))

	say('Constants:', pre = 1)
	print(str(constants))

	say('Initializations:', pre = 1)
	print(str(inits))

	### define model
	say('nimbleCode():', level = 2)
	say('For the biomass component, we assume a log-linear model model (log(biomass) ~ environment * coefficients.', post = 2)
	model_code <- nimbleCode({
	  
		# phenotype: priors for relationship of mean and variance to environment
		biomass_beta_mu[1] ~ dnorm(0, sd = 10) # intercept... not regularized
		for (i in 2:n_biomass_terms_mu) {
			biomass_beta_mu[i] ~ ddexp(0, rate = 1) # regularization toward 0
			# biomass_beta_mu[i] ~ dnorm(0, sd = 10) # broad prior
		}

		biomass_beta_sigma[1] ~ dnorm(0, sd = 10) # intercept... not regularized
		for (i in 2:n_biomass_terms_sigma) {
			biomass_beta_sigma[i] ~ ddexp(0, rate = 1) # regularization toward 0
			# biomass_beta_sigma[i] ~ dnorm(0, sd = 10) # broad prior
		}

		# likelihood of biomass
		# mean and sd of biomass at a site are latent and functions of environment
		# plant biomasses are samples from the site-level distribution defined by the mean and sd of the given site
		for (i in 1:n_pheno_sites) {
			
			# log(biomass_site_level_mu_hat[i]) <-  inprod(biomass_beta_mu[1:n_biomass_terms_mu], biomass_by_site_x_mu[i, 1:n_biomass_terms_mu])
			biomass_site_level_mu_hat[i] <-  inprod(biomass_beta_mu[1:n_biomass_terms_mu], biomass_by_site_x_mu[i, 1:n_biomass_terms_mu])

			# log(biomass_site_level_sigma_hat[i]) <-  inprod(biomass_beta_sigma[1:n_biomass_terms_sigma], biomass_by_site_x_sigma[i, 1:n_biomass_terms_sigma])
			biomass_site_level_sigma_hat[i] <-  inprod(biomass_beta_sigma[1:n_biomass_terms_sigma], biomass_by_site_x_sigma[i, 1:n_biomass_terms_sigma])
			
			biomass_site_level_tau_hat[i] <- 1 / (biomass_site_level_sigma_hat[i] ^ 2)

		}

		for (i in 1:n_biomass) {

			y_biomass[i] ~ dlnorm(
				meanlog = biomass_site_level_mu_hat[biomass_site_index[i]],
				taulog = biomass_site_level_tau_hat[biomass_site_index[i]]
			)
	
		}
			
		# posterior samplers for predictions to counties in status quo and future
		for (i in 1:n_counties) {

			# biomass mean and sigma: status quo
			biomass_county_mu_sq[i] <- inprod(biomass_beta_mu[1:n_biomass_terms_mu], pheno_x_counties_mu_sq[i, 1:n_biomass_terms_mu])

			biomass_county_sigma_sq[i] <- inprod(biomass_beta_sigma[1:n_biomass_terms_sigma], pheno_x_counties_sigma_sq[i, 1:n_biomass_terms_sigma])

			# biomass mu: future
			biomass_county_mu_ssp245_2041_2070[i] <- inprod(biomass_beta_mu[1:n_biomass_terms_mu], pheno_x_counties_mu_ssp245_2041_2070[i, 1:n_biomass_terms_mu])
			
			biomass_county_mu_ssp245_2071_2100[i] <- inprod(biomass_beta_mu[1:n_biomass_terms_mu], pheno_x_counties_mu_ssp245_2071_2100[i, 1:n_biomass_terms_mu])
			
			biomass_county_mu_ssp370_2041_2070[i] <- inprod(biomass_beta_mu[1:n_biomass_terms_mu], pheno_x_counties_mu_ssp370_2041_2070[i, 1:n_biomass_terms_mu])

			biomass_county_mu_ssp370_2071_2100[i] <- inprod(biomass_beta_mu[1:n_biomass_terms_mu], pheno_x_counties_mu_ssp370_2071_2100[i, 1:n_biomass_terms_mu])

			# biomass sigma: future
			biomass_county_sigma_ssp245_2041_2070[i] <- inprod(biomass_beta_sigma[1:n_biomass_terms_sigma], pheno_x_counties_sigma_ssp245_2041_2070[i, 1:n_biomass_terms_sigma])
			
			biomass_county_sigma_ssp245_2071_2100[i] <- inprod(biomass_beta_sigma[1:n_biomass_terms_sigma], pheno_x_counties_sigma_ssp245_2071_2100[i, 1:n_biomass_terms_sigma])
			
			biomass_county_sigma_ssp370_2041_2070[i] <- inprod(biomass_beta_sigma[1:n_biomass_terms_sigma], pheno_x_counties_sigma_ssp370_2041_2070[i, 1:n_biomass_terms_sigma])

			biomass_county_sigma_ssp370_2071_2100[i] <- inprod(biomass_beta_sigma[1:n_biomass_terms_sigma], pheno_x_counties_sigma_ssp370_2071_2100[i, 1:n_biomass_terms_sigma])

		}

		# posterior predictive sampler for response curves for biomass: site-level mu
		for (j in 1:n_biomass_predictors_mu) {

			for (i in 1:n_response_curve_values) {
					
				pheno_response_curves_biomass_mu[i, j] <-
					inprod(biomass_beta_mu[1:n_biomass_terms_mu], biomass_x_response_curve_mu[i, 1:n_biomass_terms_mu, j])
				
			}

		}

		# posterior predictive sampler for response curves for biomass: site-level sigma
		for (j in 1:n_biomass_predictors_sigma) {

			for (i in 1:n_response_curve_values) {
					
				pheno_response_curves_biomass_sigma[i, j] <-
					inprod(biomass_beta_sigma[1:n_biomass_terms_sigma], biomass_x_response_curve_sigma[i, 1:n_biomass_terms_sigma, j])
				
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
	calc <- model$calculate()
	say('model$calculate(): ', calc)

	say('configureMCMC():', level = 2)
	monitors <- c(
		'biomass_beta_mu', 'biomass_beta_sigma',
		'biomass_site_level_mu_hat', 'biomass_site_level_sigma_hat',

		'pheno_response_curves_biomass_mu', 'pheno_response_curves_biomass_sigma',

		'biomass_county_mu_sq', 'biomass_county_sigma_sq',

		'biomass_county_mu_ssp245_2041_2070', 'biomass_county_mu_ssp245_2071_2100', 'biomass_county_mu_ssp370_2041_2070', 'biomass_county_mu_ssp370_2071_2100',

		'biomass_county_sigma_ssp245_2041_2070', 'biomass_county_sigma_ssp245_2071_2100', 'biomass_county_sigma_ssp370_2041_2070', 'biomass_county_sigma_ssp370_2071_2100'
	)
	
	conf <- configureMCMC(
		model,
		monitors = monitors,
		print = TRUE,
		enableWAIC = waic
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

	chains <- readRDS(paste0(out_dir, '/chains.rds'))
	mcmc <- chains$samples

	cols <- c(
		paste0('biomass_beta_mu[', 1:n_biomass_terms_mu, ']'),
		paste0('biomass_beta_sigma[', 1:n_biomass_terms_sigma, ']'),
		paste0('biomass_site_level_mu_hat[', 1:n_pheno_sites, ']'),
		paste0('biomass_site_level_sigma_hat[', 1:n_pheno_sites, ']')
	)
	for (i in 1:nchains) {
		mcmc[[i]] <- mcmc[[i]][ , cols]
	}

	# graphing trace and density plots for all betas: mu
	pars <- paste0('biomass_beta_mu[', 1:n_biomass_terms_mu, ']')
	file <- paste0(out_dir, '/beta_biomass_mu_trace.png')
	ggsave(mcmc_trace(mcmc, pars = pars), file = file, width = 10, height = 8, dpi = 450, bg = 'white')

	file <- paste0(out_dir, '/beta_biomass_mu_density.png')
	ggsave(mcmc_dens_overlay(mcmc, pars = pars), file = file, width = 10, height = 8, dpi = 450, bg = 'white')

	# graphing trace and density plots for all betas: sigma
	pars <- paste0('biomass_beta_sigma[', 1:n_biomass_terms_sigma, ']')
	file <- paste0(out_dir, '/beta_biomass_sigma_trace.png')
	ggsave(mcmc_trace(mcmc, pars = pars), file = file, width = 10, height = 8, dpi = 450, bg = 'white')

	file <- paste0(out_dir, '/beta_biomass_sigma_density.png')
	ggsave(mcmc_dens_overlay(mcmc, pars = pars), file = file, width = 10, height = 8, dpi = 450, bg = 'white')

	# graphing trace and density plots for all site-level biomass: mu
	pars <- paste0('biomass_site_level_mu_hat[', 1:n_pheno_sites, ']')
	file <- paste0(out_dir, '/site_level_biomass_mu_trace.png')
	ggsave(mcmc_trace(mcmc, pars = pars), file = file, width = 10, height = 8, dpi = 450, bg = 'white')

	file <- paste0(out_dir, '/site_level_biomass_mu_density.png')
	ggsave(mcmc_dens_overlay(mcmc, pars = pars), file = file, width = 10, height = 8, dpi = 450, bg = 'white')

	# graphing trace and density plots for all site-level biomass: sigma
	pars <- paste0('biomass_site_level_sigma_hat[', 1:n_pheno_sites, ']')
	file <- paste0(out_dir, '/site_level_biomass_sigma_trace.png')
	ggsave(mcmc_trace(mcmc, pars = pars), file = file, width = 10, height = 8, dpi = 450, bg = 'white')

	file <- paste0(out_dir, '/site_level_biomass_sigma_density.png')
	ggsave(mcmc_dens_overlay(mcmc, pars = pars), file = file, width = 10, height = 8, dpi = 450, bg = 'white')

	# Gelman-Rubin statistic
	rhats <- tryCatch(
		gelman.diag(mcmc, autoburnin = FALSE, multivariate = TRUE),
		error = function(cond) FALSE
	)
	
	sink(paste0(out_dir, '/convergence.txt'), split = TRUE)
	say('GELMAN-RUBIN STATISTICS')
	say(date(), post = 2)
	print(rhats)
	sink()

say('##############################', pre = 1)
say('### map of current biomass ###')
say('##############################')

	# user-defined


	chains <- readRDS(paste0(out_dir, '/chains.rds'))

	# subset chain summary to just the lambdas
	summary <- chains$summary$all.chains

	nam <- vect(paste0(drive, '/Research Data/GADM/Version 4.1/High Res North America Level 1 sans Great Lakes SpatVector WGS84.gpkg'))
	nam <- project(nam, ag_vect_sq)

	resp_types <- c('mu', 'sigma')
	for (resp_type in resp_types) {

		# just counties with data
		which_resp <- grepl(rownames(summary), pattern = paste0('biomass_county_', resp_type, '_sq'))
		resp <- summary[which_resp, ]

		ag_vect_sq$resp_mean <- resp[ , 'Mean']
		ag_vect_sq$resp_0.05ci <- resp[ , '95%CI_low']
		ag_vect_sq$resp_0.95ci <- resp[ , '95%CI_upp']
		ag_vect_sq$resp_ci <- ag_vect_sq$resp_0.95ci - ag_vect_sq$resp_0.05ci

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

		if (resp_type == 'mu') {
			legend_title <- 'Biomass (g)'
		} else if (resp_type == 'sigma') {
			legend_title <- ' SD of\nbiomass (g)'
		}

		resp_limits <- range(c(ag_vect_sq$resp_mean, site_data_with_biomass_vect$Biomass))

		map <- ggplot() +
			layer_spatial(ag_vect_sq, aes(fill = exp(resp_mean)), color = NA) +
			layer_spatial(nam, color = 'gray30', fill = NA, linewidth = 0.1) +
			layer_spatial(cents_with_ag, pch = 16, color = alpha('gray20', 0.5), size = 0.35) +
			layer_spatial(site_data_with_biomass_vect, aes(fill = Biomass), pch = 21, size = 4) +
			scale_fill_continuous(
				name = legend_title,
				# low = 'lightgreen',
				# high = 'darkgreen',
				low = 'darkorange',
				high = 'darkorange4',
				limits = resp_limits,
				trans = 'log10'
			) +
			xlim(extent[1], extent[2]) + ylim(extent[3], extent[4]) +
			ggtitle(
				bquote('Present-day distribution of ' * italic('Andropogon gerardi') * ' biomass ' * italic(.(resp_type))),
				subtitle = '1991-2020 | biomass-only model') +
			theme(
				plot.title = element_text(size = 16),
				plot.subtitle = element_text(size = 14)
			)

		names(ag_vect_sq)[names(ag_vect_sq) == 'resp_mean'] <- paste0('biomass_', resp_type, '_mean')
		names(ag_vect_sq)[names(ag_vect_sq) == 'resp_0.05ci'] <- paste0('biomass_', resp_type, '_0.05ci')
		names(ag_vect_sq)[names(ag_vect_sq) == 'resp_0.95ci'] <- paste0('biomass_', resp_type, '_0.95ci')
		names(ag_vect_sq)[names(ag_vect_sq) == 'resp_ci'] <- paste0('biomass_', resp_type, '_ci')

		ggsave(plot = map, filename = paste0(out_dir, '/biomass_', resp_type, '_status_quo.png'), width = 12, height = 10, dpi = 600)

	}

	writeVector(ag_vect_sq, paste0(out_dir, '/biomass_status_quo.gpkg'), overwrite = TRUE)

say('###################################')
say('### maps of future distribution ###')
say('###################################')

	# user-defined

	chains <- readRDS(paste0(out_dir, '/chains.rds'))

	futs <- c(
		'ssp245_2041_2070',
		'ssp245_2071_2100',
		'ssp370_2041_2070',
		'ssp370_2071_2100'
	)

	nam <- vect(paste0(drive, '/Research Data/GADM/Version 4.1/High Res North America Level 1 sans Great Lakes SpatVector WGS84.gpkg'))
	nam <- project(nam, ag_vect_sq)

	# subset chain summary to just the lambdas
	summary <- chains$summary$all.chains

	resp_types <- c('mu', 'sigma')
	for (fut in futs) {

		for (resp_type in resp_types) {

			which_resp <- grepl(rownames(summary), pattern = paste0('biomass_county_', resp_type, '_', fut))
			resp <- summary[which_resp, ]

			this_ag_vect <- get(paste0('ag_vect_', fut))

			# get just counties with data
			this_ag_vect$resp_mean <- resp[ , 'Mean']
			this_ag_vect$resp_0.05ci <- resp[ , '95%CI_low']
			this_ag_vect$resp_0.95ci <- resp[ , '95%CI_upp']
			this_ag_vect$resp_ci <- this_ag_vect$resp_0.95ci - this_ag_vect$resp_0.05ci

			pretty_title <- paste0('SSP ', substr(fut, 4, 6), ' ', substr(fut, 8, 11), '-', substr(fut, 13, 16), ' | biomass-only model')

			if (resp_type == 'mu') {
				legend_title <- 'Biomass (g)'
			} else if (resp_type == 'sigma') {
				legend_title <- ' SD of\nbiomass (g)'
			}

			resp_limits <- range(c(this_ag_vect$resp_mean, site_data_with_biomass_vect$Biomass))

			map <- ggplot() +
				layer_spatial(this_ag_vect, aes(fill = exp(resp_mean)), color = NA) +
				layer_spatial(nam, color = 'gray30', fill = NA, linewidth = 0.1) +
				layer_spatial(cents_with_ag, pch = 16, color = alpha('gray20', 0.5), size = 0.35) +
				layer_spatial(site_data_with_biomass_vect, aes(fill = Biomass), pch = 21, size = 4) +
				scale_fill_continuous(
					name = legend_title,
					low = 'darkorange',
					high = 'darkorange4',
					limits = resp_limits,
					trans = 'log10'
				) +
				xlim(extent[1], extent[2]) + ylim(extent[3], extent[4]) +
					ggtitle(
						bquote('Future equilibrial distribution of ' * italic('Andropogon gerardi') * ' biomass ' * italic(.(resp_type))),
						subtitle = pretty_title
					) +
				theme(
					plot.title = element_text(size = 16),
					plot.subtitle = element_text(size = 14)
				)

			names(ag_vect_sq)[names(ag_vect_sq) == 'resp_mean'] <- paste0('biomass_', resp_type, '_mean')
			names(ag_vect_sq)[names(ag_vect_sq) == 'resp_0.05ci'] <- paste0('biomass_', resp_type, '_0.05ci')
			names(ag_vect_sq)[names(ag_vect_sq) == 'resp_0.95ci'] <- paste0('biomass_', resp_type, '_0.95ci')
			names(ag_vect_sq)[names(ag_vect_sq) == 'resp_ci'] <- paste0('biomass_', resp_type, '_ci')

			ggsave(plot = map, filename = paste0(out_dir, '/biomass_', resp_type, '_', fut, '.png'), width = 12, height = 10, dpi = 600)

		} # next resp_type

		writeVector(this_ag_vect, paste0(out_dir, '/biomass_', fut, '.gpkg'), overwrite = TRUE)

	} # next future

say('###########################')
say('### parameter estimates ###')
say('###########################')

	chains <- readRDS(paste0(out_dir, '/chains.rds'))

	# subset chain summary to just the lambdas
	summary <- chains$summary$all.chains

	resp_types <- c('mu', 'sigma')
	for (resp_type in resp_types) {
		
		form <- get(paste0('pheno_form_', resp_type))

		model_term <- attr(terms(form), 'term.labels')
		model_term_nice <- model_term
		model_term_nice <- gsub(model_term_nice, pattern = 'I\\(', replacement = '')
		model_term_nice <- gsub(model_term_nice, pattern = '\\)', replacement = '')
		model_term_nice <- c('Intercept', model_term_nice)

		pars <- paste0('biomass_beta_', resp_type, '[', 1:(1 + length(model_term)), ']')

		coeff_chains <- chains$samples
		for (i in 1:nchains) {
			coeff_chains[[i]] <- coeff_chains[[i]][ , pars]
			colnames(coeff_chains[[i]]) <- model_term_nice
		}

		caterpillars <- mcmc_intervals(coeff_chains, pars = model_term_nice)

		# plot posterior of coefficient estimates
		caterpillars <- caterpillars +
			xlab('Estimated value') +
			ggtitle(paste0('Posterior distribution of coefficients for ', resp_type, ' of biomass-only model')) +
			theme(plot.background = element_rect(fill = 'white'))

		ggsave(caterpillars, filename = paste0(out_dir, '/coefficients_biomass_', resp_type, '.png'), width = 10, height = 12, dpi = 300)

	}

say('#######################')
say('### response curves ###')
say('#######################')

	chains <- readRDS(paste0(out_dir, '/chains.rds'))

	resp_type <- c('mu', 'sigma')
	for (resp_type in resp_types) {

		if (resp_type == 'mu') {
			n_predictors <- n_biomass_predictors_mu
			y_lab <- bquote('Biomass ' * italic(.(resp_type)) * ' (g)')
		} else if (resp_type == 'sigma') {
			n_predictors <- n_biomass_predictors_sigma
			y_lab <- bquote('Biomass ' * italic(.(resp_type)) * ' (g)')
		}

		# make a plot for how biomass mu and sigma responds to each predictor
		responses <- list()
		for (i in 1:n_predictors) {

			pred <- predictor_names[i]

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

			nice_title <- bquote('Biomass ' * italic(.(resp_type)) * ' versus ' * .(nice_title))

			# unscale predictor value
			if (resp_type == 'mu') {
				x <- biomass_x_response_curve_mu[ , pred, pred]
			} else if (resp_type == 'sigma') {
				x <- biomass_x_response_curve_sigma[ , pred, pred]
			}

			x <- x * biomass_x_scales[[pred]] + biomass_x_centers[[pred]]

			# create data frame with SDM response
			pars <- paste0('pheno_response_curves_biomass_mu[', 1:n_response_curve_values, ', ', i, ']')
			response_mean <- chains$summary$all.chains[pars, 'Mean']
			response_lower <- chains$summary$all.chains[pars, '95%CI_low']
			response_upper <- chains$summary$all.chains[pars, '95%CI_upp']

			if (resp_type == 'mu') {
				response_mean <- exp(response_mean)
				response_lower <- exp(response_lower)
				response_upper <- exp(response_upper)
			}

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

			max_val <- max(-Inf, df_upper$response[!is.infinite(df_upper$response)])

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
					plot.title = element_text(size = 10)
				)

			responses[[pred]] <- response

		} # next predictor

		for (i in 1:n_predictors) {
			responses[[i]] <- responses[[i]] + ylim(0, max_val)
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

		ggsave(plot = responses, filename = paste0(out_dir, '/response_curves_biomass_', resp_type, '.png'), width = width, height = height, dpi = 600)

	}

say(date())
say('FINIS!', deco = '+', level = 1)

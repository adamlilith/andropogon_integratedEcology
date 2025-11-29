### MODELING ANDROPOGON GERARDI DISTRIBUTION, PHENOTYPE, PHYSIOLOGY, GENOTYPE, and ASSOCIATED MICROBIAL COMMUNITIES
### Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2023-12
###
### This model estimates the distribution of site-level mean Andropogon gerardi ramet biomass. The distribution of biomasses among ramets at a site follows a zero-inflated gamma distribution with individual ramet biomasses drawn from this distribution. The site-level mean biomass, standard deviation in biomass, and probability of a zero value are each functions of one or more soil/climate predictors (possibly with higher-order terms).
###
### source('C:/Kaji/R/andropogon_integratedEcology/sdm_pdm/sdm_pdm_04d_model_biomass_zig~normal_heteroscedastic_fx_of_environment.r')

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

	# trial <- TRUE # TRUE for testing
	trial <- FALSE # TRUE for testing

	# do cross-validation?
	crossvalidate <- FALSE
	# crossvalidate <- TRUE

	### formula for how aspects of species responds to environment

	# formula_biomass_mu <- ~ 1 + bio12 # response of biomass to environment
	# preds_filename <- 'bio12'

	formula_biomass_mu <- ~ 1 + bio12 + I(bio12^2) # response of biomass to environment
	preds_filename <- 'bio12^2'

	# formula_biomass_mu <- ~ 1 + bio1 # response of biomass to environment
	# preds_filename <- 'bio1'

	# formula_biomass_mu <- ~ 1 + bio1 + I(bio1^2) # response of biomass to environment
	# preds_filename <- 'bio1^2'

	# formula_biomass_mu <- ~ 1 + bio1 + bio12 # response of biomass to environment
	# preds_filename <- 'bio1_bio12'

	# formula_biomass_mu <- ~ 1 + bio1 + bio12 + bio1:bio12 # response of biomass to environment
	# preds_filename <- 'bio1_x_bio12'

	# formula_biomass_mu <- ~ 1 + bio1 + bio12 + I(bio1^2) # response of biomass to environment
	# preds_filename <- 'bio1^2_bio12'

	# formula_biomass_mu <- ~ 1 + bio1 + bio12 + I(bio12^2) # response of biomass to environment
	# preds_filename <- 'bio1_bio12^2'

	out_dir <- paste0('./outputs_loretta/integrated_sdm_pdm/models_biomass/', ifelse(trial, 'TRIAL_', ''), '[biomass_zig~normal_heteroscedastic_', preds_filename, ']/')
	previous_dir_hetero <- paste0('./outputs_loretta/integrated_sdm_pdm/models_biomass/', ifelse(trial, 'TRIAL_', ''), '[biomass_gamma~normal_heteroscedastic_', preds_filename, ']/')
	previous_dir_zig <- paste0('./outputs_loretta/integrated_sdm_pdm/models_biomass/', ifelse(trial, 'TRIAL_', ''), '[biomass_zig~normal_homoscedastic_', preds_filename, ']/')

	# calib <- TRUE # use just counties with non-NA Poaceae for calibration region
	calib <- FALSE # use all of North America for calibration region

	homoscedastic <- FALSE # SPECIFIC TO THIS SCRIPT--SHOULD NOT BE CHANGED
	zero_inflated <- TRUE # SPECIFIC TO THIS SCRIPT--SHOULD NOT BE CHANGED

	formula_biomass_sigma <- formula_biomass_mu
	formula_pzero <- formula_biomass_mu

	if (!trial) {

		### MCMC settings
		niter <- 1600000
		nburnin <- niter / 2
		thin <- 800
		nchains <- 4
		waic <- TRUE


	} else {
		
		# ### MCMC settings FOR TESTING
		# niter <- 2000
		# nburnin <- 1000
		# thin <- 1
		# nchains <- 2

		niter <- 10000
		nburnin <- 5000
		thin <- 5
		nchains <- 2

	}

#############
### model ###
#############

	if (!trial) if (file.exists(out_dir)) stop('Output folder already exists.')
	dirCreate(out_dir)

	sink(paste0(out_dir, '/runtime_log.txt'), split = TRUE)
	say('MODELING BIOMASS')
	say('Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | ', date(), post = 1)

	### data collation
	##################
	say('data collation', post = 1)

	say('This model estimates the distribution of site-level mean Andropogon gerardi ramet biomass. The distribution of biomasses among ramets at a site follows a zero-inflated gamma distribution with individual ramet biomasses drawn from this distribution. The site-level mean biomass, standard deviation in biomass, and probability of a zero value are each functions of one or more soil/climate predictors (possibly with higher-order terms).', breaks = 60, post = 1)

	say('MCMC settings:', level = 2)
	say('trial ........................ ', trial)
	say('niter ........................ ', niter)
	say('nburnin ...................... ', nburnin)
	say('thin ......................... ', thin)
	say('nchains ...................... ', nchains)
	say('formula_biomass_mu ........... ', paste(as.character(formula_biomass_mu), collapse = ' '))
	say('formula_biomass_sigma ........ ', paste(as.character(formula_biomass_sigma), collapse = ' '))
	say('formula_pzero ........ ', paste(as.character(formula_pzero), collapse = ' '))
	say('homoscedastic ................ ', homoscedastic)
	say('zero_inflated ................ ', zero_inflated)
	say('calib ........................ ', calib, post = 2)

	say('out_dir:')
	say(out_dir)

	formulae <- list(
		formula_biomass_mu = formula_biomass_mu,
		formula_biomass_sigma = formula_biomass_sigma,
		formula_pzero = formula_pzero
	)

	########################s
	### data preparation ###
	########################

	data_biomass_mu <- prepare_biomass(formula_biomass = formula_biomass_mu, n_response_curve_values = n_response_curve_values, calib = calib)
	data_occs <- prepare_occurrences(formula_occs = ~ 1, formula_occs_bias = ~ 1, n_response_curve_values = n_response_curve_values, psa_quant = psa_quant, calib = calib)
	# Calculate mean biomass by site

	#########################
	### inputs for nimble ###
	#########################

	say('Inputs:', level = 2)
	data <- list(
		y_biomass = data_biomass_mu$y_biomass		# biomass of individual plants
	)

	constants <- list(
		
		# sampled sites
		n_pheno_sites = data_biomass_mu$n_pheno_sites, # number of phenotype sample sites

		### biomass
		x_by_site_biomass = data_biomass_mu$x_by_site_biomass, # MM with covariates for biomass (scaled)
		n_biomass = data_biomass_mu$n_biomass, # number of biomass observations
		site_index_biomass = data_biomass_mu$site_index_biomass, # index of sampled site for each row in biomass data
		n_terms_biomass_mu = data_biomass_mu$n_terms_biomass, # number of terms in formula for biomass model (including intercept)
		n_terms_biomass_sigma = data_biomass_mu$n_terms_biomass, # number of terms in formula for biomass model (including intercept)

		n_covariates_biomass = data_biomass_mu$n_covariates_biomass,
		resp_curves_x_biomass = data_biomass_mu$resp_curves_x_biomass, # response curve array for biomass

		counties_x_biomass_sq = data_biomass_mu$counties_x_biomass_sq,
		counties_x_biomass_ssp245_2041_2070 = data_biomass_mu$counties_x_biomass_ssp245_2041_2070,
		counties_x_biomass_ssp245_2071_2100 = data_biomass_mu$counties_x_biomass_ssp245_2071_2100,
		counties_x_biomass_ssp370_2041_2070 = data_biomass_mu$counties_x_biomass_ssp370_2041_2070,
		counties_x_biomass_ssp370_2071_2100 = data_biomass_mu$counties_x_biomass_ssp370_2071_2100,

		counties_x_biomass_thirties = data_biomass_mu$counties_x_biomass_thirties,

		# counties
		n_counties = data_biomass_mu$n_counties, # number of counties in the dataset
		n_counties_20th_cent = data_biomass_mu$n_counties_20th_cent, # number of counties in the dataset

		# response curves (general)
		n_response_curve_values = n_response_curve_values # number of values in response curve array

	)

	constants <- c(constants, constants_shared_biomass)

	# initializations
	previous_chains <- readRDS(paste0(previous_dir_hetero, '/chains.rds'))

	beta_biomass_mu_inits <- hammer_extract(previous_chains, 'beta_biomass_mu', j = TRUE, stat = 'mean')
	beta_biomass_sigma_inits <- hammer_extract(previous_chains, 'beta_biomass_sigma', j = TRUE, stat = 'mean')

	sigma_biomass_within_sites_init <- hammer_extract(previous_chains, 'sigma_biomass_within_sites', stat = 'mean')
	sigma_biomass_within_sites_log_init <- log(sigma_biomass_within_sites_init)

	rm(previous_chains)

	previous_chains <- readRDS(paste0(previous_dir_zig, '/chains.rds'))

	beta_biomass_pzero_inits <- hammer_extract(previous_chains, 'beta_biomass_pzero', j = TRUE, stat = 'mean')

	rm(previous_chains)

	inits <- list(

		site_biomass_mu_log = rep(2, data_biomass_mu$n_pheno_sites),
		sigma_biomass_within_sites_log = sigma_biomass_within_sites_log_init,
		y_biomass_sim = data_biomass_mu$y_biomass, # simulated values for biomass (for DHARMa residuals)
		beta_biomass_mu = beta_biomass_mu_inits,
		beta_biomass_sigma = beta_biomass_sigma_inits,
		beta_biomass_pzero = beta_biomass_pzero_inits,

		log_biomass_county_mu_sq = rep(0, data_biomass_mu$n_counties),
		log_biomass_county_mu_ssp245_2041_2070 = rep(0, data_biomass_mu$n_counties),
		log_biomass_county_mu_ssp245_2071_2100 = rep(0, data_biomass_mu$n_counties),
		log_biomass_county_mu_ssp370_2041_2070 = rep(0, data_biomass_mu$n_counties),
		log_biomass_county_mu_ssp370_2071_2100 = rep(0, data_biomass_mu$n_counties),

		log_biomass_county_mu_thirties = rep(0, data_biomass_mu$n_counties_20th_cent)

	)

	say('Data:')
	print(str(data))

	say('Constants:', pre = 1)
	print(str(constants))

	say('Initializations:', pre = 1)
	print(str(inits))

	### define model
	say('nimbleCode():', level = 2)

	### model
	#########
	model_code <- nimbleCode({
	
		# BIOMASS: priors for relationship of site-level mean biomass to environment
		for (i in 1:n_terms_biomass_mu) {
			beta_biomass_mu[i] ~ dnorm(0, sd = beta_biomass_mu_prior_dnorm_sd) # broad prior
		}

		# BIOMASS: priors for relationship of site-level s.d. of biomass to environment
		for (i in 1:n_terms_biomass_sigma) {
			beta_biomass_sigma[i] ~ dnorm(0, sd = beta_biomass_sigma_prior_dnorm_sd) # broad prior
		}

		# BIOMASS: priors for relationship of site-level zero values
		for (i in 1:n_terms_biomass_mu) {
			beta_biomass_pzero[i] ~ dnorm(0, sd = beta_biomass_mu_prior_dnorm_sd) # broad prior
		}

		# prior for sd of individual plant biomass on lognormal (~ half-Cauchy), ==> vague
		sigma_biomass_within_sites_log ~ dnorm(0, sd = sigma_biomass_within_sites_log_prior_sd)
		sigma_biomass_within_sites <- exp(sigma_biomass_within_sites_log)

		# BIOMASS: parameters of biomass distribution are latent and functions of environment
		# individual plant biomasses are samples from the site-level distribution defined by the site-level distribution (next chunk after this one)
		for (i in 1:n_pheno_sites) {

			# relationship of biomass to the environment
			site_biomass_mu_log[i] ~ dnorm(site_biomass_mu_mean_log[i], sd = sigma_biomass_among_sites[i])
			site_biomass_mu_mean_log[i] <- inprod(beta_biomass_mu[1:n_terms_biomass_mu], x_by_site_biomass[i, 1:n_terms_biomass_mu])
			
			# site-level mean biomass
			mu_biomass_site[i] <- exp(site_biomass_mu_log[i])

			# site-level probability of zero biomass
			logit(pzero[i]) <- inprod(beta_biomass_pzero[1:n_terms_biomass_mu], x_by_site_biomass[i, 1:n_terms_biomass_mu])

			# site-level s.d. of biomass
			sigma_biomass_among_sites_log[i] <- inprod(beta_biomass_sigma[1:n_terms_biomass_sigma], x_by_site_biomass[i, 1:n_terms_biomass_sigma])
			sigma_biomass_among_sites[i] <- exp(sigma_biomass_among_sites_log[i])

			# moment matching to get dgamma() parameters
			shape_biomass[i] <- mu_biomass_site[i]^2 / sigma_biomass_within_sites^2
			rate_biomass[i] <- mu_biomass_site[i] / sigma_biomass_within_sites^2

		}

		# BIOMASS: likelihood of individual plants
		for (i in 1:n_biomass) {

			# likelihood
			y_biomass[i] ~ dzigamma(shape = shape_biomass[site_index_biomass[i]], rate = rate_biomass[site_index_biomass[i]], pzero = pzero[site_index_biomass[i]])

			# simulated values for unconditional DHARMa residuals
			y_biomass_sim[i] ~ dzigamma(shape = shape_biomass[site_index_biomass[i]], rate = rate_biomass[site_index_biomass[i]], pzero = pzero[site_index_biomass[i]])

			log_lik_biomass[i] <- dzigamma(y_biomass[i], shape = shape_biomass[site_index_biomass[i]], rate = rate_biomass[site_index_biomass[i]], pzero = pzero[site_index_biomass[i]])
	
		}

		log_lik <- sum(log_lik_biomass[1:n_biomass])

		# BIOMASS: posterior samplers for predictions of to counties in status quo and future
		for (i in 1:n_counties) {

			# biomass: status quo
			mu_biomass_county_sq[i] <- exp(log_biomass_county_mu_sq[i])
			log_biomass_county_mu_sq[i] ~ dnorm(log_biomass_county_mu_unscaled_sq[i], sd = sigma_biomass_among_sites_sq[i])
			log_biomass_county_mu_unscaled_sq[i] <- inprod(beta_biomass_mu[1:n_terms_biomass_mu], counties_x_biomass_sq[i, 1:n_terms_biomass_mu])

			sigma_biomass_among_sites_log_sq[i] <- inprod(beta_biomass_sigma[1:n_terms_biomass_sigma], counties_x_biomass_sq[i, 1:n_terms_biomass_sigma])
			sigma_biomass_among_sites_sq[i] <- exp(sigma_biomass_among_sites_log_sq[i])

			logit(pzero_biomass_county_sq[i]) <-
				inprod(beta_biomass_pzero[1:n_terms_biomass_mu], counties_x_biomass_sq[i, 1:n_terms_biomass_mu])

			# biomass: ssp245_2041_2070
			mu_biomass_county_ssp245_2041_2070[i] <- exp(log_biomass_county_mu_ssp245_2041_2070[i])
			log_biomass_county_mu_ssp245_2041_2070[i] ~ dnorm(log_biomass_county_mu_unscaled_ssp245_2041_2070[i], sd = sigma_biomass_among_sites_ssp245_2041_2070[i])
			log_biomass_county_mu_unscaled_ssp245_2041_2070[i] <-
				inprod(beta_biomass_mu[1:n_terms_biomass_mu], counties_x_biomass_ssp245_2041_2070[i, 1:n_terms_biomass_mu])
			
			sigma_biomass_among_sites_log_ssp245_2041_2070[i] <- inprod(beta_biomass_sigma[1:n_terms_biomass_sigma], counties_x_biomass_ssp245_2041_2070[i, 1:n_terms_biomass_sigma])
			sigma_biomass_among_sites_ssp245_2041_2070[i] <- exp(sigma_biomass_among_sites_log_ssp245_2041_2070[i])

			logit(pzero_biomass_county_ssp245_2041_2070[i]) <-
				inprod(beta_biomass_pzero[1:n_terms_biomass_mu], counties_x_biomass_ssp245_2041_2070[i, 1:n_terms_biomass_mu])

			# biomass: ssp245_2071_2100
			mu_biomass_county_ssp245_2071_2100[i] <- exp(log_biomass_county_mu_ssp245_2071_2100[i])
			log_biomass_county_mu_ssp245_2071_2100[i] ~ dnorm(log_biomass_county_mu_unscaled_ssp245_2071_2100[i], sd = sigma_biomass_among_sites_ssp245_2071_2100[i])
			log_biomass_county_mu_unscaled_ssp245_2071_2100[i] <-
				inprod(beta_biomass_mu[1:n_terms_biomass_mu], counties_x_biomass_ssp245_2071_2100[i, 1:n_terms_biomass_mu])
			
			sigma_biomass_among_sites_log_ssp245_2071_2100[i] <- inprod(beta_biomass_sigma[1:n_terms_biomass_sigma], counties_x_biomass_ssp245_2071_2100[i, 1:n_terms_biomass_sigma])
			sigma_biomass_among_sites_ssp245_2071_2100[i] <- exp(sigma_biomass_among_sites_log_ssp245_2071_2100[i])

			logit(pzero_biomass_county_ssp245_2071_2100[i]) <-
				inprod(beta_biomass_pzero[1:n_terms_biomass_mu], counties_x_biomass_ssp245_2071_2100[i, 1:n_terms_biomass_mu])

			# biomass: ssp370_2041_2070
			mu_biomass_county_ssp370_2041_2070[i] <- exp(log_biomass_county_mu_ssp370_2041_2070[i])
			log_biomass_county_mu_ssp370_2041_2070[i] ~ dnorm(log_biomass_county_mu_unscaled_ssp370_2041_2070[i], sd = sigma_biomass_among_sites_ssp370_2041_2070[i])
			log_biomass_county_mu_unscaled_ssp370_2041_2070[i] <-
				inprod(beta_biomass_mu[1:n_terms_biomass_mu], counties_x_biomass_ssp370_2041_2070[i, 1:n_terms_biomass_mu])

			sigma_biomass_among_sites_log_ssp370_2041_2070[i] <- inprod(beta_biomass_sigma[1:n_terms_biomass_sigma], counties_x_biomass_ssp370_2041_2070[i, 1:n_terms_biomass_sigma])
			sigma_biomass_among_sites_ssp370_2041_2070[i] <- exp(sigma_biomass_among_sites_log_ssp370_2041_2070[i])

			logit(pzero_biomass_county_ssp370_2041_2070[i]) <-
				inprod(beta_biomass_pzero[1:n_terms_biomass_mu], counties_x_biomass_ssp370_2041_2070[i, 1:n_terms_biomass_mu])

			# biomass: ssp370_2071_2100
			mu_biomass_county_ssp370_2071_2100[i] <- exp(log_biomass_county_mu_ssp370_2071_2100[i])
			log_biomass_county_mu_ssp370_2071_2100[i] ~ dnorm(log_biomass_county_mu_unscaled_ssp370_2071_2100[i], sd = sigma_biomass_among_sites_ssp370_2071_2100[i])
			log_biomass_county_mu_unscaled_ssp370_2071_2100[i] <-
				inprod(beta_biomass_mu[1:n_terms_biomass_mu], counties_x_biomass_ssp370_2071_2100[i, 1:n_terms_biomass_mu])

			sigma_biomass_among_sites_log_ssp370_2071_2100[i] <- inprod(beta_biomass_sigma[1:n_terms_biomass_sigma], counties_x_biomass_ssp370_2071_2100[i, 1:n_terms_biomass_sigma])
			sigma_biomass_among_sites_ssp370_2071_2100[i] <- exp(sigma_biomass_among_sites_log_ssp370_2071_2100[i])

			logit(pzero_biomass_county_ssp370_2071_2100[i]) <-
				inprod(beta_biomass_pzero[1:n_terms_biomass_mu], counties_x_biomass_ssp370_2071_2100[i, 1:n_terms_biomass_mu])

		}

		# BIOMASS: posterior samplers for predictions of to counties in 20th century
		for (i in 1:n_counties_20th_cent) {

			# thirties
			mu_biomass_county_thirties[i] <- exp(log_biomass_county_mu_thirties[i])
			log_biomass_county_mu_thirties[i] ~ dnorm(log_biomass_county_mu_unscaled_thirties[i], sd = sigma_biomass_among_sites_thirties[i])
			log_biomass_county_mu_unscaled_thirties[i] <-
				inprod(beta_biomass_mu[1:n_terms_biomass_mu], counties_x_biomass_thirties[i, 1:n_terms_biomass_mu])

			sigma_biomass_among_sites_log_thirties[i] <- inprod(beta_biomass_sigma[1:n_terms_biomass_sigma], counties_x_biomass_thirties[i, 1:n_terms_biomass_sigma])
			sigma_biomass_among_sites_thirties[i] <- exp(sigma_biomass_among_sites_log_thirties[i])

			logit(pzero_biomass_county_thirties[i]) <-
				inprod(beta_biomass_pzero[1:n_terms_biomass_mu], counties_x_biomass_thirties[i, 1:n_terms_biomass_mu])

		}

	})

	if (data_biomass_mu$n_covariates_biomass == 1) {

		### univariate
		response_curve_code <- nimbleCode({

			# BIOMASS: posterior predictive sampler for response curves: site-level mean
			# NB we assume just one predictor for biomass, so the response curve "x" is a matrix, not an array
			for (i in 1:n_response_curve_values) {
					
				log_response_curves_biomass_mu[i] <-
					inprod(beta_biomass_mu[1:n_terms_biomass_mu], resp_curves_x_biomass[i, 1:n_terms_biomass_mu])
				response_curves_biomass_mu[i] <- exp(log_response_curves_biomass_mu[i])
				
				log_response_curves_biomass_sigma[i] <-
					inprod(beta_biomass_sigma[1:n_terms_biomass_sigma], resp_curves_x_biomass[i, 1:n_terms_biomass_sigma])
				response_curves_biomass_sigma[i] <- exp(log_response_curves_biomass_sigma[i])
				
				logit(response_curves_biomass_pzero[i]) <-
					inprod(beta_biomass_pzero[1:n_terms_biomass_mu], resp_curves_x_biomass[i, 1:n_terms_biomass_mu])

			}

		})
	
	} else {
	
		### multivariate
		response_curve_code <- nimbleCode({

			# BIOMASS: posterior predictive sampler for response curves: site-level mean
			# NB we assume just one predictor for biomass, so the response curve "x" is a matrix, not an array
			for (i in 1:n_covariates_biomass) {
				
				for (j in 1:n_response_curve_values) {
						
					log_response_curves_biomass_mu[j, i] <-
						inprod(beta_biomass_mu[1:n_terms_biomass_mu], resp_curves_x_biomass[j, 1:n_terms_biomass_mu, i])
					response_curves_biomass_mu[j, i] <- exp(log_response_curves_biomass_mu[j, i])
					
					log_response_curves_biomass_sigma[j, i] <-
						inprod(beta_biomass_sigma[1:n_terms_biomass_sigma], resp_curves_x_biomass[j, 1:n_terms_biomass_sigma, i])
					response_curves_biomass_sigma[j, i] <- exp(log_response_curves_biomass_sigma[j, i])

					logit(response_curves_biomass_pzero[j, i]) <-
						inprod(beta_biomass_pzero[1:n_terms_biomass_mu], resp_curves_x_biomass[j, 1:n_terms_biomass_mu, i])

				}
			}

		})
	
	}

	model_code <- glueNimbleCode(model_code, response_curve_code)

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
	if (is.na(calc) | is.infinite(calc)) stop('Likelihood is incalculable.')

	say('configureMCMC():', level = 2)

	# monitors for coefficients that have no indexing
	monitors_coeffs_not_indexed <- c(
		'sigma_biomass_within_sites'
	)

	# coefficients that have bracketed indexing
	monitors_coeffs_single_index <- c(
		'beta_biomass_mu', 'beta_biomass_sigma', 'beta_biomass_pzero'
	)

	monitors_coeffs_double_index <- c(
	)

	monitors_derived_not_indexed <- c(
		'log_lik'
	)

	monitors_derived_single_index <- c(
		'mu_biomass_site', 'sigma_biomass_among_sites'
	)
	monitors_derived_double_index <- c(
	)

	monitors_geog_nam <- c(
		'mu_biomass_county_sq', 'mu_biomass_county_ssp245_2041_2070', 	
		'mu_biomass_county_ssp245_2071_2100', 'mu_biomass_county_ssp370_2041_2070', 'mu_biomass_county_ssp370_2071_2100',
		'pzero_biomass_county_sq', 'pzero_biomass_county_ssp245_2041_2070', 	
		'pzero_biomass_county_ssp245_2071_2100', 'pzero_biomass_county_ssp370_2041_2070', 'pzero_biomass_county_ssp370_2071_2100'
	)

	monitors_geog_conus <- c(
		'mu_biomass_county_thirties'
	)

	monitors_dharma <- c(
		'y_biomass_sim'
	)

	monitors_resp_curves <- c(
		'response_curves_biomass_mu', 'response_curves_biomass_sigma', 'response_curves_biomass_pzero'
	)

	monitors <- c(monitors_coeffs_not_indexed, monitors_coeffs_single_index, monitors_coeffs_double_index, monitors_derived_not_indexed, monitors_derived_single_index, monitors_derived_double_index, monitors_geog_nam, monitors_geog_conus, monitors_dharma, monitors_resp_curves)

	conf <- configureMCMC(
		model,
		monitors = monitors,
		print = TRUE,
		enableWAIC = TRUE
	)

	# add no U-turn sampler (Hamiltonian Monte Carlo)
	vars <- c('sigma_biomass_within_sites_log', 'beta_biomass_mu', 'beta_biomass_sigma', 'beta_biomass_pzero')
	conf$addSampler(target = vars, type = 'NUTS')
	say('NUTS sampler added to all continuous priors.')

	# vars <- c('sigma_biomass_within_sites_log', 'beta_biomass_mu', 'beta_biomass_sigma')
	# conf$removeSamplers(vars)
	# conf$addSampler(target = vars, type = 'AF_slice')
	# say('AF slice sampler added to ', paste(vars, collapse = ' & '), '.')

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
		WAIC = TRUE,
		perChainWAIC = FALSE
	)

	saveRDS(chains, paste0(out_dir, '/chains.rds'))

	say('session info', level = 2)
	print(sessionInfo())

	say(date(), pre = 1)
	sink()

######################################################
### post-modeling analysis for BIOMASS-only models ###
######################################################

	# burn spatial predictions into spatial vectors
	pred_vect_nam <- data_occs$ag_vect_sq

	for (var in monitors_geog_nam) {

		preds <- hammer_extract(chains, param = var, j = TRUE, stat = 'mean')
		pred_vect_nam[[var]] <- preds
		names(pred_vect_nam)[ncol(pred_vect_nam)] <- var

	}

	writeVector(pred_vect_nam, paste0(out_dir, '/prediction_vector_nam.gpkg'), overwrite = TRUE)

	pred_vect_conus <- data_occs$counties_thirties
	pred_vect_conus <- pred_vect_conus[ , c('country', 'state_province', 'county')]

	for (var in monitors_geog_conus) {

		preds <- hammer_extract(chains, param = var, j = TRUE, stat = 'mean')
		pred_vect_conus[[var]] <- preds
		names(pred_vect_conus)[ncol(pred_vect_conus)] <- var

	}

	writeVector(pred_vect_conus, paste0(out_dir, '/prediction_vector_conus.gpkg'), overwrite = TRUE)

	### post-modeling analysis of BIOMASS
	trait <- 'biomass'
	descrip <- paste0(trait, ': gamma~normal heteroscedastic')

	workflow_postmodeling_generic(facet = trait, formulae = formulae, descrip = descrip, out_dir = out_dir)
	workflow_postmodeling_biomass(descrip = descrip, formula_biomass_mu = formula_biomass_mu, homoscedastic = homoscedastic, zero_inflated = zero_inflated, crossvalidate = crossvalidate, pred_vect_nam = pred_vect_nam, pred_vect_conus = pred_vect_conus, out_dir = out_dir)

say(date())
say('FINIS!', deco = '+', level = 1)

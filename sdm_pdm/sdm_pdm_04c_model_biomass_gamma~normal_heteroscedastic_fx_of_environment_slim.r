### MODELING ANDROPOGON GERARDI DISTRIBUTION, PHENOTYPE, PHYSIOLOGY, GENOTYPE, and ASSOCIATED MICROBIAL COMMUNITIES
### Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2023-12
###
### This model estimates the distribution of site-level mean Andropogon gerardi ramet biomass. The distribution of biomasses among ramets at a site follows a gamma distribution with individual ramet biomasses drawn from this distribution. The site-level mean biomass is a function of one or more soil/climate predictors (possibly with higher-order terms).
###
### source('C:/Kaji/R/andropogon_integratedEcology/sdm_pdm/sdm_pdm_04c_model_biomass_gamma~normal_heteroscedastic_fx_of_environment_slim.r')

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
	# crossvalidate <- FALSE
	crossvalidate <- TRUE

	### formula for how aspects of species responds to environment

	formula_biomass_mu <- ~ 1 + bio12 # response of biomass to environment
	preds_filename <- 'bio12'

	# formula_biomass_mu <- ~ 1 + bio12 + I(bio12^2) # response of biomass to environment
	# preds_filename <- 'bio12^2'

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

	# formula_biomass_mu <- ~ 1 + bio1 + bio12 + I(bio1^2) + bio1:bio12 # response of biomass to environment
	# preds_filename <- 'bio1^2_x_bio12'

	out_dir <- paste0('./outputs_loretta/integrated_sdm_pdm/models_biomass/', ifelse(trial, 'TRIAL_', ''), '[biomass_gamma~normal_heteroscedastic_', preds_filename, ']/')
	previous_dir <- paste0('./outputs_loretta/integrated_sdm_pdm/models_biomass/[biomass_gamma~normal_homoscedastic_', preds_filename, ']/')

	# calib <- TRUE # use just counties with non-NA Poaceae for calibration region
	calib <- FALSE # use all of North America for calibration region

	formula_biomass_sigma <- formula_biomass_mu # SPECIFIC TO THIS SCRIPT--SHOULD NOT BE CHANGED
	formula_biomass_pzero <- NULL # SPECIFIC TO THIS SCRIPT--SHOULD NOT BE CHANGED

	if (!trial) {

		### MCMC settings
		# niter <- 1600000
		# nburnin <- niter / 2
		# thin <- 800
		# nchains <- 4
		# waic <- TRUE

		niter <- 400000
		nburnin <- 200000
		thin <- 200
		nchains <- 4

	} else {
		
		## MCMC settings FOR TESTING
		niter <- 1100
		nburnin <- 100
		thin <- 10
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

	say('This model estimates the distribution of site-level mean Andropogon gerardi ramet biomass. The distribution of biomasses among ramets at a site follows a gamma distribution with individual ramet biomasses drawn from this distribution. The site-level mean biomass is a function of or more soil/climate predictors (possibly with higher-order terms).', breaks = 60, post = 1)

	say('MCMC settings:', level = 2)
	say('trial ........................ ', trial)
	say('niter ........................ ', niter)
	say('nburnin ...................... ', nburnin)
	say('thin ......................... ', thin)
	say('nchains ...................... ', nchains)
	say('formula_biomass_mu ........... ', paste(as.character(formula_biomass_mu), collapse = ' '))
	say('formula_biomass_sigma ........ ', paste(as.character(formula_biomass_sigma), collapse = ' '))
	say('formula_biomass_pzero ........ ', paste(as.character(formula_biomass_pzero), collapse = ' '))
	say('homoscedastic ................ ', is.null(formula_biomass_sigma))
	say('zero_inflated ................ ', is.null(formula_biomass_pzero))
	say('calib ........................ ', calib, post = 2)

	say('out_dir:')
	say(out_dir)

	formulae <- list(
		formula_biomass_mu = formula_biomass_mu,
		formula_biomass_sigma = formula_biomass_sigma
	)

	########################s
	### data preparation ###
	########################

	data_biomass_mu <- prepare_biomass(formula_biomass = formula_biomass_mu, n_response_curve_values = n_response_curve_values, calib = calib)

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

		# response curves (general)
		n_response_curve_values = n_response_curve_values # number of values in response curve array

	)

	constants <- c(constants, constants_shared_biomass)

	# initializations
	previous_chains <- readRDS(paste0(previous_dir, '/chains.rds'))

	beta_biomass_mu_inits <- hammer_extract(previous_chains, 'beta_biomass_mu', j = TRUE, stat = 'mean')
	beta_biomass_sigma_inits <- rep(0, length(beta_biomass_mu_inits))

	sigma_biomass_within_sites_init <- hammer_extract(previous_chains, 'sigma_biomass_within_sites', stat = 'mean')
	sigma_biomass_within_sites_log_init <- log(sigma_biomass_within_sites_init)

	rm(previous_chains)

	inits <- list(

		site_biomass_mu_log = rep(2, data_biomass_mu$n_pheno_sites),
		sigma_biomass_within_sites_log = sigma_biomass_within_sites_log_init,
		y_biomass_sim = data_biomass_mu$y_biomass, # simulated values for biomass (for DHARMa residuals)
		beta_biomass_mu = beta_biomass_mu_inits,
		beta_biomass_sigma = beta_biomass_sigma_inits

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
			y_biomass[i] ~ dgamma(shape = shape_biomass[site_index_biomass[i]], rate = rate_biomass[site_index_biomass[i]])

			# simulated values for unconditional DHARMa residuals
			y_biomass_sim[i] ~ dgamma(shape = shape_biomass[site_index_biomass[i]], rate = rate_biomass[site_index_biomass[i]])

			log_lik_biomass[i] <- dgamma(y_biomass[i], shape = shape_biomass[site_index_biomass[i]], rate = rate_biomass[site_index_biomass[i]])
	
		}

		log_lik <- sum(log_lik_biomass[1:n_biomass])

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
		'beta_biomass_mu', 'beta_biomass_sigma'
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

	monitors_dharma <- c(
		'y_biomass_sim'
	)

	monitors_resp_curves <- c(
		'response_curves_biomass_mu', 'response_curves_biomass_sigma'
	)

	monitors <- c(monitors_coeffs_not_indexed, monitors_coeffs_single_index, monitors_coeffs_double_index, monitors_derived_not_indexed, monitors_derived_single_index, monitors_derived_double_index, monitors_dharma, monitors_resp_curves)

	conf <- configureMCMC(
		model,
		monitors = monitors,
		print = TRUE,
		enableWAIC = TRUE
	)

	# add no U-turn sampler (Hamiltonian Monte Carlo)
	vars <- c('sigma_biomass_within_sites_log', 'beta_biomass_mu', 'beta_biomass_sigma')
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


	### post-modeling analysis of BIOMASS
	trait <- 'biomass'
	descrip <- paste0(trait, ': gamma~normal heteroscedastic')

	workflow_postmodeling_generic(facet = trait, formulae = formulae, descrip = descrip, out_dir = out_dir)
	workflow_postmodeling_biomass(chains = chains, descrip = descrip, formula_biomass_mu = formula_biomass_mu, formula_biomass_sigma = formula_biomass_sigma, formula_biomass_pzero = formula_biomass_pzero, crossvalidate = crossvalidate, out_dir = out_dir)

say(date())
say('FINIS!', deco = '+', level = 1)

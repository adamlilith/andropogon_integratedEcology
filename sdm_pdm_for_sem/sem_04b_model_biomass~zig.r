### MODELING ANDROPOGON GERARDI DISTRIBUTION, PHENOTYPE, PHYSIOLOGY, GENOTYPE, and ASSOCIATED MICROBIAL COMMUNITIES
### Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2023-12
###
### This model estimates the distribution of site-level mean Andropogon gerardi ramet biomass. The distribution of biomasses among ramets at a site follows a zero-inflated gamma distribution with individual ramet biomasses drawn from this distribution. The site-level mean biomass is a function of soil/climate. The probability of a zero value is a function of the same covariate(s) and has the same functional form as the submodel of the mean.
###
### source('C:/Kaji/R/andropogon_integratedEcology/sdm_pdm_for_sem/sem_04b_model_biomass~zig.r')

#############
### setup ###
#############

	rm(list = ls())

	drive <- 'C:/Kaji/'

	setwd(paste0(drive, '/Research/Andropogon/Andropogon'))
	source(paste0(drive, '/R/andropogon_integratedEcology/sdm_pdm_for_sem/sem_00_shared_functions_and_variables.r'))

###########################
### user-defined values ###
###########################

	# trial <- TRUE # TRUE for testing
	trial <- FALSE # TRUE for testing

	# do cross-validation?
	# crossvalidate <- FALSE
	crossvalidate <- TRUE

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

	# formula_biomass_mu <- ~ 1 + bio12 + sand # response of biomass to environment
	# preds_filename <- 'bio12_sand'

	# formula_biomass_mu <- ~ 1 + bio12 + sand + I(bio12^2) # response of biomass to environment
	# preds_filename <- 'bio12^2_sand'

	# formula_biomass_mu <- ~ 1 + bio12 + sand + I(sand^2) # response of biomass to environment
	# preds_filename <- 'bio12_sand^2'

	# formula_biomass_mu <- ~ 1 + bio12 + sand + bio12:sand # response of biomass to environment
	# preds_filename <- 'bio12_x_sand'

	out_dir <- paste0('./outputs_loretta/sdm_pdm_for_sem/', ifelse(trial, 'TRIAL_', ''), 'model_biomass~zig~', preds_filename, '_wide_pzero_prior/')

	# calib <- TRUE # use just counties with non-NA Poaceae for calibration region
	calib <- FALSE # use all of North America for calibration region

	zero_inflated <- TRUE # SPECIFIC TO THIS SCRIPT--SHOULD NOT BE CHANGED

	if (!trial) {

		niter <- 400000
		nburnin <- 20000
		nchains <- 4

	} else {
		
		### MCMC settings FOR TESTING
		niter <- 1100
		nburnin <- 100
		nchains <- 2

	}
	thin <- (niter - nburnin) / 1000

	# do not change
	formula_pzero <- formula_biomass_mu

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

	say('This model estimates the distribution of site-level mean Andropogon gerardi ramet biomass. The distribution of biomasses among ramets at a site follows a zero-inflated gamma distribution with individual ramet biomasses drawn from this distribution. The site-level mean biomass is a function of soil/climate. The probability of a zero value is a function of the same covariate(s) and has the same functional form as the submodel of the mean.', breaks = 60, post = 1)

	say('MCMC settings:', level = 2)
	say('trial ........................ ', trial)
	say('niter ........................ ', niter)
	say('nburnin ...................... ', nburnin)
	say('thin ......................... ', thin)
	say('nchains ...................... ', nchains)
	say('formula_biomass_mu ........... ', paste(as.character(formula_biomass_mu), collapse = ' '))
	say('formula_pzero ................ ', paste(as.character(formula_pzero), collapse = ' '))
	say('zero_inflated ................ ', zero_inflated)
	say('calib ........................ ', calib, post = 2)

	say('out_dir:')
	say(out_dir)

	formulae <- list(
		formula_biomass_mu = formula_biomass_mu,
		formula_pzero = formula_pzero
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
		y = data_biomass_mu$y		# biomass of individual plants
	)

	constants <- list(
		
		# sampled sites
		n_pheno_sites = data_biomass_mu$n_pheno_sites, # number of phenotype sample sites

		### biomass
		x_by_site_biomass = data_biomass_mu$x_by_site_biomass, # MM with covariates for biomass (scaled)
		n_biomass = data_biomass_mu$n_biomass, # number of biomass observations
		site_index_biomass = data_biomass_mu$site_index_biomass, # index of sampled site for each row in biomass data
		n_terms_biomass_mu = data_biomass_mu$n_terms_biomass, # number of terms in formula for biomass model (including intercept)

		n_covariates_biomass = data_biomass_mu$n_covariates_biomass,
		resp_curves_x_biomass = data_biomass_mu$resp_curves_x_biomass, # response curve array for biomass

		# response curves (general)
		n_response_curve_values = n_response_curve_values # number of values in response curve array

	)

	constants <- c(constants, constants_shared_biomass)

	beta_biomass_mu_inits <- rep(0, constants$n_terms_biomass_mu)
	beta_pzero_inits <- rep(-1, constants$n_terms_biomass_mu)

	inits <- list(

		sigma_log = 1,
		y_sim = data_biomass_mu$y, # simulated values for biomass (for DHARMa residuals)
		beta_mu = beta_biomass_mu_inits,
		beta_pzero = beta_pzero_inits

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
			beta_mu[i] ~ dnorm(0, sd = beta_prior) # broad prior
		}

		# BIOMASS: priors for relationship of site-level zero values
		for (i in 1:n_terms_biomass_mu) {
			# beta_pzero[i] ~ dnorm(0, sd = beta_prior) # broad prior
			beta_pzero[i] ~ dnorm(0, sd = 200) # broad prior
		}

		# prior for sd of individual plant biomass on lognormal (~ half-Cauchy), ==> vague
		sigma_log ~ dnorm(0, sd = sigma_log_prior_sd)
		sigma <- exp(sigma_log)

		# BIOMASS: parameters of biomass distribution are latent and functions of environment
		# individual plant biomasses are samples from the site-level distribution (next chunk after this one)
		for (i in 1:n_pheno_sites) {

		  # site-level mean biomass
			mu[i] <- exp(lambda[i])

		  # relationship of biomass to the environment
			lambda[i] <- inprod(
				beta_mu[1:n_terms_biomass_mu],
				x_by_site_biomass[i, 1:n_terms_biomass_mu]
			)

			logit(pzero[i]) <- inprod(
				beta_pzero[1:n_terms_biomass_mu],
				x_by_site_biomass[i, 1:n_terms_biomass_mu]
			)

			# moment matching to get zigamma() parameters
			shape[i] <- mu[i]^2 / sigma^2
			rate[i] <- mu[i] / sigma^2

		}

		# BIOMASS: likelihood of individual plants at a site
		for (i in 1:n_biomass) {

			# likelihood
			y[i] ~ dzigamma(shape = shape[site_index_biomass[i]], rate = rate[site_index_biomass[i]], pzero = pzero[site_index_biomass[i]])

			# simulated values for unconditional DHARMa residuals
			y_sim[i] ~ dzigamma(shape = shape[site_index_biomass[i]], rate = rate[site_index_biomass[i]], pzero = pzero[site_index_biomass[i]])

		}

	})

	if (data_biomass_mu$n_covariates_biomass == 1) {

		### univariate
		response_curve_code <- nimbleCode({

			# BIOMASS: posterior predictive sampler for response curves: site-level mean
			# NB we assume just one predictor for biomass, so the response curve "x" is a matrix, not an array
			for (i in 1:n_response_curve_values) {
					
				log_response_curves_biomass_mu[i] <-
					inprod(beta_mu[1:n_terms_biomass_mu], resp_curves_x_biomass[i, 1:n_terms_biomass_mu])
				response_curves_mu[i] <- exp(log_response_curves_biomass_mu[i])

				logit(response_curves_pzero[i]) <-
					inprod(beta_pzero[1:n_terms_biomass_mu], resp_curves_x_biomass[i, 1:n_terms_biomass_mu])
				
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
						inprod(beta_mu[1:n_terms_biomass_mu], resp_curves_x_biomass[j, 1:n_terms_biomass_mu, i])
					response_curves_mu[j, i] <- exp(log_response_curves_biomass_mu[j, i])

					logit(response_curves_pzero[j, i]) <-
						inprod(beta_pzero[1:n_terms_biomass_mu], resp_curves_x_biomass[j, 1:n_terms_biomass_mu, i])
					
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
		'sigma'
	)

	# coefficients that have bracketed indexing
	monitors_coeffs_single_index <- c(
		'beta_mu', 'beta_pzero'
	)

	monitors_coeffs_double_index <- c(
	)

	monitors_derived_not_indexed <- c(
	)

	monitors_derived_single_index <- c(
		'mu'
	)
	monitors_derived_double_index <- c(
	)

	monitors_dharma <- c(
		'y_sim'
	)

	monitors_resp_curves <- c(
		'response_curves_mu', 'response_curves_pzero'
	)

	monitors <- c(monitors_coeffs_not_indexed, monitors_coeffs_single_index, monitors_coeffs_double_index, monitors_derived_not_indexed, monitors_derived_single_index, monitors_derived_double_index, monitors_dharma, monitors_resp_curves)

	conf <- configureMCMC(
		model,
		monitors = monitors,
		print = TRUE,
		enableWAIC = TRUE
	)

	# add no U-turn sampler (Hamiltonian Monte Carlo)
	vars <- c(monitors_coeffs_not_indexed, monitors_coeffs_single_index, monitors_coeffs_double_index)
	conf$addSampler(target = vars, type = 'NUTS')
	say('NUTS sampler added to all continuous parameters.')

	# vars <- c(monitors_coeffs_not_indexed, monitors_coeffs_single_index, monitors_coeffs_double_index)
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
	descrip <- paste0(trait, ': zero-inflated gamma')

	workflow_postmodeling_generic(facet = trait, formulae = formulae, descrip = descrip, out_dir = out_dir)
	workflow_postmodeling_biomass(chains = chains, descrip = descrip, formula_biomass_mu = formula_biomass_mu, formula_pzero = formula_pzero, crossvalidate = crossvalidate, out_dir = out_dir)

say(date())
say('FINIS!', deco = '+', level = 1)

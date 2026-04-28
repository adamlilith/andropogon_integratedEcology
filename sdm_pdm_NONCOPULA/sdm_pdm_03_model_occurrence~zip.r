### MODELING ANDROPOGON GERARDI DISTRIBUTION, PHENOTYPE, PHYSIOLOGY, GENOTYPE, and ASSOCIATED MICROBIAL COMMUNITIES
### Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2023-12
###
### This script constructs a species distribution model Andropogon gerardi where "county" is the observational unit. It assumes (latent) abundance follows a zero-inflated Poisson distribution, and the observed number of AG is a binomial distribution where the probability of observing AG is a function of county area and the number of Poaceae recorded in the county (including AG). The expected abundance is a log function of environmental predictors (climate, soil, etc.). The probability of (inflated) zero is a function of environmental covariates. The model is run using nimble. The model does not use a copula.
###
### source('C:/Kaji/R/andropogon_integratedEcology/sdm_pdm/sdm_pdm_03_model_occurrence~zip.r')
### 
#############
### setup ###
#############

	rm(list = ls())

	drive <- 'C:/Kaji/'

	setwd(paste0(drive, '/Research/Andropogon/Andropogon'))
	source('C:/Kaji/R/andropogon_integratedEcology/sdm_pdm/sdm_pdm_00_shared_functions_and_variables.r')

###########################
### user-defined values ###
###########################

	# trial <- TRUE # TRUE for testing
	trial <- FALSE # TRUE for testing

	# calib <- TRUE # use just counties with non-NA Poaceae for calibration region
	calib <- FALSE # use all of North America for calibration region

	# do cross-validation?
	do_crossvalidation <- TRUE
	# do_crossvalidation <- FALSE

	# use log of BIOs 12-14 and 16-19?
	# log_precip <- FALSE
	log_precip <- TRUE

	### formula for how aspects of species responds to environment

	formula_occs <- ~ 1 + bio1 + bio12 + bio15 + I(bio1^2) + I(bio12^2) + I(bio15^2) # response of occurrence to climate and soil
	preds_filename <- 'bio1^2_bio12^2_bio15^2'

	# formula_occs <- ~ 1 + bio1 + bio12 + bio15 + sand + I(bio1^2) + I(bio12^2) + I(bio15^2) + I(sand^2) # response of occurrence to climate and soil
	# preds_filename <- 'bio1^2_bio12^2_bio15^2_sand^2'

	# formula_psi <- ~ 1 + bio1 + bio12 + bio15 + sand + I(bio1^2) + I(bio12^2) + I(bio15^2) + I(sand^2)
	# zip_filename <- 'bio1^2_bio12^2_bio15^2_sand^2'

	# formula_psi <- ~ 1 + bio1 + bio12 + bio15 + ph + I(bio1^2) + I(bio12^2) + I(bio15^2) + I(ph^2)
	# zip_filename <- 'bio1^2_bio12^2_bio15^2_ph^2'

	formula_psi <- ~ 1 + bio1 + bio12 + bio15 + I(bio1^2) + I(bio12^2) + I(bio15^2)
	zip_filename <- 'bio1^2_bio12^2_bio15^2'

	formula_occs_bias <- ~ 1 + area_km2_log10 + n_poaceae_log10p1 # sampling bias for AG records
	bias_filename <- 'area_poaceae'
	
	# formula_occs_bias <- ~ 1 + area_km2_log10 # sampling bias for AG records
	# bias_filename <- 'area'
	
	# formula_occs_bias <- ~ 1 + area_km2_log10 # sampling bias for AG records
	# bias_filename <- 'poaceae'
	
	# formula_occs_bias <- ~ 1 # sampling bias for AG records
	# bias_filename <- '1'

	### output folder and bias formula
	out_dir <- paste0('./outputs_loretta/integrated_sdm_pdm/models_occurrence_noncopula/', ifelse(trial, 'TRIAL_', ''), '[occs_zip[psi~', zip_filename, ']=normal~', preds_filename, '_[bias~', bias_filename, ']]', ifelse(log_precip, '_log_precip', ''), '/')

	if (!trial) {

		### MCMC settings
		# need these settings to achieve independent samples using bio1^2, bio12^2, bio15^
		# niter <- 1600000
		niter <- 800000
		nchains <- 4

	} else {
		
		### MCMC settings FOR TESTING
		niter <- 2000
		nchains <- 2

	}
	nburnin <- niter / 2
	thin <- (niter - nburnin) / 1000

#############
### model ###
#############

	if (!trial) if (file.exists(out_dir)) stop('Output folder already exists.')
	dirCreate(out_dir)

	sink(paste0(out_dir, '/runtime_log.txt'), split = TRUE)
	say('MODELING OCCURRENCE')
	say('Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | ', date(), post = 1)

	### data collation
	##################
	say('data collation', post = 1)

	say('This script constructs a species distribution model Andropogon gerardi where "county" is the observational unit. It assumes (latent) abundance follows a zero-inflated Poisson distribution, and the observed number of AG is a binomial distribution where the probability of observing AG is a function of county area and the number of Poaceae recorded in the county (including AG). The expected abundance is a log function of environmental predictors (climate, soil, etc.). The probability of (inflated) zero is a function of environmental covariates. The model is run using nimble. The model does not use a copula.', breaks = 60, post = 1)

	say('MCMC settings:', level = 2)
	say('trial ........................ ', trial)
	say('niter ........................ ', niter)
	say('nburnin ...................... ', nburnin)
	say('thin ......................... ', thin)
	say('nchains ...................... ', nchains)
	say('log_precip ................... ', log_precip)
	say('formula_occs ................. ', paste(as.character(formula_occs), collapse = ' '))
	say('formula_psi .................. ', paste(as.character(formula_psi), collapse = ' '))
	say('formula_occs_bias ............ ', paste(as.character(formula_occs_bias), collapse = ' '))

	say('out_dir')
	say(out_dir, post = 2)

	formulae <- list(
		formula_occs = formula_occs,
		formula_psi = formula_psi,
		formula_occs_bias = formula_occs_bias
	)
	saveRDS(formula, paste0(out_dir, '/formulae.rds'))

	########################s
	### data preparation ###
	########################

	data_occs <- prepare_occurrence_data(formula_occs = formula_occs, formula_occs_bias = formula_occs_bias, log_precip = log_precip, n_response_curve_values = n_response_curve_values, psa_quant = psa_quant, calib = calib)

	data_occs_psi <- prepare_occurrence_data(formula_occs = formula_psi, formula_occs_bias = formula_occs_bias, log_precip = log_precip, n_response_curve_values = n_response_curve_values, psa_quant = psa_quant, calib = calib)

	#########################
	### inputs for nimble ###
	#########################

	say('Inputs:', level = 2)
	data <- list(
		y_n_ag = data_occs$y_n_ag				# number of AG observations in each county
	)

	constants <- list(
		
		# y_n_ag_min = data_occs$y_n_ag,

		### occurrences
		n_counties_occs_calib = data_occs$n_counties_occs_calib, # number of counties in calibration region
		n_terms_occs = data_occs$n_terms_occs, # number of terms in formula for occurrence model (including intercept)
		n_terms_psi = data_occs_psi$n_terms_occs, # number of terms in sampling bias model

		w_occs_bias = data_occs$w_occs_bias, # model matrix of sampling bias of AG observed occurrences
		n_terms_occs_bias = data_occs$n_terms_occs_bias, # number of terms in sampling bias model

		n_covariates_occs_bias = data_occs$n_covariates_occs_bias, # number of covariates in formula for occurrence model

		counties_x_occs_calib_sq = data_occs$counties_x_occs_sq,
		counties_x_occs_psi_calib_sq = data_occs_psi$counties_x_occs_sq

	)

	constants <- c(constants, constants_shared_occs, constants_shared_psi)

	N_inits_calib <- data_occs$y_n_ag * 2
	N_inits_all_counties <- 2 * (1 + data_occs$ag_vect_sq$n_andropogon_gerardi)

	# initial values for occ ~ f(env)

	prelim_data <- cbind(data_occs$w_occs_bias, data_occs$counties_x_occs_sq[ , 2:ncol(data_occs$counties_x_occs_sq)])
	prelim_model <- glm.fit(prelim_data, y = data_occs$y_n_ag, family = poisson(log))

	alpha_occs_inits <- prelim_model$coefficients[c('(Intercept)', data_occs$terms_occs_bias)]
	beta_occs_inits <- prelim_model$coefficients[c('(Intercept)', data_occs$terms_occs)]

	prelim_data <- cbind(data_occs_psi$w_occs_bias, data_occs_psi$counties_x_occs_sq[ , 2:ncol(data_occs_psi$counties_x_occs_sq)])
	y <- as.numeric(data_occs$y_n_ag > 0)
	prelim_model <- glm.fit(prelim_data, y = y, family = binomial())
	beta_occs_psi_inits <- prelim_model$coefficients[c('(Intercept)', data_occs$terms_occs)]

	inits <- list(

		lambda_sigma = 1,
		log_lambda_sigma = log(1),

		z_county = as.numeric(N_inits_calib > 0), # AG present?
		alpha_occs = alpha_occs_inits, # intercept, area, # of Poaceae
		beta_occs = beta_occs_inits, # occurrence ~ environment coefficients (including intercept)
		beta_psi = beta_occs_psi_inits, # psi ~ environment coefficients (including intercept)

		N = N_inits_calib # number of latent AG in calibration counties
		
	)

	say('Data:')
	print(str(data))

	say('Constants:', pre = 1)
	print(str(constants))

	say('Initializations:', pre = 1)
	print(str(inits))

	### define model
	say('nimbleCode():', level = 2)

	model_code <- nimbleCode({
	
		# OCCURRENCE: weakly regularized or regularized priors for relationship to environment
		# ddexp(): rate = 0.7675 sets |beta| < 3 90% of time and <5 97.8% of the time, and <10 99.95% of the time
		# ddexp(): rate = 0.9986 sets |beta| < 3 95% of time
		beta_occs[1] ~ dnorm(0, sd = beta_occs_prior_dnorm_sd_1)
		for (i in 2:n_terms_occs) {
			# beta_occs[i] ~ ddexp(0, rate = beta_occs_prior_ddexp_rate)
			beta_occs[i] ~ dnorm(0, sd = beta_occs_prior_dnorm_sd)
		}

		beta_psi[1] ~ dnorm(0, sd = beta_psi_prior_dnorm_sd_1)
		for (i in 2:n_terms_psi) {
			# beta_psi[i] ~ ddexp(0, rate = beta_occs_prior_ddexp_rate)
			beta_psi[i] ~ dnorm(0, sd = beta_psi_prior_dnorm_sd)
		}

		# OCCURRENCE: priors for sampling bias
		alpha_occs[1] ~ dnorm(0, sd = alpha_occs_prior_dnorm_sd_1)

		# OCCURRENCE: prior for spread of normal distribution of lambda
		log(lambda_sigma) ~ dnorm(0, sd = lambda_sigma_prior_sd) # half-Cauchy

		# OCCURRENCE: likelihood
		for (i in 1:n_counties_occs_calib) {
			
			### actual abundance (latent--unobserved)
			N[i] ~ dTruncPseudoZIP(lambda = lambda_mu_sq[i], z = z_county[i])

			### observed number of AG and sampling bias
			# logit(p[i]) <- inprod(alpha_occs[1:n_terms_occs_bias], w_occs_bias[i, 1:n_terms_occs_bias])
			y_n_ag[i] ~ dbinom(size = N[i], prob = p[i])

			# relationship between expected (latent) abundance and environment assuming NORMAL distribution
			log(lambda_mu_sq[i]) <- inprod(beta_occs[1:n_terms_occs], counties_x_occs_calib_sq[i, 1:n_terms_occs])

			# (inflated) probability of zero abundance
			logit(psi[i]) <- inprod(beta_psi[1:n_terms_psi], counties_x_occs_psi_calib_sq[i, 1:n_terms_psi])
			z_county[i] ~ dbern(psi[i])

		}


	})

	### NO bias covariate
	if (data_occs$n_covariates_occs_bias == 0) {

		bias_code <- nimbleCode({

			# PROBABILITY OF OBSERVATION
			for (i in 1:n_counties_occs_calib) {
				logit(p[i]) <- alpha_occs[1]
			}

		})

		model_code <- glueNimbleCode(model_code, bias_code)

	} else if (data_occs$n_covariates_occs_bias >= 1) {
	### ONE bias covariate

		bias_code <- nimbleCode({

			# PROBABILITY OF OBSERVATION
			for (i in 2:n_terms_occs_bias) {
				alpha_occs[i] ~ dnorm(0, sd = alpha_occs_prior_dnorm_sd)
			}

			for (i in 1:n_counties_occs_calib) {
				logit(p[i]) <- inprod(alpha_occs[1:n_terms_occs_bias], w_occs_bias[i, 1:n_terms_occs_bias])
			}

		})

		model_code <- glueNimbleCode(model_code, bias_code)

	}

	print(model_code)

	say('nimbleModel():', level = 2)
	model <- nimbleModel(
		code = model_code, # our model
		constants = constants, # constants
		data = data, # data
		inits = inits, # initialization values
		check = TRUE, # any errors?
		calculate = FALSE,
		# buildDerivs = TRUE # need for Hamiltonian Monte Carlo
		buildDerivs = FALSE # need for Hamiltonian Monte Carlo
	)

	say('initializeInfo() and $calculate():', level = 2)
	check_nodes(model)
	model$initializeInfo()
	calc <- model$calculate()
	say('model$calculate(): ', calc)
	if (is.na(calc) || is.infinite(calc)) stop('Impossible likelihood.')
	
	say('configureMCMC():', level = 2)

	monitors_coeffs_not_indexed <- 'lambda_sigma'
	monitors_coeffs_single_index <- c('beta_occs', 'beta_psi', 'alpha_occs')
	monitors_coeffs_double_index <- c()

	# monitors_derived_not_indexed <- c('log_lik')
	monitors_derived_not_indexed <- c()
	monitors_derived_single_index <- c() # c('N', 'p', 'z_county', 'psi')
	monitors_derived_double_index <- c()

	monitors <- c(monitors_coeffs_not_indexed, monitors_coeffs_single_index, monitors_coeffs_double_index, monitors_derived_not_indexed, monitors_derived_single_index, monitors_derived_double_index)

	conf <- configureMCMC(
		model,
		monitors = monitors,
		print = TRUE,
		enableWAIC = TRUE
	)

	# # add no U-turn sampler (Hamiltonian Monte Carlo)
	# vars <- c(monitors_coeffs_not_indexed, monitors_coeffs_single_index, monitors_coeffs_double_index)
	# conf$addSampler(target = vars, type = 'NUTS')
	# say('NUTS sampler added to ', paste(vars, collapse = ' & '), '.')

	# # RW block samplers for correlated parameters
	# conf$removeSamplers('beta_occs[1]')
	# conf$removeSamplers('beta_occs_vs_biomass')
	# conf$addSampler(target = c('beta_occs[1]', 'beta_occs_vs_biomass[1]', 'beta_occs_vs_biomass[2]'), type = 'RW_block')
	# say('RW_block sampler added to beta_occs[1] and beta_occs_vs_biomass[1:2].')

	# # AF slice sampler
	# vars <- c('alpha_occs', 'beta_occs')
	# if (!homoscedastic) vars <- c(vars, 'beta_occs_sigma')
	# if (zero_inflated) vars <- c(vars, 'beta_psi')
	# for (var in vars) {
	# 	conf$removeSamplers(var)
	# }
	# conf$addSampler(target = vars, type = 'AF_slice')
	# say('AF_slice sampler added to ', paste(vars, collapse = ' & '), '.')

	### compile/build/run model/save MCMC
	build <- buildMCMC(conf)

	say('Compiling ', date())
	compiled <- compileNimble(model, build, showCompilerOutput = FALSE)

	say('Sampling ', date())
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
	say(date())

	saveRDS(chains, paste0(out_dir, '/chains.rds'))

	say('session info', level = 2)
	print(sessionInfo())

	say(date(), pre = 1)
	sink()

say('#################################################')
say('### post-modeling diagnostics and predictions ###')
say('#################################################')

	descrip <- 'occurrence ~ ZIP = exp(normal(env)))'
	workflow_postmodeling_generic(facet = 'occurrence', formulae = formulae, descrip = descrip, out_dir = out_dir)

	workflow_postmodeling_occurrence(formula_occs = formula_occs, formula_psi = formula_psi, formula_occs_bias = formula_occs_bias, out_dir = out_dir)
	
	if (do_crossvalidation) workflow_postmodeling_occurrence_crossvalidation(formula_occs = formula_occs, formula_occs_bias = formula_occs_bias, formula_psi = formula_psi, constants = constants, out_dir = out_dir)

say(date())
say('FINIS!', deco = '+', level = 1)

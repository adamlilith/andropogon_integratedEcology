### MODELING ANDROPOGON GERARDI DISTRIBUTION, PHENOTYPE, PHYSIOLOGY, GENOTYPE, and ASSOCIATED MICROBIAL COMMUNITIES
### Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2023-12
###
### This script constructs a species distribution model Andropogon gerardi where "county" is the observational unit. It assumes (latent) abundance follows a hurdle (zero-inflated) Poisson distribution, and the observed number of AG is a binomial distribution where the probability of observing AG is a function of county area and the number of Poaceae recorded in the county (including AG). The expected abundance is a log function of environmental predictors (climate, soil, etc.). The probability of (inflated) zero is a function of environmental covariates. The model is run using nimble.
###
### source('C:/Kaji/R/andropogon_integratedEcology/sdm_pdm/sdm_pdm_03b_model_occurrence~hurdle_bias_covariates.r')
### 
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

	# calib <- TRUE # use just counties with non-NA Poaceae for calibration region
	calib <- FALSE # use all of North America for calibration region

	# do cross-validation?
	do_crossvalidation <- TRUE
	# do_crossvalidation <- FALSE

	### formula for how aspects of species responds to environment

	# formula_occs <- ~ 1 + bio1 + bio12_log10p1 + bio15 + I(bio1^2) + I(bio12_log10p1^2) + I(bio15^2)
	# filename_occs <- 'bio1^2_log(bio12)^2_bio15^2'

	formula_occs <- ~ 1 + bio1 + bio12_log10p1 + bio15 + I(bio1^2) + I(bio15^2)
	filename_occs <- 'bio1^2_log(bio12)_bio15^2'

	# formula_occs <- ~ 1 + bio1 + bio12 + bio15 + I(bio1^2) + I(bio12^2) + I(bio15^2)
	# filename_occs <- 'bio1^2_bio12^2_bio15^2'

	# formula_occs <- ~ 1 + bio1 + bio12_log10p1 + bio15 + sand + I(bio1^2) + I(bio12_log10p1^2) + I(bio15^2) + I(sand^2)
	# filename_occs <- 'bio1^2_log(bio12)^2_bio15^2_sand^2'

	# formula_occs <- ~ 1 + bio1 + bio12 + bio15 + sand + I(bio1^2) + I(bio12^2) + I(bio15^2) + I(sand^2)
	# filename_occs <- 'bio1^2_bio12^2_bio15^2_sand^2'

	# formula_bias <- ~ 1 + area_km2_log10 + n_poaceae_log10p1 # sampling bias for AG records
	# filename_bias <- 'area_poaceae'
	
	# formula_bias <- ~ 1 + area_km2_log10 # sampling bias for AG records
	# filename_bias <- 'area'
	
	# formula_bias <- ~ 1 + n_poaceae_log10p1 # sampling bias for AG records
	# filename_bias <- 'poaceae'
	
	formula_bias <- ~ 1 # sampling bias for AG records
	filename_bias <- '1'

	# do not change
	formula_psi <- formula_occs
	filename_psi <- filename_occs

	### output folder and bias formula
	out_dir <- paste0('./outputs_loretta/integrated_sdm_pdm/models_occurrence/', ifelse(trial, 'TRIAL_', ''), '[occs~hurdlepoisson(', filename_occs, ')]_[bias~logit(', filename_bias, ')]')

	if (!trial) {

		### MCMC settings
		niter <- 1200000 # need ~1M to achieve ESS ~1000 for beta_occs params when >1 bias parameter
		nchains <- 4

	} else {
		
		### MCMC settings FOR TESTING
		niter <- 2000
		nchains <- 1

	}
	nburnin <- niter / 2
	thin <- if ((niter - nburnin) / 1000 < 1) { 1 } else { (niter - nburnin) / 1000 }

	seed <- 1
	set.seed(seed)

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

	say('This script constructs a species distribution model Andropogon gerardi where "county" is the observational unit. It assumes (latent) abundance follows a hurdle (zero-inflated) Poisson distribution, and the observed number of AG is a binomial distribution where the probability of observing AG is a function of county area and the number of Poaceae recorded in the county (including AG). The expected abundance is a log function of environmental predictors (climate, soil, etc.). The probability of (inflated) zero is a function of environmental covariates. The model is run using nimble.', breaks = 60, post = 1)

	say('MCMC settings:', level = 2)
	say('trial ........................ ', trial)
	say('seed ......................... ', seed)
	say('niter ........................ ', niter)
	say('nburnin ...................... ', nburnin)
	say('thin ......................... ', thin)
	say('nchains ...................... ', nchains)
	say('formula_occs ................. ', paste(as.character(formula_occs), collapse = ' '))
	say('formula_psi .................. ', paste(as.character(formula_psi), collapse = ' '))
	say('formula_bias ................. ', paste(as.character(formula_bias), collapse = ' '), post = 2)

	say('out_dir')
	say(out_dir, post = 2)

	formulae <- list(
		formula_occs = formula_occs,
		formula_psi = formula_psi,
		formula_bias = formula_bias
	)
	saveRDS(formulae, paste0(out_dir, '/formulae.rds'))

	########################s
	### data preparation ###
	########################

	data_occs <- prepare_occurrence_data(formula_occs = formula_occs, formula_bias = formula_bias, n_response_curve_values = n_response_curve_values, psa_quant = psa_quant, calib = calib)

	data_psi <- prepare_occurrence_data(formula_occs = formula_psi, formula_bias = formula_bias, n_response_curve_values = n_response_curve_values, psa_quant = psa_quant, calib = calib)

	#########################
	### inputs for nimble ###
	#########################

	say('Inputs:', level = 2)
	data <- list(
		y_n_ag = data_occs$y_n_ag				# number of AG observations in each county
	)

	constants <- list(
		
		n_counties_occs_calib = data_occs$n_counties_occs_calib, # number of counties in calibration region
		n_terms_occs = data_occs$n_terms, # number of terms in formula for occurrence model (including intercept)
		n_terms_psi = data_psi$n_terms, # number of terms in sampling bias model

		w_bias = data_occs$w_bias, # model matrix of sampling bias of AG observed occurrences
		n_terms_bias = data_occs$n_terms_bias, # number of terms in sampling bias model

		counties_x_occs_calib_sq = data_occs$counties_x_sq,
		counties_x_psi_calib_sq = data_psi$counties_x_sq

	)

	constants <- c(constants, constants_shared_occs, constants_shared_psi)

	### initialization
	N_inits_calib <- data_occs$y_n_ag * 2

	prelim_data <- cbind(data_occs$w_bias, data_occs$counties_x_sq[ , 2:ncol(data_occs$counties_x_sq)])
	prelim_model <- glm.fit(prelim_data, y = data_occs$y_n_ag, family = poisson(log))

	alpha_occs_inits <- prelim_model$coefficients[c('(Intercept)', data_occs$terms_bias)]
	beta_occs_inits <- prelim_model$coefficients[c('(Intercept)', data_occs$terms)]

	prelim_data <- cbind(data_psi$w_bias, data_psi$counties_x_sq[ , 2:ncol(data_psi$counties_x_sq)])
	y <- as.numeric(data_occs$y_n_ag > 0)
	prelim_model <- glm.fit(prelim_data, y = y, family = binomial())
	beta_psi_inits <- prelim_model$coefficients[c('(Intercept)', data_occs$terms)]

	inits <- list(

		y_n_ag_sim = data_occs$y_n_ag, # simulated values of observed number of AG (for DHARMa residuals)
		alpha_occs = alpha_occs_inits, # intercept, area, # of Poaceae
		beta_occs = beta_occs_inits, # occurrence ~ environment coefficients (including intercept)
		beta_psi = beta_psi_inits, # psi ~ environment coefficients (including intercept)

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
	
		# OCCURRENCE: likelihood
		for (i in 1:n_counties_occs_calib) {
			
			# ### actual abundance (latent--unobserved)
			# N[i] ~ dHurdlePoisson(lambda = lambda[i], psi = psi[i])

			# ### observed number of AG and sampling bias
			# y_n_ag[i] ~ dbinom(size = N[i], prob = p[i])

			# joint Poisson-binomial likelihood
			y_n_ag[i] ~ dHurdlePoissonBinomial(lambda = lambda[i], psi = psi[i], p = p[i])

			# # simulate observations for DHARMa residuals
			# y_n_ag_sim[i] ~ dbinom(size = N[i], prob = p[i])
			y_n_ag_sim[i] ~ dHurdlePoissonBinomial(lambda = lambda[i], psi = psi[i], p = p[i])

			# relationship between expected (latent) abundance and environment assuming NORMAL distribution
			log(lambda[i]) <- 
				inprod(beta_occs[1:n_terms_occs], counties_x_occs_calib_sq[i, 1:n_terms_occs])

			# (inflated) probability of zero abundance
			logit(psi[i]) <- inprod(beta_psi[1:n_terms_psi], counties_x_psi_calib_sq[i, 1:n_terms_psi])

		}

		# log_lik <- sum(log_lik_y[1:n_counties_occs_calib])

	})

	## NO bias covariate
	if (data_occs$n_covariates_bias == 0) {

		bias_code <- nimbleCode({

			# OBSERVATION OF OCCURRENCE: likelihood
			logit(p[1:n_counties_occs_calib]) <- alpha_occs[1]
			# for (i in 1:n_counties_occs_calib) {
				# logit(p[i]) <- alpha_occs[1]
			# }

		})

	} else if (data_occs$n_covariates_bias >= 1) {

		bias_code <- nimbleCode({

			for (i in 2:n_terms_bias) {
				alpha_occs[i] ~ dnorm(0, sd = alpha_occs_prior_dnorm_sd)
			}

			# OBSERVATION OF OCCURRENCE: likelihood
			for (i in 1:n_counties_occs_calib) {
				logit(p[i]) <- inprod(alpha_occs[1:n_terms_bias], w_bias[i, 1:n_terms_bias])
			}

		})

	}

	model_code <- glueNimbleCode(
		model_code,
		model_code_beta_occs_alpha_occs_1_priors,
		bias_code,
		model_code_beta_psi_priors
	)

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
	model$initializeInfo()
	calc <- model$calculate()
	check_nodes(model)
	say('model$calculate(): ', calc)
	if (is.na(calc) || is.infinite(calc)) stop('Impossible likelihood.')
	
	say('configureMCMC():', level = 2)

	monitors_coeffs_not_indexed <- c()
	monitors_coeffs_single_index <- c('beta_occs', 'beta_psi', 'alpha_occs')
	monitors_coeffs_double_index <- c()

	monitors_derived_not_indexed <- c()
	monitors_derived_single_index <- c()
	monitors_derived_double_index <- c()

	monitors_dharma <- c(
		'y_n_ag_sim'
	)

	monitors_debug <- c(
		# 'N', 'lambda', 'psi', 'p'
	)

	monitors <- c(monitors_coeffs_not_indexed, monitors_coeffs_single_index, monitors_coeffs_double_index, monitors_derived_not_indexed, monitors_derived_single_index, monitors_derived_double_index, monitors_dharma, monitors_debug)

	conf <- configureMCMC(
		model,
		monitors = monitors,
		print = FALSE,
		enableWAIC = TRUE
	)

	# # add no U-turn sampler (Hamiltonian Monte Carlo)
	# vars <- c(monitors_coeffs_not_indexed, monitors_coeffs_single_index, monitors_coeffs_double_index)
	# conf$addSampler(target = vars, type = 'NUTS')
	# say('NUTS sampler added to ', paste(vars, collapse = ' & '), '.')

	# # AF slice sampler
	# vars <- c('beta_occs', 'beta_psi')
	# if (data_occs$n_covariates_bias >= 1) vars <- c(vars, 'alpha_occs')
	# for (var in vars) {
	# 	conf$removeSamplers(var) 
	# 	conf$addSampler(target = var, type = 'AF_slice')
	# 	say('AF_slice sampler added to ', var, '.')
	# }

	vars <- c('beta_occs', 'beta_psi', 'alpha_occs')
	for (var in vars) conf$removeSamplers(var)
	conf$addSampler(target = vars, type = 'AF_slice')
	say('AF_slice sampler added to ', paste(vars, collapse = ' & '), '.')

	print(conf)

	### compile/build/run model/save MCMC
	build <- buildMCMC(conf)

	say('Compiling ', date())
	compiled <- compileNimble(model, build, showCompilerOutput = TRUE)

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

	say('PRIORS', level = 1)

	say('constants_shared_occs', level = 2)
	print(constants_shared_occs)

	say('constants_shared_psi', level = 2)
	print(constants_shared_psi)

	say('session info', level = 1)
	print(sessionInfo())

	say(date(), pre = 1)
	sink()

say('#################################################')
say('### post-modeling diagnostics and predictions ###')
say('#################################################')

	descrip <- paste0('occs ~ hurdlePoisson(exp(', filename_occs, ')) | bias ~ ', filename_bias, ' | psi ~ ', filename_psi)
	workflow_postmodeling_generic(facet = 'occurrence', formulae = formulae, descrip = descrip, out_dir = out_dir)

	workflow_postmodeling_occurrence(formula_occs = formula_occs, formula_psi = formula_psi, formula_bias = formula_bias, out_dir = out_dir)
	
	if (do_crossvalidation) workflow_postmodeling_occurrence_crossvalidation(formula_occs = formula_occs, formula_bias = formula_bias, formula_psi = formula_psi, constants = constants, out_dir = out_dir)

say(date())
say('FINIS!', deco = '+', level = 1)

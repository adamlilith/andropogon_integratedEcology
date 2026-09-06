### MODELING ANDROPOGON GERARDI DISTRIBUTION, PHENOTYPE, PHYSIOLOGY, GENOTYPE, and ASSOCIATED MICROBIAL COMMUNITIES
### Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2023-12
###
### 'This script constructs a species distribution model Andropogon gerardi where "county" is the observational unit. It assumes relative abundance follows an overdispersed hurdle (zero-inflated) Poisson distribution, the expected value of which is an exponentiated draw from a normal distribution with a mean a function of environmental covariates and an offset (number of Poaceae specimens that are not A. gerardi). To obviate issues with taking the log of zero, 1 is added to the offset. However, this can create an upward bias in coefficients since some AG can be found at sites with no other species found (implying zero search effort, so the abundance of AG must be due to highly suitable climate). To correct for this, a binary dummy variable is used for cases where AG is >0 and number of other Poaceae is 0 (dummy variable = 1 in these cases). The interpretation is then that a site with no other Poaceae has an abundance of AG that is equivalent to exp(bias_coeff) non-AG plants. This correction can be paired with a tightly-regularized prior for the intercept. The probability of (inflated) zero is a function of environmental covariates. The model is run using nimble.'
###
### OR
###
### This script constructs a species distribution model Andropogon gerardi where "county" is the observational unit. It assumes relative abundance follows an overdispersed hurdle (zero-inflated) Poisson distribution, the expected value of which is an exponentiated draw from a normal distribution with a mean a function of environmental covariates and an offset (number of Poaceae specimens that are not A. gerardi). To obviate issues with taking the log of zero, 1 is added to the offset. However, this can create an upward bias in coefficients since some AG can be found at sites with no other species found (implying zero search effort, so the abundance of AG must be due to highly suitable climate). To correct for this, a tightly-regularized prior is used for the intercept. The probability of (inflated) zero is a function of environmental covariates. The model is run using nimble.
###
### OR
###
### This script constructs a species distribution model Andropogon gerardi where "county" is the observational unit. It assumes relative abundance follows an overdispersed hurdle (zero-inflated) Poisson distribution, the expected value of which is an exponentiated draw from a normal distribution with a mean a function of environmental covariates. The probability of (inflated) zero is a function of environmental covariates. The model is run using nimble.
###
### source('C:/Kaji/R/andropogon_integratedEcology/sdm_pdm/sdm_pdm_03d_model_occurrence~overdispersed_hurdle_offset.r')
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

	chain <- 4
	set.seed(chain)

	# trial <- TRUE # TRUE for testing
	trial <- FALSE # TRUE for testing

	# what kind of bias correction?
	# correction <- 'no_offset' # do not use Poaceae as offset
	correction <- 'offset_with_zero_correction' # use log(Poaceae) offset, and add 1/exp(1) + alpha * absent
	# correction <- 'offset_without_zero_correction' # use log(Poaceae) offset, but do not add zero correction

	# calib <- TRUE # use just counties with non-NA Poaceae for calibration region
	calib <- FALSE # use all of North America for calibration region

	# do cross-validation?
	# do_crossvalidation <- TRUE
	do_crossvalidation <- FALSE

	### formula for how aspects of species responds to environment

	formula_occs <- ~ 1 + bio1 + bio12_log10p1 + bio15 + I(bio1^2) + I(bio12_log10p1^2) + I(bio15^2)
	filename_occs <- 'bio1^2_log(bio12)^2_bio15^2'

	# formula_occs <- ~ 1 + bio1 + bio12 + bio15 + I(bio1^2) + I(bio12^2) + I(bio15^2)
	# filename_occs <- 'bio1^2_bio12^2_bio15^2'

	# do not change
	formula_psi <- formula_occs
	filename_psi <- filename_occs

	### output folder and bias formula
	if (correction == 'offset_with_zero_correction') {
		out_dir <- paste0('./outputs_loretta/integrated_sdm_pdm/models_occurrence/', ifelse(trial, 'TRIAL_', ''), '[occs~overdispersed_hurdlepoisson(', filename_occs, ')]_[bias~poaceae_offset_plus_exp_neg_1_WITH_zero_correction_covariate]_REDO_chain_', chain)
	} else if (correction == 'offset_without_zero_correction') {
		out_dir <- paste0('./outputs_loretta/integrated_sdm_pdm/models_occurrence/', ifelse(trial, 'TRIAL_', ''), '[occs~overdispersed_hurdlepoisson(', filename_occs, ')]_[bias~poaceae_offset_plus_exp_neg_1_SANS_zero_correction_covariate]')
	} else if (correction == 'no_offset') {
		out_dir <- paste0('./outputs_loretta/integrated_sdm_pdm/models_occurrence/', ifelse(trial, 'TRIAL_', ''), '[occs~overdispersed_hurdlepoisson(', filename_occs, ')]_[bias~no_offset]')
	}

	if (!trial) {

		### MCMC settings
		niter <- 3200000 # need ~1M to achieve ESS ~1000 for all coefficients
		# nchains <- 4
		nchains <- 1

	} else {
		
		### MCMC settings FOR TESTING
		niter <- 2000
		nchains <- 2

	}
	nburnin <- niter / 2
	thin <- if ((niter - nburnin) / 1000 < 1) { 1 } else { (niter - nburnin) / 1000 }

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

	if (correction == 'offset_with_zero_correction') {
		say('This script constructs a species distribution model Andropogon gerardi where "county" is the observational unit. It assumes relative abundance follows an overdispersed hurdle (zero-inflated) Poisson distribution, the expected value of which is an exponentiated draw from a normal distribution with a mean a function of environmental covariates and an offset (number of Poaceae specimens that are not A. gerardi). To obviate issues with taking the log of zero, 1 is added to the offset. However, this can create an upward bias in coefficients since some AG can be found at sites with no other species found (implying zero search effort, so the abundance of AG must be due to highly suitable climate). To correct for this, a binary dummy variable is used for cases where AG is >0 and number of other Poaceae is 0 (dummy variable = 1 in these cases). The interpretation is then that a site with no other Poaceae has an abundance of AG that is equivalent to exp(bias_coeff) non-AG plants. This correction can be paired with a tightly-regularized prior for the intercept. The probability of (inflated) zero is a function of environmental covariates. The model is run using nimble.', breaks = 120, post = 1)
	} else if (correction == 'offset_without_zero_correction') {
		say('This script constructs a species distribution model Andropogon gerardi where "county" is the observational unit. It assumes relative abundance follows an overdispersed hurdle (zero-inflated) Poisson distribution, the expected value of which is an exponentiated draw from a normal distribution with a mean a function of environmental covariates and an offset (number of Poaceae specimens that are not A. gerardi). To obviate issues with taking the log of zero, 1 is added to the offset. However, this can create an upward bias in coefficients since some AG can be found at sites with no other species found (implying zero search effort, so the abundance of AG must be due to highly suitable climate). To correct for this, a tightly-regularized prior is used for the intercept. The probability of (inflated) zero is a function of environmental covariates. The model is run using nimble.', breaks = 120, post = 1)
	} else if (correction == 'no_correction') {
		say('This script constructs a species distribution model Andropogon gerardi where "county" is the observational unit. It assumes relative abundance follows an overdispersed hurdle (zero-inflated) Poisson distribution, the expected value of which is an exponentiated draw from a normal distribution with a mean a function of environmental covariatesThe model is run using nimble.', breaks = 120, post = 1)
	}

	say('MCMC settings:', level = 2)
	say('trial ........................ ', trial)
	say('seed ......................... ', chain)
	say('chain ........................ ', chain)
	say('niter ........................ ', niter)
	say('nburnin ...................... ', nburnin)
	say('thin ......................... ', thin)
	say('nchains ...................... ', nchains)
	say('formula_occs ................. ', paste(as.character(formula_occs), collapse = ' '))
	say('formula_psi .................. ', paste(as.character(formula_psi), collapse = ' '))
	if (correction == 'offset_with_zero_correction') {
		say('bias_offset .................. poaceae PLUS 1/e with zero correction covariate')
	} else if (correction == 'offset_without_zero_correction') {
		say('bias_offset .................. poaceae PLUS 1/e WITHOUT zero correction covariate')
	} else if (correction == 'no_offset') {
		say('bias_offset .................. no offset')
	}

	say('out_dir')
	say(out_dir, post = 2)

	formulae <- list(
		formula_occs = formula_occs,
		formula_psi = formula_psi
	)
	saveRDS(formulae, paste0(out_dir, '/formulae.rds'))

	########################
	### data preparation ###
	########################

	data_occs <- prepare_occurrence_data(formula_occs = formula_occs, formula_bias = ~1, n_response_curve_values = n_response_curve_values, psa_quant = psa_quant, calib = calib)

	data_psi <- prepare_occurrence_data(formula_occs = formula_psi, formula_bias = ~1, n_response_curve_values = n_response_curve_values, psa_quant = psa_quant, calib = calib)

	#########################
	### inputs for nimble ###
	#########################

	say('Inputs:', level = 2)
	data <- list(
		y_n_ag = data_occs$y_n_ag				# number of AG observations in each county
	)

	bias_offset <- ceiling(data_occs$ag_vect_sq$n_poaceae) - data_occs$ag_vect_sq$n_andropogon_gerardi
	# log_bias_offset <- log1p(bias_offset)
	log_bias_offset <- log(bias_offset + 1 / exp(1))
	zero_bias_correction <- as.numeric(bias_offset == 0)

	constants <- list(
		
		n_counties_occs_calib = data_occs$n_counties_occs_calib, # number of counties in calibration region
		n_terms_occs = data_occs$n_terms, # number of terms in formula for occurrence model (including intercept)
		n_terms_psi = data_psi$n_terms, # number of terms in sampling bias model

		bias_offset = log_bias_offset,

		counties_x_occs_calib_sq = data_occs$counties_x_sq,
		counties_x_psi_calib_sq = data_psi$counties_x_sq

	)
	if (correction == 'offset_with_zero_correction') constants$zero_bias_correction <- zero_bias_correction

	constants <- c(constants, constants_shared_occs, constants_shared_psi)

	### initialization
	beta_occs_inits <- rep(0, data_occs$n_terms)
	beta_psi_inits <- rep(1, data_psi$n_terms)

	inits <- list(

		y_n_ag_sim = data_occs$y_n_ag, # simulated values of observed number of AG (for DHARMa residuals)
		beta_occs = beta_occs_inits, # occurrence ~ environment coefficients (including intercept)
		beta_psi = beta_psi_inits, # psi ~ environment coefficients (including intercept)

		lambda_star_prime = rep(0, data_occs$n_counties_occs_calib), # expected log abundance of AG with bias and overdispersion for each county

		log_sigma_occs = 1
	
	)
	if (correction == 'offset_with_zero_correction') inits$alpha_occs <- 0

	say('Data:')
	print(str(data))

	say('Constants:', pre = 1)
	print(str(constants))

	say('Initializations:', pre = 1)
	print(str(inits))

	### define model
	say('nimbleCode():', level = 2)

	model_code_base <- nimbleCode({

		# expected density
		lambda[1:n_counties_occs_calib] <- (counties_x_occs_calib_sq[1:n_counties_occs_calib, 1:n_terms_occs] %*% beta_occs[1:n_terms_occs])[ , 1]
		
		log(sigma_occs) ~ dnorm(0, sd = 2.5)

		# zero abundance component
		psi_star[1:n_counties_occs_calib] <- (counties_x_psi_calib_sq[1:n_counties_occs_calib, 1:n_terms_psi] %*% beta_psi[1:n_terms_psi])[ , 1]
		logit(psi[1:n_counties_occs_calib]) <- psi_star[1:n_counties_occs_calib]

		# OCCURRENCE: likelihood
		for (i in 1:n_counties_occs_calib) {
			
			lambda_star_prime[i] ~ dnorm(mean = lambda_star[i], sd = sigma_occs)
			exp_lambda_star_prime[i] <- exp(lambda_star_prime[i])

			### likelihood of observed number of AG
			y_n_ag[i] ~ dHurdlePoisson(lambda = exp_lambda_star_prime[i], psi = psi[i])

			# simulate observations for DHARMa residuals
			y_n_ag_sim[i] ~ dHurdlePoisson(lambda = exp_lambda_star_prime[i], psi = psi[i])

		}

	})

	if (correction == 'offset_with_zero_correction') {

		code_correction <- nimbleCode({

			# version 1: expected density with offset with correction for zero Poaceae
			lambda_star[1:n_counties_occs_calib] <-
				lambda[1:n_counties_occs_calib] +
				bias_offset[1:n_counties_occs_calib] +
				alpha_occs * zero_bias_correction[1:n_counties_occs_calib]
		
		})

	} else if (correction == 'offset_without_zero_correction') {

		code_correction <- nimbleCode({

			# version 2: expected density with offset with NO correction for zero Poaceae
			lambda_star[1:n_counties_occs_calib] <-
				lambda[1:n_counties_occs_calib] +
				bias_offset[1:n_counties_occs_calib]
		
		})


	} else if (correction == 'no_offset') {

		code_correction <- nimbleCode({

			# version 3: expected density with bias with NO offset
			lambda_star[1:n_counties_occs_calib] <-
				lambda[1:n_counties_occs_calib]
		
		})

	}

	if (correction == 'offset_with_zero_correction') {

		model_code <- glueNimbleCode(
			model_code_base,
			code_correction,
			model_code_beta_occs_priors,
			model_code_beta_psi_priors,
			model_code_alpha_occs_offset_correction_priors
		)

	} else {

		model_code <- glueNimbleCode(
			model_code_base,
			code_correction,
			model_code_beta_occs_priors,
			model_code_beta_psi_priors
		)
	}

	print(model_code)
	say('')

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

	if (correction == 'offset_with_zero_correction') {
		monitors_coeffs_not_indexed <- c('sigma_occs', 'alpha_occs')
	} else if (correction %in% c('offset_without_zero_correction', 'no_offset')) {
		monitors_coeffs_not_indexed <- c('sigma_occs')
	}
	monitors_coeffs_single_index <- c('beta_occs', 'beta_psi')
	monitors_coeffs_double_index <- c()

	monitors_derived_not_indexed <- c()
	monitors_derived_single_index <- c()
	monitors_derived_double_index <- c()

	monitors_dharma <- c(
		'y_n_ag_sim'
	)

	monitors_debug <- c(
		# 'lambda', 'psi'
	)

	monitors <- c(monitors_coeffs_not_indexed, monitors_coeffs_single_index, monitors_coeffs_double_index, monitors_derived_not_indexed, monitors_derived_single_index, monitors_derived_double_index, monitors_dharma, monitors_debug)

	conf <- configureMCMC(
		model,
		monitors = monitors,
		print = FALSE,
		enableWAIC = TRUE
	)

	vars <- c('beta_occs', 'beta_psi')
	for (var in vars) conf$removeSamplers(var)
	conf$addSampler(target = vars, type = 'AF_slice')
	say('Added AF_slice sampler to ', paste(vars, collapse = ' & '), '.')

	print(conf)

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

	say('PRIORS', level = 1)

	say('constants_shared_occs', level = 2)
	print(constants_shared_occs)

	say('constants_shared_psi', level = 2)
	print(constants_shared_psi)

	say('session info', level = 1)
	print(sessionInfo())

	say(date(), pre = 1)
	sink()

# say('#################################################')
# say('### post-modeling diagnostics and predictions ###')
# say('#################################################')

# 	descrip <- paste0('occs ~ OD hurdlePoisson(exp(', filename_occs, ')) | bias ~ ', ifelse(correction == 'offset_with_zero_correction', 'offset_with_zero_correction', ifelse(correction == 'offset_without_zero_correction', 'offset_without_zero_correction', 'no_offset')), ' | psi ~ ', filename_psi)
# 	workflow_postmodeling_generic(facet = 'occurrence', formulae = formulae, descrip = descrip, out_dir = out_dir)

# 	workflow_postmodeling_occurrence(
# 		formula_occs = formula_occs,
# 		formula_psi = formula_psi,
# 		formula_bias = ~ 1,
# 		out_dir = out_dir,
# 		overdispersed = TRUE
# 	)
	
# 	if (do_crossvalidation) workflow_postmodeling_occurrence_crossvalidation(formula_occs = formula_occs, formula_bias = ~ 1, formula_psi = formula_psi, constants = constants, out_dir = out_dir)

say(date())
say('FINIS!', deco = '+', level = 1)

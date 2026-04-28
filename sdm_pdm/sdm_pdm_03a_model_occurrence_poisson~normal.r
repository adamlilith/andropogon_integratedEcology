### MODELING ANDROPOGON GERARDI DISTRIBUTION, PHENOTYPE, PHYSIOLOGY, GENOTYPE, and ASSOCIATED MICROBIAL COMMUNITIES
### Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2023-12
###
### This script constructs a species distribution model Andropogon gerardi where "county" is the observational unit. It assumes (latent) abundance follows a Poisson distribution, and the observed number of AG is a binomial distribution where the probability of observing AG is a function of county area and the number of Poaceae recorded in the county (including AG). The expected abundance drawn from a normal distribution where the mean value is given by a function of environmental predictors (climate, soil, etc.). The model is run using nimble.
###
### source('C:/Kaji/R/andropogon_integratedEcology/sdm_pdm/sdm_pdm_03a_model_occurrence_poisson~normal.r')
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

	# formula_occs <- ~ 1 + bio1 + bio12 + bio15 + I(bio1^2) +  I(bio12^2) + I(bio15^2) # response of occurrence to climate and soil
	# preds_filename <- 'bio1^2_bio12^2_bio15^2'

	formula_occs <- ~ 1 + bio1 + bio12 + bio15 + I(bio1^2) + sand + ph + I(bio12^2) + I(bio15^2) + I(sand^2)# response of occurrence to climate and soil
	preds_filename <- 'bio1^2_bio12^2_bio15^2_sand^2'

	# formula_occs <- ~ 1 + bio1 + bio12 + bio15 + I(bio1^2) + sand + ph + I(bio12^2) + I(bio15^2) + I(ph^2)# response of occurrence to climate and soil
	# preds_filename <- 'bio1^2_bio12^2_bio15^2_ph^2'

	# formula_bias <- ~ 1 + n_poaceae_log10p1 # sampling bias for AG records
	# bias_filename <- 'poaceae'

	# formula_bias <- ~ 1 + area_km2_log10 # sampling bias for AG records
	# bias_filename <- 'area'

	# formula_bias <- ~ 1 + area_km2_log10 + n_poaceae_log10p1 # sampling bias for AG records
	# bias_filename <- 'area_poaceae'

	formula_bias <- ~ 1 # sampling bias for AG records
	bias_filename <- '1'

	### output folder and bias formula
	out_dir <- paste0('./outputs_loretta/integrated_sdm_pdm/models_occurrence/', ifelse(trial, 'TRIAL_', ''), '[occs_poisson~normal_homoscedastic~', preds_filename, '_[bias~', bias_filename, ']]/')

	if (!trial) {

		### MCMC settings
		# need these settings to achieve independent samples using bio1^2, bio12^2, bio15^
		niter <- 1600000
		nburnin <- niter / 2
		thin <- 800
		nchains <- 4
		waic <- TRUE

	} else {
		
		### MCMC settings FOR TESTING
		niter <- 1100
		nburnin <- 100
		thin <- 1
		nchains <- 2

	}

	# DO NOT CHANGE--SPECIFIC TO THIS SCRIPT
	formula_psi <- NULL

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

	say('This script constructs a species distribution model Andropogon gerardi where "county" is the observational unit. It assumes (latent) abundance follows a Poisson distribution, and the observed number of AG is a binomial distribution where the probability of observing AG is a function of county area and the number of Poaceae recorded in the county (including AG). The expected abundance drawn from a normal distribution where the mean value is given by a function of environmental predictors (climate, soil, etc.). The model is run using nimble.', breaks = 60, post = 1)

	say('MCMC settings:', level = 2)
	say('trial ........................ ', trial)
	say('niter ........................ ', niter)
	say('nburnin ...................... ', nburnin)
	say('thin ......................... ', thin)
	say('nchains ...................... ', nchains)
	say('formula_occs ................. ', paste(as.character(formula_occs), collapse = ' '))
	say('formula_psi ................ ', paste(as.character(formula_psi), collapse = ' '))
	say('formula_bias ............ ', paste(as.character(formula_bias), collapse = ' '))

	say('out_dir')
	say(out_dir, post = 2)

	formulae <- list(
		formula_occs = formula_occs,
		formula_psi = formula_psi,
		formula_bias = formula_bias
	)
	saveRDS(formula, paste0(out_dir, '/formulae.rds'))

	########################s
	### data preparation ###
	########################

	data_occs <- prepare_occurrence_data(formula_occs = formula_occs, formula_bias = formula_bias, n_response_curve_values = n_response_curve_values, psa_quant = psa_quant, calib = calib)

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
		n_terms_occs = data_occs$n_terms, # number of terms in formula for occurrence model (including intercept)
		w_occs_bias = data_occs$w_bias, # model matrix of sampling bias of AG observed occurrences
		n_terms_occs_bias = data_occs$n_terms_bias, # number of terms in sampling bias model

		n_covariates_occs = data_occs$n_covariates_occs, # number of covariates in formula for occurrence model
		n_covariates_occs_bias = data_occs$n_covariates_bias, # number of covariates in formula for occurrence model
		resp_curves_x_occs = data_occs$resp_curves_x, # response curve array for occurrences vs environment
		resp_curves_w_occs = data_occs$resp_curves_w, # response curve array for occurrences vs environment

		counties_x_occs_calib_sq = data_occs$counties_x_sq,

		# response curves (general)
		n_response_curve_values = n_response_curve_values # number of values in response curve array

	)

	constants <- c(constants, constants_shared_occs)

	N_inits_calib <- data_occs$y_n_ag * 2
	N_inits_all_counties <- 2 * (1 + data_occs$ag_vect_sq$n_andropogon_gerardi)

	prelim_data <- cbind(data_occs$w_bias, data_occs$counties_x_sq[ , 2:ncol(data_occs$counties_x_sq)])
	prelim_model <- glm.fit(prelim_data, y = data_occs$y_n_ag, family = poisson(log))

	alpha_occs_inits <- prelim_model$coefficients[c('(Intercept)', data_occs$terms_bias)]
	beta_occs_inits <- prelim_model$coefficients[c('(Intercept)', data_occs$terms)]

	response_curves_occs_mu_inits <- matrix(2, nrow = n_response_curve_values, ncol = data_occs$n_covariates_occs)

	inits <- list(

		lambda_sigma = 1,
		log_lambda_sigma = log(1),

		y_n_ag_sim = data_occs$y_n_ag, # simulated values of observed number of AG (for DHARMa residuals)
		log_lambda_mu_sq = rep(1, data_occs$n_counties_occs_calib), # expected value of number of AG
		alpha_occs = alpha_occs_inits, # intercept, area, # of Poaceae
		beta_occs = beta_occs_inits, # occurrence ~ environment coefficients (including intercept)

		N = N_inits_calib, # number of latent AG in calibration counties

		log_lambda_resp_curves_mu = response_curves_occs_mu_inits,
		response_curves_occs_mu = response_curves_occs_mu_inits
		
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
	
		# OCCURRENCE: prior for spread of normal distribution of lambda
		log(lambda_sigma) ~ dnorm(0, sd = lambda_sigma_prior_sd) # half-Cauchy

		# OCCURRENCE: likelihood
		for (i in 1:n_counties_occs_calib) {
			
			### actual abundance (latent--unobserved)
			N[i] ~ dpois(lambda_mu_sq[i])

			### observed number of AG and sampling bias
			# logit(p[i]) <- inprod(alpha_occs[1:n_terms_occs_bias], w_occs_bias[i, 1:n_terms_occs_bias])
			y_n_ag[i] ~ dbinom(prob = p[i], size = N[i])

			# simulate observations for DHARMa residuals
			y_n_ag_sim[i] ~ dbinom(prob = p[i], size = N[i])

			# # relationship between expected (latent) abundance and environment assuming NORMAL distribution
			# log(lambda_mu_sq[i]) ~ dnorm(phi_occs[i], sd = lambda_sigma)
			# phi_occs[i] <- inprod(beta_occs[1:n_terms_occs], counties_x_occs_calib_sq[i, 1:n_terms_occs])

			# relationship between expected (latent) abundance and environment
			log(lambda_mu_sq[i]) <- inprod(beta_occs[1:n_terms_occs], counties_x_occs_calib_sq[i, 1:n_terms_occs])

			# # likelihood
			# log_lik_y[i] <- dbinom(y_n_ag[i], prob = p[i], size = N[i], log = 1)

		}

		# log_lik <- sum(log_lik_y[1:n_counties_occs_calib])

		# OCCURRENCE: posterior predictive sampler for ENVIRONMENTAL response curves
		# We're assuming occurrence responds to two or more environmental predictors, so the response curve "x" is an array with one "page" per predictor and output is a matrix with one column per predictor
		for (i in 1:n_covariates_occs) {

			for (j in 1:n_response_curve_values) {
				
				response_curves_occs_mu[j, i] ~ dpois(lambda_resp_curves_mu[j, i])
				log(lambda_resp_curves_mu[j, i]) ~ dnorm(phi_lambda_resp_curves_mu[j, i], sd = lambda_sigma)
				phi_lambda_resp_curves_mu[j, i] <-
					inprod(beta_occs[1:n_terms_occs], resp_curves_x_occs[j, 1:n_terms_occs, i])

			}

		}

	})

	### NO bias covariate
	if (data_occs$n_covariates_bias == 0) {

		bias_response_curve_code <- nimbleCode({

			# OCCURRENCE: likelihood
			for (i in 1:n_counties_occs_calib) {
				logit(p[i]) <- alpha_occs[1]
			}

			# OCCURRENCE: posterior predictive sampler for BIAS response curves
			logit(response_curves_occs_bias) <- alpha_occs[1]

		})

		model_code <- glueNimbleCode(model_code, bias_response_curve_code)

	} else if (data_occs$n_covariates_bias == 1) {
	### ONE bias covariate

		bias_response_curve_code <- nimbleCode({

			for (i in 2:n_terms_occs_bias) {
				alpha_occs[i] ~ dnorm(0, sd = alpha_occs_prior_dnorm_sd)
			}

			# OCCURRENCE: likelihood
			for (i in 1:n_counties_occs_calib) {
				logit(p[i]) <- inprod(alpha_occs[1:n_terms_occs_bias], w_occs_bias[i, 1:n_terms_occs_bias])
			}

			# OCCURRENCE: posterior predictive sampler for BIAS response curves
			for (j in 1:n_response_curve_values) {
				
				logit(response_curves_occs_bias[j]) <-
					inprod(alpha_occs[1:n_terms_occs_bias], resp_curves_w_occs[j, 1:n_terms_occs_bias])

			}

		})

		model_code <- glueNimbleCode(model_code, bias_response_curve_code)

	### MORE THAN ONE bias covariate
	} else if (data_occs$n_covariates_bias > 1) {

		bias_response_curve_code <- nimbleCode({

			for (i in 2:n_terms_occs_bias) {
				alpha_occs[i] ~ dnorm(0, sd = alpha_occs_prior_dnorm_sd)
			}

			# OCCURRENCE: likelihood
			for (i in 1:n_counties_occs_calib) {
				logit(p[i]) <- inprod(alpha_occs[1:n_terms_occs_bias], w_occs_bias[i, 1:n_terms_occs_bias])
			}

			# OCCURRENCE: posterior predictive sampler for BIAS response curves
			for (i in 1:n_covariates_occs_bias) {
				
				for (j in 1:n_response_curve_values) {
					
					logit(response_curves_occs_bias[j, i]) <-
						inprod(alpha_occs[1:n_terms_occs_bias], resp_curves_w_occs[j, 1:n_terms_occs_bias, i])

				}

			}

		})

		model_code <- glueNimbleCode(model_code, bias_response_curve_code)

	}

	model_code <- glueNimbleCode(model_code, model_code_beta_occs_alpha_occs_1_priors)

	print(model_code)

	say('nimbleModel():', level = 2)
	model <- nimbleModel(
		code = model_code, # our model
		constants = constants, # constants
		data = data, # data
		inits = inits, # initialization values
		check = TRUE, # any errors?
		calculate = FALSE,
		buildDerivs = TRUE # need for NUTS
		# buildDerivs = FALSE
	)

	say('initializeInfo() and $calculate():', level = 2)
	model$initializeInfo()
	calc <- model$calculate()
	say('model$calculate(): ', calc)
	if (is.na(calc) || is.infinite(calc)) stop('Impossible likelihood.')

	say('configureMCMC():', level = 2)

	monitors_coeffs_not_indexed <- 'lambda_sigma'
	monitors_coeffs_single_index <- c('beta_occs', 'alpha_occs')
	monitors_coeffs_double_index <- c()

	monitors_derived_not_indexed <- c()
	monitors_derived_single_index <- c()
	monitors_derived_double_index <- c()

	monitors_dharma <- c(
		'y_n_ag_sim'
	)

	monitors_resp_curves <- c(
		'response_curves_occs_mu'
	)
	if (data_occs$n_covariates_bias > 0) monitors_resp_curves <- c(monitors_resp_curves, 'response_curves_occs_bias')

	monitors <- c(monitors_coeffs_not_indexed, monitors_coeffs_single_index, monitors_coeffs_double_index, monitors_derived_not_indexed, monitors_derived_single_index, monitors_derived_double_index, monitors_dharma, monitors_resp_curves)

	conf <- configureMCMC(
		model,
		monitors = monitors,
		print = TRUE,
		enableWAIC = TRUE
	)

	# add no U-turn sampler (Hamiltonian Monte Carlo)
	vars <- c('log_lambda_sigma', 'alpha_occs')
	for (var in vars) {
		conf$removeSamplers(var)
		conf$addSampler(target = var, type = 'NUTS')
		say('NUTS sampler added to ', var, '.')
	}

	# AF slice sampler
	var <- 'beta_occs'
	conf$removeSamplers(var)
	conf$addSampler(target = var, type = 'AF_slice')
	say('AF_slice sampler added to ', paste(vars, collapse = ' & '), '.')

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

	descrip <- 'occurrence: Poisson ~ normal homoscedastic'
	workflow_postmodeling_generic(facet = 'occurrence', formulae = formulae, descrip = descrip, out_dir = out_dir)
	
	workflow_postmodeling_occurrence(formula_occs = formula_occs, formula_psi = formula_psi, formula_bias = formula_bias, out_dir = out_dir)
	
	if (do_crossvalidation) workflow_postmodeling_occurrence_crossvalidation(formula_occs = formula_occs, formula_bias = formula_bias, formula_occs_sigma = formula_occs_sigma, formula_psi = formula_psi, constants = constants, out_dir = out_dir)


say(date())
say('FINIS!', deco = '+', level = 1)

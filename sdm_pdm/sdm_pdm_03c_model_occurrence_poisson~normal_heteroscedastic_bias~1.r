### MODELING ANDROPOGON GERARDI DISTRIBUTION, PHENOTYPE, PHYSIOLOGY, GENOTYPE, and ASSOCIATED MICROBIAL COMMUNITIES
### Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2023-12
###
### This script constructs a species distribution model Andropogon gerardi where "county" is the observational unit. It assumes (latent) abundance follows a Poisson distribution, and the observed number of AG is a binomial distribution where the probability of observing AG is an estimated constant. The expected abundance drawn from a normal distribution where the mean value and the standard deviation are functions of environmental predictors (climate, soil, etc.). The model is run using nimble.
###
### source('C:/Kaji/R/andropogon_integratedEcology/sdm_pdm/sdm_pdm_03c_model_occurrence_poisson~normal_heteroscedastic_bias~1.r')
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
	# do_crossvalidation <- TRUE
	do_crossvalidation <- FALSE

	homoscedastic <- FALSE # SPECIFIC TO THIS SCRIPT--SHOULD NOT BE CHANGED
	zero_inflated <- FALSE # SPECIFIC TO THIS SCRIPT--SHOULD NOT BE CHANGED

	### formula for how aspects of species responds to environment

	# formula_occs <- ~ 1 + bio1 + bio12 + bio15 + ph + I(bio1^2) + I(bio12^2) + I(bio15^2) + I(ph^2) # response of occurrence to climate and soil
	# preds_filename <- 'bio1^2_bio12^2_bio15^2_ph^2'

	# formula_occs <- ~ 1 + bio1 + bio12 + bio15 + sand + I(bio1^2) + I(bio12^2) + I(bio15^2) + I(sand^2) # response of occurrence to climate and soil
	# preds_filename <- 'bio1^2_bio12^2_bio15^2_sand^2'

	formula_occs <- ~ 1 + bio1 + bio12 + bio15 + I(bio1^2) + I(bio12^2) + I(bio15^2) # response of occurrence to climate and soil
	preds_filename <- 'bio1^2_bio12^2_bio15^2'

	### output folder and bias formula
	out_dir <- paste0('./outputs_loretta/integrated_sdm_pdm/models_occurrence/[occs_poisson~normal_heteroscedastic_', preds_filename, '_[bias~1]]', ifelse(trial, '_TRIAL', ''), '/')
	previous_dir <- paste0('./outputs_loretta/integrated_sdm_pdm/models_occurrence/[occs_poisson~normal_homoscedastic_', preds_filename, '_[bias~1]]', ifelse(trial, '_TRIAL', ''), '/')
	formula_occs_bias <- ~ 1 # sampling bias for AG records

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
		nchains <- 4

	}

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

	say('This script constructs a species distribution model Andropogon gerardi where "county" is the observational unit. It assumes (latent) abundance follows a Poisson distribution, and the observed number of AG is a binomial distribution where the probability of observing AG is an estimated constant. The expected abundance drawn from a normal distribution where the mean value and the standard deviation are functions of environmental predictors (climate, soil, etc.). The model is run using nimble.', breaks = 60, post = 1)

	say('MCMC settings:', level = 2)
	say('trial ........................ ', trial)
	say('niter ........................ ', niter)
	say('nburnin ...................... ', nburnin)
	say('thin ......................... ', thin)
	say('nchains ...................... ', nchains)
	say('formula_occs ................. ', paste(as.character(formula_occs), collapse = ' '))
	say('formula_occs_bias ............ ', paste(as.character(formula_occs_bias), collapse = ' '))
	say('homoscedastic ................ ', homoscedastic)
	say('zero_inflated ................ ', zero_inflated, post = 2)

	say('out_dir')
	say(out_dir, post = 2)

	formulae <- list(
		formula_occs = formula_occs,
		formula_occs_bias = formula_occs_bias
	)
	saveRDS(formula, paste0(out_dir, '/formulae.rds'))

	########################s
	### data preparation ###
	########################

	data_occs <- prepare_occurrences(formula_occs = formula_occs, formula_occs_bias = formula_occs_bias, n_response_curve_values = n_response_curve_values, psa_quant = psa_quant, calib = calib)

	# needed for plots
	data_biomass_mu <- prepare_biomass(formula = ~ 1, n_response_curve_values = n_response_curve_values, calib = calib)
	data_traits <- prepare_nonbiomass_traits(trait = 'height', formula = ~ 1, n_response_curve_values = n_response_curve_values, calib = calib)

	#########################
	### inputs for nimble ###
	#########################

	say('Inputs:', level = 2)
	data <- list(
		y_n_ag = data_occs$y_n_ag				# number of AG observations in each county
	)

	constants <- list(
		
		### occurrences
		n_counties_occs_calib = data_occs$n_counties_occs_calib, # number of counties in calibration region
		n_terms_occs = data_occs$n_terms_occs, # number of terms in formula for occurrence model (including intercept)

		n_covariates_occs = data_occs$n_covariates_occs, # number of covariates in formula for occurrence model
		resp_curves_x_occs = data_occs$resp_curves_x_occs, # response curve array for occurrences vs environment

		counties_x_occs_sq = data_occs$counties_x_occs_sq,
		counties_x_occs_calib_sq = data_occs$counties_x_occs_sq,
		counties_x_occs_ssp245_2041_2070 = data_occs$counties_x_occs_ssp245_2041_2070,
		counties_x_occs_ssp245_2071_2100 = data_occs$counties_x_occs_ssp245_2071_2100,
		counties_x_occs_ssp370_2041_2070 = data_occs$counties_x_occs_ssp370_2041_2070,
		counties_x_occs_ssp370_2071_2100 = data_occs$counties_x_occs_ssp370_2071_2100,

		counties_x_occs_thirties = data_occs$counties_x_occs_thirties,
		counties_x_occs_fifties = data_occs$counties_x_occs_fifties,

		# counties
		n_counties = data_occs$n_counties, # number of counties in the dataset
		n_counties_20th_cent = data_occs$n_counties_20th_cent, # number of counties in the dataset

		# response curves (general)
		n_response_curve_values = n_response_curve_values # number of values in response curve array

	)

	constants <- c(constants, constants_shared_occs)

	N_inits_calib <- data_occs$y_n_ag * 2
	N_inits_all_counties <- data_occs$ag_vect_sq$n_andropogon_gerardi * 2

	previous_chains <- readRDS(paste0(previous_dir, '/chains.rds'))
	alpha_occs_inits <- hammer_extract(previous_chains, 'alpha_occs', stat = 'mean')
	beta_occs_mu_inits <- hammer_extract(previous_chains, 'beta_occs_mu', j = TRUE, stat = 'mean')
	beta_occs_sigma_inits <- hammer_extract(previous_chains, 'beta_occs_mu', j = TRUE, stat = 'mean')
	rm(previous_chains)

	response_curves_occs_mu_inits <- matrix(2, nrow = n_response_curve_values, ncol = data_occs$n_covariates_occs)

	inits <- list(

		y_n_ag_sim = data_occs$y_n_ag, # simulated values of observed number of AG (for DHARMa residuals)
		log_lambda_mu_sq = rep(1, data_occs$n_counties_occs_calib), # expected value of number of AG
		alpha_occs = -3, # intercept, area, # of Poaceae
		beta_occs_mu = beta_occs_mu_inits, # occurrence ~ environment coefficients (including intercept)
		beta_occs_sigma = beta_occs_sigma_inits, # s.d. as function of environment

		N = N_inits_calib, # number of latent AG in calibration counties
		N_ag_county_sq = N_inits_all_counties, # number of latent AG in all counties
		N_ag_county_ssp245_2041_2070 = N_inits_all_counties,
		N_ag_county_ssp245_2071_2100 = N_inits_all_counties,
		N_ag_county_ssp370_2041_2070 = N_inits_all_counties,
		N_ag_county_ssp370_2071_2100 = N_inits_all_counties,

		log_county_lambda_sq = rep(1, data_occs$n_counties), # number of latent AG in all counties
		log_county_lambda_ssp245_2041_2070 = rep(1, data_occs$n_counties),
		log_county_lambda_ssp245_2071_2100 = rep(1, data_occs$n_counties),
		log_county_lambda_ssp370_2041_2070 = rep(1, data_occs$n_counties),
		log_county_lambda_ssp370_2071_2100 = rep(1, data_occs$n_counties),

		N_ag_county_thirties = rep(10, data_occs$n_counties_20th_cent),
		N_ag_county_fifties = rep(10, data_occs$n_counties_20th_cent),

		log_county_lambda_thirties = rep(log(10), data_occs$n_counties_20th_cent),
		log_county_lambda_fifties = rep(log(10), data_occs$n_counties_20th_cent),

		response_curves_occs_mu = response_curves_occs_mu_inits,
		log_lambda_resp_curves_mu = response_curves_occs_mu_inits,

		log_lik_y = rep(1, data_occs$n_counties_occs_calib)
		
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
		beta_occs_mu[1] ~ dnorm(0, sd = beta_occs_mu_prior_dnorm_sd_1)
		for (i in 2:n_terms_occs) {
			beta_occs_mu[i] ~ ddexp(0, rate = beta_occs_mu_prior_ddexp_rate)
		}

		# OCCURRENCE: priors for standard deviation of normal
		beta_occs_sigma[1] ~ dnorm(0, sd = beta_occs_sigma_prior_dnorm_sd_1)
		for (i in 2:n_terms_occs) {
			beta_occs_sigma[i] ~ ddexp(0, rate = beta_occs_sigma_prior_ddexp_rate)
		}

		# OCCURRENCE: priors for sampling bias
		alpha_occs ~ dnorm(0, sd = alpha_occs_mu_prior_dnorm_sd_1)
		logit(p) <- alpha_occs

		# OCCURRENCE: likelihood
		for (i in 1:n_counties_occs_calib) {
			
			### actual abundance (latent--unobserved)
			N[i] ~ dpois(lambda_mu_sq[i])

			### observed number of AG and sampling bias
			y_n_ag[i] ~ dbinom(prob = p, size = N[i])

			# simulate observations for DHARMa residuals
			y_n_ag_sim[i] ~ dbinom(prob = p, size = N[i])

			# relationship between expected (latent) abundance and environment assuming NORMAL distribution
			log(lambda_mu_sq[i]) ~ dnorm(phi_mu_sq[i], sd = lambda_sigma[i])
			phi_mu_sq[i] <- inprod(beta_occs_mu[1:n_terms_occs], counties_x_occs_calib_sq[i, 1:n_terms_occs])

			lambda_sigma[i] <- exp(inprod(beta_occs_sigma[1:n_terms_occs], counties_x_occs_calib_sq[i, 1:n_terms_occs]))

			# likelihood
			log_lik_y[i] <- dbinom(y_n_ag[i], prob = p, size = N[i], log = 1)

		}

		log_lik <- sum(log_lik_y[1:n_counties_occs_calib])

		# OCCURRENCE: posterior samplers for geographic predictions
		for (i in 1:n_counties) {

			# sq (status quo)
			N_ag_county_sq[i] ~ dpois(county_lambda_sq[i])
			log(county_lambda_sq[i]) ~ dnorm(phi_county_lambda_mu_sq[i], sd = county_lambda_sigma_sq[i])
			phi_county_lambda_mu_sq[i] <-
				inprod(beta_occs_mu[1:n_terms_occs], counties_x_occs_sq[i, 1:n_terms_occs])
			county_lambda_sigma_sq[i] <-
				exp(inprod(beta_occs_sigma[1:n_terms_occs], counties_x_occs_sq[i, 1:n_terms_occs]))

			# ssp245_2041_2070
			N_ag_county_ssp245_2041_2070[i] ~ dpois(county_lambda_ssp245_2041_2070[i])
			log(county_lambda_ssp245_2041_2070[i]) ~ dnorm(phi_county_lambda_mu_ssp245_2041_2070[i], sd = county_lambda_sigma_ssp245_2041_2070[i])
			phi_county_lambda_mu_ssp245_2041_2070[i] <-
				inprod(beta_occs_mu[1:n_terms_occs], counties_x_occs_ssp245_2041_2070[i, 1:n_terms_occs])
			county_lambda_sigma_ssp245_2041_2070[i] <-
				exp(inprod(beta_occs_sigma[1:n_terms_occs], counties_x_occs_ssp245_2041_2070[i, 1:n_terms_occs]))

			# ssp245_2071_2100
			N_ag_county_ssp245_2071_2100[i] ~ dpois(county_lambda_ssp245_2071_2100[i])
			log(county_lambda_ssp245_2071_2100[i]) ~ dnorm(phi_county_lambda_mu_ssp245_2071_2100[i], sd = county_lambda_sigma_ssp245_2071_2100[i])
			phi_county_lambda_mu_ssp245_2071_2100[i] <-
				inprod(beta_occs_mu[1:n_terms_occs], counties_x_occs_ssp245_2071_2100[i, 1:n_terms_occs])
			county_lambda_sigma_ssp245_2071_2100[i] <-
				exp(inprod(beta_occs_sigma[1:n_terms_occs], counties_x_occs_ssp245_2071_2100[i, 1:n_terms_occs]))

			# ssp370_2041_2070
			N_ag_county_ssp370_2041_2070[i] ~ dpois(county_lambda_ssp370_2041_2070[i])
			log(county_lambda_ssp370_2041_2070[i]) ~ dnorm(phi_county_lambda_mu_ssp370_2041_2070[i], sd = county_lambda_sigma_ssp370_2041_2070[i])
			phi_county_lambda_mu_ssp370_2041_2070[i] <-
				inprod(beta_occs_mu[1:n_terms_occs], counties_x_occs_ssp370_2041_2070[i, 1:n_terms_occs])
			county_lambda_sigma_ssp370_2041_2070[i] <-
				exp(inprod(beta_occs_sigma[1:n_terms_occs], counties_x_occs_ssp370_2041_2070[i, 1:n_terms_occs]))			

			# ssp370_2071_2100
			N_ag_county_ssp370_2071_2100[i] ~ dpois(county_lambda_ssp370_2071_2100[i])
			log(county_lambda_ssp370_2071_2100[i]) ~ dnorm(phi_county_lambda_mu_ssp370_2071_2100[i], sd = county_lambda_sigma_ssp370_2071_2100[i])
			phi_county_lambda_mu_ssp370_2071_2100[i] <-
				inprod(beta_occs_mu[1:n_terms_occs], counties_x_occs_ssp370_2071_2100[i, 1:n_terms_occs])
			county_lambda_sigma_ssp370_2071_2100[i] <-
				exp(inprod(beta_occs_sigma[1:n_terms_occs], counties_x_occs_ssp370_2071_2100[i, 1:n_terms_occs]))

		}

		# OCCURRENCE: posterior samplers for geographic predictions to 20th century time periods
		for (i in 1:n_counties_20th_cent) {

			# thirties
			N_ag_county_thirties[i] ~ dpois(county_lambda_thirties[i])
			log(county_lambda_thirties[i]) ~ dnorm(phi_county_lambda_mu_thirties[i], sd = county_lambda_sigma_thirties[i])
			phi_county_lambda_mu_thirties[i] <-
				inprod(beta_occs_mu[1:n_terms_occs], counties_x_occs_thirties[i, 1:n_terms_occs])
			county_lambda_sigma_thirties[i] <-
				exp(inprod(beta_occs_sigma[1:n_terms_occs], counties_x_occs_thirties[i, 1:n_terms_occs]))

			# fifties
			N_ag_county_fifties[i] ~ dpois(county_lambda_fifties[i])
			log(county_lambda_fifties[i]) ~ dnorm(phi_county_lambda_mu_fifties[i], sd = county_lambda_sigma_fifties[i])
			phi_county_lambda_mu_fifties[i] <-
				inprod(beta_occs_mu[1:n_terms_occs], counties_x_occs_fifties[i, 1:n_terms_occs])
			county_lambda_sigma_fifties[i] <-
				exp(inprod(beta_occs_sigma[1:n_terms_occs], counties_x_occs_fifties[i, 1:n_terms_occs]))

		}

		# OCCURRENCE: posterior predictive sampler for ENVIRONMENTAL response curves
		# We're assuming occurrence responds to two or more environmental predictors, so the response curve "x" is an array with one "page" per predictor and output is a matrix with one column per predictor
		for (i in 1:n_covariates_occs) {
			
			for (j in 1:n_response_curve_values) {
				
				response_curves_occs_mu[j, i] ~ dpois(lambda_resp_curves_mu[j, i])
				log(lambda_resp_curves_mu[j, i]) ~ dnorm(phi_lambda_resp_curves_mu[j, i], sd = response_curves_occs_sigma[j, i])
				phi_lambda_resp_curves_mu[j, i] <-
					inprod(beta_occs_mu[1:n_terms_occs], resp_curves_x_occs[j, 1:n_terms_occs, i])
				response_curves_occs_sigma[j, i] <-
					exp(inprod(beta_occs_sigma[1:n_terms_occs], resp_curves_x_occs[j, 1:n_terms_occs, i]))

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
		# buildDerivs = TRUE # need for Hamiltonian Monte Carlo
		buildDerivs = FALSE # need for Hamiltonian Monte Carlo
	)

	say('initializeInfo() and $calculate():', level = 2)
	model$initializeInfo()
	calc <- model$calculate()
	say('model$calculate(): ', calc)
	if (is.na(calc) || is.infinite(calc)) stop('Impossible likelihood.')
	# say('model$simulate()')
	# model$simulate()
	# say('model$calculate(): ', model$calculate())

	say('configureMCMC():', level = 2)

	monitors_coeffs_not_indexed <- c('alpha_occs', 'p')
	monitors_coeffs_single_index <- c('beta_occs_mu', 'beta_occs_sigma')
	monitors_coeffs_double_index <- c()

	monitors_derived_not_indexed <- c('log_lik')
	monitors_derived_single_index <- c()
	monitors_derived_double_index <- c()

	monitors_geog_nam <- c(
		'N_ag_county_sq', 'N_ag_county_ssp245_2041_2070', 'N_ag_county_ssp245_2071_2100', 'N_ag_county_ssp370_2041_2070', 
		'N_ag_county_ssp370_2071_2100'
	)

	monitors_geog_conus <- c(
		'N_ag_county_thirties', 'N_ag_county_fifties'
	)

	monitors_dharma <- c(
		'y_n_ag_sim', 'lambda_mu_sq'
	)

	monitors_resp_curves <- c(
		'response_curves_occs_mu', 'response_curves_occs_sigma'
	)

	monitors <- c(monitors_coeffs_not_indexed, monitors_coeffs_single_index, monitors_coeffs_double_index, monitors_derived_not_indexed, monitors_derived_single_index, monitors_derived_double_index, monitors_geog_nam, monitors_geog_conus, monitors_dharma, monitors_resp_curves)

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

	# # AF slice sampler
	# vars <- c(monitors_coeffs_not_indexed, monitors_coeffs_single_index, monitors_coeffs_double_index)
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

	saveRDS(chains, paste0(out_dir, '/chains.rds'))

	say('session info', level = 2)
	print(sessionInfo())

	say(date(), pre = 1)
	sink()

say('###################################################')
say('### collate all predictions into one SpatVector ###')
say('###################################################')

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

say('#################################################')
say('### post-modeling diagnostics and predictions ###')
say('#################################################')

	descrip <- 'occurrence: Poisson ~ normal homoscedastic'

	workflow_postmodeling_generic(facet = 'occurrence', formulae = formulae, descrip = descrip, out_dir = out_dir)

	workflow_postmodeling_occurrence(homoscedastic = homoscedastic, zero_inflated = zero_inflated, formula_occs = formula_occs, formula_occs_bias = formula_occs_bias, pred_vect_nam = pred_vect_nam, out_dir = out_dir)
	
	if (do_crossvalidation) workflow_postmodeling_occurrence_crossvalidation(homoscedastic = homoscedastic, zero_inflated = zero_inflated, formula_occs = formula_occs, formula_occs_bias = formula_occs_bias, constants = constants, out_dir = out_dir)

say(date())
say('FINIS!', deco = '+', level = 1)

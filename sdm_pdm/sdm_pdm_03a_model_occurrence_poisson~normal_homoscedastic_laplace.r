### MODELING ANDROPOGON GERARDI DISTRIBUTION, PHENOTYPE, PHYSIOLOGY, GENOTYPE, and ASSOCIATED MICROBIAL COMMUNITIES
### Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2023-12
###
### This script constructs a species distribution model Andropogon gerardi where "county" is the observational unit. It assumes (latent) abundance follows a Poisson distribution, and the observed number of AG is a binomial distribution where the probability of observing AG is a function of county area and the number of Poaceae recorded in the county (including AG). The expected abundance drawn from a normal distribution where the mean value is given by a function of environmental predictors (climate, soil, etc.). The model is run using nimble with Laplacian approximation.
###
### source('C:/Kaji/R/andropogon_integratedEcology/sdm_pdm/sdm_pdm_03a_model_occurrence_poisson~normal_homoscedastic_laplace.r')
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

	homoscedastic <- TRUE # SPECIFIC TO THIS SCRIPT--SHOULD NOT BE CHANGED
	zero_inflated <- FALSE # SPECIFIC TO THIS SCRIPT--SHOULD NOT BE CHANGED

	### formula for how aspects of species responds to environment
	formula_occs <- ~ 1 + bio1 + bio12 + bio15 + I(bio1^2) + I(bio12^2) + I(bio15^2) # response of occurrence to climate and soil
	# formula_occs <- ~ 1 + bio1 + bio12 + bio15 + ph + I(bio1^2) + I(bio12^2) + I(bio15^2) + I(ph^2)# response of occurrence to climate and soil

	### output folder and bias formula

	# out_dir <- paste0('./outputs_loretta/integrated_sdm_pdm/models_occurrence/[occs_poisson~normal_homoscedastic_bio1^2_bio12^2_bio15^2[bias~quad_ia]]', ifelse(trial, '_TRIAL', ''), '_laplace/')
	# formula_occs_bias <- ~ 1 + area_km2_log10 + n_poaceae_log10p1 + I(area_km2_log10^2) + I(n_poaceae_log10p1^2) + area_km2_log10:n_poaceae_log10p1 # sampling bias for AG records

	# out_dir <- paste0('./outputs_loretta/integrated_sdm_pdm/models_occurrence/[occs_poisson~normal_homoscedastic_bio1^2_bio12^2_bio15^2[bias~quad]]', ifelse(trial, '_TRIAL', ''), '_laplace/')
	# formula_occs_bias <- ~ 1 + area_km2_log10 + n_poaceae_log10p1 + I(area_km2_log10^2) + I(n_poaceae_log10p1^2) # sampling bias for AG records

	out_dir <- paste0('./outputs_loretta/integrated_sdm_pdm/models_occurrence/[occs_poisson~normal_homoscedastic_bio1^2_bio12^2_bio15^2[bias~ia]]', ifelse(trial, '_TRIAL', ''), '_laplace/')
	formula_occs_bias <- ~ 1 + area_km2_log10 + n_poaceae_log10p1 + area_km2_log10:n_poaceae_log10p1 # sampling bias for AG records

	# out_dir <- paste0('./outputs_loretta/integrated_sdm_pdm/models_occurrence/[occs_poisson~normal_homoscedastic_bio1^2_bio12^2_bio15^2[bias~linear]]', ifelse(trial, '_TRIAL', ''), '_laplace/')
	# formula_occs_bias <- ~ 1 + area_km2_log10 + n_poaceae_log10p1 # sampling bias for AG records

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
	data_traits <- prepare_nonbiomass_traits(trait = 'height', formula_traits = ~ bio1, n_response_curve_values = n_response_curve_values, calib = calib)

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
		w_occs_bias = data_occs$w_occs_bias, # model matrix of sampling bias of AG observed occurrences
		n_terms_occs_bias = data_occs$n_terms_occs_bias, # number of terms in sampling bias model

		n_covariates_occs = data_occs$n_covariates_occs, # number of covariates in formula for occurrence model
		n_covariates_occs_bias = data_occs$n_covariates_occs_bias, # number of covariates in formula for occurrence model
		resp_curves_x_occs = data_occs$resp_curves_x_occs, # response curve array for occurrences vs environment
		resp_curves_w_occs = data_occs$resp_curves_w_occs, # response curve array for occurrences vs environment

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

	N_inits_calib <- data_occs$y_n_ag + 1

	prelim_data <- cbind(data_occs$w_occs_bias, data_occs$counties_x_occs_sq[ , 2:ncol(data_occs$counties_x_occs_sq)])
	prelim_model <- glm.fit(prelim_data, y = data_occs$y_n_ag, family = poisson(log))

	alpha_occs_inits <- prelim_model$coefficients[c('(Intercept)', data_occs$terms_occs_bias)]
	beta_occs_inits <- prelim_model$coefficients[c('(Intercept)', data_occs$terms_occs)]

	inits <- list(

		# lambda_sigma = 1,
		log_lambda_sigma = log(1),

		log_lambda_mu_sq = data_occs$y_n_ag, # expected value of number of AG
		alpha_occs = alpha_occs_inits, # intercept, area, # of Poaceae
		beta_occs_mu = beta_occs_inits, # occurrence ~ environment coefficients (including intercept)

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
		beta_occs_mu[1] ~ dnorm(0, sd = beta_occs_mu_prior_dnorm_sd_1)
		for (i in 2:n_terms_occs) {
			beta_occs_mu[i] ~ ddexp(0, rate = beta_occs_mu_prior_ddexp_rate)
		}

		# OCCURRENCE: priors for sampling bias
		alpha_occs[1] ~ dnorm(0, sd = alpha_occs_mu_prior_dnorm_sd_1)
		for (i in 2:n_terms_occs_bias) {
			alpha_occs[i] ~ dnorm(0, sd = alpha_occs_mu_prior_dnorm_sd)
		}

		# OCCURRENCE: prior for spread of normal distribution of lambda
		log(lambda_sigma) ~ dnorm(0, sd = lambda_sigma_prior_sd) # half-Cauchy

		# OCCURRENCE: likelihood
		for (i in 1:n_counties_occs_calib) {
			
			### actual abundance (latent--unobserved)
			N[i] ~ dpois(lambda_mu_sq[i])

			### observed number of AG and sampling bias
			logit(p[i]) <- inprod(alpha_occs[1:n_terms_occs_bias], w_occs_bias[i, 1:n_terms_occs_bias])
			y_n_ag[i] ~ dbinom(prob = p[i], size = N[i])

			# relationship between expected (latent) abundance and environment assuming NORMAL distribution
			log(lambda_mu_sq[i]) ~ dnorm(phi_mu_sq[i], sd = lambda_sigma)
			phi_mu_sq[i] <- inprod(beta_occs_mu[1:n_terms_occs], counties_x_occs_calib_sq[i, 1:n_terms_occs])

		}

	})

	print(model_code)

	say('nimbleModel():', level = 2)
	model <- nimbleModel(
		code = model_code, # our model
		constants = constants, # constants
		data = data, # data
		inits = inits,
		check = TRUE, # any errors?
		calculate = FALSE,
		buildDerivs = TRUE # need for Hamiltonian Monte Carlo
	)

	say('initializeInfo() and $calculate():', level = 2)
	model$initializeInfo()
	calc <- model$calculate()
	say('model$calculate(): ', calc)
	if (is.na(calc) | is.infinite(calc)) stop('Invalid likelihood.')

	say('Compiling ', date())
	compiled <- compileNimble(model, showCompilerOutput = FALSE)

	say('Laplace ', date())
	# also try 'CG', 'L-BFGS-B', 'SANN'.
	build <- buildLaplace(model, control = list(trace = TRUE, maxit = 1000))
	# laplace <- buildLaplace(model, control = list(optimMethod = 'BFGS', trace = TRUE, maxit = 1000))

	say('Compiling Laplace ', date())
	compiled_laplace <- compileNimble(build, project = model)

	laplace <- runLaplace(compiled_laplace)

	saveRDS(laplace, paste0(out_dir, '/laplace.rds'))

#####################
### post-modeling ###
#####################

	descrip <- 'occurrence: Poisson ~ normal homoscedastic'

	workflow_postmodeling_generic_laplace(
		facet = 'occurrence',
		laplace = laplace,
		formulae = formulae,
		descrip = descrip,
		homoscedastic = homoscedastic,
		zero_inflated = zero_inflated,
		out_dir = out_dir
	 )

	workflow_postmodeling_occurrence_laplace(
		laplace = laplace,
		homoscedastic = homoscedastic,
		zero_inflated = zero_inflated,
		formula_occs = formula_occs,
		formula_occs_bias = formula_occs_bias,
		data_occs = data_occs,
		out_dir = out_dir
	)

say(date())
say('FINIS!', deco = '+', level = 1)
sink()

### MODELING ANDROPOGON GERARDI DISTRIBUTION, PHENOTYPE, PHYSIOLOGY, GENOTYPE, and ASSOCIATED MICROBIAL COMMUNITIES
### Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2023-12
###
### This model estimates the distribution of site-level mean Andropogon gerardi ramet biomass. The distribution of biomasses among ramets at a site follows a zero-inflated gamma distribution with individual ramet biomasses drawn from this distribution. The site-level mean biomass is a function of soil/climate. The mean is assumed to respond to one or more climatic predictors (possibly with higher-order terms). The probability of a zero value is also a function of environmental covariates and has the same functional form (but with different coefficients) as the mean.
###
### source('C:/Kaji/R/andropogon_integratedEcology/sdm_pdm/sdm_pdm_06b_model_biomass_site_ziln-plant_ziln_homoscedastic_fx_of_environment.r')

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

	trial <- TRUE # TRUE for testing
	# trial <- FALSE # TRUE for testing

	# calib <- TRUE # use just counties with non-NA Poaceae for calibration region
	calib <- FALSE # use all of North America for calibration region

	homoscedastic <- TRUE # SPECIFIC TO THIS SCRIPT--SHOULD NOT BE CHANGED
	zero_inflated <- TRUE # SPECIFIC TO THIS SCRIPT--SHOULD NOT BE CHANGED

	# formula for how aspects of species responds to environment
	# formula_biomass_mu <- ~ 1 + bio12 + I(bio12^2)# response of biomass to environment
	formula_biomass_mu <- ~ 1 + bio12 # response of biomass to environment

	out_dir <- paste0('./outputs_loretta/integrated_sdm_pdm/models_biomass/[biomass_zigamma~normal_homoscedastic_bio12]', ifelse(trial, '_TRIAL', ''), '/')

	if (!trial) {

		### MCMC settings
		niter <- 1600000
		nburnin <- niter / 2
		thin <- 800
		nchains <- 4

	} else {
		
		### MCMC settings FOR TESTING
		niter <- 1100
		nburnin <- 100
		thin <- 1
		nchains <- 2

		# niter <- 200000
		# nburnin <- niter / 2
		# thin <- 100
		# nchains <- 4

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

	say('This model estimates the distribution of site-level mean Andropogon gerardi ramet biomass. The distribution of biomasses among ramets at a site follows a zero-inflated gamma distribution with individual ramet biomasses drawn from this distribution. The site-level mean biomass is a function of soil/climate. The mean is assumed to respond to one or more climatic predictors (possibly with higher-order terms). The probability of a zero value is also a function of environmental covariates and has the same functional form (but with different coefficients) as the mean.', breaks = 60, post = 1)

	say('MCMC settings:', level = 2)
	say('trial ........................ ', trial)
	say('niter ........................ ', niter)
	say('nburnin ...................... ', nburnin)
	say('thin ......................... ', thin)
	say('nchains ...................... ', nchains)
	say('formula_biomass_mu ........... ', paste(as.character(formula_biomass_mu), collapse = ' '))
	say('homoscedastic ................ ', homoscedastic)
	say('zero_inflated ................ ', zero_inflated, post = 2)

	say('out_dir:')
	say(out_dir)

	formulae <- list(
		formula_biomass_mu = formula_biomass_mu
	)

	########################s
	### data preparation ###
	########################

	data_biomass_mu <- prepare_biomass(formula_biomass = formula_biomass_mu, n_response_curve_values = n_response_curve_values, calib = calib)

	data_occs <- prepare_occurrences(formula_occs = ~ 1, formula_occs_bias = ~ 1, n_response_curve_values = n_response_curve_values, psa_quant = psa_quant, calib = calib)
	
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

		n_covariates_biomass = data_biomass_mu$n_covariates_biomass,
		resp_curves_x_biomass = data_biomass_mu$resp_curves_x_biomass, # response curve array for biomass

		counties_x_biomass_sq = data_biomass_mu$counties_x_biomass_sq,
		counties_x_biomass_ssp245_2041_2070 = data_biomass_mu$counties_x_biomass_ssp245_2041_2070,
		counties_x_biomass_ssp245_2071_2100 = data_biomass_mu$counties_x_biomass_ssp245_2071_2100,
		counties_x_biomass_ssp370_2041_2070 = data_biomass_mu$counties_x_biomass_ssp370_2041_2070,
		counties_x_biomass_ssp370_2071_2100 = data_biomass_mu$counties_x_biomass_ssp370_2071_2100,

		counties_x_biomass_thirties = data_biomass_mu$counties_x_biomass_thirties,
		counties_x_biomass_fifties = data_biomass_mu$counties_x_biomass_fifties,

		# counties
		n_counties = data_biomass_mu$n_counties, # number of counties in the dataset
		n_counties_20th_cent = data_biomass_mu$n_counties_20th_cent, # number of counties in the dataset

		# response curves (general)
		n_response_curve_values = n_response_curve_values # number of values in response curve array

	)

	constants <- c(constants, constants_shared_biomass)

	beta_biomass_mu_inits <- rep(0, constants$n_terms_biomass_mu)
	inits <- list(

		beta_biomass_mu = beta_biomass_mu_inits,
		beta_biomass_pzero = beta_biomass_mu_inits,

		sigma_biomass_within_sites_log = 0,
		sigma_biomass_among_sites_log = 0,

		mu_biomass_site = rep(log(mean(data_biomass_mu$y_biomass)), data_biomass_mu$n_pheno_sites),

		y_biomass_sim = data_biomass_mu$y_biomass # simulated values for biomass (for DHARMa residuals)

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
	# biomass: gamma distribution with mean a functions of ONE environmental predictor... if more, then need an n-dimensional response array
	model_code <- nimbleCode({
	
		# BIOMASS: priors for relationship of mean to environment
		for (i in 1:n_terms_biomass_mu) {
			beta_biomass_mu[i] ~ dnorm(0, sd = beta_biomass_mu_prior_dnorm_sd) # broad prior
		}

		# BIOMASS: priors for probability of zero value across a site
		for (i in 1:n_terms_biomass_mu) {
			beta_biomass_pzero[i] ~ dnorm(0, sd = beta_biomass_pzero_prior_dnorm_sd) # broad prior
		}

		# prior for sd of mean of biomass at a site on lognormal (~ half-Cauchy), ==> vague
		sigma_biomass_among_sites_log ~ dnorm(0, sd = sigma_biomass_among_sites_log_prior_sd)
		sigma_biomass_among_sites <- exp(sigma_biomass_among_sites_log)

		# prior for sd of individual plant biomass on lognormal (~ half-Cauchy), ==> vague
		sigma_biomass_within_sites_log ~ dnorm(0, sd = sigma_biomass_within_sites_log_prior_sd)
		sigma_biomass_within_sites <- exp(sigma_biomass_within_sites_log)

		# BIOMASS: parameters of biomass distribution are latent and functions of environment
		# individual plant biomasses are samples from the site-level distribution defined by the site-level distribution (next chunk after this one)
		for (i in 1:n_pheno_sites) {

			# relationship of site-level mean biomass to the environment
			mu_biomass_site[i] ~ dziln(meanlog = site_biomass_mu_mean_log[i], sdlog = sigma_biomass_among_sites, pzero = pzero_biomass_site[i])
			site_biomass_mu_mean_log[i] <- inprod(beta_biomass_mu[1:n_terms_biomass_mu], x_by_site_biomass[i, 1:n_terms_biomass_mu])
			
			# probability of 0
			logit(pzero_biomass_site[i]) <- inprod(beta_biomass_pzero[1:n_terms_biomass_mu], x_by_site_biomass[i, 1:n_terms_biomass_mu])

			# # moment matching to get dgamma() parameters for distribution of biomasses of individual plants
			# shape_biomass[i] <- mu_biomass_site[i]^2 / site_biomass_plant_sigma_log^2
			# rate_biomass[i] <- mu_biomass_site[i] / site_biomass_plant_sigma_log^2

		}

		# BIOMASS: likelihood of individual plants
		for (i in 1:n_biomass) {

			# ### likelihood using zero-inflated gamma
			# y_biomass[i] ~ dzig(shape = shape_biomass[site_index_biomass[i]], rate = rate_biomass[site_index_biomass[i]], pzero = pzero_biomass_site[site_index_biomass[i]])

			# # simulated values for unconditional DHARMa residuals
			# y_biomass_sim[i] ~ dzig(shape = shape_biomass[site_index_biomass[i]], rate = rate_biomass[site_index_biomass[i]], pzero = pzero_biomass_site[site_index_biomass[i]])

			# log_lik_biomass[i] <- dzig(y_biomass[i], shape = shape_biomass[site_index_biomass[i]], rate = rate_biomass[site_index_biomass[i]], pzero = pzero_biomass_site[site_index_biomass[i]])

			### likelihood using zero-inflated lognormal
			y_biomass[i] ~ dziln(meanlog = mu_biomass_site[site_index_biomass[i]], sdlog = sigma_biomass_within_sites, pzero = pzero_biomass_site[site_index_biomass[i]])

			# simulated values for unconditional DHARMa residuals
			y_biomass_sim[i] ~ dziln(meanlog = mu_biomass_site[site_index_biomass[i]], sdlog = sigma_biomass_within_sites, pzero = pzero_biomass_site[site_index_biomass[i]])

			log_lik_biomass[i] <- dziln(y_biomass[i], meanlog = mu_biomass_site[site_index_biomass[i]], sdlog = sigma_biomass_within_sites, pzero = pzero_biomass_site[site_index_biomass[i]])

		}

		log_lik <- sum(log_lik_biomass[1:n_biomass])

		# BIOMASS: posterior samplers for predictions of to counties in status quo and future
		for (i in 1:n_counties) {

			# biomass: status quo
			log_biomass_county_mu_sq[i] <-
				inprod(beta_biomass_mu[1:n_terms_biomass_mu], counties_x_biomass_sq[i, 1:n_terms_biomass_mu])
			mu_biomass_county_sq[i] <- exp(log_biomass_county_mu_sq[i])

			logit(pzero_biomass_county_sq[i]) <- inprod(beta_biomass_pzero[1:n_terms_biomass_mu], counties_x_biomass_sq[i, 1:n_terms_biomass_mu])


			# biomass: ssp245_2041_2070
			log_biomass_county_mu_ssp245_2041_2070[i] <-
				inprod(beta_biomass_mu[1:n_terms_biomass_mu], counties_x_biomass_ssp245_2041_2070[i, 1:n_terms_biomass_mu])
			mu_biomass_county_ssp245_2041_2070[i] <- exp(log_biomass_county_mu_ssp245_2041_2070[i])
			
			logit(pzero_biomass_county_ssp245_2041_2070[i]) <- inprod(beta_biomass_pzero[1:n_terms_biomass_mu], counties_x_biomass_ssp245_2041_2070[i, 1:n_terms_biomass_mu])

			# biomass: ssp245_2071_2100
			log_biomass_county_mu_ssp245_2071_2100[i] <-
				inprod(beta_biomass_mu[1:n_terms_biomass_mu], counties_x_biomass_ssp245_2071_2100[i, 1:n_terms_biomass_mu])
			mu_biomass_county_ssp245_2071_2100[i] <- exp(log_biomass_county_mu_ssp245_2071_2100[i])

			logit(pzero_biomass_county_ssp245_2071_2100[i]) <- inprod(beta_biomass_pzero[1:n_terms_biomass_mu], counties_x_biomass_ssp245_2071_2100[i, 1:n_terms_biomass_mu])

			# biomass: ssp370_2041_2070
			log_biomass_county_mu_ssp370_2041_2070[i] <-
				inprod(beta_biomass_mu[1:n_terms_biomass_mu], counties_x_biomass_ssp370_2041_2070[i, 1:n_terms_biomass_mu])
			mu_biomass_county_ssp370_2041_2070[i] <- exp(log_biomass_county_mu_ssp370_2041_2070[i])

			logit(pzero_biomass_county_ssp370_2041_2070[i]) <- inprod(beta_biomass_pzero[1:n_terms_biomass_mu], counties_x_biomass_ssp370_2041_2070[i, 1:n_terms_biomass_mu])

			# biomass: ssp370_2071_2100
			log_biomass_county_mu_ssp370_2071_2100[i] <-
				inprod(beta_biomass_mu[1:n_terms_biomass_mu], counties_x_biomass_ssp370_2071_2100[i, 1:n_terms_biomass_mu])
			mu_biomass_county_ssp370_2071_2100[i] <- exp(log_biomass_county_mu_ssp370_2071_2100[i])

			logit(pzero_biomass_county_ssp370_2071_2100[i]) <- inprod(beta_biomass_pzero[1:n_terms_biomass_mu], counties_x_biomass_ssp370_2071_2100[i, 1:n_terms_biomass_mu])

		}

		# BIOMASS: posterior samplers for predictions of to counties in 20th century
		for (i in 1:n_counties_20th_cent) {

			# thirties
			log_biomass_county_mu_thirties[i] <-
				inprod(beta_biomass_mu[1:n_terms_biomass_mu], counties_x_biomass_thirties[i, 1:n_terms_biomass_mu])
			biomass_mu_county_thirties[i] <- exp(log_biomass_county_mu_thirties[i])

			logit(pzero_biomass_county_thirties[i]) <- inprod(beta_biomass_pzero[1:n_terms_biomass_mu], counties_x_biomass_thirties[i, 1:n_terms_biomass_mu])

			# fifties
			log_biomass_county_mu_fifties[i] <-
				inprod(beta_biomass_mu[1:n_terms_biomass_mu], counties_x_biomass_fifties[i, 1:n_terms_biomass_mu])
			biomass_mu_county_fifties[i] <- exp(log_biomass_county_mu_fifties[i])

			logit(pzero_biomass_county_fifties[i]) <- inprod(beta_biomass_pzero[1:n_terms_biomass_mu], counties_x_biomass_fifties[i, 1:n_terms_biomass_mu])

		}

	})

	if (data_biomass_mu$n_covariates_biomass == 1) {

		### univariate model
		response_curve_code <- nimbleCode({

			# # BIOMASS: posterior predictive sampler for response curves: site-level mean
			# # NB we assume just one predictor for biomass, so the response curve "x" is a matrix, not an array
			# for (i in 1:n_response_curve_values) {
					
			# 	log_response_curves_biomass_mu[i] <-
			# 		inprod(beta_biomass_mu[1:n_terms_biomass_mu], resp_curves_x_biomass[i, 1:n_terms_biomass_mu])
			# 	response_curves_biomass_mu[i] <- exp(log_response_curves_biomass_mu[i])
				
			# 	logit(response_curves_biomass_pzero[i]) <- inprod(beta_biomass_pzero[1:n_terms_biomass_mu], resp_curves_x_biomass[i, 1:n_terms_biomass_mu])

			# }

		})
	
	} else {
	
		### multivariate model
		response_curve_code <- nimbleCode({

			# BIOMASS: posterior predictive sampler for response curves: site-level mean
			# NB we assume just one predictor for biomass, so the response curve "x" is a matrix, not an array
			for (i in 1:n_covariates_biomass) {
				
				for (j in 1:n_response_curve_values) {
						
					log_response_curves_biomass_mu[j, i] <-
						inprod(beta_biomass_mu[1:n_terms_biomass_mu], resp_curves_x_biomass[j, 1:n_terms_biomass_mu, i])
					response_curves_biomass_mu[j, i] <- exp(log_response_curves_biomass_mu[j, i])
					
					logit(response_curves_biomass_pzero[j, i]) <- inprod(beta_biomass_pzero[1:n_terms_biomass_mu], resp_curves_x_biomass[j, 1:n_terms_biomass_mu, i])

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
		'sigma_biomass_among_sites', 'sigma_biomass_within_sites'
	)

	# coefficients that have bracketed indexing
	monitors_coeffs_single_index <- c(
		'beta_biomass_mu', 'beta_biomass_pzero'
	)

	monitors_coeffs_double_index <- c(
	)

	monitors_derived_not_indexed <- c(
		'log_lik'
	)

	monitors_derived_single_index <- c(
		'mu_biomass_site', 'pzero_biomass_site'
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
		'biomass_mu_county_thirties', 'biomass_mu_county_fifties',
		'pzero_biomass_county_thirties', 'pzero_biomass_county_fifties'
	)

	monitors_dharma <- c(
		'y_biomass_sim'
	)

	monitors_resp_curves <- c(
		'response_curves_biomass_mu', 'response_curves_biomass_pzero'
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
	# say('NUTS sampler added to all continuous parameters.')

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
	descrip <- paste0(trait, ': gamma~normal homoscedastic')

	workflow_postmodeling_generic(facet = trait, formulae = formulae, descrip = descrip, out_dir = out_dir)

	workflow_postmodeling_biomass(descrip = descrip, formula_biomass_mu = formula_biomass_mu, homoscedastic = homoscedastic, zero_inflated = zero_inflated, pred_vect_nam = pred_vect_nam, pred_vect_conus = pred_vect_conus, out_dir = out_dir)


say(date())
say('FINIS!', deco = '+', level = 1)

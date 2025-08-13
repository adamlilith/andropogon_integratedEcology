### MODELING ANDROPOGON GERARDI DISTRIBUTION, PHENOTYPE, PHYSIOLOGY, GENOTYPE, and ASSOCIATED MICROBIAL COMMUNITIES
### Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2023-12
###
### This model estimates the distribution of site-level mean Andropogon gerardi ramet biomass. The distribution of biomasses among ramets at a site follows a lognormal distribution with individual ramet biomasses drawn from this distribution. The site-level mean biomass is a draw from a lognormal distribution, the mean of which is a a draw from another normal. The mean of this second normal is a function of soil/climate, and assumed to respond to one climatic predictor (possibly with higher-order terms).
###
### source('C:/Kaji/R/andropogon_integratedEcology/sdm_pdm/sdm_pdm_06e_model_biomass_multivariate_lognormal_dnorm_fx_of_environment.r')

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

	out_dir <- paste0('./outputs_loretta/integrated_sdm_pdm/models_biomass_fx_of_environment/biomass_lognormal_dnorm_fx_of_environment_bio12_bio18', ifelse(trial, '_TRIAL', ''), '/')
	if (!trial & file.exists(out_dir)) stop('Output folder already exists.')

	# formula for how aspects of species responds to environment
	# formula_biomass_mu <- ~ 1 + bio12 + I(bio12^2) # response of biomass to environment
	formula_biomass_mu <- ~ 1 + bio12 + bio18 # response of biomass to environment
	formula_occs <- ~ 1 # not used except for prepare_occurrences()

	# ### MCMC settings
	# niter <- 240000 * 4
	# nburnin <- 40000 * 4
	# thin <- 200 * 4
	# nchains <- 4
	# waic <- TRUE

	### MCMC settings
	niter <- 240000 * 2
	nburnin <- 40000 * 2
	thin <- 200 * 2
	nchains <- 4
	waic <- TRUE

	# # ### MCMC settings FOR TESTING
	# niter <- 1200
	# nburnin <- 20
	# thin <- 1
	# nchains <- 2
	# waic <- TRUE

#############
### model ###
#############

	dirCreate(out_dir)
	sink(paste0(out_dir, '/runtime_log.txt'), split = TRUE)
	say('MODELING BIOMASS')
	say('Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | ', date(), post = 1)

	### data collation
	##################
	say('data collation', post = 1)

	say('This model estimates the distribution of site-level mean Andropogon gerardi ramet biomass. The distribution of biomasses among ramets at a site follows a lognormal distribution with individual ramet biomasses drawn from this distribution. The site-level mean biomass is a draw from a lognormal distribution, the mean of which is a a draw from another normal. The mean of this second normal is a function of soil/climate, and assumed to respond to one climatic predictor (possibly with higher-order terms).', breaks = 60, post = 1)

	say('MCMC settings:', level = 2)
	say('niter ........................ ', niter)
	say('nburnin ...................... ', nburnin)
	say('thin ......................... ', thin)
	say('nchains ...................... ', nchains)
	say('formula_biomass_mu ........... ', paste(as.character(formula_biomass_mu), collapse = ' '))
	say('trial ........................ ', trial, post = 2)

	say('out_dir')
	say(out_dir, post = 2)

	formulae <- list(
		formula_biomass_mu = formula_biomass_mu
	)

	########################s
	### data preparation ###
	########################

	data_biomass_mu <- prepare_biomass(formula = formula_biomass_mu, n_response_curve_values = n_response_curve_values, calib = calib)
	data_occs <- prepare_occurrences(formula = formula_occs, n_response_curve_values = n_response_curve_values, psa_quant = psa_quant, calib = calib)
	
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

		n_covariates_biomass = data_biomass_mu$n_covariates_biomass, # do not need if just one predictor for biomass
		resp_curves_x_biomass = data_biomass_mu$resp_curves_x_biomass, # response curve array for biomass

		counties_x_biomass_sq = data_biomass_mu$counties_x_biomass_sq,
		counties_x_biomass_ssp245_2041_2070 = data_biomass_mu$counties_x_biomass_ssp245_2041_2070,
		counties_x_biomass_ssp245_2071_2100 = data_biomass_mu$counties_x_biomass_ssp245_2071_2100,
		counties_x_biomass_ssp370_2041_2070 = data_biomass_mu$counties_x_biomass_ssp370_2041_2070,
		counties_x_biomass_ssp370_2071_2100 = data_biomass_mu$counties_x_biomass_ssp370_2071_2100,

		# counties
		n_counties = data_occs$n_counties, # number of counties in the dataset

		# response curves (general)
		n_response_curve_values = n_response_curve_values # number of values in response curve array

	)

	beta_biomass_mu_inits <- rep(0, constants$n_terms_biomass_mu)
	observed_site_mean_biomass <- data_biomass_mu$raw_data_biomass[ , .(mean_biomass = mean(Biomass)), by = SITE][['mean_biomass']]

	inits <- list(

		y_biomass_sim = data_biomass_mu$y_biomass, # simulated values for biomass (for DHARMa residuals)
		mu_biomass_site = observed_site_mean_biomass, # mean site-level biomass

		site_biomass_plant_sigma_log = sd(data_biomass_mu$y_biomass),
		site_biomass_mean_sigma = sd(data_biomass_mu$y_biomass),

		beta_biomass_mu = beta_biomass_mu_inits

	)

	say('Data:')
	print(str(data))
	assignList(data)

	say('Constants:', pre = 1)
	print(str(constants))
	assignList(constants)

	say('Initializations:', pre = 1)
	print(str(inits))

	### define model
	say('nimbleCode():', level = 2)

	### model
	#########
	# biomass: gamma distribution with mean a functions of ONE environmental predictor... if more, then need an n-dimensional response array
	model_code <- nimbleCode({
	
		# BIOMASS: priors for relationship of mean and variance to environment
		for (i in 1:n_terms_biomass_mu) {
			beta_biomass_mu[i] ~ dnorm(0, sd = 10) # broad prior
		}

		# prior for sd of site-level biomass on lognormal (~ half-Cauchy), ==> vague
		site_biomass_plant_sigma_log ~ dhalfflat()
		site_biomass_mean_sigma ~ dhalfflat()

		# BIOMASS: parameters of biomass distribution are latent and functions of environment
		# individual plant biomasses are samples from the site-level distribution (next chunk after this one)
		for (i in 1:n_pheno_sites) {

			# relationship of biomass to the environment
			phi_site_biomass_mu[i] <- inprod(beta_biomass_mu[1:n_terms_biomass_mu], x_by_site_biomass[i, 1:n_terms_biomass_mu])
			mu_biomass_site[i] ~ dnorm(phi_site_biomass_mu[i], sd = site_biomass_plant_sigma_log)

		}

		# BIOMASS: likelihood of individual plants
		for (i in 1:n_biomass) {

			# likelihood
			y_biomass[i] ~ dlnorm(meanlog = mu_biomass_site[site_index_biomass[i]], taulog = 1 / site_biomass_mean_sigma)

			# simulated values for unconditional DHARMa residuals
			y_biomass_sim[i] ~ dlnorm(meanlog = mu_biomass_site[site_index_biomass[i]], taulog = 1 / site_biomass_mean_sigma)
	
			# likelihood
			log_lik_y[i] <- dlnorm(y_biomass[i], meanlog = mu_biomass_site[site_index_biomass[i]], sdlog = site_biomass_mean_sigma, log = 1)

		}

		log_lik <- sum(log_lik_y[1:n_biomass])

		# BIOMASS: posterior samplers for predictions of to counties in status quo and future
		for (i in 1:n_counties) {

			# biomass: status quo
			mu_biomass_county_sq[i] <-
				exp(inprod(beta_biomass_mu[1:n_terms_biomass_mu], counties_x_biomass_sq[i, 1:n_terms_biomass_mu]))

			# biomass: future
			mu_biomass_county_ssp245_2041_2070[i] <-
				exp(inprod(beta_biomass_mu[1:n_terms_biomass_mu], counties_x_biomass_ssp245_2041_2070[i, 1:n_terms_biomass_mu]))
			
			mu_biomass_county_ssp245_2071_2100[i] <-
				exp(inprod(beta_biomass_mu[1:n_terms_biomass_mu], counties_x_biomass_ssp245_2071_2100[i, 1:n_terms_biomass_mu]))
			
			mu_biomass_county_ssp370_2041_2070[i] <-
				exp(inprod(beta_biomass_mu[1:n_terms_biomass_mu], counties_x_biomass_ssp370_2041_2070[i, 1:n_terms_biomass_mu]))

			mu_biomass_county_ssp370_2071_2100[i] <-
				exp(inprod(beta_biomass_mu[1:n_terms_biomass_mu], counties_x_biomass_ssp370_2071_2100[i, 1:n_terms_biomass_mu]))

		}

		# BIOMASS: posterior predictive sampler for response curves: site-level mean
		# NB we assume >1 predictor for biomass, so the response curve "x" is a 3D array, not a matrix, with one "page" per covariate
		for (j in 1:n_covariates_biomass) {
			for (i in 1:n_response_curve_values) {
					
				response_curves_biomass_mu[i, j] <-
					exp(inprod(beta_biomass_mu[1:n_terms_biomass_mu], resp_curves_x_biomass[i, 1:n_terms_biomass_mu, j]))
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
		calculate = FALSE
		# buildDerivs = TRUE # need for Hamiltonian Monte Carlo
	)

	say('$initializeInfo() and $calculate():', level = 2)
	model$initializeInfo()
	calc <- model$calculate()
	say('model$calculate(): ', calc)
	if (is.na(calc)) stop('Likelihood is incalculable.')

	say('configureMCMC():', level = 2)

	# monitors for coefficients that have no indexing
	monitors_coeffs_not_indexed <- c(
		'site_biomass_plant_sigma_log', 'site_biomass_mean_sigma'
	)

	# coefficients that have bracketed indexing
	monitors_coeffs_single_index <- c(
		'beta_biomass_mu'
	)

	monitors_coeffs_double_index <- c()

	monitors_derived_not_indexed <- c(
		'log_lik'
	)

	monitors_derived_indexed <- c(
		'mu_biomass_site'
	)

	monitors_geog <- c(
		'mu_biomass_county_sq', 'mu_biomass_county_ssp245_2041_2070', 	
		'mu_biomass_county_ssp245_2071_2100', 'mu_biomass_county_ssp370_2041_2070', 'mu_biomass_county_ssp370_2071_2100'
	)

	monitors_dharma <- c(
		'y_biomass_sim'
	)

	monitors_resp_curves <- c(
		'response_curves_biomass_mu'
	)

	monitors <- c(monitors_coeffs_not_indexed, monitors_coeffs_single_index, monitors_coeffs_double_index, monitors_derived_not_indexed, monitors_derived_indexed, monitors_geog, monitors_dharma, monitors_resp_curves)

	conf <- configureMCMC(
		model,
		monitors = monitors,
		print = TRUE,
		enableWAIC = waic
	)

	# # add no U-turn sampler (Hamiltonian Monte Carlo)
	# conf$addSampler(target = monitors_coeffs, type = 'NUTS')
	# say('NUTS sampler added to all continuous parameters.')

	# conf$removeSamplers('beta_occs_vs_biomass')
	# conf$addSampler(target = 'beta_occs_vs_biomass', type = 'AF_slice')
	# say('AF slice sampler added to beta_occs_vs_biomass.')

	conf$removeSamplers('beta_biomass_mu')
	conf$addSampler(target = c('beta_biomass_mu'), type = 'AF_slice')
	say('AF slice sampler added to beta_biomass_mu_rate.')

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
		WAIC = waic,
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

	### collate all predictions into one SpatVector
	###############################################

		pred_vect_nam <- data_occs$ag_vect_sq

		for (var in monitors_geog) {
		
			preds <- hammer_extract(chains, param = var, j = TRUE, stat = 'mean')
			pred_vect_nam[[var]] <- preds
			names(pred_vect_nam)[ncol(pred_vect_nam)] <- var

		}

		writeVector(pred_vect_nam, paste0(out_dir, '/prediction_vector.gpkg'), overwrite = TRUE)

	### post-modeling analysis of BIOMASS
	descrip <- paste0('biomass: lognormal~dnorm homoscedastic')
	workflow_postmodeling_generic(facet = trait, formulae = formulae, descrip = descrip, out_dir = out_dir)
	workflow_postmodeling_biomass(homoscedastic = TRUE, pred_vect_nam = pred_vect_nam)

say(date())
say('FINIS!', deco = '+', level = 1)

### MODELING ANDROPOGON GERARDI DISTRIBUTION, PHENOTYPE, PHYSIOLOGY, GENOTYPE, and ASSOCIATED MICROBIAL COMMUNITIES
### Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2023-12
###
### This script constructs a distribution model of site-level mean biomass for Andropogon gerardi, which is in turn a covariate for a distribution model. It assumes that the mean site-level response of biomass is a function of climate (and perhaps also soil) variables, and that the standard deviation of biomass within a site is constant across sites.
###
### source('C:/Kaji/R/andropogon_integratedEcology/sdm_pdm/sdm_pdm_08_model_occurrence_fx_of_environment_biomass_homoscedastic.r')
### 
### CONTENTS ###
### setup ###
### model ###

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

	out_dir <- paste0('./outputs_loretta/integrated_sdm_pdm/models_occurrence_fx_of_biomass_environment/occurrence_as_fx_of_biomass_environment_nam_calib', ifelse(trial, '_TRIAL', ''), '/')

	# formula for how aspects of species responds to environment
	formula_occs <- ~ 1 + bio1 + bio15 + I(bio1^2) + I(bio15^2) # response of occurrence to climate
	formula_biomass <- ~ 1 + bio12 # response of biomass to environment
	formula_occs_vs_biomass <- ~ -1 + biomass + I(biomass^2) # response of occurrence to biomass

	dirCreate(out_dir)
	formula <- list(
		formula_occs = formula_occs,
		formula_biomass = formula_biomass,
		formula_occs_vs_biomass = formula_occs_vs_biomass
	)
	saveRDS(formula, paste0(out_dir, '/formula.rds'))

	### MCMC settings
	niter <- 240000
	nburnin <- 40000
	thin <- 200
	nchains <- 4
	waic <- TRUE
	k_folds <- 5

	# # ### MCMC settings FOR TESTING
	# niter <- 2200
	# nburnin <- 200
	# thin <- 2
	# nchains <- 2
	# waic <- TRUE
	# k_folds <- 2

	# # ### MCMC settings FOR TESTING
	# niter <- 22000
	# nburnin <- 2000
	# thin <- 20
	# nchains <- 2
	# waic <- TRUE
	# k_folds <- 2

	# ### MCMC settings FOR NUTS
	# niter <- 3000
	# nburnin <- 1000
	# thin <- 2
	# nchains <- 4
	# waic <- TRUE
	# k_folds <- 5

#############
### model ###
#############

	sink(paste0(out_dir, '/runtime_log.txt'), split = TRUE)
	say('MODELING OCCURRENCE AS A FUNCTION OF BIOMASS')
	say('Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | ', date(), post = 1)

	### data collation
	##################
	say('data collation', post = 1)

	say('This script constructs a distribution model Andropogon gerardi abundance, which is a function of climate/soil covariates and of biomass. The distribution of biomasses among ramets at a site follows a gamma distribution with individual ramet biomasses drawn from this distribution. The site-level mean and standard deviation of biomass is a function of soil/climate. Abundance is an N-mixture model.', breaks = 60, post = 1)

	say('MCMC settings:', level = 2)
	say('niter ........................ ', niter)
	say('nburnin ...................... ', nburnin)
	say('thin ......................... ', thin)
	say('nchains ...................... ', nchains)
	say('formula_biomass .............. ', paste(as.character(formula_biomass), collapse = ' '))
	say('formula_occs ................. ', paste(as.character(formula_occs), collapse = ' '))
	say('formula_occs_vs_biomass ...... ', paste(as.character(formula_occs_vs_biomass), collapse = ' '))
	say('trial ........................ ', trial, post = 2)

	########################s
	### data preparation ###
	########################

	data_biomass <- prepare_biomass(formula = formula_biomass, n_response_curve_values = n_response_curve_values)
	
	data_occs <- prepare_occurrences(formula = formula_occs, n_response_curve_values = n_response_curve_values, psa_quant = psa_quant, calib = calib)
	
	data_response_curves_occs_vs_biomass <- create_response_curve_array_occ_vs_biomass(formula = formula_occs_vs_biomass, data_biomass = data_biomass, n_response_curve_values = n_response_curve_values, calib = calib)

	#########################
	### inputs for nimble ###
	#########################

	say('Inputs:', level = 2)
	data <- list(
		y_biomass = data_biomass$y_biomass,		# biomass of individual plants
		y_n_ag = data_occs$y_n_ag				# number of AG observations in each county
	)

	n_terms_occs_vs_biomass <- length(attr(terms(formula_occs_vs_biomass), 'term.labels'))

	constants <- list(
		
		# sampled sites
		n_pheno_sites = data_biomass$n_pheno_sites, # number of phenotype sample sites

		### biomass
		x_by_site_biomass = data_biomass$x_by_site_biomass, # MM with covariates for biomass (scaled)
		n_biomass = data_biomass$n_biomass, # number of biomass observations
		site_index_biomass = data_biomass$site_index_biomass, # index of sampled site for each row in biomass data
		n_terms_biomass = data_biomass$n_terms_biomass + 1, # number of terms in formula for biomass model (including intercept)

		site_biomass_mean = data_biomass$x_centers_biomass, # mean biomass for scaling
		site_biomass_sd = data_biomass$x_scales_biomass, # sd of biomass for scaling

		# n_covariates_biomass = data_biomass$n_covariates_biomass, # do not need if just one predictor for biomass
		resp_curves_x_biomass = data_biomass$resp_curves_x_biomass, # response curve array for biomass

		counties_x_biomass_sq = data_biomass$counties_x_biomass_sq,
		counties_x_biomass_ssp245_2041_2070 = data_biomass$counties_x_biomass_ssp245_2041_2070,
		counties_x_biomass_ssp245_2071_2100 = data_biomass$counties_x_biomass_ssp245_2071_2100,
		counties_x_biomass_ssp370_2041_2070 = data_biomass$counties_x_biomass_ssp370_2041_2070,
		counties_x_biomass_ssp370_2071_2100 = data_biomass$counties_x_biomass_ssp370_2071_2100,

		index_biomass_county_in_calib_county = data_biomass$index_biomass_county_in_calib_county, # index of biomass county in calibration region

		### occurrences
		n_counties_occs_calib = data_occs$n_counties_occs_calib, # number of counties in calibration region
		n_terms_occs = data_occs$n_terms_occs + 1, # number of terms in formula for occurrence model (including intercept)
		area_km2_log10 = data_occs$area_km2_log10, # log10(area_km2) for each county
		n_poaceae_log10p1 = data_occs$n_poaceae_log10p1, # log10(n_poaceae + 1) for each county

		n_terms_occs_vs_biomass = n_terms_occs_vs_biomass, # number of terms in formula for occurrence as fx of biomass		
		
		n_covariates_occs = data_occs$n_covariates_occs, # number of covariates in formula for occurrence model
		resp_curves_x_occs = data_occs$resp_curves_x_occs, # response curve array for occurrences vs environment
		response_curve_x_occs_vs_biomass = data_response_curves_occs_vs_biomass$response_curves_x_occs_vs_biomass, # response curve array for occurrence vs biomass

		counties_x_occs_sq = data_occs$counties_x_occs_sq,
		counties_x_occs_ssp245_2041_2070 = data_occs$counties_x_occs_ssp245_2041_2070,
		counties_x_occs_ssp245_2071_2100 = data_occs$counties_x_occs_ssp245_2071_2100,
		counties_x_occs_ssp370_2041_2070 = data_occs$counties_x_occs_ssp370_2041_2070,
		counties_x_occs_ssp370_2071_2100 = data_occs$counties_x_occs_ssp370_2071_2100,

		# counties
		n_counties = data_occs$n_counties, # number of counties in the dataset

		# response curves (general)
		n_response_curve_values = n_response_curve_values # number of values in response curve array

	)

	# initial values from occurrence-only model
	# SDM-only formula ~ bio1 + bio12 + bio15 + bio1^2 + bio12^2 + bio15^2:
	# beta0: 2.0279181
	# beta1 (bio1): 0.6847217
	# beta2 (bio15): -0.1890767
	# beta3 (bio1^2): -1.4889148
	# beta4 (bio15^2): -1.3496611
	# alpha0: 1.4302661
	# alpha1 (area): 3.5144105
	# alpha2 (Poaceae): 3.6310264
	beta_occs_inits <- c(2, 0.7, -0.2, -1.5, -1.3) # beta0, beta1, beta2, beta3, beta4
	alpha_occs_inits <- c(1.4, 3.5, 3.6) # alpha0, alpha1, alpha2

	# initial values from biomass-only model
	# beta0_rate: 3.8275738
	# beta1_rate (bio12): 0.6434223
	# beta0_shape: 0.9504877
	# beta1_shape (bio12): 0.6610501
	beta_biomass_mu_rate_inits <- c(3.8, 0.65) # beta0_rate, beta1_rate
	beta_biomass_mu_shape_inits <- c(0.95, 0.65) # beta0_shape, beta1_shape

	N_inits_calib <- data_occs$y_n_ag * 2 + 10
	N_inits_all_counties <- data_occs$ag_vect_sq$n_andropogon_gerardi * 2 + 10

	inits <- list(

		y_biomass_sim = data_biomass$y_biomass, # simulated values for biomass (for DHARMa residuals)

		beta_biomass_mu_rate = beta_biomass_mu_rate_inits,
		beta_biomass_mu_shape = beta_biomass_mu_shape_inits,

		y_n_ag_sim = data_occs$y_n_ag, # simulated values of observed number of AG (for DHARMa residuals)
		alpha_occs = alpha_occs_inits, # intercept, area, # of Poaceae
		beta_occs = beta_occs_inits, # occurrence ~ environment coefficients (including intercept)
		beta_occs_vs_biomass = c(0, 1), # linear + quadratic

		N = N_inits_calib, # number of latent AG in calibration counties
		N_ag_county_sq = N_inits_all_counties, # number of latent AG in all counties
		N_ag_county_ssp245_2041_2070 = N_inits_all_counties,
		N_ag_county_ssp245_2071_2100 = N_inits_all_counties,
		N_ag_county_ssp370_2041_2070 = N_inits_all_counties,
		N_ag_county_ssp370_2071_2100 = N_inits_all_counties
		
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
	# biomass: gamma distribution with mean and variance as functions of ONE environmental predictor (bio12)
	# SDM: N-mixture model with Poisson distribution for abundance and binomial distribution for observed number of AG where latent abundance is a function of TWO or more environmental predictors and of biomass
	model_code <- nimbleCode({
	
		### BIOMASS
		###########

			# BIOMASS: priors for relationship of mean and variance to environment
			for (i in 1:n_terms_biomass) {
				beta_biomass_mu_rate[i] ~ dnorm(0, sd = 10) # broad prior
				beta_biomass_mu_shape[i] ~ dnorm(0, sd = 10) # broad prior
			}

			# BIOMASS: parameters of biomass distribution are latent and functions of environment
			# individual plant biomasses are samples from the site-level distribution defined by the site-level distribution (next chunk after this one)
			for (i in 1:n_pheno_sites) {

				# log of biomass and shape parameter of the gamma distribution
				log_site_biomass_mu_hat[i] <- inprod(beta_biomass_mu_rate[1:n_terms_biomass], x_by_site_biomass[i, 1:n_terms_biomass])
				log_shape_biomass[i] <- inprod(beta_biomass_mu_shape[1:n_terms_biomass], x_by_site_biomass[i, 1:n_terms_biomass])
				
				# site-level mean biomass and the shape parameter of the gamma distribution
				mu_biomass_site[i] <- exp(log_site_biomass_mu_hat[i])
				shape_biomass[i] <- exp(log_shape_biomass[i])
				
				rate_biomass[i] <- shape_biomass[i] / mu_biomass_site[i]

				# sd of biomass
				site_biomass_sigma[i] <- mu_biomass_site[i] / sqrt(rate_biomass[i])

			}

			# BIOMASS: likelihood of individual plants
			for (i in 1:n_biomass) {

				# likelihood
				y_biomass[i] ~ dgamma(shape = shape_biomass[site_index_biomass[i]], rate = rate_biomass[site_index_biomass[i]])

				# simulated values for unconditional DHARMa residuals
				y_biomass_sim[i] ~ dgamma(shape = shape_biomass[site_index_biomass[i]], rate = rate_biomass[site_index_biomass[i]])
		
			}

			# BIOMASS: posterior samplers for predictions of to counties in status quo and future
			for (i in 1:n_counties) {

				# biomass: status quo
				log_biomass_county_mu_sq[i] <-
					inprod(beta_biomass_mu_rate[1:n_terms_biomass], counties_x_biomass_sq[i, 1:n_terms_biomass])
				mu_biomass_county_sq[i] <- exp(log_biomass_county_mu_sq[i])

				log_biomass_county_sigma_sq[i] <-
					inprod(beta_biomass_mu_shape[1:n_terms_biomass], counties_x_biomass_sq[i, 1:n_terms_biomass])
				biomass_sigma_county_sq[i] <- exp(log_biomass_county_sigma_sq[i])

				# biomass: future
				log_biomass_county_mu_ssp245_2041_2070[i] <-
					inprod(beta_biomass_mu_rate[1:n_terms_biomass], counties_x_biomass_ssp245_2041_2070[i, 1:n_terms_biomass])
				mu_biomass_county_ssp245_2041_2070[i] <- exp(log_biomass_county_mu_ssp245_2041_2070[i])
				
				log_biomass_county_mu_ssp245_2071_2100[i] <-
					inprod(beta_biomass_mu_rate[1:n_terms_biomass], counties_x_biomass_ssp245_2071_2100[i, 1:n_terms_biomass])
				mu_biomass_county_ssp245_2071_2100[i] <- exp(log_biomass_county_mu_ssp245_2071_2100[i])
				
				log_biomass_county_mu_ssp370_2041_2070[i] <-
					inprod(beta_biomass_mu_rate[1:n_terms_biomass], counties_x_biomass_ssp370_2041_2070[i, 1:n_terms_biomass])
				mu_biomass_county_ssp370_2041_2070[i] <- exp(log_biomass_county_mu_ssp370_2041_2070[i])

				log_biomass_county_mu_ssp370_2071_2100[i] <-
					inprod(beta_biomass_mu_rate[1:n_terms_biomass], counties_x_biomass_ssp370_2071_2100[i, 1:n_terms_biomass])
				mu_biomass_county_ssp370_2071_2100[i] <- exp(log_biomass_county_mu_ssp370_2071_2100[i])

			}

			# BIOMASS: posterior predictive sampler for response curves: site-level mean
			# NB we assume just one predictor for biomass, so the response curve "x" is a matrix, not an array
			for (i in 1:n_response_curve_values) {
					
				log_response_curves_biomass_mu[i] <-
					inprod(beta_biomass_mu_rate[1:n_terms_biomass], resp_curves_x_biomass[i, 1:n_terms_biomass])
				
				response_curves_biomass_mu[i] <- exp(log_response_curves_biomass_mu[i])
				
				log_response_curves_biomass_sigma[i] <-
					inprod(beta_biomass_mu_shape[1:n_terms_biomass], resp_curves_x_biomass[i, 1:n_terms_biomass])
				
				response_curves_biomass_sigma[i] <- exp(log_response_curves_biomass_sigma[i])
				
			}
	
		### OCCURRENCE
		##############

			# OCCURRENCE: priors for relationship to environment
			for (i in 1:n_terms_occs) {
				beta_occs[i] ~ dnorm(0, sd = 10) # broad prior
			}

			# OCCURRENCE: priors for relationship to biomass
			for (i in 1:2) {
				beta_occs_vs_biomass[i] ~ dnorm(0, sd = 10)
			}

			# OCCURRENCE: priors for sampling bias
			for (i in 1:3) {
				alpha_occs[i] ~ dnorm(0, sd = 10) # 
			}

			# OCCURRENCE: likelihood
			for (i in 1:n_counties_occs_calib) {
				
				### actual abundance (latent--unobserved)
				N[i] ~ dpois(lambda_sq[i])

				### observed number of AG
				y_n_ag[i] ~ dbin(prob = p[i], size = N[i])
				y_n_ag_sim[i] ~ dbin(prob = p[i], size = N[i]) # simulate values for DHARMa residuals

				# sampling bias
				logit(p[i]) <- alpha_occs[1] + alpha_occs[2] * area_km2_log10[i] + alpha_occs[3] * n_poaceae_log10p1[i]

				# scaled, site-level mean biomass for this county
				biomass_county_mu_sq_calib_scaled[i] <-
					(mu_biomass_county_sq[index_biomass_county_in_calib_county[i]] - site_biomass_mean) / site_biomass_sd

				# relationship between expected (latent) abundance and environment and biomass
				log(lambda_sq[i]) <- inprod(beta_occs[1:n_terms_occs], counties_x_occs_sq[i, 1:n_terms_occs]) +
					beta_occs_vs_biomass[1] * biomass_county_mu_sq_calib_scaled[i] +
					beta_occs_vs_biomass[2] * biomass_county_mu_sq_calib_scaled[i]^2

			}

			# OCCURRENCE: posterior samplers for geographic predictions
			for (i in 1:n_counties) {

				N_ag_county_sq[i] ~ dpois(county_lambda_sq[i])
				log(county_lambda_sq[i]) <-
					inprod(beta_occs[1:n_terms_occs], counties_x_occs_sq[i, 1:n_terms_occs]) +
					beta_occs_vs_biomass[1] * biomass_county_mu_sq_scaled[i] +
					beta_occs_vs_biomass[2] * biomass_county_mu_sq_scaled[i]^2
				biomass_county_mu_sq_scaled[i] <-
					(mu_biomass_county_sq[i] - site_biomass_mean) / site_biomass_sd

				N_ag_county_ssp245_2041_2070[i] ~ dpois(county_lambda_ssp245_2041_2070[i])
				log(county_lambda_ssp245_2041_2070[i]) <-
					inprod(beta_occs[1:n_terms_occs], counties_x_occs_ssp245_2041_2070[i, 1:n_terms_occs]) +
					beta_occs_vs_biomass[1] * biomass_county_mu_ssp245_2041_2070_scaled[i] +
					beta_occs_vs_biomass[2] * biomass_county_mu_ssp245_2041_2070_scaled[i]^2
				biomass_county_mu_ssp245_2041_2070_scaled[i] <-
					(mu_biomass_county_ssp245_2041_2070[i] - site_biomass_mean) / site_biomass_sd

				N_ag_county_ssp245_2071_2100[i] ~ dpois(county_lambda_ssp245_2071_2100[i])
				log(county_lambda_ssp245_2071_2100[i]) <-
					inprod(beta_occs[1:n_terms_occs], counties_x_occs_ssp245_2071_2100[i, 1:n_terms_occs]) +
					beta_occs_vs_biomass[1] * biomass_county_mu_ssp245_2071_2100_scaled[i] +
					beta_occs_vs_biomass[2] * biomass_county_mu_ssp245_2071_2100_scaled[i]^2
				biomass_county_mu_ssp245_2071_2100_scaled[i] <-
					(mu_biomass_county_ssp245_2071_2100[i] - site_biomass_mean) / site_biomass_sd

				N_ag_county_ssp370_2041_2070[i] ~ dpois(county_lambda_ssp370_2041_2070[i])
				log(county_lambda_ssp370_2041_2070[i]) <-
					inprod(beta_occs[1:n_terms_occs], counties_x_occs_ssp370_2041_2070[i, 1:n_terms_occs]) +
					beta_occs_vs_biomass[1] * biomass_county_mu_ssp370_2041_2070_scaled[i] +
					beta_occs_vs_biomass[2] * biomass_county_mu_ssp370_2041_2070_scaled[i]^2
				biomass_county_mu_ssp370_2041_2070_scaled[i] <-
					(mu_biomass_county_ssp370_2041_2070[i] - site_biomass_mean) / site_biomass_sd

				N_ag_county_ssp370_2071_2100[i] ~ dpois(county_lambda_ssp370_2071_2100[i])
				log(county_lambda_ssp370_2071_2100[i]) <-
					inprod(beta_occs[1:n_terms_occs], counties_x_occs_ssp370_2071_2100[i, 1:n_terms_occs]) +
					beta_occs_vs_biomass[1] * biomass_county_mu_ssp370_2071_2100_scaled[i] +
					beta_occs_vs_biomass[2] * biomass_county_mu_ssp370_2071_2100_scaled[i]^2
				biomass_county_mu_ssp370_2071_2100_scaled[i] <-
					(mu_biomass_county_ssp370_2071_2100[i] - site_biomass_mean) / site_biomass_sd

			}

			# OCCURRENCE: posterior predictive sampler for ENVIRONMENTAL response curves
			# NB we're assuming occurrence responds to just one environmental predictor, so the response curve "x" is a matrix, not an array
			for (i in 1:n_covariates_occs) {
				
				for (j in 1:n_response_curve_values) {
					
					log(response_curves_occs_mu[j, i]) <-
						inprod(beta_occs[1:n_terms_occs], resp_curves_x_occs[j, 1:n_terms_occs, i])

					# exp_response_curves_occ_lambda_vs_env[j, i] <- exp(response_curves_occs_mu[j, i])

				}

			}

			# OCCURRENCE: posterior predictive sampler for BIOMASS response curves
			for (i in 1:n_response_curve_values) {
				
				log(response_curve_occs_vs_biomass[i]) <-
					inprod(beta_occs_vs_biomass[1:2], response_curve_x_occs_vs_biomass[i, 1:n_terms_occs_vs_biomass])

				# exp_response_curve_occs_vs_biomass[i] <- exp(response_curve_occs_vs_biomass[i])

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

	say('configureMCMC():', level = 2)

	monitors_coeffs <- c(
		'beta_biomass_mu_rate', 'beta_biomass_mu_shape',
		'beta_occs', 'beta_occs_vs_biomass', 'alpha_occs'
	)

	monitors_derived <- c(
		'mu_biomass_site', 'site_biomass_sigma'
	)

	monitors_geog <- c(
		'mu_biomass_county_sq', 'biomass_sigma_county_sq', 'mu_biomass_county_ssp245_2041_2070', 	
		'mu_biomass_county_ssp245_2071_2100', 'mu_biomass_county_ssp370_2041_2070', 'mu_biomass_county_ssp370_2071_2100',
		
		'N_ag_county_sq', 'N_ag_county_ssp245_2041_2070', 'N_ag_county_ssp245_2071_2100', 'N_ag_county_ssp370_2041_2070', 'N_ag_county_ssp370_2071_2100'
		# 'county_lambda_sq', 'county_lambda_ssp245_2041_2070', 'county_lambda_ssp245_2071_2100', 'county_lambda_ssp370_2041_2070', 'county_lambda_ssp370_2071_2100'
	)

	monitors_dharma <- c(
		'y_biomass_sim', 
		'y_n_ag_sim', 'lambda_sq'
	)

	monitors_resp_curves <- c(
		'response_curves_biomass_mu', 'response_curves_biomass_sigma',
		'response_curves_occs_mu',
		'response_curve_occs_vs_biomass'
	)

	monitors <- c(monitors_coeffs, monitors_derived, monitors_geog, monitors_dharma, monitors_resp_curves)

	conf <- configureMCMC(
		model,
		monitors = monitors,
		print = TRUE,
		enableWAIC = waic
	)

	# # add no U-turn sampler (Hamiltonian Monte Carlo)
	# conf$addSampler(target = monitors_coeffs, type = 'NUTS')
	# say('NUTS sampler added to all continuous parameters.')

	# # slice samplers for correlated parameters
	# conf$removeSamplers('beta_occs[1]')
	# conf$removeSamplers('beta_occs_vs_biomass')
	# conf$addSampler(target = 'beta_occs[1]', type = 'slice')
	# conf$addSampler(target = 'beta_occs_vs_biomass[1]', type = 'slice')
	# conf$addSampler(target = 'beta_occs_vs_biomass[2]', type = 'slice')
	# say('Slice sampler added to beta_occs[1] and beta_occs_vs_biomass[1:2].')

	# # RW block samplers for correlated parameters
	# conf$removeSamplers('beta_occs[1]')
	# conf$removeSamplers('beta_occs_vs_biomass')
	# conf$addSampler(target = c('beta_occs[1]', 'beta_occs_vs_biomass[1]', 'beta_occs_vs_biomass[2]'), type = 'RW_block')
	# say('RW_block sampler added to beta_occs[1] and beta_occs_vs_biomass[1:2].')

	AF_slice samplers for correlated parameters
	conf$removeSamplers('beta_occs[1]', 'beta_occs_vs_biomass[1]', 'beta_occs_vs_biomass[2]')
	conf$addSampler(target = c('beta_occs[1]', 'beta_occs_vs_biomass[1]', 'beta_occs_vs_biomass[2]'), type = 'AF_slice')
	say('AF_slice sampler added to beta_occs[1] and beta_occs_vs_biomass[1:2].')

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

say('#########################')
say('### model convergence ###')
say('#########################')

	# chains <- readRDS(paste0(out_dir, '/chains.rds'))
	mcmc <- chains$samples
	cols <- c(
		paste0('beta_biomass_mu_rate[', 1:n_terms_biomass, ']'),
		paste0('beta_biomass_mu_shape[', 1:n_terms_biomass, ']'),
		paste0('beta_occs[', 1:n_terms_occs, ']'),
		paste0('beta_occs_vs_biomass[', 1:2, ']'),
		paste0('alpha_occs[', 1:3, ']')
	)
	for (i in 1:nchains) {
		mcmc[[i]] <- mcmc[[i]][ , cols]
	}

	ggs_mcmc <- ggs(mcmc)

	### trace and density plots
	###########################

	# BIOMASS: graphing trace and density plots for all betas
	pars <- paste0('beta_biomass_mu_rate')
	file <- paste0(out_dir, '/beta_rate_biomass_mu_trace.png')
	ggsave(ggs_traceplot(ggs_mcmc, family = pars), file = file, width = 10, height = 8, dpi = 450, bg = 'white')

	file <- paste0(out_dir, '/beta_rate_biomass_mu_density.png')
	ggsave(ggs_density(ggs_mcmc, family = pars, rug = TRUE, hpd = TRUE), file = file, width = 6, height = 8, dpi = 450, bg = 'white')

	# BIOMASS: graphing trace and density plots for all betas
	pars <- paste0('beta_biomass_mu_shape')
	file <- paste0(out_dir, '/beta_shape_biomass_mu_trace.png')
	ggsave(ggs_traceplot(ggs_mcmc, family = pars), file = file, width = 10, height = 8, dpi = 450, bg = 'white')

	file <- paste0(out_dir, '/beta_shape_biomass_mu_density.png')
	ggsave(ggs_density(ggs_mcmc, family = pars, rug = TRUE, hpd = TRUE), file = file, width = 6, height = 8, dpi = 450, bg = 'white')

	# BIOMASS: correlations between samples of biomass shape and rate
	pars <- c('beta_biomass_mu_shape')
	file <- paste0(out_dir, '/beta_shape_biomass_correlations.png')
	ggsave(ggs_crosscorrelation(ggs_mcmc, family = pars), file = file, width = 8, height = 8, dpi = 450, bg = 'white')

	pars <- c('beta_biomass_mu_rate')
	file <- paste0(out_dir, '/beta_rate_biomass_correlations.png')
	ggsave(ggs_crosscorrelation(ggs_mcmc, family = pars), file = file, width = 8, height = 8, dpi = 450, bg = 'white')

	# OCCURRENCE: graphing trace and density plots for all betas
	pars <- paste0('beta_occs')
	file <- paste0(out_dir, '/beta_occs_trace.png')
	ggsave(ggs_traceplot(ggs_mcmc, family = pars), file = file, width = 10, height = 8, dpi = 450, bg = 'white')

	file <- paste0(out_dir, '/beta_occs_density.png')
	ggsave(ggs_density(ggs_mcmc, family = pars, rug = TRUE, hpd = TRUE), file = file, width = 6, height = 8, dpi = 450, bg = 'white')

	# OCCURRENCE: graphing trace and density plots for all alphas
	pars <- paste0('alpha_occs')
	file <- paste0(out_dir, '/alpha_occs_trace.png')
	ggsave(ggs_traceplot(ggs_mcmc, family = pars), file = file, width = 10, height = 8, dpi = 450, bg = 'white')

	file <- paste0(out_dir, '/alpha_occs_density.png')
	ggsave(ggs_density(ggs_mcmc, family = pars, rug = TRUE, hpd = TRUE), file = file, width = 6, height = 8, dpi = 450, bg = 'white')

	### Gelman-Rubin convergence statistic
	######################################
	mcmc_coeffs <- chains$samples
	cols <- c(
		paste0('beta_biomass_mu_rate[', 1:n_terms_biomass, ']'),
		paste0('beta_biomass_mu_shape[', 1:n_terms_biomass, ']'),
		paste0('beta_occs[', 1:n_terms_occs, ']'),
		paste0('beta_occs_vs_biomass[', 1:2, ']'),
		paste0('alpha_occs[', 1:3, ']')
	)
	for (i in 1:nchains) {
		mcmc_coeffs[[i]] <- mcmc_coeffs[[i]][ , cols]
	}

	rhats <- tryCatch(
		gelman.diag(mcmc_coeffs, autoburnin = FALSE, multivariate = TRUE),
		error = function(cond) FALSE
	)
	
	sink(paste0(out_dir, '/convergence.txt'), split = TRUE)
	say('GELMAN-RUBIN STATISTICS')
	say(date(), post = 2)
	print(rhats)
	sink()

say('#############################')
say('### effective sample size ###')
say('#############################')

	mcmc <- chains$samples
	cols <- c(
		paste0('beta_biomass_mu_rate[', 1:n_terms_biomass, ']'),
		paste0('beta_biomass_mu_shape[', 1:n_terms_biomass, ']'),
		paste0('mu_biomass_site[', 1:data_biomass$n_pheno_sites, ']'),
		paste0('site_biomass_sigma[', 1:data_biomass$n_pheno_sites, ']'),
		paste0('beta_occs[', 1:n_terms_occs, ']'),
		paste0('beta_occs_vs_biomass[', 1:2, ']'),
		paste0('alpha_occs[', 1:3, ']')
	)
	for (i in 1:nchains) {
		mcmc[[i]] <- mcmc[[i]][ , cols]
	}

	sink(paste0(out_dir, '/effective_sample_size.txt'), split = TRUE)
		say('EFFECTIVE SAMPLE SIZE', post = 2)
		print(effectiveSize(mcmc))
		say('')
	sink()

say('########################################')
say('### model fit and residuals analysis ###')
say('########################################')

	### WAIC
	########
	sink(paste0(out_dir, '/waic.txt'), split = TRUE)
	say('WAIC')
	say(date(), post = 2)
	print(chains$WAIC)
	sink()

	### BIOMASS: plant-level residuals analysis
	###########################################

	# NB this uses the mean predicted value of a site as a plant-level prediction
	sims_by_plant <- hammer_subset(chains, param = 'y_biomass_sim', j = TRUE)
	sims_by_plant <- hammer_rbind(sims_by_plant)
	sims_by_plant <- t(sims_by_plant)

	estimates <- hammer_subset(chains, param = 'mu_biomass_site', j = TRUE)
	estimates <- hammer_rbind(estimates)

	site_counts <- data_biomass$raw_data_biomass[ , .N, by = SITE]
	fit_by_site <- apply(estimates, 2, median)
	fits <- numeric()
	for (i in seq_along(fit_by_site)) {
		fits <- c(fits, rep(fit_by_site[i], site_counts$N[i]))
	}

	dharma <- createDHARMa(simulatedResponse = sims_by_plant, observedResponse = data_biomass$y_biomass, fittedPredictedResponse = fits, integerResponse = FALSE)

	file <- paste0(out_dir, '/dharma_by_plant_biomass.png')
	png(file, width = 1200, height = 800)
		plot(dharma)
	dev.off()
	
	file <- paste0(out_dir, '/dharma_by_plant_biomass_residuals.png')
	png(file, width = 1400, height = 1000)
		hist(dharma$scaledResiduals, main = 'DHARMa residuals for biomass by plant', xlab = 'Scaled residuals', breaks = 30)
	dev.off()

	### BIOMASS: plant-level residuals vs covariates
	################################################

	resids <- dharma$scaledResiduals
	resids_vs_covariates <- list()
	for (i in seq_along(data_biomass$terms_biomass)) {

		this_x <- data.frame(
			x = data_biomass$raw_data_biomass[[data_biomass$terms_biomass[i]]],
			y = resids
		)

		resids_vs_covariates[[i]] <- ggplot(this_x, aes(x = x, y = y)) +
			geom_point() +
			geom_smooth(method = loess, formula = y ~ x, se = FALSE) +
			xlab(data_biomass$terms_biomass[i]) +
			ylab('Scaled residuals') +
			ggtitle('Biomass-only model')

	}

	if (data_biomass$n_covariates == 1) {
		nrow <- 1
		width <- 10
		height <- 8
	} else if (data_biomass$n_covariates <= 4) {
		nrow <- 1
		width <- 14
		height <- 4
	} else {
		nrow <- 2
		width <- 16
		height <- 10	
	}

	resids_vs_covariates <- plot_grid(plotlist = resids_vs_covariates, nrow = nrow)
	ggsave(plot = resids_vs_covariates, filename = paste0(out_dir, '/dharma_residuals_vs_covariates_biomass.png'), width = width, height = height, dpi = 600)

	### OCCURRENCE: DHARMa residuals
	################################

	sims <- hammer_subset(chains, param = 'y_n_ag_sim', j = TRUE)
	sims <- hammer_rbind(sims)
	sims <- sims[ , data_occs$ag_vect_sq$focal_region]
	sims <- t(sims)

	fits <- hammer_extract(chains, param = 'lambda_sq', j = TRUE, stat = 'mean')
	fits <- fits[data_occs$ag_vect_sq$focal_region]

	observed_y <- data_occs$y_n_ag[data_occs$ag_vect_sq$focal_region]
	
	dharma <- createDHARMa(simulatedResponse = sims, observedResponse = observed_y, fittedPredictedResponse = fits, integerResponse = TRUE)

	file <- paste0(out_dir, '/dharma_n_ag.png')
	png(file, width = 1200, height = 800)
		plot(dharma)
	dev.off()
	
	file <- paste0(out_dir, '/dharma_lambda_residuals_n_ag.png')
	png(file, width = 1200, height = 800)
		hist(dharma$scaledResiduals, main = 'DHARMa residuals for number of observed AG (y_n_ag)', xlab = 'Scaled residuals', breaks = 30)
	dev.off()

say('###########################')
say('### parameter estimates ###')
say('###########################')

	# chains <- readRDS(paste0(out_dir, '/chains.rds'))
	vars <- c(monitors_coeffs, monitors_derived)

	caterpillars <- list()
	for (i in seq_along(vars)) {

		var <- vars[i]

		mcmc <- hammer_subset(chains, param = var, j = TRUE)
		mcmc <- hammer_samples(mcmc)

		ggs_mcmc <- ggs(mcmc)

		caterpillars[[i]] <- ggs_caterpillar(ggs_mcmc, family = var) +
			xlab('Estimated value') +
			ggtitle(var) +
			theme(
				plot.title = element_text(size = 12),
				axis.title.y = element_blank()
			)

	
	}

	combo <- plot_grid(plotlist = caterpillars, ncol = 3, align = 'v', axis = 'l')
	ggsave(combo, filename = paste0(out_dir, '/coefficients.png'), width = 12, height = 8, dpi = 300, bg = 'white')

say('########################################################')
say('### response curves: BIOMASS site-level mu and sigma ###')
say('########################################################')

	# chains <- readRDS(paste0(out_dir, '/chains.rds'))

	n_covariates <- data_biomass$n_covariates_biomass
	covariates <- data_biomass$terms_biomass
	resp_curves_x <- data_biomass$resp_curves_x_biomass
	resp_curves_unscaled <- data_biomass$resp_curve_x_biomass_unscaled
	centers <- data_biomass$x_centers_biomass
	scales <- data_biomass$x_scales_biomass
	y_lab <- bquote('Biomass ' * ' (g)')

	# make a plot for how biomass mu and sigma respond to each predictor
	responses <- list()
	for (response_type in c('mu', 'sigma')) {

		max_val <- -Inf
		responses[[length(responses) + 1]] <- list()
		for (i in 1:n_covariates) {

			pred <- covariates[i]

			if (pred == 'aridity') {
				nice_title <- 'Aridity ((temp. + 10) / (precip. / 1000))'
				nice_axis <- 'Aridity'
			} else if (pred == 'bio1') {
				nice_title <- 'Mean Annual Temperature (BIO01)'
				nice_axis <- 'Mean annual temperature (°C)'
			} else if (pred == 'bio5') {
				nice_title <- 'Temperature of the Hottest Month (BIO05)'
				nice_axis <- 'Temperature of the hottest month (°C)'
			} else if (pred == 'bio6') {
				nice_title <- 'Temperature of Coldest Month (BIO06)'
				nice_axis <- 'Temperature of the coldest month (°C)'
			} else if (pred == 'bio7') {
				nice_title <- 'Temperature Annual Range (BIO07)'
				nice_axis <- 'Temperature annual range (°C)'
			} else if (pred == 'bio12') {
				nice_title <- 'Total Annual Precipitation (BIO12)'
				nice_axis <- 'Total annual precipitation (mm)'
			} else if (pred == 'bio15') {
				nice_title <- 'Precipitation Seasonality (BIO15)'
				nice_axis <- 'Precipitation seasonality'
			} else if (pred == 'bio18') {
				nice_title <- 'Precipitation of Warmest Quarter (BIO18)'
				nice_axis <- 'Precipitation of warmest quarter (mm)'
			} else if (pred == 'ph') {
				nice_title <- 'Soil pH'
				nice_axis <- 'pH'
			} else if (pred == 'sand') {
				nice_title <- 'Soil Proportion Sand'
				nice_axis <- 'Proportion sand'
			} else if (pred == 'silt') {
				nice_title <- 'Soil Proportion Silt'
				nice_axis <- 'Proportion silt'
			} else if (pred == 'soc') {
				nice_title <- 'Soil Organic Matter'
				nice_axis <- 'Soil organic matter'
			}

			if (response_type == 'mu') {
				response_type_nice <- 'mean (μ)'
			} else if (response_type == 'sigma') {
				response_type_nice <- 'S.D. (σ)'
			}

			nice_title <- bquote('Biomass ' * ' versus ' * .(nice_title) * ' ' * .(response_type_nice))

			# unscaled predictor value
			x <- resp_curves_unscaled[ , pred]

			# create data frame with SDM response
			if (n_covariates == 1) {
				pars <- paste0('response_curves_biomass_', response_type, '[', 1:n_response_curve_values, ']')
			} else {
				pars <- paste0('response_curves_biomass_', response_type, '[', 1:n_response_curve_values, ', ', i, ']')
			}
			
			param <- paste0('response_curves_biomass_', response_type)
			response_mean <- hammer_extract(chains, param = param, j = 1:n_response_curve_values, stat = 'mean')
			response_lower <- hammer_extract(chains, param = param, j = 1:n_response_curve_values, stat = 'lower')
			response_upper <- hammer_extract(chains, param = param, j = 1:n_response_curve_values, stat = 'upper')

			# data frames to hold predictions in long format
			df_mean <- data.frame(
				x = x,
				response = response_mean
			)

			df_lower <- data.frame(
				x = x,
				response = response_lower
			)

			df_upper <- data.frame(
				x = x,
				response = response_upper
			)

			this_df_ci <- data.frame(
				x = c(x, rev(x)),
				response = c(df_upper$response, rev(df_lower$response))
			)

			response <- ggplot() +
				geom_polygon(
					data = this_df_ci,
					mapping = aes(x = x, y = response),
					color = NA,
					fill = alpha('blue', 0.1)
				) +
				geom_line(
					data = df_mean,
					mapping = aes(x = x, y = response)
				) +
				xlab(nice_axis) +
				ylab(paste0('Site-level biomass ', response_type_nice)) +
				ggtitle(nice_title) +
				theme(
					legend.position = 'none',
					plot.title = element_text(size = 10)
				)

				data_TEMP <- data_biomass$raw_data_biomass
				data_TEMP$x <- data_TEMP[[pred]]
				data_TEMP$y <- data_TEMP$Biomass
				data_TEMP$site <- data_TEMP$SITE

				max_val <- max(
					max_val,
					quantile(df_upper$response[!is.infinite(df_upper$response)], 0.5)
				)

				response <- response + geom_point(
					data = data_TEMP,
					mapping = aes(x = x, y = y, color = site)
				)

			responses[[length(responses)]][[i]] <- response

		} # next predictor

		for (i in 1:n_covariates) {
			responses[[length(responses)]][[i]] <- responses[[length(responses)]][[i]] + coord_cartesian(ylim = c(0, max_val)) 
		}

	} # next response type

	if (n_covariates == 1) {
		nrow <- 1
		width <- 10
		height <- 8
	} else if (n_covariates < 4) {
		nrow <- 1
		width <- 14
		height <- 4
	} else {
		nrow <- 2
		width <- 16
		height <- 10	
	}
	responses <- unlist(responses, recursive = FALSE)

	responses <- plot_grid(plotlist = responses, nrow = nrow)
	ggsave(plot = responses, filename = paste0(out_dir, '/response_curves_biomass.png'), width = width, height = height, dpi = 600)

say('##################################################')
say('### response curves: OCCURRENCE vs environment ###')
say('##################################################')

	# chains <- readRDS(paste0(out_dir, '/chains.rds'))

	n_covariates <- data_occs$n_covariates_occs
	covariates <- data_occs$terms_occs
	resp_curves_x <- data_occs$resp_curves_x_occs
	resp_curves_unscaled <- data_occs$resp_curves_x_occs_unscaled
	centers <- data_occs$x_centers_occs
	scales <- data_occs$x_scales_occs
	y_lab <- bquote('Expected abundance')

	# make a plot for how biomass mu and sigma respond to each predictor
	responses <- list()
	max_val <- -Inf
	for (i in 1:n_covariates) {

		pred <- covariates[i]

		if (pred == 'aridity') {
			nice_title <- 'Aridity ((temp. + 10) / (precip. / 1000))'
			nice_axis <- 'Aridity'
		} else if (pred == 'bio1') {
			nice_title <- 'Mean Annual Temperature (BIO01)'
			nice_axis <- 'Mean annual temperature (°C)'
		} else if (pred == 'bio5') {
			nice_title <- 'Temperature of the Hottest Month (BIO05)'
			nice_axis <- 'Temperature of the hottest month (°C)'
		} else if (pred == 'bio6') {
			nice_title <- 'Temperature of Coldest Month (BIO06)'
			nice_axis <- 'Temperature of the coldest month (°C)'
		} else if (pred == 'bio7') {
			nice_title <- 'Temperature Annual Range (BIO07)'
			nice_axis <- 'Temperature annual range (°C)'
		} else if (pred == 'bio12') {
			nice_title <- 'Total Annual Precipitation (BIO12)'
			nice_axis <- 'Total annual precipitation (mm)'
		} else if (pred == 'bio15') {
			nice_title <- 'Precipitation Seasonality (BIO15)'
			nice_axis <- 'Precipitation seasonality'
		} else if (pred == 'bio18') {
			nice_title <- 'Precipitation of Warmest Quarter (BIO18)'
			nice_axis <- 'Precipitation of warmest quarter (mm)'
		} else if (pred == 'ph') {
			nice_title <- 'Soil pH'
			nice_axis <- 'pH'
		} else if (pred == 'sand') {
			nice_title <- 'Soil Proportion Sand'
			nice_axis <- 'Proportion sand'
		} else if (pred == 'silt') {
			nice_title <- 'Soil Proportion Silt'
			nice_axis <- 'Proportion silt'
		} else if (pred == 'soc') {
			nice_title <- 'Soil Organic Matter'
			nice_axis <- 'Soil organic matter'
		}

		nice_title <- bquote('Abundance' * ' versus ' * .(nice_title))

		# unscaled predictor value
		x <- resp_curves_unscaled[ , pred]
		
		param <-'response_curves_occs_mu'
		response_mean <- hammer_extract(chains, param = param, j = 1:n_response_curve_values, k = i, stat = 'mean')
		response_lower <- hammer_extract(chains, param = param, j = 1:n_response_curve_values, k = i, stat = 'lower')
		response_upper <- hammer_extract(chains, param = param, j = 1:n_response_curve_values, k = i, stat = 'upper')

		# data frames to hold predictions in long format
		df_mean <- data.frame(
			x = x,
			response = response_mean
		)

		df_lower <- data.frame(
			x = x,
			response = response_lower
		)

		df_upper <- data.frame(
			x = x,
			response = response_upper
		)

		this_df_ci <- data.frame(
			x = c(x, rev(x)),
			response = c(df_upper$response, rev(df_lower$response))
		)

		response <- ggplot() +
			geom_polygon(
				data = this_df_ci,
				mapping = aes(x = x, y = response),
				color = NA,
				fill = alpha('blue', 0.1)
			) +
			geom_line(
				data = df_mean,
				mapping = aes(x = x, y = response)
			) +
			xlab(nice_axis) +
			ylab(paste0('Expected abundance')) +
			ggtitle(nice_title) +
			theme(
				legend.position = 'none',
				plot.title = element_text(size = 9)
			)

			max_val <- max(
				max_val,
				quantile(df_upper$response[!is.infinite(df_upper$response)], 0.5)
			)

		responses[[i]] <- response

	} # next predictor

	for (i in 1:n_covariates) {
		responses[[i]] <- responses[[i]] + coord_cartesian(ylim = c(0, max_val)) 
	}

	if (n_covariates == 1) {
		nrow <- 1
		width <- 10
		height <- 8
	} else if (n_covariates <= 2) {
		nrow <- 1
		width <- 8
		height <- 4
	} else if (n_covariates <= 3) {
		nrow <- 1
		width <- 14
		height <- 4
	} else {
		nrow <- 2
		width <- 16
		height <- 10	
	}

	responses <- plot_grid(plotlist = responses, nrow = nrow)
	ggsave(plot = responses, filename = paste0(out_dir, '/response_curves_abundance.png'), width = width, height = height, dpi = 600)

say('##############################################')
say('### response curves: OCCURRENCE vs BIOMASS ###')
say('##############################################')

	# chains <- readRDS(paste0(out_dir, '/chains.rds'))
	
	y_lab <- bquote('Expected abundance')
	nice_title <- bquote('Abundance versus biomass')

	# unscaled predictor value... note, we need to do this the same
	x <- data_response_curves_occs_vs_biomass$unscaled_biomass

	param <-'response_curve_occs_vs_biomass'
	response_mean <- hammer_extract(chains, param = param, j = 1:n_response_curve_values, stat = 'mean')
	response_lower <- hammer_extract(chains, param = param, j = 1:n_response_curve_values, stat = 'lower')
	response_upper <- hammer_extract(chains, param = param, j = 1:n_response_curve_values, stat = 'upper')

	# data frames to hold predictions in long format
	df_mean <- data.frame(
		x = x,
		response = response_mean
	)

	df_lower <- data.frame(
		x = x,
		response = response_lower
	)

	df_upper <- data.frame(
		x = x,
		response = response_upper
	)

	this_df_ci <- data.frame(
		x = c(x, rev(x)),
		response = c(df_upper$response, rev(df_lower$response))
	)

	max_val <- max(df_mean$response, quantile(df_upper$response[!is.infinite(df_upper$response)], 0.5))

	response <- ggplot() +
		geom_polygon(
			data = this_df_ci,
			mapping = aes(x = x, y = response),
			color = NA,
			fill = alpha('blue', 0.1)
		) +
		geom_line(
			data = df_mean,
			mapping = aes(x = x, y = response)
		) +
		coord_cartesian(ylim = c(0, max_val)) +
		xlab('Mean (μ) site-level biomass (g)') +
		ylab(paste0('Expected abundance (λ)')) +
		ggtitle(nice_title) +
		theme(
			legend.position = 'none',
			plot.title = element_text(size = 10)
		)

	ggsave(plot = response, filename = paste0(out_dir, '/response_curves_abundance_vs_biomass.png'), width = 12, height = 8, dpi = 600)

say('#################################', pre = 1)
say('### map of residuals: BIOMASS ###')
say('#################################')

	# Make map of standardized residuals (predicted - observed) / (observed_mean)

	# # user-defined
	# chains <- readRDS(paste0(out_dir, '/chains.rds'))

	# calculate residuals
	resp <- hammer_extract(chains, param = 'mu_biomass_site', j = TRUE)

	mean_biomass_by_site <- data_biomass$raw_data_biomass[ , .(mean_biomass = mean(Biomass, na.rm = TRUE)), by = SITE]
	stand_resid <- (mean_biomass_by_site$mean_biomass - resp) / mean_biomass_by_site$mean_biomass

	site_vect_biomass_resid <- data_biomass$site_vect_biomass
	site_vect_biomass_resid$stand_resid <- stand_resid

	# administrative boundaries
	nam <- vect('./data_from_gadm/gadm_4pt1_level_1_north_america_sans_alaska_lambert.gpkg')
	nam <- project(nam, site_vect_biomass_resid)

	# extent
	extent <- buffer(site_vect_biomass_resid, width = 200 * 1000)
	extent <- ext(extent)
	extent <- as.vector(extent)

	legend_title <- 'Biomass\nmean-\nstandardized\nresidual'

	map <- ggplot() +
		layer_spatial(nam, color = 'gray30', fill = 'white', linewidth = 0.3) +
		layer_spatial(site_vect_biomass_resid, aes(fill = stand_resid), pch = 21, size = 6) +
		scale_fill_gradient2(
			name = legend_title,
			low = 'red',
			mid = 'beige',
			high = 'blue',
			midpoint = 0
		) +
		xlim(extent[1], extent[2]) + ylim(extent[3], extent[4]) +
		ggtitle(
			bquote('Residuals for ' * italic('Andropogon gerardi') * ' biomass '),
			subtitle = '1991-2020 | biomass-only model') +
		theme(
			plot.title = element_text(size = 16),
			plot.subtitle = element_text(size = 14)
		)

	ggsave(plot = map, filename = paste0(out_dir, '/map_residuals_biomass_by_site.png'), width = 12, height = 9, dpi = 600)

say('###################################################')
say('### collate all predictions into one SpatVector ###')
say('###################################################')

	pred_vect_nam <- data_occs$ag_vect_sq

	for (var in monitors_geog) {
	
		preds <- hammer_extract(chains, param = var, j = TRUE, stat = 'mean')
		pred_vect_nam[[var]] <- preds
		names(pred_vect_nam)[ncol(pred_vect_nam)] <- var

	}

	writeVector(pred_vect_nam, paste0(out_dir, '/prediction_vector.gpkg'), overwrite = TRUE)

say('############################', pre = 1)
say('### current map: BIOMASS ###')
say('############################')

	map_biomass_mu_sq <- map_biomass(
		out_dir = out_dir,
		pred_vect_nam = pred_vect_nam,
		response_var = 'mu_biomass_county_sq',
		response_var_type = 'mu',
		data_occs = data_occs,
		data_biomass = data_biomass,
		title = bquote('Present-day distribution of mean ' * italic('Andropogon gerardi') * ' biomass '),
		subtitle = 'Abundance ~ f(environment, biomass ~ f(environment))  | 1991-2020',
		ag_core_quant = 0.95
	)

	map_biomass_sigma_sq <- map_biomass(
		out_dir = out_dir,
		pred_vect_nam = pred_vect_nam,
		response_var = 'biomass_sigma_county_sq',
		response_var_type = 'sigma',
		data_occs = data_occs,
		data_biomass = data_biomass,
		title = bquote('Present-day distribution of standard deviation of ' * italic('Andropogon gerardi') * ' biomass '),
		subtitle = 'Abundance ~ f(environment, biomass ~ f(environment))  | 1991-2020',
		ag_core_quant = 0.95
	)

say('###########################', pre = 1)
say('### future map: BIOMASS ###')
say('###########################')

	for (fut in futs) {

		say(fut)

		subtitle <- paste0('Abundance ~ f(environment, biomass ~ f(environment))  | SSP ', substr(fut, 4, 6), ' ', substr(fut, 8, 11), '-', substr(fut, 13, 16))

		map_biomass_fut <- map_biomass(
			out_dir = out_dir,
			pred_vect_nam = pred_vect_nam,
			response_var = response_var,
			response_var_type = 'mu',
			data_occs = data_occs,
			data_biomass = data_biomass,
			title = bquote('Future distribution of mean ' * italic('Andropogon gerardi') * ' biomass '),
			subtitle = subtitle,
			ag_core_quant = 0.95
		)

	}

say('##################################', pre = 1)
say('### future change map: BIOMASS ###')
say('##################################')

	maps_biomass_change <- list()
	for (fut in futs) {

		say(fut)

		title <- bquote('Change in ' * italic('Andropogon gerardi') * ' biomass ')
		subtitle <- paste0('Abundance ~ f(environment, biomass ~ f(environment))  | SSP ', substr(fut, 4, 6), ' ', substr(fut, 8, 11), '-', substr(fut, 13, 16))

		maps_biomass_change[[length(maps_biomass_change) + 1]] <- map_biomass_change(
			out_dir = out_dir,
			pred_vect_nam = pred_vect_nam,
			response_var = response_var,
			data_occs = data_occs,
			data_biomass = data_biomass,
			title = title,
			subtitle = subtitle,
			ag_core_quant = 0.95
		)

	}

say('###############################', pre = 1)
say('### current map: OCCURRENCE ###')
say('###############################')

	map_occ_sq <- map_occurrence(
		pred_vect_nam = pred_vect_nam,
		response_var = 'N_ag_county_sq',
		data_occs = data_occs,
		data_biomass = data_biomass,
		title = bquote('Present-day distribution of ' * italic('Andropogon gerardi') * ' abundance'),
		subtitle = 'Abundance ~ f(environment, biomass ~ f(environment))  | 1961-2020',
		ag_core_quant = ag_core_quant
	)

	ggsave(plot = map, filename = paste0(out_dir, '/map_abundance_present.png'), width = 12, height = 10, dpi = 600)

say('###############################', pre = 1)
say('### future maps: OCCURRENCE ###')
say('###############################')

	maps_occs_fut <- list()
	for (fut in futs) {

		say(fut)

		title <- bquote('Future distribution of ' * italic('Andropogon gerardi') * ' abundance')
		subtitle <- paste0('Abundance ~ f(environment, biomass ~ f(environment))  | SSP ', substr(fut, 4, 6), ' ', substr(fut, 8, 11), '-', substr(fut, 13, 16))
		response_var <- paste0('N_ag_county_', fut)

		maps_occs_fut[[length(maps_occs_fut) + 1]] <- map_occurrence(
			pred_vect_nam = pred_vect_nam,
			response_var = response_var,
			data_occs = data_occs,
			data_biomass = data_biomass,
			title = title,
			subtitle = subtitle,
			ag_core_quant = ag_core_quant
		)

		ggsave(plot = maps_occs_fut[[length(maps_occs_fut)]], filename = paste0(out_dir, '/map_abundance_', fut, '.png'), width = 12, height = 10, dpi = 600)

	}

say('######################################', pre = 1)
say('### future change maps: OCCURRENCE ###')
say('######################################')

	for (fut in futs) {

		say(fut)

		response_var <- paste0('N_ag_county_', fut)

		map <- map_occurrence_change(
			pred_vect_nam = pred_vect_nam,
			response_var = response_var,
			data_occs = data_occs,
			data_biomass = data_biomass,
			title = bquote('Change in ' * italic('Andropogon gerardi') * ' abundance'),
			subtitle = paste0('Abundance ~ f(environment, biomass ~ f(environment))  | SSP ', substr(fut, 4, 6), ' ', substr(fut, 8, 11), '-', substr(fut, 13, 16)),
			ag_core_quant = 0.95
		)

		ggsave(plot = map, filename = paste0(out_dir, '/map_abundance_change_', fut, '.png'), width = 12, height = 10, dpi = 600)

	}

say('########################')
say('### cross validation ###')
say('########################')

	if (!trial) {

		sink(paste0(out_dir, '/cross_validation.txt'), split = TRUE)
		say('GEO-FOLD CROSS VALIDATION')
		say(date(), post = 2)
		say('Number of geo-folds: ', k_folds, post = 2)

		# CV for BIOMASS and OCCURRENCE
		folds_fx <- folds_for_biomass_and_occurrences # change according to response data type we're using

		cv <- runCrossValidate(
			MCMCconfiguration = conf,
			k = k_folds, # universal setting
			foldFunction = folds_fx,
			lossFunction = 'MSE',
			MCMCcontrol = list(niter = niter, nburnin = nburnin),
			returnSamples = FALSE,
			nCores = 1,
			nBootReps = 200,
			silent = FALSE
		)

		say('CV statistics for BIOMASS and OCCURRENCE:', post = 2)
		print(cv)

		# CV for BIOMASS
		folds_fx <- folds_for_biomass # change according to response data type we're using

		cv <- runCrossValidate(
			MCMCconfiguration = conf,
			k = k_folds, # universal setting
			foldFunction = folds_fx,
			lossFunction = 'MSE',
			MCMCcontrol = list(niter = niter, nburnin = nburnin),
			returnSamples = FALSE,
			nCores = 1,
			nBootReps = 200,
			silent = TRUE
		)

		say('CV statistics for BIOMASS:', post = 2)
		print(cv)

		# CV for OCCURRENCE
		folds_fx <- folds_for_occurrence # change according to response data type we're using

		cv <- runCrossValidate(
			MCMCconfiguration = conf,
			k = k_folds, # universal setting
			foldFunction = folds_fx,
			lossFunction = 'MSE',
			MCMCcontrol = list(niter = niter, nburnin = nburnin),
			returnSamples = FALSE,
			nCores = 1,
			nBootReps = 200,
			silent = TRUE
		)

		say('CV statistics for OCCURRENCE:', post = 2)
		print(cv)

		sink()

	}

say(date())
say('FINIS!', deco = '+', level = 1)

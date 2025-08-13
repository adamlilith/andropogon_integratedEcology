### MODELING ANDROPOGON GERARDI DISTRIBUTION, PHENOTYPE, PHYSIOLOGY, GENOTYPE, and ASSOCIATED MICROBIAL COMMUNITIES
### Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2023-12
###
### This script constructs an integrated biomass/species distribution model for Andropogon gerardi. The model is written for nimble.
### The SDM assumes an N-mixture model where observations are number of herbarium records of AG per county, and the quantity of interest is the (latent) number of AG per county. Number of observations is assumed to follow a binomial distribution from the latent abundance, with a probability that is a function of county area and number of Poaceae records in the county. Mean latent abundance is assumed to follow a lognormal distribution where the mean is a function of environmental covariates.
### The biomass component is based on field observations at 26 sites at which the aboveground biomass of several plants were measured. The model assumes that the distribution of biomass *at a site* is drawn from a lognormal distribution, where the mean is a function of environmental covariates. Each ramet's biomass is a draw from the site-level lognormal distribution. The precision of the lognormal is assumed to be the same across all sites and is estimated from an uninformed prior (i.e., not a function of covariates).
### The integration assumes that the "normal" component of the lognormal distribution of the SDM and "normal" component of the lognormal of biomass are drawn from a multivariate normal distribution. Off-diagonal elements of the VCV matrix are estimated.
###
###
### source('C:/Kaji/R/andropogon_integratedEcology/sdm_pdm/sdm_pdm_09_model_[integration_mvnorm]_[occurrence]_[biomass_univariate_lognormal].r')
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

	trial <- TRUE # TRUE for testing
	# trial <- FALSE # TRUE for testing

	# calib <- TRUE # use just counties with non-NA Poaceae for calibration region
	calib <- FALSE # use all of North America for calibration region

	out_dir <- paste0('./outputs_loretta/integrated_sdm_pdm/models_integrated_occurrence_biomass/integrated_[occurrence_bio1^2_bio12^2_bio15^2_[bias_***]]_[biomass_lognormal_bio12]', ifelse(trial, '_TRIAL', ''), '/')

	# formula for how aspects of species responds to environment
	formula_occs <- ~ 1 + bio1 + bio12 + bio15 + I(bio1^2) + I(bio12^2) + I(bio15^2) # response of occurrence to climate and soil
	# formula_occs_bias <- ~ 1 + area_km2_log10 + n_poaceae_log10p1 + I(area_km2_log10^2) + I(n_poaceae_log10p1^2) + area_km2_log10:n_poaceae_log10p1 # sampling bias for AG records
	formula_occs_bias <- ~ 1 + area_km2_log10 + n_poaceae_log10p1 # sampling bias for AG records
	formula_biomass_mu <- ~ 1 + bio12

	dirCreate(out_dir)
	formula <- list(
		formula_occs = formula_occs,
		formula_occs_bias = formula_occs_bias,
		formula_biomass_mu = formula_biomass_mu
	)
	saveRDS(formula, paste0(out_dir, '/formula.rds'))

	# # ### MCMC settings
	# # need these settings to achieve independent samples using bio1^2, bio12^2, bio15^
	# niter <- 240000 * 200
	# nburnin <- 40000 * 200
	# thin <- 200 * 200
	# nchains <- 4
	# waic <- TRUE

	# ### MCMC settings
	# niter <- 240000 * 4
	# nburnin <- 40000 * 4
	# thin <- 200 * 4
	# nchains <- 4
	# waic <- TRUE

	# ### MCMC settings FOR TUNING
	# niter <- 44000
	# nburnin <- 4000
	# thin <- 40
	# nchains <- 4
	# waic <- FALSE

	# ### MCMC settings FOR TESTING
	# niter <- 22000
	# nburnin <- 2000
	# thin <- 20
	# nchains <- 2
	# waic <- TRUE

	### MCMC settings FOR TESTING
	niter <- 30000
	nburnin <- 10000
	thin <- 1
	nchains <- 2
	waic <- TRUE

#############
### model ###
#############

	sink(paste0(out_dir, '/runtime_log.txt'), split = TRUE)
	say('MODELING OCCURRENCE')
	say('Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | ', date(), post = 1)

	### data collation
	##################
	say('data collation', post = 1)

	say('This script constructs an integrated biomass/species distribution model for Andropogon gerardi. The model is written for nimble.', breaks = 80)
	say('The SDM assumes an N-mixture model where observations are number of herbarium records of AG per county, and the quantity of interest is the (latent) number of AG per county. Number of observations is assumed to follow a binomial distribution from the latent abundance, with a probability that is a function of county area and number of Poaceae records in the county. Mean latent abundance is assumed to follow a lognormal distribution where the mean is a function of environmental covariates.', breaks = 80)
	say('The biomass component is based on field observations at 26 sites at which the aboveground biomass of several plants were measured. The model assumes that the distribution of biomass *at a site* is drawn from a lognormal distribution, where the mean is a function of environmental covariates. Each ramet\'s biomass is a draw from the site-level lognormal distribution. The precision of the lognormal is assumed to be the same across all sites and is estimated from an uninformed prior (i.e., not a function of covariates).', breaks = 80)
	say('The integration assumes that the "normal" component of the lognormal distribution of the SDM and "normal" component of the lognormal of biomass are drawn from a multivariate normal distribution. Off-diagonal elements of the VCV matrix are estimated. The VCV is estimated using Cholesky decomposition.', breaks = 80)


	say('MCMC settings:', level = 2)
	say('niter ........................ ', niter)
	say('nburnin ...................... ', nburnin)
	say('thin ......................... ', thin)
	say('nchains ...................... ', nchains)
	say('formula_occs ................. ', paste(as.character(formula_occs), collapse = ' '))
	say('formula_occs_bias ............ ', paste(as.character(formula_occs_bias), collapse = ' '))
	say('formula_biomass_mu .............. ', paste(as.character(formula_biomass_mu), collapse = ' '))
	say('trial ........................ ', trial, post = 2)

	say('out_dir')
	say(out_dir, post = 2)

	formulae <- list(
		formula_occs = formula_occs,
		formula_occs_bias = formula_occs_bias,
		formula_biomass_mu = formula_biomass_mu
	)

	########################s
	### data preparation ###
	########################

	data_biomass_mu <- prepare_biomass(formula = formula_biomass_mu, n_response_curve_values = n_response_curve_values, calib = calib)
	data_occs <- prepare_occurrences(formula_occs = formula_occs, formula_occs_bias = formula_occs_bias, n_response_curve_values = n_response_curve_values, psa_quant = psa_quant, calib = calib)

	#########################
	### inputs for nimble ###
	#########################

	say('Inputs:', level = 2)
	data <- list(
		y_n_ag = data_occs$y_n_ag,				# number of AG observations in each county
		y_biomass = data_biomass_mu$y_biomass		# biomass of individual plants

	)

	constants <- list(

		n_responses = 2, # total number of response variables (occurrence = 1, biomass = 1, etc.)

		### occurrences
		n_terms_occs = data_occs$n_terms_occs, # number of terms in formula for occurrence model (including intercept)
		w_occs_bias = data_occs$w_occs_bias, # model matrix of sampling bias of AG observed occurrences
		n_terms_occs_bias = data_occs$n_terms_occs_bias, # number of terms in sampling bias model

		n_covariates_occs = data_occs$n_covariates_occs, # number of covariates in formula for occurrence model
		resp_curves_x_occs = data_occs$resp_curves_x_occs, # response curve array for occurrences vs environment

		counties_x_occs_sq = data_occs$counties_x_occs_sq,
		counties_x_occs_ssp245_2041_2070 = data_occs$counties_x_occs_ssp245_2041_2070,
		counties_x_occs_ssp245_2071_2100 = data_occs$counties_x_occs_ssp245_2071_2100,
		counties_x_occs_ssp370_2041_2070 = data_occs$counties_x_occs_ssp370_2041_2070,
		counties_x_occs_ssp370_2071_2100 = data_occs$counties_x_occs_ssp370_2071_2100,

		### biomass
		x_by_site_biomass = data_biomass_mu$x_by_site_biomass, # MM with covariates for biomass (scaled)
		n_biomass = data_biomass_mu$n_biomass, # number of biomass observations
		site_index_biomass = data_biomass_mu$site_index_biomass, # index of sampled site for each row in biomass data
		n_terms_biomass_mu = data_biomass_mu$n_terms_biomass, # number of terms in formula for biomass model (including intercept)

		# n_covariates_biomass = data_biomass_mu$n_covariates_biomass, # do not need if just one predictor for biomass
		resp_curves_x_biomass = data_biomass_mu$resp_curves_x_biomass, # response curve array for biomass

		counties_x_biomass_sq = data_biomass_mu$counties_x_biomass_sq,
		counties_x_biomass_ssp245_2041_2070 = data_biomass_mu$counties_x_biomass_ssp245_2041_2070,
		counties_x_biomass_ssp245_2071_2100 = data_biomass_mu$counties_x_biomass_ssp245_2071_2100,
		counties_x_biomass_ssp370_2041_2070 = data_biomass_mu$counties_x_biomass_ssp370_2041_2070,
		counties_x_biomass_ssp370_2071_2100 = data_biomass_mu$counties_x_biomass_ssp370_2071_2100,

		### general

		# counties
		n_counties = data_occs$n_counties, # number of counties in the dataset

		# sites sampled for phenotyping
		n_pheno_sites = data_biomass_mu$n_pheno_sites, # number of phenotype sample sites

		# number of values to estimate for response curves
		n_response_curve_values = n_response_curve_values # number of values in response curve array


	)

	### initializations for SDM

	N_inits_calib <- data_occs$y_n_ag * 2
	N_inits_all_counties <- data_occs$ag_vect_sq$n_andropogon_gerardi * 2

	prelim_data <- cbind(data_occs$w_occs_bias, data_occs$counties_x_occs_sq[ , 2:ncol(data_occs$counties_x_occs_sq)])
	prelim_model <- glm.fit(prelim_data, y = data_occs$y_n_ag, family = poisson(log))

	alpha_occs_inits <- prelim_model$coefficients[c('(Intercept)', 'area_km2_log10', 'n_poaceae_log10p1')]
	beta_occs_inits <- prelim_model$coefficients[c('(Intercept)', data_occs$terms_occs)]

	### initializations for biomass model
	raw_data_biomass <- data_biomass_mu$raw_data_biomass
	preds <- data_biomass_mu$covariates
	raw_data_biomass <- raw_data_biomass[ , c('Biomass', ..preds)]
	y <- raw_data_biomass$Biomass
	x <- scale(raw_data_biomass[ , ..preds])
	x <- as.data.table(x)
	mm <- model.matrix(formula_biomass_mu, x)
	biomass_model_prelim <- glm.fit(x = mm, y = y, family = Gamma(link = 'log'))
	beta_biomass_mu_inits <- biomass_model_prelim$coefficients

	observed_site_mean_biomass <- data_biomass_mu$raw_data_biomass[ , .(mean_biomass = mean(Biomass)), by = SITE][['mean_biomass']]


	inits <- list(

		### occurrence
		y_n_ag_sim = data_occs$y_n_ag, # simulated values of observed number of AG (for DHARMa residuals)
		alpha_occs = alpha_occs_inits, # intercept, area, # of Poaceae
		beta_occs_mu = beta_occs_inits, # occurrence ~ environment coefficients (including intercept)
		beta_occs_mu_reg_log = 1, # regularization penalty for occurrence betas

		N = N_inits_calib, # number of latent AG in calibration counties
		N_ag_county_sq = N_inits_calib,
		N_ag_county_ssp245_2041_2070 = N_inits_all_counties,
		N_ag_county_ssp245_2071_2100 = N_inits_all_counties,
		N_ag_county_ssp370_2041_2070 = N_inits_all_counties,
		N_ag_county_ssp370_2071_2100 = N_inits_all_counties,

		Phi_county_ssp245_2041_2070 = matrix(0, ncol = 2, nrow = data_occs$n_counties),
		Phi_county_ssp245_2071_2100 = matrix(0, ncol = 2, nrow = data_occs$n_counties),
		Phi_county_ssp370_2041_2070 = matrix(0, ncol = 2, nrow = data_occs$n_counties),
		Phi_county_ssp370_2071_2100 = matrix(0, ncol = 2, nrow = data_occs$n_counties),

		### biomass
		y_biomass_sim = data_biomass_mu$y_biomass, # simulated values for biomass (for DHARMa residuals)
		mu_biomass_site = observed_site_mean_biomass, # mean site-level biomass

		site_biomass_plant_sigma_log = log(sd(data_biomass_mu$y_biomass)),
		sigma_biomass_among_sites = log(sd(data_biomass_mu$y_biomass)),

		beta_biomass_mu = beta_biomass_mu_inits,
		beta_biomass_mu_reg = 1, # regularization penalty for biomass betas

		### integration
		L = matrix(c(1, 0.5, 0.5, 1), nrow = 2, byrow = TRUE),
		eta = 1

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
	model_code <- nimbleCode({

		### INTEGRATION
		###############

		# LKJ prior on correlation matrix between abundance and traits
		L[1:n_responses, 1:n_responses] ~ dlkj_corr_cholesky(eta = eta, p = n_responses)

		# MVN precision matrix
		Omega[1:n_responses, 1:n_responses] <- L[1:n_responses, 1:n_responses] %*% t(L[1:n_responses, 1:n_responses])

		# prior for hyperparameter of Cholesky
		eta ~ dgamma(2, 1)

		for (i in 1:n_counties) {

			### OCCURRENCE

			# occurrence: actual abundance (latent--unobserved)
			N[i] ~ dpois(lambda_mu_sq[i])
			N_ag_county_sq[i] ~ dpois(lambda_mu_sq[i])

			# occurrence: observed number of AG and sampling bias
			logit(p[i]) <- inprod(alpha_occs[1:n_terms_occs_bias], w_occs_bias[i, 1:n_terms_occs_bias])
			y_n_ag[i] ~ dbin(prob = p[i], size = N[i])

			# occurrence: simulate observations of AG occurrences for DHARMa residuals
			y_n_ag_sim[i] ~ dbin(prob = p[i], size = N[i])

			# occurrence: likelihood
			log_lik_occs_y[i] <- dbinom(y_n_ag[i], prob = p[i], size = N[i], log = 1)

			# occurrence: relationship between expected (latent) abundance and environment assuming NORMAL distribution
			log(lambda_mu_sq[i]) <- Phi[i, 1]
			phi_occs_mu_county_sq[i] <- inprod(beta_occs_mu[1:n_terms_occs], counties_x_occs_sq[i, 1:n_terms_occs])

			### BIOMASSS

			# biomass: site-level mean for county
			phi_biomass_mu_county_sq[i] <- inprod(beta_biomass_mu[1:n_terms_biomass_mu], counties_x_biomass_sq[i, 1:n_terms_biomass_mu])
			log(mu_biomass_county_sq[i]) <- Phi[i, 2]

			### INTEGRATION

			# means
			mvnorm_mu[i, 1:n_responses] <- c(phi_occs_mu_county_sq[i], phi_biomass_mu_county_sq[i])

			# likelihood
			Phi[i, 1:n_responses] ~ dmnorm(mean = mvnorm_mu[i, 1:n_responses], cholesky = Omega[1:n_responses, 1:n_responses], prec_param = 0)

		}

		# INTEGRATED: posterior samplers for geographic predictions in future
		for (i in 1:n_counties) {

			### posterior predictive nodes for ssp245_2041_2070
			Phi_county_ssp245_2041_2070[i, 1:n_responses] ~ dmnorm(mean = mvnorm_mu_ssp245_2041_2070[i, 1:n_responses], cholesky = Omega[1:n_responses, 1:n_responses], prec_param = 0)

			mvnorm_mu_ssp245_2041_2070[i, 1:n_responses] <-
				c(phi_occs_mu_county_ssp245_2041_2070[i], phi_biomass_mu_county_ssp245_2041_2070[i])

			phi_occs_mu_county_ssp245_2041_2070[i] <- inprod(beta_occs_mu[1:n_terms_occs], counties_x_occs_ssp245_2041_2070[i, 1:n_terms_occs])

			phi_biomass_mu_county_ssp245_2041_2070[i] <- inprod(beta_biomass_mu[1:n_terms_biomass_mu], counties_x_biomass_ssp245_2041_2070[i, 1:n_terms_biomass_mu])

			log(lambda_county_ssp245_2041_2070[i]) <- Phi_county_ssp245_2041_2070[i, 1]
			log(mu_biomass_county_ssp245_2041_2070[i]) <- Phi_county_ssp245_2041_2070[i, 2]

			N_ag_county_ssp245_2041_2070[i] ~ dpois(lambda_county_ssp245_2041_2070[i])

			### posterior predictive nodes for ssp245_2071_2100
			Phi_county_ssp245_2071_2100[i, 1:n_responses] ~ dmnorm(mean = mvnorm_mu_ssp245_2071_2100[i, 1:n_responses], cholesky = Omega[1:n_responses, 1:n_responses], prec_param = 0)

			mvnorm_mu_ssp245_2071_2100[i, 1:n_responses] <-
				c(phi_occs_mu_county_ssp245_2071_2100[i], phi_biomass_mu_county_ssp245_2071_2100[i])

			phi_occs_mu_county_ssp245_2071_2100[i] <- inprod(beta_occs_mu[1:n_terms_occs], counties_x_occs_ssp245_2071_2100[i, 1:n_terms_occs])

			phi_biomass_mu_county_ssp245_2071_2100[i] <- inprod(beta_biomass_mu[1:n_terms_biomass_mu], counties_x_biomass_ssp245_2071_2100[i, 1:n_terms_biomass_mu])

			log(lambda_county_ssp245_2071_2100[i]) <- Phi_county_ssp245_2071_2100[i, 1]
			log(mu_biomass_county_ssp245_2071_2100[i]) <- Phi_county_ssp245_2071_2100[i, 2]

			N_ag_county_ssp245_2071_2100[i] ~ dpois(lambda_county_ssp245_2071_2100[i])

			### posterior predictive nodes for ssp370_2041_2070
			Phi_county_ssp370_2041_2070[i, 1:n_responses] ~ dmnorm(mean = mvnorm_mu_ssp370_2041_2070[i, 1:n_responses], cholesky = Omega[1:n_responses, 1:n_responses], prec_param = 0)

			mvnorm_mu_ssp370_2041_2070[i, 1:n_responses] <-
				c(phi_occs_mu_county_ssp370_2041_2070[i], phi_biomass_mu_county_ssp370_2041_2070[i])

			phi_occs_mu_county_ssp370_2041_2070[i] <- inprod(beta_occs_mu[1:n_terms_occs], counties_x_occs_ssp370_2041_2070[i, 1:n_terms_occs])

			phi_biomass_mu_county_ssp370_2041_2070[i] <- inprod(beta_biomass_mu[1:n_terms_biomass_mu], counties_x_biomass_ssp370_2041_2070[i, 1:n_terms_biomass_mu])

			log(lambda_county_ssp370_2041_2070[i]) <- Phi_county_ssp370_2041_2070[i, 1]
			log(mu_biomass_county_ssp370_2041_2070[i]) <- Phi_county_ssp370_2041_2070[i, 2]

			N_ag_county_ssp370_2041_2070[i] ~ dpois(lambda_county_ssp370_2041_2070[i])

			### posterior predictive nodes for ssp370_2071_2100
			Phi_county_ssp370_2071_2100[i, 1:n_responses] ~ dmnorm(mean = mvnorm_mu_ssp370_2071_2100[i, 1:n_responses], cholesky = Omega[1:n_responses, 1:n_responses], prec_param = 0)

			mvnorm_mu_ssp370_2071_2100[i, 1:n_responses] <-
				c(phi_occs_mu_county_ssp370_2071_2100[i], phi_biomass_mu_county_ssp370_2071_2100[i])

			phi_occs_mu_county_ssp370_2071_2100[i] <- inprod(beta_occs_mu[1:n_terms_occs], counties_x_occs_ssp370_2071_2100[i, 1:n_terms_occs])

			phi_biomass_mu_county_ssp370_2071_2100[i] <- inprod(beta_biomass_mu[1:n_terms_biomass_mu], counties_x_biomass_ssp370_2071_2100[i, 1:n_terms_biomass_mu])

			log(lambda_county_ssp370_2071_2100[i]) <- Phi_county_ssp370_2071_2100[i, 1]
			log(mu_biomass_county_ssp370_2071_2100[i]) <- Phi_county_ssp370_2071_2100[i, 2]

			N_ag_county_ssp370_2071_2100[i] ~ dpois(lambda_county_ssp370_2071_2100[i])

		}

	### OCCURRENCE
	##############

		# OCCURRENCE: weakly regularized or regularized priors for relationship to environment
		# ddexp(): rate = 0.7675 allows |beta| >= 3 10% of time, 0.9986 allows |beta| >= 3 5% of time
		
		beta_occs_mu_reg_log ~ dnorm(0, sd = 2.5)
		beta_occs_mu_reg <- exp(beta_occs_mu_reg_log)
		beta_occs_mu[1] ~ dnorm(0, sd = 20)
		for (i in 2:n_terms_occs) {
			# beta_occs_mu[i] ~ dnorm(0, sd = 20)
			beta_occs_mu[i] ~ ddexp(0, rate = beta_occs_mu_reg)
		}

		# OCCURRENCE: priors for sampling bias
		alpha_occs[1] ~ dnorm(0, sd = 100)
		for (i in 2:n_terms_occs_bias) {
			alpha_occs[i] ~ dnorm(0, sd = 10)
			# alpha_occs[i] ~ ddexp(0, rate = 0.7675)
		}

		log_lik_occs <- sum(log_lik_occs_y[1:n_counties])

		# OCCURRENCE: posterior predictive sampler for ENVIRONMENTAL response curves
		# We're assuming occurrence responds to two or more environmental predictors, so the response curve "x" is an array with one "page" per predictor and output is a matrix with one column per predictor

		for (i in 1:n_covariates_occs) {

			for (j in 1:n_response_curve_values) {

				log(response_curves_occs_mu[j, i]) <-
					inprod(beta_occs_mu[1:n_terms_occs], resp_curves_x_occs[j, 1:n_terms_occs, i])

				# exp_response_curves_occ_lambda_vs_env[j, i] <- exp(response_curves_occs_mu[j, i])

			}

		}

	### BIOMASS
	###########

		# BIOMASS: priors for relationship of mean and variance to environment
		beta_biomass_mu_reg_log ~ dnorm(0, sd = 2.5)
		beta_biomass_mu_reg <- exp(beta_biomass_mu_reg_log)
		beta_biomass_mu[1] ~ dnorm(0, sd = 20)
		for (i in 2:n_terms_biomass_mu) {
			# beta_biomass_mu[i] ~ dnorm(0, sd = 10) # broad prior
			beta_biomass_mu[i] ~ ddexp(0, rate = beta_biomass_mu_reg)
		}

		# site_biomass_plant_sigma_log ~ dhalfflat()
		# site_biomass_mean_sigma ~ dhalfflat()

		# site_biomass_plant_sigma_log ~ dunif(0, 100)
		# site_biomass_mean_sigma ~ dunif(0, 20)

		site_biomass_plant_sigma_log ~ dnorm(0, sd = 10) # broad half-Cauchy
		site_biomass_plant_sigma_log <- exp(site_biomass_plant_sigma_log)

		sigma_biomass_among_sites ~ dnorm(0, sd = 2.5) # narrow half-Cauchy
		site_biomass_mean_sigma <- exp(sigma_biomass_among_sites)

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
			log_lik_biomass_y[i] <- dlnorm(y_biomass[i], meanlog = mu_biomass_site[site_index_biomass[i]], sdlog = site_biomass_mean_sigma, log = 1)

		}

		log_lik_biomass <- sum(log_lik_biomass_y[1:n_biomass])

		# BIOMASS: posterior predictive sampler for response curves: site-level mean
		# NB we assume 1 predictor for biomass, so the response curve "x" is a matrix, not an array (as with >1 predictor)
		for (i in 1:n_response_curve_values) {

			response_curves_biomass_mu[i] <-
				exp(inprod(beta_biomass_mu[1:n_terms_biomass_mu], resp_curves_x_biomass[i, 1:n_terms_biomass_mu]))
		}

		### overall
		###########

		log_lik <- log_lik_occs + log_lik_biomass

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
	say('model$calculate(): ', model$calculate())
	model$simulate()
	say('model$calculate(): ', model$calculate())

	say('configureMCMC():', level = 2)

	monitors_coeffs_not_indexed <- c(
		'site_biomass_plant_sigma_log', 'site_biomass_mean_sigma',
		'beta_biomass_mu_reg', 'beta_occs_mu_reg',
		'eta'
	)

	# monitors for coefficients that have no indexing
	monitors_coeffs_single_index <- c(
		'beta_occs_mu', 'alpha_occs',
		'beta_biomass_mu'
	)

	monitors_coeffs_double_index <- c(
		'L', 'Omega'
	)

	monitors_derived_not_indexed <- c(
		'log_lik',
		'log_lik_occs', 'log_lik_biomass'
	)

	monitors_derived_indexed <- c(
		'mu_biomass_site'
	)

	monitors_geog <- c(
		'N_ag_county_sq', 'mu_biomass_county_sq',
		'N_ag_county_ssp245_2041_2070', 'mu_biomass_county_ssp245_2041_2070',
		'N_ag_county_ssp245_2071_2100', 'mu_biomass_county_ssp245_2071_2100',
		'N_ag_county_ssp370_2041_2070', 'mu_biomass_county_ssp370_2041_2070',
		'N_ag_county_ssp370_2071_2100', 'mu_biomass_county_ssp370_2071_2100'
	)

	monitors_dharma <- c(
		'y_n_ag_sim',
		'y_biomass_sim'
	)

	monitors_resp_curves <- c(
		'response_curves_occs_mu',
		'response_curves_biomass_mu'
	)

	monitors <- c(monitors_coeffs_not_indexed, monitors_coeffs_single_index, monitors_coeffs_double_index, monitors_derived_not_indexed, monitors_derived_indexed, monitors_geog, monitors_dharma, monitors_resp_curves)
	# monitors <- c(monitors_coeffs_not_indexed, monitors_coeffs_single_index, monitors_coeffs_double_index, monitors_derived_not_indexed, monitors_derived_indexed)

	conf <- configureMCMC(
		model,
		monitors = monitors,
		print = TRUE,
		enableWAIC = waic
	)

	# add no U-turn sampler (Hamiltonian Monte Carlo)
	# nut_samplers <- c(monitors_coeffs_not_indexed, monitors_coeffs_single_index, monitors_coeffs_double_index)
	# conf$addSampler(target = nut_samplers, type = 'NUTS')
	# say('NUTS sampler added to all continuous parameters.')

	# # slice samplers for correlated parameters
	# conf$removeSamplers('beta_occs_mu[1]')
	# conf$removeSamplers('beta_occs_vs_biomass')
	# conf$addSampler(target = 'beta_occs_mu[1]', type = 'slice')
	# conf$addSampler(target = 'beta_occs_vs_biomass[1]', type = 'slice')
	# conf$addSampler(target = 'beta_occs_vs_biomass[2]', type = 'slice')
	# say('Slice sampler added to beta_occs_mu[1] and beta_occs_vs_biomass[1:2].')

	# # RW block samplers for correlated parameters
	# conf$removeSamplers('beta_occs_mu[1]')
	# conf$removeSamplers('beta_occs_vs_biomass')
	# conf$addSampler(target = c('beta_occs_mu[1]', 'beta_occs_vs_biomass[1]', 'beta_occs_vs_biomass[2]'), type = 'RW_block')
	# say('RW_block sampler added to beta_occs_mu[1] and beta_occs_vs_biomass[1:2].')

	# # AF_slice samplers for correlated parameters
	# conf$removeSamplers('alpha_occs')
	# conf$removeSamplers('beta_occs_mu')
	# conf$addSampler(target = c('alpha_occs', 'beta_occs_mu'), type = 'AF_slice')
	# say('AF_slice sampler added to alpha_occs and beta_occs_mu.')

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

say('#########################')
say('### model diagnostics ###')
say('#########################')

	workflow_postmodel_generic(facet = 'occurrence & biomass', formulae = formulae, descrip = 'occurrence: Poisson lognormal homoscedastic; biomass: gamma homoscedastic', out_dir = out_dir)
	workflow_postmodel_occurrence(homoscedastic = TRUE, pred_vect_nam = pred_vect_nam, out_dir = out_dir)
	workflow_postmodel_biomass(homoscedastic = TRUE, pred_vect_nam = pred_vect_nam, out_dir = out_dir)


say(date())
say('FINIS!', deco = '+', level = 1)

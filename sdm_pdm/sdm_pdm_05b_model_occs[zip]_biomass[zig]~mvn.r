### MODELING ANDROPOGON GERARDI DISTRIBUTION, PHENOTYPE, PHYSIOLOGY, GENOTYPE, and ASSOCIATED MICROBIAL COMMUNITIES
### Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2023-12
###
### This is a joint distribution-trait model for Andropogon gerardi. For the occurrence model, "county" is the observational unit but for biomass "site" is the observational unit. It assumes (latent) abundance follows a zero-inflated Poisson distribution, and the observed number of AG is a binomial distribution where the probability of observing AG is optionally a function of county area and the number of Poaceae recorded in the county (including AG). The expected abundance drawn from a normal distribution where the mean value is given by a function of environmental predictors (climate, soil, etc.). The biomass model assumes that the site-level mean biomass of plants is a function of environmental covariates, and that the biomass of plant biomass at a site is drawn from a zero-inflated gamma distribution with this environmentally-determined mean. The occurrence and biomass components communicate via a latent multivariate-normal parameter Phi (estimated at the county and site levels), which provides the link to the occurrence- or biomass-environment relationship. The normal draws of Phi are appropriately transformed to represent the mean (expected) value of site- or county-level abundance and biomass. The probability of zero inflation is a function of environmental variables. The model is run using nimble.
###
### source('C:/Kaji/R/andropogon_integratedEcology/sdm_pdm/sdm_pdm_05b_model_occs[zip]_biomass[zig]~mvn.r')
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
	crossvalidate <- TRUE
	# crossvalidate <- FALSE

	# log BIOs 12-14 and 16-19?
	# log_precip <- FALSE
	log_precip <- TRUE

	### formula for how aspects of species responds to environment

	formula_occs <- ~ 1 + bio1 + bio12 + bio15 + I(bio1^2) + I(bio12^2) + I(bio15^2) # response of occurrence to climate and soil
	occs_filename <- 'bio1^2_bio12^2_bio15^2'

	formula_occs_bias <- ~ 1 # sampling bias for AG records
	bias_filename <- '1'

	formula_biomass <- ~ 1 + bio12 # response of biomass to environment
	biomass_filename <- 'bio12'

	formula_psi <- ~ 1 + bio1 + bio12 + bio15 + I(bio1^2) + I(bio12^2) + I(bio15^2)
	psi_filename <- 'bio1^2_bio12^2_bio15^2'

	### output folder and bias formula
	out_dir <- paste0('./outputs_loretta/integrated_sdm_pdm/models_occs_biomass/', ifelse(trial, 'TRIAL_', ''), '[occs~zip[', occs_filename, '_bias~', bias_filename, ']]_[biomass~', biomass_filename, ']_[psi_', psi_filename, ']', ifelse(log_precip, '_log_precip', ''))

	if (!trial) {

		### MCMC settings
		# need these settings to achieve independent samples using bio1^2, bio12^2, bio15^
		niter <- 1600000
		nchains <- 1
		# nchains <- 4
		waic <- TRUE

	} else {
		
		### MCMC settings FOR TESTING
		niter <- 4000
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

	say('This is a joint distribution-trait model for Andropogon gerardi. For the occurrence model, "county" is the observational unit but for biomass "site" is the observational unit. It assumes (latent) abundance follows a zero-inflated Poisson distribution, and the observed number of AG is a binomial distribution where the probability of observing AG is optionally a function of county area and the number of Poaceae recorded in the county (including AG). The expected abundance drawn from a normal distribution where the mean value is given by a function of environmental predictors (climate, soil, etc.). The biomass model assumes that the site-level mean biomass of plants is a function of environmental covariates, and that the biomass of plant biomass at a site is drawn from a zero-inflated gamma distribution with this environmentally-determined mean. The occurrence and biomass components communicate via a latent multivariate-normal parameter Phi (estimated at the county and site levels), which provides the link to the occurrence- or biomass-environment relationship. The normal draws of Phi are appropriately transformed to represent the mean (expected) value of site- or county-level abundance and biomass. The probability of zero inflation is a function of environmental variables. The model is run using nimble.', breaks = 60, post = 1)

	say('MCMC settings:', level = 2)
	say('trial ........................ ', trial)
	say('niter ........................ ', niter)
	say('nburnin ...................... ', nburnin)
	say('thin ......................... ', thin)
	say('nchains ...................... ', nchains)
	say('formula_occs ................. ', paste(as.character(formula_occs), collapse = ' '))
	say('formula_occs_bias ............ ', paste(as.character(formula_occs_bias), collapse = ' '))
	say('formula_biomass ........... ', paste(as.character(formula_biomass), collapse = ' '))
	say('formula_psi ................ ', paste(as.character(formula_psi), collapse = ' '))

	say('out_dir')
	say(out_dir, post = 2)

	formulae <- list(
		formula_occs = formula_occs,
		formula_occs_bias = formula_occs_bias,
		formula_biomass = formula_biomass,
		formula_psi = formula_psi
	)
	saveRDS(formula, paste0(out_dir, '/formulae.rds'))

	########################s
	### data preparation ###
	########################

	# data for occurrences at counties
	data_occs <- prepare_occurrence_data(formula_occs = formula_occs, formula_occs_bias = formula_occs_bias, log_precip = log_precip, n_response_curve_values = n_response_curve_values, psa_quant = psa_quant, calib = calib)

	# data for biomass at sites
	data_biomass <- prepare_biomass_data(formula_biomass = formula_biomass, log_precip = log_precip, n_response_curve_values = n_response_curve_values, calib = calib)

	# data for occurrences using site-level environment
	data_occs_as_sites <- prepare_biomass_data(formula_biomass = formula_occs, log_precip = log_precip, n_response_curve_values = n_response_curve_values, calib = calib)

	# data for biomass using county-level environment
	data_biomass_as_counties <- prepare_occurrence_data(formula_occs = formula_biomass, formula_occs_bias = ~ 1, log_precip = log_precip, n_response_curve_values = n_response_curve_values, psa_quant = psa_quant, calib = calib)

	# data for psi at counties
	data_county_psi <- prepare_occurrence_data(formula_occs = formula_psi, formula_occs_bias = ~1, log_precip = log_precip, n_response_curve_values = n_response_curve_values, psa_quant = psa_quant, calib = calib)

	# data for psi at sites
	data_site_psi <- prepare_biomass_data(formula_biomass = formula_psi, log_precip = log_precip, n_response_curve_values = n_response_curve_values, calib = calib)

	#########################
	### inputs for nimble ###
	#########################

	say('Inputs:', level = 2)
	data <- list(
		y_n_ag = data_occs$y_n_ag,					# number of AG observations in each county
		y_biomass = data_biomass$y_biomass		# biomass of individual plants
	)

	constants <- list(
		
		### integration
		n_facets = 2, # number of facets being modeled

		### occurrences
		n_counties_occs_calib = data_occs$n_counties_occs_calib, # number of counties in calibration region
		n_terms_occs = data_occs$n_terms_occs, # number of terms in formula for occurrence model (including intercept)
		n_terms_occs_bias = data_occs$n_terms_occs_bias, # number of terms in sampling bias model
		w_occs_bias = data_occs$w_occs_bias, # model matrix of sampling bias of AG observed occurrences
		
		n_covariates_occs = data_occs$n_covariates_occs, # number of covariates in formula for occurrence model
		n_covariates_occs_bias = data_occs$n_covariates_occs_bias, # number of covariates in formula for occurrence model
		resp_curves_x_occs = data_occs$resp_curves_x_occs, # response curve array for occurrences vs environment
		resp_curves_w_occs = data_occs$resp_curves_w_occs, # response curve array for occurrences vs environment

		x_by_county_occs = data_occs$counties_x_occs_sq,
		x_by_site_occs = data_occs_as_sites$x_by_site_biomass, # MM with covariates for biomass (scaled)

		### biomass
		n_pheno_sites = data_biomass$n_pheno_sites, # number of phenotype sample sites
		x_by_county_biomass = data_biomass_as_counties$counties_x_occs_sq,

		x_by_site_biomass = data_biomass$x_by_site_biomass, # MM with covariates for biomass (scaled)
		n_biomass = data_biomass$n_biomass, # number of biomass observations
		site_index_biomass = data_biomass$site_index_biomass, # index of sampled site for each row in biomass data
		n_terms_biomass = data_biomass$n_terms_biomass, # number of terms in formula for biomass model (including intercept)

		n_covariates_biomass = data_biomass$n_covariates_biomass,
		resp_curves_x_biomass = data_biomass$resp_curves_x_biomass, # response curve array for biomass

		# response curves (general)
		n_response_curve_values = n_response_curve_values, # number of values in response curve array

		# zero inflation
		n_terms_psi = data_county_psi$n_terms_occs, # number of terms in sampling bias model
		n_covariates_psi = data_county_psi$n_covariates_occs, # number of terms in sampling bias model
		x_by_county_psi = data_county_psi$counties_x_occs_sq,
		x_by_site_psi = data_site_psi$x_by_site_biomass,
		resp_curves_x_psi = data_county_psi$resp_curves_x_occs

	)

	constants <- c(constants, constants_shared_occs, constants_shared_biomass, constants_shared_psi)

	# occurrences initializations
	file <- paste0('./outputs_loretta/integrated_sdm_pdm/models_occurrence/[occs_zip[psi~', psi_filename, ']~normal~', occs_filename, '_[bias~1]]', ifelse(log_precip, '_log_precip', ''), '/chains.rds')

	chains_occs <- readRDS(file)
	
	alpha_occs_inits <- mc_extract(chains_occs, 'alpha_occs', j = TRUE)
	beta_occs_mu_inits <- mc_extract(chains_occs, 'beta_occs', j = TRUE)
	beta_psi_inits <- mc_extract(chains_occs, 'beta_psi', j = TRUE)
	
	lambda_mu_sq_inits <- mc_extract(chains_occs, 'lambda_mu_sq', j = TRUE)
	lambda_mu_sq_inits[lambda_mu_sq_inits > 10] <- 10

	lambda_sigma_inits <- mc_extract(chains_occs, 'lambda_sigma')
	log_lambda_sigma_inits <- log(lambda_sigma_inits)

	rm(chains_occs)

	response_curves_occs_mu_inits <- matrix(2, nrow = n_response_curve_values, ncol = data_occs$n_covariates_occs)

	N_inits_calib <- data_occs$y_n_ag * 2
	z_county_inits <- as.numeric(N_inits_calib > 0)

	# biomass initializations
	chains_biomass <- readRDS(paste0('./outputs_loretta/integrated_sdm_pdm/models_biomass/[biomass_zig~normal_',biomass_filename,  ']', ifelse(log_precip, '_log_precip', ''), '/chains.rds'))
	beta_biomass_inits <- mc_extract(chains_biomass, 'beta_biomass', j = TRUE)
	mu_biomass_site_inits <- mc_extract(chains_biomass, 'mu_biomass_site', j = TRUE)
	sigma_biomass_among_sites_inits <- mc_extract(chains_biomass, 'sigma_biomass_among_sites')
	log_sigma_biomass_among_sites_inits <- log(sigma_biomass_among_sites_inits)
	sigma_biomass_within_sites_inits <- mc_extract(chains_biomass, 'sigma_biomass_within_sites')
	log_sigma_biomass_within_sites_inits <- log(sigma_biomass_within_sites_inits)

	rm(chains_biomass)

	lambda_resp_curves_z_inits <- matrix(1, nrow = n_response_curve_values, ncol = data_occs$n_covariates_occs)

	inits <- list(

		### occurrences
		y_n_ag_sim = data_occs$y_n_ag, # simulated values of observed number of AG (for DHARMa residuals)
		lambda_mu_sq = lambda_mu_sq_inits, # expected value of number of AG
		alpha_occs = alpha_occs_inits, # intercept, area, # of Poaceae
		beta_occs = beta_occs_mu_inits, # occurrence ~ environment coefficients (including intercept)

		N = N_inits_calib, # number of latent AG in calibration counties

		log_lambda_resp_curves_mu = response_curves_occs_mu_inits,
		response_curves_occs_mu = response_curves_occs_mu_inits,
		lambda_resp_curves_z = lambda_resp_curves_z_inits,

		### biomass
		mu_biomass_site = mu_biomass_site_inits,
		sigma_biomass_within_sites_log = log_sigma_biomass_within_sites_inits,
		y_biomass_sim = data_biomass$y_biomass, # simulated values for biomass (for DHARMa residuals)
		beta_biomass = beta_biomass_inits,

		### integration
		eta = 1,
		U_star = diag(1, nrow = 2, ncol = 2),
		log_sigmas = c(log_lambda_sigma_inits, log_sigma_biomass_among_sites_inits),
		Phi_county = matrix(1, nrow = data_occs$n_counties, ncol = 2),
		Phi_site = matrix(1, nrow = data_biomass$n_pheno_sites, ncol = 2),

		### zero inflation
		beta_psi = beta_psi_inits,
		z_site = rep(1, data_biomass$n_pheno_sites), # present/absent at each site
		z_county = z_county_inits # present/absent in county

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

		# INTEGRATION: LJK prior for standard deviations and correlations between latent occurrence and biomass processes
		eta ~ dgamma(2, 1)
		U_star[1:n_facets, 1:n_facets] ~ dlkj_corr_cholesky(eta = eta, p = n_facets)
		U[1:n_facets, 1:n_facets] <- uppertri_mult_diag(
			U_star[1:n_facets, 1:n_facets],
			sigmas[1:n_facets]
		)
	  
		# INTEGRATION: standard deviations of latent occurrence and biomass
		for (i in 1:n_facets) {
			log(sigmas[i]) ~ dnorm(0, sd = lambda_sigma_prior_sd) # half-Cauchy
		}

		# correlation matrix
		correlation[1:n_facets, 1:n_facets] <- t(U_star[1:n_facets, 1:n_facets]) %*% U_star[1:n_facets, 1:n_facets]

		# INTEGRATION: county-level
		# NB We have measurements of biomass at the site level, but not county. We thus assume that the estimates of biomass using county-level covariates are indicative of virtual sites that have the same environments as counties. Biomass is still estimated at the site level.
		for (i in 1:n_counties_occs_calib) {

			phis_county[i, 1:n_facets] <- c(phi_occs_mu_sq[i], phi_biomass_mu_sq[i])
			Phi_county[i, 1:n_facets] ~ dmnorm(phis_county[i, 1:n_facets], cholesky = U[1:n_facets, 1:n_facets], prec_param = 0)

			phi_occs_mu_sq[i] <- inprod(beta_occs[1:n_terms_occs], x_by_county_occs[i, 1:n_terms_occs])
			phi_biomass_mu_sq[i] <- inprod(beta_biomass[1:n_terms_biomass], x_by_county_biomass[i, 1:n_terms_biomass])

		}

		# INTEGRATION: site-level calculations
		# This is how the phenotypic variables communicate with abundance. Abundance is naturally estimated at the county level, so we assume that the site-level environmental covariates act like faux counties for abundance, but are indicative of the sites for phenotypic variables.
		for (i in 1:n_pheno_sites) {

			phis_site[i, 1:n_facets] <- c(phi_occs_mu_site[i], phi_biomass_mu_site[i])
			Phi_site[i, 1:n_facets] ~ dmnorm(phis_site[i, 1:n_facets], cholesky = U[1:n_facets, 1:n_facets], prec_param = 0)

			phi_occs_mu_site[i] <- inprod(beta_occs[1:n_terms_occs], x_by_site_occs[i, 1:n_terms_occs])
			phi_biomass_mu_site[i] <- inprod(beta_biomass[1:n_terms_biomass], x_by_site_biomass[i, 1:n_terms_biomass])

		}

		# OCCURRENCE: weakly regularized or regularized priors for relationship to environment
		beta_occs[1] ~ dnorm(0, sd = beta_occs_prior_dnorm_sd_1)
		for (i in 2:n_terms_occs) {
			beta_occs[i] ~ ddexp(0, rate = beta_occs_prior_ddexp_rate)
		}

		# OCCURRENCE: priors for sampling bias
		alpha_occs[1] ~ dnorm(0, sd = alpha_occs_prior_dnorm_sd_1)

		# BIOMASS: priors for probability of presence
		beta_psi[1] ~ dnorm(0, sd = beta_psi_prior_dnorm_sd_1)
		for (i in 2:n_terms_psi) {
			beta_psi[i] ~ ddexp(0, rate = beta_psi_prior_dnorm_sd)
		}

		# OCCURRENCE: likelihood
		for (i in 1:n_counties_occs_calib) {
			
			### actual abundance (latent--unobserved)
			N[i] ~ dTruncPseudoZIP(lambda_mu_sq[i], z = z_county[i])

			### observed number of AG and sampling bias
			y_n_ag[i] ~ dbinom(prob = p[i], size = N[i])

			# simulate observations for DHARMa residuals
			y_n_ag_sim[i] ~ dbinom(prob = p[i], size = N[i])

			# relationship between expected (latent) abundance and environment assuming MV NORMAL distribution
			log(lambda_mu_sq[i]) <- exp(Phi_county[i, 1])

			# probability of zero inflation
			logit(psi_county[i]) <- inprod(beta_psi[1:n_terms_psi], x_by_county_psi[i, 1:n_terms_psi])
			z_county[i] ~ dbern(psi_county[i])

			# likelihood
			log_lik_occs_y[i] <- dbinom(y_n_ag[i], prob = p[i], size = N[i], log = 1)

		}

		log_lik_occs <- sum(log_lik_occs_y[1:n_counties_occs_calib])

		# OCCURRENCE: posterior predictive sampler for ENVIRONMENTAL response curves
		# We're assuming occurrence responds to two or more environmental predictors, so the response curve "x" is an array with one "page" per predictor and output is a matrix with one column per predictor
		for (i in 1:n_covariates_occs) {

			for (j in 1:n_response_curve_values) {

				# NB respond curve is conditional on presence (z := 1)				
				response_curves_occs_mu[j, i] ~ dTruncPseudoZIP(lambda_resp_curves_mu[j, i], z = lambda_resp_curves_z[j, i])
				log(lambda_resp_curves_mu[j, i]) ~ dnorm(phi_lambda_resp_curves_mu[j, i], sd = U[1, 1])
				phi_lambda_resp_curves_mu[j, i] <-
					inprod(beta_occs[1:n_terms_occs], resp_curves_x_occs[j, 1:n_terms_occs, i])
				lambda_resp_curves_z[j, i] ~ dbern(response_curves_psi[j, i])
				logit(response_curves_psi[j, i]) <- inprod(beta_psi[1:n_terms_psi], resp_curves_x_psi[j, 1:n_terms_psi, i])

			}

		}

		# BIOMASS: priors for relationship of site-level mean biomass to environment
		for (i in 1:n_terms_biomass) {
			beta_biomass[i] ~ dnorm(0, sd = beta_biomass_prior_dnorm_sd) # broad prior
		}

		# prior for sd of individual plant biomass on lognormal (~ half-Cauchy), ==> vague
		sigma_biomass_within_sites_log ~ dnorm(0, sd = sigma_biomass_within_sites_log_prior_sd)
		sigma_biomass_within_sites <- exp(sigma_biomass_within_sites_log)

		# BIOMASS: parameters of biomass distribution are latent and functions of environment
		# individual plant biomasses are samples from the site-level distribution defined by the site-level distribution (next chunk after this one)
		for (i in 1:n_pheno_sites) {

			# relationship of biomass to the environment
			
			# site-level mean biomass
			mu_biomass_site[i] <- exp(Phi_site[i, 2])

			# moment matching to get dgamma() parameters
			shape_biomass[i] <- mu_biomass_site[i]^2 / sigma_biomass_within_sites^2
			rate_biomass[i] <- mu_biomass_site[i] / sigma_biomass_within_sites^2
	
			# probability of zero inflation
			logit(psi_site[i]) <- inprod(beta_psi[1:n_terms_psi], x_by_site_psi[i, 1:n_terms_psi])

			# presence/absence at the site level... NB somewhat of a cheat bc we could have absence at county level but presence in site in county
			z_site[i] ~ dbern(psi_site[i])

		}


		# BIOMASS: likelihood of individual plants
		for (i in 1:n_biomass) {

			# likelihood
			y_biomass[i] ~ dZIG(shape = shape_biomass[site_index_biomass[i]], rate = rate_biomass[site_index_biomass[i]], z = z_site[site_index_biomass[i]])

			# simulated values for unconditional DHARMa residuals
			y_biomass_sim[i] ~ dZIG(shape = shape_biomass[site_index_biomass[i]], rate = rate_biomass[site_index_biomass[i]], z = z_site[site_index_biomass[i]])

			log_lik_biomass[i] <- dZIG(y_biomass[i], shape = shape_biomass[site_index_biomass[i]], rate = rate_biomass[site_index_biomass[i]], z = z_site[site_index_biomass[i]])
	
		}

		log_lik_biomass <- sum(log_lik_biomass[1:n_biomass])
		log_lik <- log_lik_occs + log_lik_biomass

	})

	### OCCURRENCE: bias response curves
	# OCCURRENCE RESPONSE CURVES: NO bias covariate
	if (data_occs$n_covariates_occs_bias == 0) {

		bias_response_curve_code <- nimbleCode({

			# OCCURRENCE: likelihood
			for (i in 1:n_counties_occs_calib) {
				logit(p[i]) <- alpha_occs[1]
			}

			# OCCURRENCE: posterior predictive sampler for BIAS response curves
			logit(response_curves_occs_bias) <- alpha_occs[1]

		})

		model_code <- glueNimbleCode(model_code, bias_response_curve_code)

	} else if (data_occs$n_covariates_occs_bias == 1) {
	# OCCURRENCE RESPONSE CURVES: ONE bias covariate

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

	# OCCURRENCE RESPONSE CURVES: MORE THAN ONE bias covariate
	} else if (data_occs$n_covariates_occs_bias > 1) {

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

	### BIOMASS: response curves
	if (data_biomass$n_covariates_biomass == 1) {

		### univariate
		response_curve_code <- nimbleCode({

			# BIOMASS: posterior predictive sampler for response curves: site-level mean
			# NB we assume just one predictor for biomass, so the response curve "x" is a matrix, not an array
			# Also assume presence
			for (i in 1:n_response_curve_values) {
					
				log_response_curves_biomass_mu[i] <-
					inprod(beta_biomass[1:n_terms_biomass], resp_curves_x_biomass[i, 1:n_terms_biomass])
				response_curves_biomass_mu[i] <- exp(log_response_curves_biomass_mu[i])
				
			}

		})

		model_code <- glueNimbleCode(model_code, response_curve_code)
	
	} else {
	
		### multivariate
		response_curve_code <- nimbleCode({

			# BIOMASS: posterior predictive sampler for response curves: site-level mean
			# NB we assume just one predictor for biomass, so the response curve "x" is a matrix, not an array
			# Also assume presence
			for (i in 1:n_covariates_biomass) {
				
				for (j in 1:n_response_curve_values) {
						
					log_response_curves_biomass_mu[j, i] <-
						inprod(beta_biomass[1:n_terms_biomass], resp_curves_x_biomass[j, 1:n_terms_biomass, i])
					response_curves_biomass_mu[j, i] <- exp(log_response_curves_biomass_mu[j, i])
					
				}
			}

		})

		model_code <- glueNimbleCode(model_code, response_curve_code)
	
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
		buildDerivs = TRUE # need for Hamiltonian Monte Carlo
		# buildDerivs = FALSE # need for Hamiltonian Monte Carlo
	)

	say('initializeInfo() and $calculate():', level = 2)
	model$initializeInfo()
	calc <- model$calculate()
	say('model$calculate(): ', calc)
	if (is.na(calc) || is.infinite(calc)) stop('Impossible likelihood.')

	say('configureMCMC():', level = 2)

	monitors_coeffs_not_indexed <- c('eta', 'sigma_biomass_within_sites', 'correlation[1, 2]')
	monitors_coeffs_single_index <- c('beta_occs', 'alpha_occs', 'beta_biomass', 'beta_psi', 'sigmas')
	monitors_coeffs_double_index <- c()

	monitors_derived_not_indexed <- c('log_lik')
	monitors_derived_single_index <- c('mu_biomass_site', 'z_site')
	monitors_derived_double_index <- c('U')

	monitors_dharma <- c(
		'lambda_mu_sq', 'y_n_ag_sim', 'y_biomass_sim'
	)

	monitors_resp_curves <- c(
		'response_curves_occs_mu', 'response_curves_biomass_mu', 'response_curves_psi'
	)
	if (data_occs$n_covariates_occs_bias > 0) monitors_resp_curves <- c(monitors_resp_curves, 'response_curves_occs_bias')

	monitors <- c(monitors_coeffs_not_indexed, monitors_coeffs_single_index, monitors_coeffs_double_index, monitors_derived_not_indexed, monitors_derived_single_index, monitors_derived_double_index, monitors_dharma, monitors_resp_curves)

	conf <- configureMCMC(
		model,
		monitors = monitors,
		print = TRUE,
		enableWAIC = TRUE
	)

	# add no U-turn sampler (Hamiltonian Monte Carlo)
	vars <- 'eta'
	conf$addSampler(target = vars, type = 'NUTS')
	say('NUTS sampler added to ', paste(vars, collapse = ' & '), '.')

	vars <- 'sigmas'
	conf$addSampler(target = vars, type = 'NUTS')
	say('NUTS sampler added to ', paste(vars, collapse = ' & '), '.')

	vars <- 'sigma_biomass_within_sites_log'
	conf$addSampler(target = vars, type = 'NUTS')
	say('NUTS sampler added to ', paste(vars, collapse = ' & '), '.')

	vars <- 'alpha_occs'
	conf$addSampler(target = vars, type = 'NUTS')
	say('NUTS sampler added to ', paste(vars, collapse = ' & '), '.')

	var <- 'beta_occs'
	conf$removeSamplers(var)
	conf$addSampler(target = var, type = 'AF_slice')
	say('AF_slice sampler added to ', paste(var, collapse = ' & '), '.')

	var <- 'beta_biomass'	
	conf$removeSamplers(var)
	conf$addSampler(target = var, type = 'AF_slice')
	say('AF_slice sampler added to ', paste(var, collapse = ' & '), '.')

	var <- 'beta_psi'	
	conf$removeSamplers(var)
	conf$addSampler(target = var, type = 'AF_slice')
	say('AF_slice sampler added to ', paste(var, collapse = ' & '), '.')

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

say('#################################################')
say('### post-modeling diagnostics and predictions ###')
say('#################################################')

	descrip <- 'occurrence ~ ZIP(MVN), biomass ~ ZIG(MVN)'
	workflow_postmodeling_generic(facet = 'occurrence + biomass', formulae = formulae, descrip = descrip, out_dir = out_dir)
	
	resp_distrib <- c(
		occs = 'ZIP',
		biomass = 'ZIG'
	)

	transform <- c(
		occs = 'exponential',
		biomass = 'exponential'
	)

	workflow_postmodeling_occurrence_biomass(
		formula_occs = formula_occs,
		formula_occs_bias = formula_occs_bias,
		formula_biomass = formula_biomass,
		formula_psi = formula_psi,
		resp_distrib = resp_distrib,
		transform = transform,
		log_precip = log_precip,
		out_dir = out_dir
	)

	if (crossvalidate) {
	
		workflow_postmodeling_occurrence_biomass_crossvalidation(
			chains = chains,
			constants = constants,
			inits = inits,
			formula_occs = formula_occs,
			formula_occs_bias = formula_occs_bias,
			formula_biomass = formula_biomass,
			formula_psi = formula_psi,
			out_dir = out_dir
		)
	
	}


say(date())
say('FINIS!', deco = '+', level = 1)

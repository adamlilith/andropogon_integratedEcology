## MODELING ANDROPOGON GERARDI DISTRIBUTION, PHENOTYPE, PHYSIOLOGY, GENOTYPE, and ASSOCIATED MICROBIAL COMMUNITIES
### Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2023-12
###
### This model estimates the occurrence (present/absent), relative abundance, and site-level mean values traits of Andropogon gerardi. Traits include aboveground individual biomass, height, canopy diameter, and others. All of these "facets" integrate through a shared multivariate normal distribution that accounts for correlations among facets. Mean values (or transforms of the means) are drawn from the MVN, then used as parameters for facet-specific distributions (e.g., zero-inflated Poisson for abundance, zero-inflated lognormal for biomass, etc.). The probability of occurrence is modeled as a function of the environment, and is reflected by hurdle models for each facet. Ergo, if the species is predicted to be absent in a location, all of its facet values will be forced to 0, as well.
###
### source('C:/Kaji/R/andropogon_integratedEcology/sdm_pdm/sdm_pdm_08b_model_fully_integrated~mvn.r')
### 
#############
### setup ###
#############

# if (TRUE) {
if (FALSE) {

	rm(list = ls())

	drive <- 'C:/Kaji/'

	setwd(paste0(drive, '/Research/Andropogon/Andropogon'))
	source(paste0(drive, '/R/andropogon_integratedEcology/sdm_pdm/sdm_pdm_00_shared_functions_and_variables.r'))

##########################
## user-defined values ###
##########################

	trial <- TRUE # TRUE for testing
	# trial <- FALSE # TRUE for testing

	# calib <- TRUE # use just counties with non-NA Poaceae for calibration region
	calib <- FALSE # use all of North America for calibration region

	# do cross-validation?
	crossvalidate <- TRUE
	# crossvalidate <- FALSE

	# log BIOs 12-14 and 16-19?
	# log_precip_occs <- FALSE

	### formula for how aspects of species responds to environment

	formula_occs <- ~ 1 + bio1 + bio12 + bio15 + I(bio1^2) + I(bio12^2) + I(bio15^2) # response of occurrence to climate and soil
	formula_occs_bias <- ~ 1 # sampling bias for AG records
	occs_filename <- 'bio1^2_bio12^2_bio15^2'
	log_precip_occs <- TRUE

	formula_psi <- ~ 1 + bio1 + bio12 + bio15 + I(bio1^2) + I(bio12^2) + I(bio15^2)
	psi_filename <- 'bio1^2_bio12^2_bio15^2'

	formula_biomass <- ~ 1 + bio12 # response of biomass to environment
	biomass_filename <- 'bio12'
	log_precip_biomass <- TRUE

	nonbiomasss_facets <- list(
		height = list(
			formula = ~ 1 + bio12,
			filename = 'bio12',
			resp_distrib = 'ZILN',
			transform = 'identity',
			log_precip = TRUE
		),
		blade_width = list(
			formula = ~ 1 + bio1 + bio12,
			filename = 'bio1_bio12',
			resp_distrib = 'ZILN',
			transform = 'exponential',
			log_precip = TRUE
		)
	)

	nonbiomasss_facet_names <- names(nonbiomasss_facets)
	nonbiomasss_facet_names_short <- paste(substr(nonbiomasss_facet_names, 1, 3), collapse = '_')

	for (f in seq_along(nonbiomasss_facets)) {

		facet <- names(nonbiomasss_facets)[f]
		form <- nonbiomasss_facets[[facet]]$formula
		terms <- attr(terms(form), 'term.labels')
		nonbiomasss_facets[[facet]]$n_terms <- length(terms) + 1

	}

	### output folder and bias formula
	out_dir <- paste0('./outputs_loretta/integrated_sdm_pdm/models_integrated/', ifelse(trial, 'TRIAL_', ''), '[occs~zip]_[biomass~ziln]_[', nonbiomasss_facet_names_short, ']_[psi]')

	if (!trial) {

		### MCMC settings
		# need these settings to achieve independent samples using bio1^2, bio12^2, bio15^
		niter <- 1600000
		# nchains <- 1
		nchains <- 4

	} else {
		
		### MCMC settings FOR TESTING
		niter <- 2000
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
	say('MODELING OCCURRENCE, BIOMASS, and NON-BIOMASS FACETS')
	say('Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | ', date(), post = 1)

	### data collation
	##################
	say('data collation', post = 1)

	say('This model estimates the occurrence (present/absent), relative abundance, and site-level mean values traits of Andropogon gerardi. Traits include aboveground individual biomass, height, canopy diameter, and others. All of these "facets" integrate through a shared multivariate normal distribution that accounts for correlations among facets. Mean values (or transforms of the means) are drawn from the MVN, then used as parameters for facet-specific distributions (e.g., zero-inflated Poisson for abundance, zero-inflated lognormal for biomass, etc.). The probability of occurrence is modeled as a function of the environment, and is reflected by hurdle models for each facet. Ergo, if the species is predicted to be absent in a location, all of its facet values will be forced to 0, as well.', breaks = 60, post = 1)

	say('MCMC settings:', level = 2)
	say('trial ........................ ', trial)
	say('niter ........................ ', niter)
	say('nburnin ...................... ', nburnin)
	say('thin ......................... ', thin)
	say('nchains ...................... ', nchains)
	say('formula_occs ................. ', paste(as.character(formula_occs), collapse = ' '))
	say('formula_occs_bias ............ ', paste(as.character(formula_occs_bias), collapse = ' '))
	say('formula_biomass .............. ', paste(as.character(formula_biomass), collapse = ' '))
	for (f in seq_along(nonbiomasss_facets)) {
		say('facet ', f, ' ...................... ', names(nonbiomasss_facets)[f])
		say('formula ', f, ' ....................... ', paste(as.character(nonbiomasss_facets[[f]]$formula), collapse = ' '))
		say('transform ', f, ' ..................... ', nonbiomasss_facets[[f]]$transform)
	}
	say('formula_psi .................. ', paste(as.character(formula_psi), collapse = ' '))

	say('out_dir')
	say(out_dir, post = 2)

	formulae <- list(
		formula_occs = formula_occs,
		formula_occs_bias = formula_occs_bias,
		formula_psi = formula_psi,
		formula_biomass = formula_biomass,
		nonbiomasss_facets = nonbiomasss_facets
	)
	saveRDS(formulae, paste0(out_dir, '/formulae.rds'))

	########################s
	### data preparation ###
	########################

	# data for occurrences at counties
	data_occs_counties <- prepare_occurrence_data(formula_occs = formula_occs, formula_occs_bias = formula_occs_bias, log_precip = log_precip_occs, n_response_curve_values = n_response_curve_values, psa_quant = psa_quant, calib = calib)

	# data for occurrences using site-level environment
	data_occs_sites <- prepare_biomass_data(formula_biomass = formula_occs, log_precip = log_precip_occs, n_response_curve_values = n_response_curve_values, calib = calib)

	# data for biomass at sites
	data_biomass_sites <- prepare_biomass_data(formula_biomass = formula_biomass, log_precip = log_precip_biomass, n_response_curve_values = n_response_curve_values, calib = calib)

	# data for biomass using county-level environment
	data_biomass_counties <- prepare_occurrence_data(formula_occs = formula_biomass, formula_occs_bias = ~ 1, log_precip = log_precip_biomass, n_response_curve_values = n_response_curve_values, psa_quant = psa_quant, calib = calib)

	# data for facets at sites and counties
	data_nonbiomass_sites <- data_nonbiomass_counties <- list()
	for (f in seq_along(nonbiomasss_facets)) {

		facet <- names(nonbiomasss_facets)[f]
		formula_facet <- nonbiomasss_facets[[facet]]$formula
		log_precip <- nonbiomasss_facets[[facet]]$log_precip

		data_nonbiomass_sites[[f]] <- prepare_nonbiomass_data(facet = facet, formula_facet = formula_facet, log_precip = log_precip, n_response_curve_values = n_response_curve_values, calib = calib)

		data_nonbiomass_counties[[f]] <- prepare_occurrence_data(formula_occs = formula_facet, formula_occs_bias = ~ 1, log_precip = log_precip, n_response_curve_values = n_response_curve_values, psa_quant = psa_quant, calib = calib)

	}

	names(data_nonbiomass_sites) <- names(data_nonbiomass_counties) <- names(nonbiomasss_facets)

	# data for psi at counties
	data_psi_counties <- prepare_occurrence_data(formula_occs = formula_psi, formula_occs_bias = ~1, log_precip = log_precip_occs, n_response_curve_values = n_response_curve_values, psa_quant = psa_quant, calib = calib)

	# data for psi at sites
	data_psi_sites <- prepare_biomass_data(formula_biomass = formula_psi, log_precip = log_precip_occs, n_response_curve_values = n_response_curve_values, calib = calib)

	#########################
	### inputs for nimble ###
	#########################

	say('Inputs:', level = 2)
	data <- list(
		y_n_ag = data_occs_counties$y_n_ag,				# number of AG observations in each county
		y_biomass = data_biomass_sites$y_biomass		# biomass of individual plants
	)

	for (f in seq_along(nonbiomasss_facets)) {

		data$DUMMY <- data_nonbiomass_sites[[f]]$y_facet
		names(data)[length(data)] <- paste0('y_facet_', f)

	}

	n_nonbiomasss_facets <- length(nonbiomasss_facets)
	n_facets <- n_nonbiomasss_facets + 2

	constants <- list(
		
		### integration
		n_facets = n_facets,
		n_nonbiomass_facets = n_nonbiomasss_facets,

		### occurrences
		n_counties_occs_calib = data_occs_counties$n_counties_occs_calib, # number of counties in calibration region
		n_terms_occs = data_occs_counties$n_terms_occs, # number of terms in formula for occurrence model (including intercept)
		n_terms_occs_bias = data_occs_counties$n_terms_occs_bias, # number of terms in sampling bias model
		w_occs_bias = data_occs_counties$w_occs_bias, # model matrix of sampling bias of AG observed occurrences
		
		# n_covariates_occs = data_occs_counties$n_covariates_occs, # number of covariates in formula for occurrence model
		# n_covariates_occs_bias = data_occs_counties$n_covariates_occs_bias, # number of covariates in formula for occurrence model
		# resp_curves_x_occs = data_occs_counties$resp_curves_x_occs, # response curve array for occurrences vs environment
		# resp_curves_w_occs = data_occs_counties$resp_curves_w_occs, # response curve array for occurrences vs environment

		x_by_county_occs = data_occs_counties$counties_x_occs_sq,
		x_by_site_occs = data_occs_sites$x_by_site_biomass, # MM with covariates for biomass (scaled)

		### sites (generic)
		n_pheno_sites = data_biomass_sites$n_pheno_sites, # number of phenotype sample sites

		### biomass
		x_by_county_biomass = data_biomass_counties$counties_x_occs_sq,

		x_by_site_biomass = data_biomass_sites$x_by_site_biomass, # MM with covariates for biomass (scaled)
		n_biomass = data_biomass_sites$n_biomass, # number of biomass observations
		site_index_biomass = data_biomass_sites$site_index_biomass, # index of sampled site for each row in biomass data
		n_terms_biomass = data_biomass_sites$n_terms_biomass, # number of terms in formula for biomass model (including intercept)

		# n_covariates_biomass = data_biomass_sites$n_covariates_biomass,
		# resp_curves_x_biomass = data_biomass_sites$resp_curves_x_biomass, # response curve array for biomass

		# facets (generic)
		n_non_biomass = data_nonbiomass_sites[[1]]$n_plants,
		site_index_facet = data_nonbiomass_sites[[1]]$site_index_facet, # index of sampled site for each row in facet data

		# # response curves (general)
		# n_response_curve_values = n_response_curve_values, # number of values in response curve array

		# zero inflation
		n_terms_psi = data_psi_counties$n_terms_occs, # number of terms in sampling bias model
		# n_covariates_psi = data_psi_counties$n_covariates_occs, # number of terms in sampling bias model
		x_by_county_psi = data_psi_counties$counties_x_occs_sq,
		x_by_site_psi = data_psi_sites$x_by_site_biomass
		# resp_curves_x_psi = data_psi_counties$resp_curves_x_occs

	)

	for (f in seq_along(nonbiomasss_facets)) {

		facet <- names(nonbiomasss_facets)[f]
		constants_facet <- list(

			x_by_site_facet = data_nonbiomass_sites[[facet]]$x_by_site_facet, # MM with covariates
			x_by_county_facet = data_nonbiomass_counties[[facet]]$counties_x_occs_sq, # MM with covariates
			n_terms_facet = nonbiomasss_facets[[facet]]$n_terms # includes intercept

		)
		names(constants_facet) <- paste0(names(constants_facet), '_', f)
		constants <- c(constants, constants_facet)

	}

	constants <- c(constants, constants_shared_occs, constants_shared_biomass, constants_shared_psi, constants_shared_facet)

	# occurrences initializations
	file <- paste0('./outputs_loretta/integrated_sdm_pdm/models_occurrence/[occs_zip[psi~', psi_filename, ']~normal~', occs_filename, '_[bias~1]]', ifelse(log_precip, '_log_precip', ''), '/chains.rds')
	chains_inits <- readRDS(file)
	
	alpha_occs_inits <- mc_extract(chains_inits, 'alpha_occs', j = TRUE)
	beta_occs_mu_inits <- mc_extract(chains_inits, 'beta_occs', j = TRUE)
	beta_psi_inits <- mc_extract(chains_inits, 'beta_psi', j = TRUE)
	
	lambda_mu_sq_inits <- mc_extract(chains_inits, 'lambda_mu_sq', j = TRUE)
	lambda_mu_sq_inits[lambda_mu_sq_inits > 10] <- 10

	lambda_sigma_inits <- mc_extract(chains_inits, 'lambda_sigma')
	log_lambda_sigma_inits <- log(lambda_sigma_inits)

	rm(chains_inits)

	N_inits_calib <- data_occs_counties$y_n_ag * 2 + 1
	z_county_inits <- as.numeric(N_inits_calib > 0)

	# biomass initializations
	file <- paste0('./outputs_loretta/integrated_sdm_pdm/models_biomass/[biomass_ziln~normal~', biomass_filename,  ']', ifelse(log_precip, '_log_precip', ''), '/chains.rds')
	chains_inits <- readRDS(file)
	beta_biomass_inits <- mc_extract(chains_inits, 'beta_biomass', j = TRUE)
	mu_biomass_site_inits <- data_biomass_sites$site_biomass_mean
	sigma_biomass_among_sites_inits <- mc_extract(chains_inits, 'sigma_biomass_among_sites')
	log_sigma_biomass_among_sites_inits <- log(sigma_biomass_among_sites_inits)
	sigma_biomass_within_sites_inits <- mc_extract(chains_inits, 'sigma_biomass_within_sites')
	log_sigma_biomass_within_sites_inits <- log(sigma_biomass_within_sites_inits)

	rm(chains_inits)

	inits <- list(

		### occurrences
		y_n_ag_sim = data_occs_counties$y_n_ag, # simulated values of observed number of AG (for DHARMa residuals)
		lambda_mu_sq = lambda_mu_sq_inits, # expected value of number of AG
		alpha_occs = alpha_occs_inits, # intercept, area, # of Poaceae
		beta_occs = beta_occs_mu_inits, # occurrence ~ environment coefficients (including intercept)

		N = N_inits_calib, # number of latent AG in calibration counties

		### biomass
		sigma_biomass_within_sites_log = log_sigma_biomass_within_sites_inits,
		y_biomass_sim = data_biomass_sites$y_biomass, # simulated values for biomass (for DHARMa residuals)
		beta_biomass = beta_biomass_inits,

		### integration
		eta = 1,
		U_star = diag(1, nrow = n_facets, ncol = n_facets),
		log_sigmas = c(log_lambda_sigma_inits, log_sigma_biomass_among_sites_inits),
		Phi_county = matrix(1, nrow = data_occs_counties$n_counties, ncol = n_facets),
		Phi_site = matrix(1, nrow = data_biomass_sites$n_pheno_sites, ncol = n_facets),

		### zero inflation
		beta_psi = beta_psi_inits,
		z_site = rep(1, data_biomass_sites$n_pheno_sites), # present/absent at each site
		z_county = z_county_inits # present/absent in county

	)

	# initialize each non-biomass facet
	for (f in seq_along(nonbiomasss_facets)) {

		facet <- names(nonbiomasss_facets)[f]

		filename <- nonbiomasss_facets[[facet]]$filename
		resp_distrib <- nonbiomasss_facets[[facet]]$resp_distrib
		this_log_precip <- nonbiomasss_facets[[facet]]$log_precip

		file <- paste0('./outputs_loretta/integrated_sdm_pdm/models_', facet, '/[', facet, '_', tolower(resp_distrib), '~normal~', filename,  ']', ifelse(log_precip, '_log_precip', ''), '/chains.rds')

		chains_inits <- readRDS(file)

		sigma_facet_within_sites_init <- mc_extract(chains_inits, 'sigma_facet_within_sites')
		sigma_facet_among_sites_init <- mc_extract(chains_inits, 'sigma_facet_among_sites')
		beta_facet_inits <- mc_extract(chains_inits, 'beta_facet', j = TRUE)

		y_facet_sim <- data_nonbiomass_sites[[facet]]$y_facet # simulated values for facet (for DHARMa residuals)

		inits_facet <- list(
			log_sigma_facet_within_sites = log(sigma_facet_within_sites_init),
			beta_facet = beta_facet_inits,
			y_facet_sim = y_facet_sim
		)

		names(inits_facet) <- paste0(names(inits_facet), '_', f)

		inits <- c(inits, inits_facet)

		inits$log_sigmas <- c(inits$log_sigmas, log(sigma_facet_among_sites_init))

	}

	say('Data:')
	print(str(data))

	say('Constants:', pre = 1)
	print(str(constants))

	say('Initializations:', pre = 1)
	print(str(inits))

	### define model
	say('nimbleCode():', level = 2)

	} # if prep

	model_code <- nimbleCode({

		# INTEGRATION: LJK prior for standard deviations and correlations between latent occurrence and biomass processes
		eta ~ dgamma(2, 1)
		U_star[1:n_facets, 1:n_facets] ~ dlkj_corr_cholesky(eta = eta, p = n_facets)
		U[1:n_facets, 1:n_facets] <- uppertri_mult_diag(
			U_star[1:n_facets, 1:n_facets],
			sigmas[1:n_facets]
		)
	  
		# INTEGRATION: standard deviations of latent occurrence and biomass
		log(sigmas[1]) ~ dnorm(0, sd = lambda_sigma_prior_sd) # half-Cauchy
		log(sigmas[2]) ~ dnorm(0, sd = lambda_sigma_prior_sd) # half-Cauchy

		for (i in 3:n_facets) {
			log(sigmas[i]) ~ dnorm(0, sd = sigma_facet_among_sites_log_prior_sd) # half-Cauchy
		}

		# correlation matrix
		correlation[1:n_facets, 1:n_facets] <- t(U_star[1:n_facets, 1:n_facets]) %*% U_star[1:n_facets, 1:n_facets]

		# INTEGRATION: county-level
		# NB We have measurements of biomass at the site level, but not county. We thus assume that the estimates of biomass using county-level covariates are indicative of virtual sites that have the same environments as counties. Biomass is still estimated at the site level.
		for (i in 1:n_counties_occs_calib) {

			phis_county[i, 1:n_facets] <- c(
				phi_occs_mu_sq[i],
				phi_biomass_mu_sq[i],
				phi_facet_mu_sq[i, 1],
				phi_facet_mu_sq[i, 2]
			)
			
			Phi_county[i, 1:n_facets] ~ dmnorm(phis_county[i, 1:n_facets], cholesky = U[1:n_facets, 1:n_facets], prec_param = 0)

			phi_occs_mu_sq[i] <- inprod(beta_occs[1:n_terms_occs], x_by_county_occs[i, 1:n_terms_occs])
			phi_biomass_mu_sq[i] <- inprod(beta_biomass[1:n_terms_biomass], x_by_county_biomass[i, 1:n_terms_biomass])

			# for (j in 1:n_nonbiomass_facets) {

				phi_facet_mu_sq[i, 1] <- inprod(beta_facet_1[1:n_terms_facet_1], x_by_county_facet_1[i, 1:n_terms_facet_1])
				phi_facet_mu_sq[i, 2] <- inprod(beta_facet_2[1:n_terms_facet_2], x_by_county_facet_2[i, 1:n_terms_facet_2])

			# }

		}

		# INTEGRATION: site-level calculations
		# This is how the phenotypic variables communicate with abundance. Abundance is naturally estimated at the county level, so we assume that the site-level environmental covariates act like faux counties for abundance, but are indicative of the sites for phenotypic variables.
		for (i in 1:n_pheno_sites) {

			phis_site[i, 1:n_facets] <- c(
				phi_occs_mu_site[i],
				phi_biomass_mu_site[i],
				phi_facet_mu_site[i, 1],
				phi_facet_mu_site[i, 2]
			)
			Phi_site[i, 1:n_facets] ~ dmnorm(phis_site[i, 1:n_facets], cholesky = U[1:n_facets, 1:n_facets], prec_param = 0)

			phi_occs_mu_site[i] <- inprod(beta_occs[1:n_terms_occs], x_by_site_occs[i, 1:n_terms_occs])
			phi_biomass_mu_site[i] <- inprod(beta_biomass[1:n_terms_biomass], x_by_site_biomass[i, 1:n_terms_biomass])

			# for (j in 1:n_nonbiomass_facets) {

				phi_facet_mu_site[i, 1] <- inprod(beta_facet_1[1:n_terms_facet_1], x_by_site_facet_1[i, 1:n_terms_facet_1])
				phi_facet_mu_site[i, 2] <- inprod(beta_facet_2[1:n_terms_facet_2], x_by_site_facet_2[i, 1:n_terms_facet_2])

			# }


		}

		# OCCURRENCE: weakly regularized or regularized priors for relationship to environment
		beta_occs[1] ~ dnorm(0, sd = beta_occs_prior_dnorm_sd_1)
		for (i in 2:n_terms_occs) {
			beta_occs[i] ~ ddexp(0, rate = beta_occs_prior_ddexp_rate)
		}

		# OCCURRENCE: priors for sampling bias
		alpha_occs[1] ~ dnorm(0, sd = alpha_occs_prior_dnorm_sd_1)

		# ZERO-INFLATION: weakly regularized or regularized priors for relationship to environment
		beta_psi[1] ~ dnorm(0, sd = beta_psi_prior_dnorm_sd_1)
		for (i in 2:n_terms_psi) {
			beta_psi[i] ~ ddexp(0, rate = beta_psi_prior_exp_rate)
		}

		# OCCURRENCE: likelihood
		for (i in 1:n_counties_occs_calib) {
			
			### actual abundance (latent--unobserved)
			N[i] ~ dTruncPseudoZIP(lambda_mu_sq[i], z = z_county[i])

			### observed number of AG and sampling bias
			y_n_ag[i] ~ dbinom(prob = p[i], size = N[i])

			# simulate observations for DHARMa residuals
			y_n_ag_sim[i] ~ dbinom(prob = p[i], size = N[i])

			# probability of observing a single AG
			logit(p[i]) <- alpha_occs[1]

			# relationship between expected (latent) abundance and environment assuming MV NORMAL distribution
			log(lambda_mu_sq[i]) <- exp(Phi_county[i, 1])

			# probability of zero inflation
			logit(psi_county[i]) <- inprod(beta_psi[1:n_terms_psi], x_by_county_psi[i, 1:n_terms_psi])
			z_county[i] ~ dbern(psi_county[i])

			# # likelihood
			# log_lik_occs_y[i] <- dbinom(y_n_ag[i], prob = p[i], size = N[i], log = 1)

		}

		# log_lik_occs <- sum(log_lik_occs_y[1:n_counties_occs_calib])

		# BIOMASS: priors for relationship of site-level mean biomass to environment
		for (i in 1:n_terms_biomass) {
			beta_biomass[i] ~ dnorm(0, sd = beta_biomass_prior_dnorm_sd) # broad prior
		}

		# NON-BIOMASS FACETS
		for (i in 1:n_terms_facet_1) {
			beta_facet_1[i] ~ dnorm(0, sd = beta_facet_prior_dnorm_sd) # broad prior
		}

		for (i in 1:n_terms_facet_2) {
			beta_facet_2[i] ~ dnorm(0, sd = beta_facet_prior_dnorm_sd) # broad prior
		}

		# prior for sd of individual plant biomass on lognormal (~ half-Cauchy), ==> vague
		sigma_biomass_within_sites_log ~ dnorm(0, sd = sigma_biomass_within_sites_log_prior_sd)
		sigma_biomass_within_sites <- exp(sigma_biomass_within_sites_log)

		# BIOMASS and NON-BIOMASS FACETS: parameters of biomass distribution are latent and functions of environment
		# individual plant traits are samples from the site-level distribution defined by the site-level distribution
		for (i in 1:n_pheno_sites) {

			# PROBABILITY OF ZERO INFLATION
			logit(psi_site[i]) <- inprod(beta_psi[1:n_terms_psi], x_by_site_psi[i, 1:n_terms_psi])

			# presence/absence at the site level... NB somewhat of a cheat bc we could have absence at county level but presence in site in county
			z_site[i] ~ dbern(psi_site[i])

			### BIOMASS
			mu_biomass_site[i] <- exp(Phi_site[i, 2])

			### FACET 1
			mu_facet_site_1[i] <- exp(Phi_site[i, 3])

			### FACET 2
			mu_facet_site_2[i] <- exp(Phi_site[i, 4])

		}

		# BIOMASS and NON-BIOMASS: likelihood of individual plants
		for (i in 1:n_biomass) {

			### BIOMASS

			# likelihood
			y_biomass[i] ~ dZILN(meanlog = Phi_site[site_index_biomass[i], 2], sdlog = sigma_biomass_within_sites, z = z_site[site_index_biomass[i]])

			# simulated values for unconditional DHARMa residuals
			y_biomass_sim[i] ~ dZILN(meanlog = Phi_site[site_index_biomass[i], 2], sdlog = sigma_biomass_within_sites, z = z_site[site_index_biomass[i]])

		}

		# FACETs prior for sd of individual plant facet value on lognormal (~ half-Cauchy), ==> vague
		log(sigma_facet_within_sites_1) ~ dnorm(0, sd = sigma_facet_within_sites_log_prior_sd)
		log(sigma_facet_within_sites_2) ~ dnorm(0, sd = sigma_facet_within_sites_log_prior_sd)

		# FACETS: likelihood of individual plant's non-biomass traits
		for (i in 1:n_non_biomass) {

			### FACET 1
			y_facet_1[i] ~ dZILN(meanlog = Phi_site[site_index_facet[i], 3], sdlog = sigma_facet_within_sites_1, z = z_site[site_index_facet[i]])

			# simulated values for unconditional DHARMa residuals
			y_facet_sim_1[i] ~ dZILN(meanlog = Phi_site[site_index_facet[i], 3], sdlog = sigma_facet_within_sites_1, z = z_site[site_index_facet[i]])

			### FACET 2
			y_facet_2[i] ~ dZILN(meanlog = Phi_site[site_index_facet[i], 4], sdlog = sigma_facet_within_sites_2, z = z_site[site_index_facet[i]])

			# simulated values for unconditional DHARMa residuals
			y_facet_sim_2[i] ~ dZILN(meanlog = Phi_site[site_index_facet[i], 4], sdlog = sigma_facet_within_sites_2, z = z_site[site_index_facet[i]])

		}

	})

	print(model_code)

	say('nimbleModel():', level = 2)
	model <- nimbleModel(
		code = model_code,
		constants = constants,
		data = data,
		inits = inits,
		check = TRUE,
		calculate = FALSE,
		# buildDerivs = TRUE # need for Hamiltonian Monte Carlo
		buildDerivs = FALSE # should be TRUE with dmnorm()
	)

	say('initializeInfo() and $calculate():', level = 2)
	model$initializeInfo()
	calc <- model$calculate()
	say('model$calculate(): ', calc)
	if (is.na(calc) || is.infinite(calc)) stop('Impossible likelihood.')

	say('configureMCMC():', level = 2)

	monitors_coeffs_not_indexed <- c('eta', 'sigma_biomass_within_sites', 'sigma_facet_within_sites_1', 'sigma_facet_within_sites_2')
	monitors_coeffs_single_index <- c('beta_occs', 'alpha_occs', 'beta_biomass', 'beta_psi', 'sigmas', 'beta_facet_1', 'beta_facet_2')
	monitors_coeffs_double_index <- c('correlation')

	monitors_derived_not_indexed <- c()
	monitors_derived_single_index <- c()
	monitors_derived_double_index <- c('U')

	monitors_dharma <- c(
		'lambda_mu_sq', 'y_n_ag_sim', 'y_biomass_sim', 'y_facet_sim_1', 'y_facet_sim_2'
	)

	monitors <- c(monitors_coeffs_not_indexed, monitors_coeffs_single_index, monitors_coeffs_double_index, monitors_derived_not_indexed, monitors_derived_single_index, monitors_derived_double_index, monitors_dharma)

	conf <- configureMCMC(
		model,
		monitors = monitors,
		print = TRUE,
		enableWAIC = TRUE
	)

	# # add no U-turn sampler (Hamiltonian Monte Carlo)
	# var <- 'eta'
	# conf$addSampler(target = var, type = 'NUTS')
	# say('NUTS sampler added to ', paste(var, collapse = ' & '), '.')

	# var <- 'sigmas'
	# conf$addSampler(target = var, type = 'NUTS')
	# say('NUTS sampler added to ', paste(var, collapse = ' & '), '.')

	# var <- 'sigma_biomass_within_sites_log'
	# conf$addSampler(target = var, type = 'NUTS')
	# say('NUTS sampler added to ', paste(var, collapse = ' & '), '.')

	# var <- 'sigma_facet_within_sites_log_1'
	# conf$addSampler(target = var, type = 'NUTS')
	# say('NUTS sampler added to ', paste(var, collapse = ' & '), '.')

	# var <- 'sigma_facet_within_sites_log_2'
	# conf$addSampler(target = var, type = 'NUTS')
	# say('NUTS sampler added to ', paste(var, collapse = ' & '), '.')

	# var <- 'alpha_occs'
	# conf$addSampler(target = var, type = 'NUTS')
	# say('NUTS sampler added to ', paste(var, collapse = ' & '), '.')

	var <- 'beta_occs'
	conf$removeSamplers(var)
	conf$addSampler(target = var, type = 'AF_slice')
	say('AF_slice sampler added to ', paste(var, collapse = ' & '), '.')

	var <- 'beta_psi'
	conf$removeSamplers(var)
	conf$addSampler(target = var, type = 'AF_slice')
	say('AF_slice sampler added to ', paste(var, collapse = ' & '), '.')

	var <- 'beta_biomass'
	conf$removeSamplers(var)
	conf$addSampler(target = var, type = 'AF_slice')
	say('AF_slice sampler added to ', paste(var, collapse = ' & '), '.')

	for (f in seq_along(nonbiomasss_facets)) {

		facet <- names(nonbiomasss_facets)[f]
		var <- paste0('beta_facet_', f)
		# if (nonbiomasss_facets[[f]]$n_terms == 2) {

			# conf$addSampler(target = var, type = 'NUTS')
			# say('NUTS sampler added to ', paste(var, collapse = ' & '), '.')

		# } else if (nonbiomasss_facets[[f]]$n_terms > 2) {
		if (nonbiomasss_facets[[f]]$n_terms > 2) {

			conf$removeSamplers(var)
			conf$addSampler(target = var, type = 'AF_slice')
			say('AF_slice sampler added to ', paste(var, collapse = ' & '), '.')

		}

	}

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

	descrip <- 'occurrence ~ ZIP(MVN), biomass ~ ZILN(MVN), facets~fx(MNV)'
	facets_nice <- paste0('occurrence + biomass + ', paste(names(nonbiomasss_facets), collapse = ' + '))
	workflow_postmodeling_generic(facet = facets_nice, formulae = formulae, descrip = descrip, out_dir = out_dir)

print(NON)

	workflow_postmodeling_fully_integrated(
		formula_occs = formula_occs,
		formula_occs_bias = formula_occs_bias,
		log_precip_occs = log_precip_occs,
		formula_psi = formula_psi,
		formula_biomass = formula_biomass,
		log_precip_biomass = log_precip_biomass,
		nonbiomasss_facets = nonbiomasss_facets,
		out_dir = out_dir
	)
	
	if (crossvalidate) {
	
		workflow_postmodeling_fully_integrated_crossvalidation(
			chains = chains,
			constants = constants,
			inits = inits,
			formula_occs = formula_occs,
			formula_occs_bias = formula_occs_bias,
			formula_psi = formula_psi,
			formula_biomass = formula_biomass,
			nonbiomasss_facets = nonbiomasss_facets,
			out_dir = out_dir
		)
	
	}

say(date())
say('FINIS!', deco = '+', level = 1)

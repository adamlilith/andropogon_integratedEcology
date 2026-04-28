## MODELING ANDROPOGON GERARDI DISTRIBUTION, PHENOTYPE, PHYSIOLOGY, GENOTYPE, and ASSOCIATED MICROBIAL COMMUNITIES
### Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2023-12
###
### This model estimates the occurrence (present/absent), relative abundance, and site-level mean values traits of Andropogon gerardi. Traits include aboveground individual biomass, height, canopy diameter, and others. All of these "facets" integrate through a shared multivariate normal distribution that accounts for correlations among facets. Mean values (or transforms of the means) are drawn from the MVN, then used as parameters for facet-specific distributions (e.g., zero-inflated Poisson for abundance, zero-inflated lognormal for biomass, etc.). The probability of occurrence is modeled as a function of the environment, and is reflected by hurdle models for each facet. Ergo, if the species is predicted to be absent in a location, all of its facet values will be forced to 0, as well. This script models TWO non-biomass facets.
###
### source('C:/Kaji/R/andropogon_integratedEcology/sdm_pdm/sdm_pdm_08a_model_fully_integrated~mvn_two_traits.r')
### 
#############
### setup ###
#############

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

	### formula for how aspects of species responds to environment

	formula_occs <- ~ 1 + bio1 + bio12_log10p1 + bio15 + I(bio1^2) + I(bio12_log10p1^2) + I(bio15^2) # response of occurrence to climate and soil
	formula_bias <- ~ 1 # sampling bias for AG records
	occs_filename <- 'bio1^2_log(bio12)^2_bio15^2'
	bias_filename <- '1'

	formula_psi <- ~ 1 + bio1 + bio12_log10p1 + bio15 + I(bio1^2) + I(bio12_log10p1^2) + I(bio15^2)
	psi_filename <- 'bio1^2_log(bio12)^2_bio15^2'

	formula_biomass <- ~ 1 + bio12_log10p1 # response of biomass to environment
	biomass_filename <- 'log(bio12)'
	resp_distrib_biomass <- 'hurdleLN'
	transform_biomass <- if (resp_distrib_biomass == 'hGamma') { 'exponential' } else if (resp_distrib_biomass == 'hurdleLN') { 'identity' }

	nonbiomass_facets <- list(
		height = list(
			formula = ~ 1 + bio12_log10p1,
			filename = 'log(bio12)',
			resp_distrib = 'hurdleLN',
			transform = 'identity'
		),
		blade_width = list(
			formula = ~ 1 + bio1 + bio12,
			filename = 'bio1_bio12',
			resp_distrib = 'hurdleLN',
			transform = 'identity'
		)
	)

	nonbiomass_facets <- nonbiomass_facets[sort(names(nonbiomass_facets))]

	nonbiomasss_facet_names <- names(nonbiomass_facets)
	nonbiomasss_facet_names_short <- paste(substr(nonbiomasss_facet_names, 1, 3), collapse = '_')

	for (f in seq_along(nonbiomass_facets)) {

		facet <- names(nonbiomass_facets)[f]
		form <- nonbiomass_facets[[facet]]$formula
		terms <- attr(terms(form), 'term.labels')
		nonbiomass_facets[[facet]]$n_terms <- length(terms) + 1

	}

	### output folder and bias formula
	out_dir <- paste0('./outputs_loretta/integrated_sdm_pdm/models_integrated/', ifelse(trial, 'TRIAL_', ''), '[occs~hurdlepoisson]_[biomass~hurdleln]_[', nonbiomasss_facet_names_short, ']')

	if (!trial) {

		### MCMC settings
		# need these settings to achieve independent samples using bio1^2, bio12^2, bio15^
		niter <- 1600000
		# nchains <- 1
		nchains <- 4

	} else {
		
		### MCMC settings FOR TESTING
		niter <- 20000
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

	say('This model estimates the occurrence (present/absent), relative abundance, and site-level mean values traits of Andropogon gerardi. Traits include aboveground individual biomass, height, canopy diameter, and others. All of these "facets" integrate through a shared multivariate normal distribution that accounts for correlations among facets. Mean values (or transforms of the means) are drawn from the MVN, then used as parameters for facet-specific distributions (e.g., zero-inflated Poisson for abundance, zero-inflated lognormal for biomass, etc.). The probability of occurrence is modeled as a function of the environment, and is reflected by hurdle models for each facet. Ergo, if the species is predicted to be absent in a location, all of its facet values will be forced to 0, as well. This script models TWO non-biomass facets.', breaks = 60, post = 1)

	say('MCMC settings:', level = 2)
	say('trial ........................ ', trial)
	say('niter ........................ ', niter)
	say('nburnin ...................... ', nburnin)
	say('thin ......................... ', thin)
	say('nchains ...................... ', nchains)
	say('formula_occs ................. ', paste(as.character(formula_occs), collapse = ' '))
	say('formula_psi .................. ', paste(as.character(formula_psi), collapse = ' '))
	say('formula_bias ................. ', paste(as.character(formula_bias), collapse = ' '))
	say('formula_biomass .............. ', paste(as.character(formula_biomass), collapse = ' '))
	for (f in seq_along(nonbiomass_facets)) {
		say('facet ', f, ' ...................... ', names(nonbiomass_facets)[f])
		say('formula ', f, ' ........................ ', paste(as.character(nonbiomass_facets[[f]]$formula), collapse = ' '))
		say('resp_distrib ', f, ' .................. ', nonbiomass_facets[[f]]$resp_distrib)
		say('transform ', f, ' ..................... ', nonbiomass_facets[[f]]$transform)
	}

	say('out_dir ......................... ')
	say(out_dir, post = 2)

	formulae <- list(
		formula_occs = formula_occs,
		formula_bias = formula_bias,
		formula_psi = formula_psi,
		formula_biomass = formula_biomass,
		nonbiomass_facets = nonbiomass_facets
	)
	saveRDS(formulae, paste0(out_dir, '/formulae.rds'))

	facet_table <- data.table(
		index = 1:(2 + length(nonbiomass_facets)),
		facet = c('occurrence', 'biomass', names(nonbiomass_facets))
	)

	fwrite(facet_table, paste0(out_dir, '/!facet_codes.csv'))

	########################s
	### data preparation ###
	########################

	say('preparing data for occurrences', level = 2)

	# data for OCCURRENCES at counties
	data_occs_counties <- prepare_occurrence_data(formula_occs = formula_occs, formula_bias = formula_bias, n_response_curve_values = n_response_curve_values, psa_quant = psa_quant, calib = calib)

	# data for OCCURRENCES using site-level environment
	data_occs_sites <- prepare_biomass_data(formula_biomass = formula_occs, n_response_curve_values = n_response_curve_values, calib = calib)

	say('preparing data for psi', level = 2)

	# data for ZERO-INFLATION at counties
	data_psi_counties <- prepare_occurrence_data(formula_occs = formula_psi, formula_bias = formula_bias, n_response_curve_values = n_response_curve_values, psa_quant = psa_quant, calib = calib)

	# data for psi at sites
	data_psi_sites <- prepare_biomass_data(formula_biomass = formula_psi, n_response_curve_values = n_response_curve_values, calib = calib)

	say('preparing data for biomass', level = 2)
	# data for BIOMASS at sites
	data_biomass_sites <- prepare_biomass_data(formula_biomass = formula_biomass, n_response_curve_values = n_response_curve_values, calib = calib)

	# data for BIOMASS using county-level environment
	data_biomass_counties <- prepare_occurrence_data(formula_occs = formula_biomass, formula_bias = ~ 1, n_response_curve_values = n_response_curve_values, psa_quant = psa_quant, calib = calib)

	# data for NON-BIOMASS FACETS at sites and counties
	data_nonbiomass_sites <- data_nonbiomass_counties <- list()
	for (f in seq_along(nonbiomass_facets)) {

		facet <- names(nonbiomass_facets)[f]
		say('preparing data for ', facet, level = 2)
		
		formula_facet <- nonbiomass_facets[[facet]]$formula

		data_nonbiomass_sites[[f]] <- prepare_nonbiomass_data(facet = facet, formula_facet = formula_facet, n_response_curve_values = n_response_curve_values, calib = calib)

		data_nonbiomass_counties[[f]] <- prepare_occurrence_data(formula_occs = formula_facet, formula_bias = ~ 1, n_response_curve_values = n_response_curve_values, psa_quant = psa_quant, calib = calib)

	}

	names(data_nonbiomass_sites) <- names(data_nonbiomass_counties) <- names(nonbiomass_facets)

	#########################
	### inputs for nimble ###
	#########################

	say('Inputs:', level = 2)
	data <- list(
		y_n_ag = data_occs_counties$y_n_ag,				# number of AG observations in each county
		y_biomass = data_biomass_sites$y_biomass		# biomass of individual plants
	)

	for (f in seq_along(nonbiomass_facets)) {

		data$DUMMY <- data_nonbiomass_sites[[f]]$y_facet
		names(data)[length(data)] <- paste0('y_facet_', f)

	}

	n_nonbiomasss_facets <- length(nonbiomass_facets)
	n_facets <- n_nonbiomasss_facets + 2

	constants <- list(
		
		### integration
		n_facets = n_facets,
		# n_nonbiomass_facets = n_nonbiomasss_facets,

		### occurrences
		n_counties_occs_calib = data_occs_counties$n_counties_occs_calib, # number of counties in calibration region
		n_terms_occs = data_occs_counties$n_terms, # number of terms in formula for occurrence model (including intercept)
		n_terms_bias = data_occs_counties$n_terms_bias, # number of terms in sampling bias model
		w_bias = data_occs_counties$w_bias, # model matrix of sampling bias of AG observed occurrences
		
		x_by_county_occs = data_occs_counties$counties_x_sq,
		x_by_site_occs = data_occs_sites$x_by_site, # MM with covariates for biomass (scaled)

		### sites (generic)
		n_pheno_sites = data_biomass_sites$n_pheno_sites, # number of phenotype sample sites

		### biomass
		x_by_county_biomass = data_biomass_counties$counties_x_sq,

		x_by_site_biomass = data_biomass_sites$x_by_site, # MM with covariates for biomass (scaled)
		n_biomass = data_biomass_sites$n_biomass, # number of biomass observations
		site_index_biomass = data_biomass_sites$site_index_biomass, # index of sampled site for each row in biomass data
		n_terms_biomass = data_biomass_sites$n_terms, # number of terms in formula for biomass model (including intercept)

		# facets (generic)
		n_non_biomass = data_nonbiomass_sites[[1]]$n_plants,
		site_index_facet = data_nonbiomass_sites[[1]]$site_index_facet, # index of sampled site for each row in facet data

		# zero inflation
		n_terms_psi = data_psi_counties$n_terms, # number of terms in sampling bias model
		x_by_county_psi = data_psi_counties$counties_x_sq,
		x_by_site_psi = data_psi_sites$x_by_site

	)

	for (f in seq_along(nonbiomass_facets)) {

		facet <- names(nonbiomass_facets)[f]
		constants_facet <- list(

			x_by_site_facet = data_nonbiomass_sites[[facet]]$x_by_site, # MM with covariates
			x_by_county_facet = data_nonbiomass_counties[[facet]]$counties_x_sq, # MM with covariates
			n_terms_facet = nonbiomass_facets[[facet]]$n_terms # includes intercept

		)
		names(constants_facet) <- paste0(names(constants_facet), '_', f)
		constants <- c(constants, constants_facet)

	}

	constants <- c(constants, constants_shared_occs, constants_shared_biomass, constants_shared_psi, constants_shared_facet)

	# occurrences initializations
	file <- paste0('./outputs_loretta/integrated_sdm_pdm/models_occurrence/[occs~hurdle_', occs_filename, ']_[bias~', bias_filename, ']_niter_160000/chains.rds')
	chains_inits <- readRDS(file)
	
	alpha_occs_inits <- mc_extract(chains_inits, 'alpha_occs', j = TRUE)
	beta_occs_mu_inits <- mc_extract(chains_inits, 'beta_occs', j = TRUE)
	beta_psi_inits <- mc_extract(chains_inits, 'beta_psi', j = TRUE)
	
	rm(chains_inits)

	N_inits_calib <- data_occs_counties$y_n_ag * 2 + 1

	# biomass initializations
	file <- paste0('./outputs_loretta/integrated_sdm_pdm/models_biomass/[biomass~hurdleln_', biomass_filename,  ']/chains.rds')
	chains_inits <- readRDS(file)
	beta_biomass_inits <- mc_extract(chains_inits, 'beta_biomass', j = TRUE)
	biomass_site_inits <- data_biomass_sites$site_mean
	sigma_biomass_within_sites_inits <- mc_extract(chains_inits, 'sigma_biomass_within_sites')
	log_sigma_biomass_within_sites_inits <- log(sigma_biomass_within_sites_inits)

	rm(chains_inits)

	inits <- list(

		### occurrences
		y_n_ag_sim = data_occs_counties$y_n_ag, # simulated values of observed number of AG (for DHARMa residuals)
		alpha_occs = alpha_occs_inits, # intercept, area, # of Poaceae
		beta_occs = beta_occs_mu_inits, # occurrence ~ environment coefficients (including intercept)

		# N = N_inits_calib, # number of latent AG in calibration counties

		### biomass
		sigma_biomass_within_sites_log = log_sigma_biomass_within_sites_inits,
		y_biomass_sim = data_biomass_sites$y_biomass, # simulated values for biomass (for DHARMa residuals)
		beta_biomass = beta_biomass_inits,

		### integration
		eta = 1,
		U_star = diag(1, nrow = n_facets, ncol = n_facets),
		log_sigmas = c(1, 1),
		Phi_county = matrix(1, nrow = data_occs_counties$n_counties, ncol = n_facets),
		Phi_site = matrix(1, nrow = data_biomass_sites$n_pheno_sites, ncol = n_facets),
		z_county = matrix(0, nrow = data_occs_counties$n_counties, ncol = n_facets),
		z_site = matrix(0, nrow = data_biomass_sites$n_pheno_sites, ncol = n_facets),

		### zero inflation
		beta_psi = beta_psi_inits

	)

	# initialize each non-biomass facet
	for (f in seq_along(nonbiomass_facets)) {

		facet <- names(nonbiomass_facets)[f]

		filename <- nonbiomass_facets[[facet]]$filename
		resp_distrib <- nonbiomass_facets[[facet]]$resp_distrib

		file <- paste0('./outputs_loretta/integrated_sdm_pdm/models_', facet, '/[', facet, '~', tolower(resp_distrib), '(', filename,  ')]/chains.rds')
		chains_inits <- readRDS(file)

		sigma_facet_within_sites_init <- mc_extract(chains_inits, 'sigma_facet_within_sites')
		beta_facet_inits <- mc_extract(chains_inits, 'beta_facet', j = TRUE)

		y_facet_sim <- data_nonbiomass_sites[[facet]]$y_facet # simulated values for facet (for DHARMa residuals)

		inits_facet <- list(
			log_sigma_facet_within_sites = log(sigma_facet_within_sites_init),
			beta_facet = beta_facet_inits,
			y_facet_sim = y_facet_sim
		)

		names(inits_facet) <- paste0(names(inits_facet), '_', f)

		inits <- c(inits, inits_facet)

		inits$log_sigmas <- c(inits$log_sigmas, 1)

	}

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
		# U[1:n_facets, 1:n_facets] <- uppertri_mult_diag(U_star[1:n_facets, 1:n_facets], sigmas[1:n_facets])
	  
		# INTEGRATION: standard deviations of latent occurrence and biomass
		log(sigmas[1]) ~ dnorm(0, sd = lambda_sigma_prior_sd) # lognormal prior
		log(sigmas[2]) ~ dnorm(0, sd = sigma_biomass_among_sites_log_prior_sd) # lognormal prior

		for (i in 3:n_facets) {
			log(sigmas[i]) ~ dnorm(0, sd = sigma_facet_among_sites_log_prior_sd) # lognormal prior
		}

		# correlation matrix
		correlation[1:n_facets, 1:n_facets] <- t(U_star[1:n_facets, 1:n_facets]) %*% U_star[1:n_facets, 1:n_facets]

		# INTEGRATION: county-level
		# NB We have measurements of biomass at the site level, but not county. We thus assume that the estimates of biomass using county-level covariates are indicative of virtual sites that have the same environments as counties. Biomass is still estimated at the site level.
		
		# latent variables for correlation
		for (i in 1:n_counties_occs_calib) {
			for (j in 1:n_facets) {
				z_county[i, j] ~ dnorm(0, tau = 1)
			}
		}

		# relationship between abundance, biomass, and facets at COUNTY level with environment
		phi_occs_county[1:n_counties_occs_calib] <- (x_by_county_occs[1:n_counties_occs_calib, 1:n_terms_occs] %*% beta_occs[1:n_terms_occs])[ , 1]
		phi_biomass_county[1:n_counties_occs_calib] <- (x_by_county_biomass[1:n_counties_occs_calib, 1:n_terms_biomass] %*% beta_biomass[1:n_terms_biomass])[, 1]

		for (i in 1:n_counties_occs_calib) {

			phi_facet_county[i, 1] <- inprod(beta_facet_1[1:n_terms_facet_1], x_by_county_facet_1[i, 1:n_terms_facet_1])
			phi_facet_county[i, 2] <- inprod(beta_facet_2[1:n_terms_facet_2], x_by_county_facet_2[i, 1:n_terms_facet_2])

			phis_county[i, 1] <- phi_occs_county[i]
			phis_county[i, 2] <- phi_biomass_county[i]
			phis_county[i, 3] <- phi_facet_county[i, 1]
			phis_county[i, 4] <- phi_facet_county[i, 2]
						
		}

		# use offset method to avoid sampling from non-centered distribution
		# stands in for: Phi_county[i, 1:n_facets] ~ dmnorm(phis_county[i, 1:n_facets], cholesky = U[1:n_facets, 1:n_facets], prec_param = 0)
		for (i in 1:n_counties_occs_calib) {
			for (j in 1:n_facets) {
				Phi_county[i, j] <- phis_county[i, j] + sigmas[j] * inprod(U_star[1:n_facets, j], z_county[i, 1:n_facets])
			}
		}

		# INTEGRATION: site-level calculations
		# This is how the phenotypic variables communicate with abundance. Abundance is naturally estimated at the county level, so we assume that the site-level environmental covariates act like faux counties for abundance, but are indicative of the sites for phenotypic variables.

		# latent variables for correlation		
		for (i in 1:n_pheno_sites) {
			for (j in 1:n_facets) {
				z_site[i, j] ~ dnorm(0, tau = 1)
			}
		}

		# relationship between abundance, biomass, and facets at SITE level with environment
		phi_occs_site[1:n_pheno_sites] <- (x_by_site_occs[1:n_pheno_sites, 1:n_terms_occs] %*% beta_occs[1:n_terms_occs])[ , 1]
		phi_biomass_site[1:n_pheno_sites] <- (x_by_site_biomass[1:n_pheno_sites, 1:n_terms_biomass] %*% beta_biomass[1:n_terms_biomass])[, 1]

		for (i in 1:n_pheno_sites) {

			phi_facet_site[i, 1] <- inprod(beta_facet_1[1:n_terms_facet_1], x_by_site_facet_1[i, 1:n_terms_facet_1])
			phi_facet_site[i, 2] <- inprod(beta_facet_2[1:n_terms_facet_2], x_by_site_facet_2[i, 1:n_terms_facet_2])

			phis_site[i, 1] <- phi_occs_site[i]
			phis_site[i, 2] <- phi_biomass_site[i]
			phis_site[i, 3] <- phi_facet_site[i, 1]
			phis_site[i, 4] <- phi_facet_site[i, 2]

		}

		# use offset method to avoid sampling from non-centered distribution
		# stands in for: Phi_site[i, 1:n_facets] ~ dmnorm(phis_site[i, 1:n_facets], cholesky = U[1:n_facets, 1:n_facets], prec_param = 0)
		for (i in 1:n_pheno_sites) {
			for (j in 1:n_facets) {
				Phi_site[i, j] <- phis_site[i, j] + sigmas[j] * inprod(U_star[1:n_facets, j], z_site[i, 1:n_facets])
			}
		}

		# OCCURRENCE: likelihood
		for (i in 1:n_counties_occs_calib) {
			
			# probability of observing a single AG
			logit(p[i]) <- alpha_occs[1]

			# relationship between expected (latent) abundance and environment assuming MV NORMAL distribution
			log(lambda[i]) <- Phi_county[i, 1]

			# probability of zero inflation
			logit(psi_county[i]) <- inprod(beta_psi[1:n_terms_psi], x_by_county_psi[i, 1:n_terms_psi])

		}

		for (i in 1:n_counties_occs_calib) {
			
			### actual abundance (latent--unobserved)
			# N[i] ~ dHurdlePoisson(lambda = lambda[i], psi = psi_county[i])

			# ### observed number of AG and sampling bias
			# y_n_ag[i] ~ dbinom(size = N[i], prob = p[i])
	
			# joint Poisson-binomial likelihood
			y_n_ag[i] ~ dHurdlePoissonBinomial(lambda = lambda[i], psi = psi_county[i], p = p[i])
	
			# simulate observations for DHARMa residuals
			# y_n_ag_sim[i] ~ dbinom(size = N[i], prob = p[i])
			y_n_ag_sim[i] ~ dHurdlePoissonBinomial(lambda = lambda[i], psi = psi_county[i], p = p[i])

		}

		# log_lik_occs <- sum(log_lik_occs_y[1:n_counties_occs_calib])

		# NON-BIOMASS FACETS
		beta_facet_1[1] ~ dnorm(0, sd = beta_facet_prior_dnorm_sd_1)
		for (i in 2:n_terms_facet_1) {
			beta_facet_1[i] ~ dnorm(0, sd = beta_facet_prior_dnorm_sd) # broad prior
		}

		beta_facet_2[1] ~ dnorm(0, sd = beta_facet_prior_dnorm_sd_1)
		for (i in 2:n_terms_facet_2) {
			beta_facet_2[i] ~ dnorm(0, sd = beta_facet_prior_dnorm_sd) # broad prior
		}

		# prior for sd of individual plant biomass on lognormal, ==> vague
		sigma_biomass_within_sites_log ~ dnorm(0, sd = sigma_biomass_within_sites_log_prior_sd)
		sigma_biomass_within_sites <- exp(sigma_biomass_within_sites_log)

		# BIOMASS and NON-BIOMASS FACETS: parameters of biomass distribution are latent and functions of environment
		# individual plant traits are samples from the site-level distribution defined by the site-level distribution
		for (i in 1:n_pheno_sites) {

			# PROBABILITY OF ZERO INFLATION
			logit(psi_site[i]) <- inprod(beta_psi[1:n_terms_psi], x_by_site_psi[i, 1:n_terms_psi])

		}

		# BIOMASS and NON-BIOMASS: likelihood of individual plants
		for (i in 1:n_biomass) {

			### BIOMASS

			# likelihood
			y_biomass[i] ~ dHLN(meanlog = Phi_site[site_index_biomass[i], 2], sdlog = sigma_biomass_within_sites, psi = psi_site[site_index_biomass[i]])

			# simulated values for unconditional DHARMa residuals
			y_biomass_sim[i] ~ dHLN(meanlog = Phi_site[site_index_biomass[i], 2], sdlog = sigma_biomass_within_sites, psi = psi_site[site_index_biomass[i]])

		}

		# FACETs prior for sd of individual plant facet value on lognormal, ==> vague
		log(sigma_facet_within_sites_1) ~ dnorm(0, sd = sigma_facet_within_sites_log_prior_sd)
		log(sigma_facet_within_sites_2) ~ dnorm(0, sd = sigma_facet_within_sites_log_prior_sd)

		# FACETS: likelihood of individual plant's non-biomass traits
		for (i in 1:n_non_biomass) {

			### FACET 1
			# y_facet_1[i] ~ dHLN(meanlog = Phi_site[site_index_facet[i], 3], sdlog = sigma_facet_within_sites_1, z = z_site[site_index_facet[i]])
			y_facet_1[i] ~ dHLN(meanlog = Phi_site[site_index_facet[i], 3], sdlog = sigma_facet_within_sites_1, psi = psi_site[site_index_facet[i]])

			# simulated values for unconditional DHARMa residuals
			# y_facet_sim_1[i] ~ dHLN(meanlog = Phi_site[site_index_facet[i], 3], sdlog = sigma_facet_within_sites_1, z = z_site[site_index_facet[i]])
			y_facet_sim_1[i] ~ dHLN(meanlog = Phi_site[site_index_facet[i], 3], sdlog = sigma_facet_within_sites_1, psi = psi_site[site_index_facet[i]])

			### FACET 2
			# y_facet_2[i] ~ dHLN(meanlog = Phi_site[site_index_facet[i], 4], sdlog = sigma_facet_within_sites_2, z = z_site[site_index_facet[i]])
			y_facet_2[i] ~ dHLN(meanlog = Phi_site[site_index_facet[i], 4], sdlog = sigma_facet_within_sites_2, psi = psi_site[site_index_facet[i]])

			# simulated values for unconditional DHARMa residuals
			# y_facet_sim_2[i] ~ dHLN(meanlog = Phi_site[site_index_facet[i], 4], sdlog = sigma_facet_within_sites_2, z = z_site[site_index_facet[i]])
			y_facet_sim_2[i] ~ dHLN(meanlog = Phi_site[site_index_facet[i], 4], sdlog = sigma_facet_within_sites_2, psi = psi_site[site_index_facet[i]])

		}

	})

	model_code <- glueNimbleCode(
		model_code,
		model_code_beta_occs_alpha_occs_1_priors,
		model_code_beta_psi_priors,
		model_code_beta_biomass_priors
	)

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
	check_nodes(model)
	say('model$calculate(): ', calc)
	if (is.na(calc) || is.infinite(calc)) stop('Impossible likelihood.')

	say('configureMCMC():', level = 2)

	monitors_coeffs_not_indexed <- c('eta', 'sigma_biomass_within_sites', 'sigma_facet_within_sites_1', 'sigma_facet_within_sites_2')
	monitors_coeffs_single_index <- c('beta_occs', 'alpha_occs', 'beta_biomass', 'beta_psi', 'sigmas', 'beta_facet_1', 'beta_facet_2')
	monitors_coeffs_double_index <- c('correlation')

	monitors_derived_not_indexed <- c()
	monitors_derived_single_index <- c()
	monitors_derived_double_index <- c('U_star')

	monitors_dharma <- c(
		'y_n_ag_sim', 'y_biomass_sim', 'y_facet_sim_1', 'y_facet_sim_2'
	)

	monitors_debug <- c(
		# 'lambda', 'p', 'Phi_county', 'Phi_site'
	)

	monitors <- c(monitors_coeffs_not_indexed, monitors_coeffs_single_index, monitors_coeffs_double_index, monitors_derived_not_indexed, monitors_derived_single_index, monitors_derived_double_index, monitors_dharma, monitors_debug)

	conf <- configureMCMC(
		model,
		monitors = monitors,
		print = FALSE,
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

	vars <- c('beta_occs', 'beta_psi', 'alpha_occs')
	for (var in vars) conf$removeSamplers(var)
	conf$addSampler(target = vars, type = 'AF_slice')
	say('AF_slice sampler added to ', paste(vars, collapse = ' & '), '.')

	# var <- 'beta_psi'
	# conf$removeSamplers(var)
	# conf$addSampler(target = var, type = 'AF_slice')
	# say('AF_slice sampler added to ', paste(var, collapse = ' & '), '.')

	# var <- 'beta_biomass'
	# conf$removeSamplers(var)
	# conf$addSampler(target = var, type = 'AF_slice')
	# say('AF_slice sampler added to ', paste(var, collapse = ' & '), '.')

	for (f in seq_along(nonbiomass_facets)) {

		facet <- names(nonbiomass_facets)[f]
		var <- paste0('beta_facet_', f)
		if (nonbiomass_facets[[f]]$n_terms > 2) {

			conf$removeSamplers(var)
			conf$addSampler(target = var, type = 'AF_slice')
			say('AF_slice sampler added to ', paste(var, collapse = ' & '), '.')

		}

	}

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

	say('session info', level = 2)
	print(sessionInfo())

	say(date(), pre = 1)
	sink()

say('#################################################')
say('### post-modeling diagnostics and predictions ###')
say('#################################################')

	descrip <- paste0('occurrence ~ hurdlePoisson(MVN), biomass ~ hurdleLN(MVN), ', paste(names(nonbiomass_facets), '~ fx(MNV)', collapse = ', '))
	facets_nice <- paste0('occurrence + biomass + ', paste(names(nonbiomass_facets), collapse = ' + '))
	workflow_postmodeling_generic(facet = facets_nice, formulae = formulae, descrip = descrip, out_dir = out_dir)

	workflow_postmodeling_fully_integrated(
		formula_occs = formula_occs,
		formula_bias = formula_bias,
		formula_psi = formula_psi,
		formula_biomass = formula_biomass,
		nonbiomass_facets = nonbiomass_facets,
		out_dir = out_dir
	)
	
	if (crossvalidate) {
	
		workflow_postmodeling_fully_integrated_crossvalidation(
			chains = chains,
			constants = constants,
			inits = inits,
			formula_occs = formula_occs,
			formula_bias = formula_bias,
			formula_psi = formula_psi,
			formula_biomass = formula_biomass,
			nonbiomass_facets = nonbiomass_facets,
			out_dir = out_dir
		)
	
	}

say(date())
say('FINIS!', deco = '+', level = 1)

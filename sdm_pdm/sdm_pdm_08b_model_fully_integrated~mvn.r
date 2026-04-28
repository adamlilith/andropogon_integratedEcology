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

	rm(list = ls())

	drive <- 'C:/Kaji/'

	setwd(paste0(drive, '/Research/Andropogon/Andropogon'))
	source(paste0(drive, '/R/andropogon_integratedEcology/sdm_pdm/sdm_pdm_00_shared_functions_and_variables.r'))

##########################
## user-defined values ###
##########################

	chain <- 6 # also used as seed

	# trial <- TRUE # TRUE for testing
	trial <- FALSE # TRUE for testing

	# calib <- TRUE # use just counties with non-NA Poaceae for calibration region
	calib <- FALSE # use all of North America for calibration region

	# do cross-validation?
	# crossvalidate <- TRUE
	crossvalidate <- FALSE

	### formula for how aspects of species responds to environment

	formula_occs <- ~ 1 + bio1 + bio12_log10p1 + bio15 + I(bio1^2) + I(bio12_log10p1^2) + I(bio15^2) # response of occurrence to climate and soil
	formula_bias <- ~ 1 # sampling bias for AG records
	occs_filename <- 'bio1^2_log(bio12)^2_bio15^2'
	file_occs <- paste0('C:/Kaji/Research/Andropogon/Andropogon/outputs_loretta/integrated_sdm_pdm/models_occurrence/[occs~hurdle_bio1^2_log(bio12)^2_bio15^2]_[bias~1]_niter_160000/chains.rds')

	formula_psi <- ~ 1 + bio1 + bio12_log10p1 + bio15 + I(bio1^2) + I(bio12_log10p1^2) + I(bio15^2)

	formula_biomass <- ~ 1 + bio12 # response of biomass to environment
	biomass_filename <- 'bio12'
	resp_distrib_biomass <- 'hurdleLN'
	transform_biomass <- if (resp_distrib_biomass == 'hGamma') { 'exponential' } else if (resp_distrib_biomass == 'hurdleLN') { 'identity' }
	file_biomass <- paste0('./outputs_loretta/integrated_sdm_pdm/models_biomass/[biomass~', tolower(resp_distrib_biomass), '_', biomass_filename, ']/chains.rds')

	# facets <- c('blade_width', 'canopy_diameter', 'cn_ratio', 'height', 'internal_co2', 'leaf_thickness', 'n_concentration', 'photosynthetic_rate', 'spad', 'stomatal_conductance', 'transpiration_rate')

	nonbiomass_facets <- list(
		blade_width = list(
			formula = ~ 1 + ph + I(ph^2),
			filename = 'ph^2',
			resp_distrib = 'hurdleLN',
			transform = 'identity'
		),
		canopy_diameter = list(
			formula = ~ 1 + bio1 + bio12 + bio1:bio12,
			filename = 'bio1_x_bio12',
			resp_distrib = 'hurdleLN',
			transform = 'identity'
		),
		cn_ratio = list(
			formula = ~ 1 + bio12,
			filename = 'bio12',
			resp_distrib = 'hurdleLN',
			transform = 'identity'
		),
		height = list(
			formula = ~ 1 + bio1 + bio12 + bio1:bio12,
			filename = 'bio1_x_bio12',
			resp_distrib = 'hurdleLN',
			transform = 'identity'
		),
		# internal_co2 = list(
			# formula = ~ 1 + bio12 + nitrogen + bio12:nitrogen,
			# filename = 'bio12_x_nitrogen',
			# resp_distrib = 'hurdleLN',
			# transform = 'identity'
		# ),
		# leaf_thickness = list(
		# 	formula = ~ 1 + bio12,
		# 	filename = 'bio12',
		# 	resp_distrib = 'hurdleLN',
		# 	transform = 'identity'
		# ),
		# n_concentration = list(
		# 	formula = ~ 1 + bio12,
		# 	filename = 'bio12',
		# 	resp_distrib = 'hurdleLN',
		# 	transform = 'identity'
		# ),
		photosynthetic_rate = list(
			formula = ~ 1 + bio1 + bio12 + bio1:bio12,
			filename = 'bio1_x_bio12',
			resp_distrib = 'hurdleLN',
			transform = 'identity'
		),
		spad = list(
			formula = ~ 1 + bio12 + insolation_2000_growing_season_kWh_per_m2,
			filename = 'bio12_srad',
			resp_distrib = 'hurdleLN',
			transform = 'identity'
		),
		stomatal_conductance = list(
			formula = ~ 1 + bio1 + bio12_log10p1 + I(bio1^2),
			filename = 'bio1^2_log(bio12)',
			resp_distrib = 'hurdleLN',
			transform = 'identity'
		),
		transpiration_rate = list(
			formula = ~ 1 + bio1 + bio12_log10p1 + bio1:bio12_log10p1,
			filename = 'bio1_x_log(bio12)',
			resp_distrib = 'hurdleLN',
			transform = 'identity'
		)
	)

	nonbiomass_facets <- nonbiomass_facets[sort(names(nonbiomass_facets))]

	for (f in seq_along(nonbiomass_facets)) {

		facet <- names(nonbiomass_facets)[f]
		form <- nonbiomass_facets[[facet]]$formula
		terms <- attr(terms(form), 'term.labels')
		nonbiomass_facets[[facet]]$n_terms <- length(terms) + 1

	}

	### output folder and bias formula
	nonbiomasss_facet_names <- names(nonbiomass_facets)
	nonbiomasss_facet_names <- gsub(nonbiomasss_facet_names, pattern = '_', replacement = '')
	nonbiomasss_facet_names_short <- paste(substr(nonbiomasss_facet_names, 1, 3), collapse = '_')

	out_dir <- paste0('./outputs_loretta/integrated_sdm_pdm/models_integrated/', ifelse(trial, 'TRIAL_', ''), '[occs~hurdlepoisson(', occs_filename, ')]_[biomass~hurdleln(', biomass_filename, ')]_[', nonbiomasss_facet_names_short, ']_chain_', chain)

	if (!trial) {

		### MCMC settings
		# need these settings to achieve independent samples using bio1^2, bio12^2, bio15^
		niter <- 400000
		# nchains <- 1
		nchains <- 1

	} else {
		
		### MCMC settings FOR TESTING
		niter <- 2000
		nchains <- 2

	}
	nburnin <- niter / 2
	thin <- (niter - nburnin) / 1000

	set.seed(chain)

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
	say('trial ........................ ', trial)
	say('seed ......................... ', chain)
	say('nburnin ...................... ', nburnin)
	say('thin ......................... ', thin)
	say('nchains ...................... ', nchains)
	say('formula_occs ................. ', paste(as.character(formula_occs), collapse = ' '))
	say('formula_psi .................. ', paste(as.character(formula_psi), collapse = ' '))
	say('formula_bias ................. ', paste(as.character(formula_bias), collapse = ' '))
	say('formula_biomass .............. ', paste(as.character(formula_biomass), collapse = ' '))
	say('n_nonbiomass_facets .......... ', length(nonbiomass_facets))
	for (f in seq_along(nonbiomass_facets)) {
		say('facet ', f, ' ...................... ', names(nonbiomass_facets)[f])
		say('formula ', f, ' ....................... ', paste(as.character(nonbiomass_facets[[f]]$formula), collapse = ' '))
		say('resp_distrib ', f, ' .................. ', nonbiomass_facets[[f]]$resp_distrib)
		say('transform ', f, ' ..................... ', nonbiomass_facets[[f]]$transform)
	}

	say('out_dir')
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

	########################
	### data preparation ###
	########################

	# data for OCCURRENCES at counties
	data_occs_counties <- prepare_occurrence_data(formula_occs = formula_occs, formula_bias = formula_bias, n_response_curve_values = n_response_curve_values, psa_quant = psa_quant, calib = calib)

	# data for OCCURRENCES using site-level environment
	data_occs_sites <- prepare_biomass_data(formula_biomass = formula_occs, n_response_curve_values = n_response_curve_values, calib = calib)

	# data for ZERO-INFLATION at counties
	data_psi_counties <- prepare_occurrence_data(formula_occs = formula_psi, formula_bias = formula_bias, n_response_curve_values = n_response_curve_values, psa_quant = psa_quant, calib = calib)

	# data for ZERO-INFLATION at sites
	data_psi_sites <- prepare_biomass_data(formula_biomass = formula_psi, n_response_curve_values = n_response_curve_values, calib = calib)

	# data for BIOMASS at sites
	data_biomass_sites <- prepare_biomass_data(formula_biomass = formula_biomass, n_response_curve_values = n_response_curve_values, calib = calib)

	# data for BIOMASS using county-level environment
	data_biomass_counties <- prepare_occurrence_data(formula_occs = formula_biomass, formula_bias = ~ 1, n_response_curve_values = n_response_curve_values, psa_quant = psa_quant, calib = calib)

	# data for NON-BIOMASS FACETS at sites and counties
	data_nonbiomass_sites <- data_nonbiomass_counties <- list()
	for (f in seq_along(nonbiomass_facets)) {

		facet <- names(nonbiomass_facets)[f]
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
		n_nonbiomass_facets = n_nonbiomasss_facets,

		### occurrences
		n_counties_occs_calib = data_occs_counties$n_counties_occs_calib, # number of counties in calibration region
		n_terms_occs = data_occs_counties$n_terms, # number of terms in formula for occurrence model (including intercept)
		n_terms_bias = data_occs_counties$n_terms_bias, # number of terms in sampling bias model
		w_bias = data_occs_counties$w_bias, # model matrix of sampling bias of AG observed occurrences
		
		# n_covariates_occs = data_occs_counties$n_covariates, # number of covariates in formula for occurrence model
		# n_covariates_bias = data_occs_counties$n_covariates_bias, # number of covariates in formula for occurrence model
		# resp_curves_x_occs = data_occs_counties$resp_curves_x, # response curve array for occurrences vs environment
		# resp_curves_w = data_occs_counties$resp_curves_w, # response curve array for occurrences vs environment

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

		# n_covariates_biomass = data_biomass_sites$n_covariates,
		# resp_curves_x_biomass = data_biomass_sites$resp_curves_x, # response curve array for biomass

		# facets (generic)
		n_non_biomass = data_nonbiomass_sites[[1]]$n_plants,
		site_index_facet = data_nonbiomass_sites[[1]]$site_index_facet, # index of sampled site for each row in facet data

		# zero inflation
		n_terms_psi = data_psi_counties$n_terms, # number of terms in sampling bias model
		# n_covariates_psi = data_psi_counties$n_covariates, # number of terms in sampling bias model
		x_by_county_psi = data_psi_counties$counties_x_sq,
		x_by_site_psi = data_psi_sites$x_by_site
		# resp_curves_x_psi = data_psi_counties$resp_curves_x

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
	chains_inits <- readRDS(file_occs)
	
	alpha_occs_inits <- mc_extract(chains_inits, 'alpha_occs', j = TRUE)
	beta_occs_mu_inits <- mc_extract(chains_inits, 'beta_occs', j = TRUE)
	beta_psi_inits <- mc_extract(chains_inits, 'beta_psi', j = TRUE)
	
	# lambda_sq_inits <- mc_extract(chains_inits, 'lambda', j = TRUE)
	# lambda_sq_inits[lambda_sq_inits > 10] <- 10

	# lambda_sigma_inits <- mc_extract(chains_inits, 'lambda_sigma')
	# log_lambda_sigma_inits <- log(lambda_sigma_inits)

	rm(chains_inits)

	N_inits_calib <- data_occs_counties$y_n_ag * 2

	# biomass initializations
	chains_inits <- readRDS(file_biomass)
	beta_biomass_inits <- mc_extract(chains_inits, 'beta_biomass', j = TRUE)
	biomass_site_inits <- data_biomass_sites$site_mean
	sigma_biomass_within_sites_inits <- mc_extract(chains_inits, 'sigma_biomass_within_sites')
	log_sigma_biomass_within_sites_inits <- log(sigma_biomass_within_sites_inits)

	rm(chains_inits)

	inits <- list(

		### occurrences
		y_n_ag_sim = data_occs_counties$y_n_ag, # simulated values of observed number of AG (for DHARMa residuals)
		# lambda = lambda_sq_inits, # expected value of number of AG
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
		log_sigmas = c(1, 1),
		Phi_county = matrix(1, nrow = data_occs_counties$n_counties, ncol = n_facets),
		Phi_site = matrix(1, nrow = data_biomass_sites$n_pheno_sites, ncol = n_facets),

		### zero inflation
		beta_psi = beta_psi_inits#,
		# z_site = rep(1, data_biomass_sites$n_pheno_sites), # present/absent at each site
		# z_county = z_county_inits # present/absent in county

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

	occs_biomass_code <- nimbleCode({

		# INTEGRATION: LJK prior for standard deviations and correlations between latent occurrence and biomass processes
		eta ~ dgamma(2, 1)
		U_star[1:n_facets, 1:n_facets] ~ dlkj_corr_cholesky(eta = eta, p = n_facets)
		U[1:n_facets, 1:n_facets] <- uppertri_mult_diag(
			U_star[1:n_facets, 1:n_facets],
			sigmas[1:n_facets]
		)
	  
		# INTEGRATION: standard deviations of latent occurrence and biomass
		log(sigmas[1]) ~ dnorm(0, sd = lambda_sigma_prior_sd) # half-Cauchy
		log(sigmas[2]) ~ dnorm(0, sd = sigma_biomass_among_sites_log_prior_sd) # half-Cauchy

		# INTEGRATION: standard deviations of latent non-biomass
		# start counter at 3 bc occurrences are indexed by 1 and biomass by 2
		for (i in 3:n_facets) {
			log(sigmas[i]) ~ dnorm(0, sd = sigma_facet_among_sites_log_prior_sd) # half-Cauchy
		}

		# correlation matrix
		correlation[1:n_facets, 1:n_facets] <- t(U_star[1:n_facets, 1:n_facets]) %*% U_star[1:n_facets, 1:n_facets]

		# OCCURRENCE: likelihood
		for (i in 1:n_counties_occs_calib) {
			
			### actual abundance (latent--unobserved)
			N[i] ~ dHurdlePoisson(lambda[i], psi = psi_county[i])

			### observed number of AG and sampling bias
			y_n_ag[i] ~ dbinom(size = N[i], prob = p[i])

			# simulate observations for DHARMa residuals
			y_n_ag_sim[i] ~ dbinom(size = N[i], prob = p[i])

			# probability of observing a single AG
			logit(p[i]) <- alpha_occs[1]

			# relationship between expected (latent) abundance and environment assuming MV NORMAL distribution
			log(lambda[i]) <- Phi_county[i, 1]

			# probability of zero inflation
			logit(psi_county[i]) <- inprod(beta_psi[1:n_terms_psi], x_by_county_psi[i, 1:n_terms_psi])
			# z_county[i] ~ dbern(psi_county[i])

			# # likelihood
			# log_lik_occs_y[i] <- dbinom(y_n_ag[i], size = N[i], prob = p[i], log = 1)

		}

		# log_lik_occs <- sum(log_lik_occs_y[1:n_counties_occs_calib])

		# prior for sd of individual plant biomass on lognormal (~ half-Cauchy), ==> vague
		sigma_biomass_within_sites_log ~ dnorm(0, sd = sigma_biomass_within_sites_log_prior_sd)
		sigma_biomass_within_sites <- exp(sigma_biomass_within_sites_log)

		# BIOMASS: parameters of biomass distribution are latent and functions of environment
		# individual plant traits are samples from the site-level distribution defined by the site-level distribution
		for (i in 1:n_pheno_sites) {

			# PROBABILITY OF ZERO INFLATION
			logit(psi_site[i]) <- inprod(beta_psi[1:n_terms_psi], x_by_site_psi[i, 1:n_terms_psi])

			# # presence/absence at the site level... NB somewhat of a cheat bc we could have absence at county level but presence in site in county
			# z_site[i] ~ dbern(psi_site[i])

			# ### BIOMASS
			# biomass_site[i] <- exp(Phi_site[i, 2])

		}

		# BIOMASS: likelihood of individual plants
		for (i in 1:n_biomass) {

			# likelihood
			# y_biomass[i] ~ dHLN(meanlog = Phi_site[site_index_biomass[i], 2], sdlog = sigma_biomass_within_sites, z = z_site[site_index_biomass[i]])
			y_biomass[i] ~ dHLN(meanlog = Phi_site[site_index_biomass[i], 2], sdlog = sigma_biomass_within_sites, psi = psi_site[site_index_biomass[i]])

			# simulated values for unconditional DHARMa residuals
			# y_biomass_sim[i] ~ dHLN(meanlog = Phi_site[site_index_biomass[i], 2], sdlog = sigma_biomass_within_sites, z = z_site[site_index_biomass[i]])
			y_biomass_sim[i] ~ dHLN(meanlog = Phi_site[site_index_biomass[i], 2], sdlog = sigma_biomass_within_sites, psi = psi_site[site_index_biomass[i]])

		}


	})

	### model code for NON-BIOMASS FACETS
	#####################################

	# Written generically so we can model code each facet in a loop, but the code is the same for each facet except for the index of the facet in the Phi_site matrix and the number of terms in the formula for that facet. We add 2 to the index of the facet in the Phi_site matrix to skip the occurrence and biomass terms.

	code_facets_raw	 <- '{

		# location of FACET XYZ in Phi_site matrix... add 2 to skip occurrence and biomass term
		index_XYZ <- XYZ + 2

		# NON-BIOMASS FACET XYZ: priors
		beta_facet_XYZ[1] ~ dnorm(0, sd = beta_facet_prior_dnorm_sd_1)
		for (i in 2:n_terms_facet_XYZ) {
			beta_facet_XYZ[i] ~ dnorm(0, sd = beta_facet_prior_dnorm_sd) # broad prior
		}

		# NON-BIOMASS FACET XYZ: parameters of biomass distribution are latent and functions of environment
		# individual plant traits are samples from the site-level distribution defined by the site-level distribution
		for (i in 1:n_pheno_sites) {
			mu_facet_site_XYZ[i] <- exp(Phi_site[i, index_XYZ])
		}

		# NON-BIOMASS FACET XYZ: prior for sd of individual plant facet value on lognormal (~ half-Cauchy), ==> vague
		log(sigma_facet_within_sites_XYZ) ~ dnorm(0, sd = sigma_facet_within_sites_log_prior_sd)

		# NON-BIOMASS FACET XYZ: likelihood of the non-biomass trait of an individual plant
		for (i in 1:n_non_biomass) {

			# likelihood
			y_facet_XYZ[i] ~ dHLN(meanlog = Phi_site[site_index_facet[i], index_XYZ], sdlog = sigma_facet_within_sites_XYZ, psi = psi_site[site_index_facet[i]])

			# simulated values for unconditional DHARMa residuals... add 2 
			y_facet_sim_XYZ[i] ~ dHLN(meanlog = Phi_site[site_index_facet[i], index_XYZ], sdlog = sigma_facet_within_sites_XYZ, psi = psi_site[site_index_facet[i]])

		}
		
		# INTEGRATION: site-level calculations
		# This is how the phenotypic variables communicate with abundance. Abundance is naturally estimated at the county level, so we assume that the site-level environmental covariates act like faux counties for abundance, but are indicative of the sites for phenotypic variables.
		for (i in 1:n_pheno_sites) {
			phi_facet_site[i, XYZ] <- inprod(beta_facet_XYZ[1:n_terms_facet_XYZ], x_by_site_facet_XYZ[i, 1:n_terms_facet_XYZ])
		}

		# INTEGRATION: county-level
		# NB We have measurements of biomass at the site level, but not county. We thus assume that the estimates of non-biomass facets (and biomass) using county-level covariates are indicative of virtual sites that have the same environments as counties. Biomass and non-biomass facets are still estimated at the site level.
		for (i in 1:n_counties_occs_calib) {
			phi_facet_county[i, XYZ] <- inprod(beta_facet_XYZ[1:n_terms_facet_XYZ], x_by_county_facet_XYZ[i, 1:n_terms_facet_XYZ])
		}

	}'

	code_nonbiomass_facets <- list()
	for (i in seq_along(nonbiomass_facets)) {

		code_facets_raw_this_facet <- gsub(pattern = 'XYZ', replacement = i, x = code_facets_raw)
		code_nonbiomass_facets[[i]] <- code_facets_raw_this_facet

	}
	code_nonbiomass_facets <- lapply(code_nonbiomass_facets, function(x) parse(text = x)[[1]])

	# code for sampling MVN for sites and counties
	code_integration <- paste0('{

		for (i in 1:n_pheno_sites) {

			phis_site[i, 1:n_facets] <- c(
				phi_occs_site[i],
				phi_biomass_site[i], ',
				paste0('\nphi_facet_site[i, ', seq_along(nonbiomass_facets), ']', collapse = ', '),
			')
			Phi_site[i, 1:n_facets] ~ dmnorm(phis_site[i, 1:n_facets], cholesky = U[1:n_facets, 1:n_facets], prec_param = 0)

			phi_occs_site[i] <- inprod(beta_occs[1:n_terms_occs], x_by_site_occs[i, 1:n_terms_occs])
			phi_biomass_site[i] <- inprod(beta_biomass[1:n_terms_biomass], x_by_site_biomass[i, 1:n_terms_biomass])

		}

		# INTEGRATION: county-level
		# NB We have measurements of biomass at the site level, but not county. We thus assume that the estimates of biomass using county-level covariates are indicative of virtual sites that have the same environments as counties. Biomass is still estimated at the site level.
		for (i in 1:n_counties_occs_calib) {

			phis_county[i, 1:n_facets] <- c(
				phi_occs_county[i],
				phi_biomass_county[i], ',
				paste0('\nphi_facet_county[i, ', seq_along(nonbiomass_facets), ']', collapse = ', '),
			')
			
			Phi_county[i, 1:n_facets] ~ dmnorm(phis_county[i, 1:n_facets], cholesky = U[1:n_facets, 1:n_facets], prec_param = 0)

			phi_occs_county[i] <- inprod(beta_occs[1:n_terms_occs], x_by_county_occs[i, 1:n_terms_occs])
			phi_biomass_county[i] <- inprod(beta_biomass[1:n_terms_biomass], x_by_county_biomass[i, 1:n_terms_biomass])

		}

	}')
	code_integration <- parse(text = code_integration)[[1]]

	model_code <- glueNimbleCode(
		code_integration,
		occs_biomass_code,
		model_code_beta_occs_alpha_occs_1_priors,
		model_code_beta_psi_priors,
		model_code_beta_biomass_priors
	)
	model_code <- Reduce(
		f = function(code_a, code_b) glueNimbleCode(code_a, code_b),
		x = code_nonbiomass_facets,
		init = model_code
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

	monitors_coeffs_not_indexed <- c('eta', 'sigma_biomass_within_sites')
	monitors_coeffs_single_index <- c('beta_occs', 'alpha_occs', 'beta_psi', 'beta_biomass', 'sigmas')
	monitors_coeffs_double_index <- c('U')

	monitors_derived_not_indexed <- c()
	monitors_derived_single_index <- c()
	monitors_derived_double_index <- c('correlation')

	for (f in seq_along(nonbiomass_facets)) {

		monitors_coeffs_not_indexed <- c(monitors_coeffs_not_indexed, paste0('sigma_facet_within_sites_', f))
		monitors_coeffs_single_index <- c(monitors_coeffs_single_index, paste0('beta_facet_', f))

		monitors_derived_single_index <- c(monitors_derived_single_index, paste0('mu_facet_site_', f))


	}

	monitors_dharma <- c(
		'y_n_ag_sim', 'y_biomass_sim', paste0('y_facet_sim_', seq_along(nonbiomass_facets))
	)

	monitors <- c(monitors_coeffs_not_indexed, monitors_coeffs_single_index, monitors_coeffs_double_index, monitors_derived_not_indexed, monitors_derived_single_index, monitors_derived_double_index, monitors_dharma)

	conf <- configureMCMC(
		model,
		monitors = monitors,
		print = TRUE,
		enableWAIC = FALSE
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

	vars <- c('beta_occs', 'beta_psi', 'beta_biomass')
	if (data_occs_counties$n_covariates_bias >= 1) vars <- c(vars, 'alpha_occs')
	for (var in vars) {
		conf$removeSamplers(var) 
		conf$addSampler(target = var, type = 'AF_slice')
		say('AF_slice sampler added to ', var, '.')
	}

	for (f in seq_along(nonbiomass_facets)) {

		facet <- names(nonbiomass_facets)[f]
		var <- paste0('beta_facet_', f)
		# if (nonbiomass_facets[[f]]$n_terms > 2) {

			conf$removeSamplers(var)
			conf$addSampler(target = var, type = 'AF_slice')
			say('AF_slice sampler added to ', paste(var, collapse = ' & '), '.')

		# }

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
		WAIC = FALSE,
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

	descrip <- paste0('occurrence ~ hurdlePoisson(MVN), biomass ~ hurdleLN(MVN)')
	for (f in seq_along(nonbiomass_facets)) {
		descrip <- paste0(descrip, ', ', names(nonbiomass_facets)[f], ' ~ ', tolower(nonbiomass_facets[[f]]$resp_distrib), '(MVN)')
	}
	facets_nice <- paste0('occurrence + biomass + ', paste(names(nonbiomass_facets), collapse = ' + '))
	workflow_postmodeling_generic(facet = facets_nice, formulae = formulae, descrip = descrip, out_dir = out_dir)
print(STOP)
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

### MODELING ANDROPOGON GERARDI DISTRIBUTION, PHENOTYPE, PHYSIOLOGY, GENOTYPE, and ASSOCIATED MICROBIAL COMMUNITIES
### Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2023-12
###
### This model estimates the distribution of site-level mean values of two non-biomass "facets" of Andropogon gerardi ramets (e.g., height, SPAD, leaf nitrogen, etc.). The distribution of facet values among ramets at a site follows a zero-inflated gamma or zero-inflated lognormal. The site-level mean of the facet is a function of soil/climate and is drawn from a multivariate normal distribution which accounts for inter-facet correlations. The probability of a zero value is a function of the same covariate(s) and has the same functional form as the submodel of the mean.
###
### source('C:/Kaji/R/andropogon_integratedEcology/sdm_pdm/sdm_pdm_07a_model_two_nonbiomass_facets_zero_inflated~mvn.r')

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

	facet_1 <- 'height'
	facet_2 <- 'canopy_diameter'

	# trial <- TRUE # TRUE for testing
	trial <- FALSE # FALSE for running for real

	# do cross-validation?
	crossvalidate <- TRUE
	# crossvalidate <- FALSE

	formula_facet_1 <- ~ 1 + bio12
	formula_facet_2 <- ~ 1 + bio12
	formula_psi <- ~ 1 + bio12

	filename_occs_facet_1 <- 'bio12'
	filename_occs_facet_2 <- 'bio12'

	# distribution of response variable
	# ZIG = zero-inflated gamma
	# 'hurdleLN' = zero-inflated lognormal
	resp_distrib_facet_1 <- 'hurdleLN'
	resp_distrib_facet_2 <- 'hurdleLN'

	### MCMC settings
	if (!trial) {

		niter <- 400000
		nchains <- 4

	} else {
		
		### MCMC settings FOR TESTING
		niter <- 2000
		nchains <- 2

	}
	nburnin <- niter / 2
	thin <- (niter - nburnin) / 1000

	# do not change
	zero_inflated <- TRUE # SPECIFIC TO THIS SCRIPT--SHOULD NOT BE CHANGED

#################
### main loop ###
#################

	# cycle across all facets, response (zero-inflated gamma or lognormal), and formulae

	transform_facet_1 <- if (resp_distrib_facet_1 == 'hGamma') { 'exponential' } else if (resp_distrib_facet_1 == 'hurdleLN') { 'identity' }
	transform_facet_2 <- if (resp_distrib_facet_2 == 'hGamma') { 'exponential' } else if (resp_distrib_facet_2 == 'hurdleLN') { 'identity' }

	out_dir <- paste0('./outputs_loretta/integrated_sdm_pdm/models_', facet_1, '_', facet_2, '/', ifelse(trial, 'TRIAL_', ''), '[', facet_1, '~', tolower(resp_distrib_facet_1), '~', filename_occs_facet_1, ']_', '[', facet_2, '~', tolower(resp_distrib_facet_2), '~', filename_occs_facet_2, ']')

	if (!trial & file.exists(out_dir)) stop('Output folder already exists.')

	dirCreate(out_dir)

	sink(paste0(out_dir, '/runtime_log.txt'), split = TRUE)
	say('MODELING ', toupper(facet_1), ' AND ', toupper(facet_2))
	say('Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | ', date(), post = 2)

	say('This model estimates the distribution of site-level mean values of two non-biomass "facets" of Andropogon gerardi ramets (e.g., height, SPAD, leaf nitrogen, etc.). The distribution of facet values among ramets at a site follows a zero-inflated gamma or zero-inflated lognormal. The site-level mean of the facet is a function of soil/climate and is drawn from a multivariate normal distribution which accounts for inter-facet correlations. The probability of a zero value is a function of the same covariate(s) and has the same functional form as the submodel of the mean.', breaks = 60, post = 1)

	say('MCMC settings:', level = 2)
	say('facet_1 ......................... ', facet_1)
	say('facet_2 ......................... ', facet_2)
	say('trial ........................... ', trial)
	say('niter ........................... ', niter)
	say('nburnin ......................... ', nburnin)
	say('thin ............................ ', thin)
	say('nchains ......................... ', nchains)
	say('formula_facet_1 ................. ', paste(as.character(formula_facet_1), collapse = ' '))
	say('formula_facet_2 ................. ', paste(as.character(formula_facet_2), collapse = ' '))
	say('formula_psi ..................... ', paste(as.character(formula_psi), collapse = ' '))
	say('resp_distrib_facet_1 ............ ', resp_distrib_facet_1)
	say('resp_distrib_facet_2 ............ ', resp_distrib_facet_2)
	say('transform_facet_1 ............... ', transform_facet_1)
	say('transform_facet_2 ............... ', transform_facet_2)
	say('zero_inflated ................... ', zero_inflated)
	say('calib ........................... ', calib, post = 2)

	say('out_dir:')
	say(out_dir)

	formulae <- list(
		formula_facet_1 = formula_facet_1,
		formula_facet_2 = formula_facet_2,
		formula_psi = formula_psi
	)
	saveRDS(formulae, paste0(out_dir, '/formulae.rds'))

	########################s
	### data preparation ###
	########################

		data_facet_1 <- prepare_nonbiomass_data(facet = facet_1, formula_facet = formula_facet_1, n_response_curve_values = n_response_curve_values, calib = calib)

		data_facet_2 <- prepare_nonbiomass_data(facet = facet_2, formula_facet = formula_facet_2, n_response_curve_values = n_response_curve_values, calib = calib)

		# using "height" bc the function needs a valid facet name
		data_psi <- prepare_nonbiomass_data(facet = 'height', formula_facet = formula_psi, n_response_curve_values = n_response_curve_values, calib = calib)

	#########################
	### inputs for nimble ###
	#########################

	say('Inputs:', level = 2)
	data <- list(
		y_facet_1 = data_facet_1$y_facet,
		y_facet_2 = data_facet_2$y_facet
	)

	constants <- list(
		
		n_facets = 2,

		# sampled sites
		n_pheno_sites = data_facet_1$n_pheno_sites, # number of phenotype sample sites

		# facet #1
		x_by_site_facet_1 = data_facet_1$x_by_site, # MM with covariates for facet (scaled)
		n_terms_facet_1 = data_facet_1$covariates, # number of terms in formula for facet model (including intercept)
		n_covariates_facet_1 = data_facet_1$n_covariates,
		resp_curves_x_facet_1 = data_facet_1$resp_curves_x, # response curve array for facet

		# facet #2
		x_by_site_facet_2 = data_facet_2$x_by_site, # MM with covariates for facet (scaled)
		n_terms_facet_2 = data_facet_2$covariates, # number of terms in formula for facet model (including intercept)
		n_covariates_facet_2 = data_facet_2$n_covariates,
		resp_curves_x_facet_2 = data_facet_2$resp_curves_x, # response curve array for facet

		# zero-inflation
		x_by_site_psi = data_psi$x_by_site,
		n_covariates_psi = data_psi$n_covariates,
		n_terms_psi = data_psi$covariates,

		n_plants = data_facet_1$n_plants, # number of facet observations
		site_index_facet = data_facet_1$site_index_facet # index of sampled site for each row in facet data

	)

	constants <- c(constants, constants_shared_facet)
	if (!is.null(formula_psi)) constants <- c(constants, constants_shared_psi)

	# get initial values from standalone models
	chains_facet_1 <- readRDS(paste0('./outputs_loretta/integrated_sdm_pdm/models_', facet_1, '/[', facet_1, '_', tolower(resp_distrib_facet_1), '~normal~', filename_occs_facet_1, ']/chains.rds'))
	beta_facet_1_inits <- mc_extract(chains_facet_1, 'beta_facet', j = TRUE)
	sigma_facet_1_within_sites_init <- exp(mc_extract(chains_facet_1, 'sigma_facet_within_sites'))
	sigma_facet_1_among_sites_init <- mc_extract(chains_facet_1, 'sigma_facet_among_sites')

	chains_facet_2 <- readRDS(paste0('./outputs_loretta/integrated_sdm_pdm/models_', facet_2, '/[', facet_2, '_', tolower(resp_distrib_facet_2), '~normal~', filename_occs_facet_2, ']/chains.rds'))
	beta_facet_2_inits <- mc_extract(chains_facet_2, 'beta_facet', j = TRUE)
	sigma_facet_2_within_sites_init <- mc_extract(chains_facet_2, 'sigma_facet_within_sites')
	sigma_facet_2_among_sites_init <- mc_extract(chains_facet_2, 'sigma_facet_among_sites')

	Phi_site_inits <- matrix(NA_real_, nrow = constants$n_pheno_sites, ncol = constants$n_facets)
	
	preds_facet_1 <- predict_nonbiomass_single_trait(chains_facet_1, x = constants$x_by_site_facet_1, resp_distrib = resp_distrib_facet_1, transform = transform_facet_1)
	preds_facet_2 <- predict_nonbiomass_single_trait(chains_facet_2, x = constants$x_by_site_facet_2, resp_distrib = resp_distrib_facet_2, transform = transform_facet_2)

	Phi_site_inits[ , 1] <- log(colMeans(preds_facet_1))
	Phi_site_inits[ , 2] <- log(colMeans(preds_facet_2))

	rm(chains_facet_1, chains_facet_2)

	beta_psi_inits <- rep(0, 1 + length(attr(terms(formula_psi), 'term.labels')))

	inits <- list(

		# integration
		eta = 1,
		U_star = diag(1, 2),
		log_sigmas = log(c(sigma_facet_1_among_sites_init, sigma_facet_2_among_sites_init)),
		Phi_site = Phi_site_inits,

		# facet #1
		log_site_facet_1_mu = log(unlist(data_facet_1$site_vect[[paste0(facet_1, '_mean')]])),
		log_sigma_facet_1_within_sites = log(sigma_facet_1_within_sites_init),
		y_facet_1_sim = data_facet_1$y_facet, # simulated values for facet (for DHARMa residuals)
		beta_facet_1 = beta_facet_1_inits,

		# facet #2
		log_site_facet_2_mu = log(unlist(data_facet_2$site_vect[[paste0(facet_2, '_mean')]])),
		log_sigma_facet_2_within_sites = log(sigma_facet_2_within_sites_init),
		y_facet_2_sim = data_facet_2$y_facet, # simulated values for facet (for DHARMa residuals)
		beta_facet_2 = beta_facet_2_inits,

		# zero-fixation
		beta_psi = beta_psi_inits#,
		# z_site = rep(1, data_facet_1$n_pheno_sites)

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
	model_code <- nimbleCode({
	
		# INTEGRATION: LJK prior for standard deviations and correlations between latent facets 1 and 2
		eta ~ dgamma(2, 1)
		U_star[1:n_facets, 1:n_facets] ~ dlkj_corr_cholesky(eta = eta, p = n_facets)
		U[1:n_facets, 1:n_facets] <- uppertri_mult_diag(
			U_star[1:n_facets, 1:n_facets],
			sigmas[1:n_facets]
		)
	  
		# INTEGRATION: standard deviations of latent facets 1 and 2
		for (i in 1:n_facets) {
			log(sigmas[i]) ~ dnorm(0, sd = sigma_facet_among_sites_log_prior_sd) # half-Cauchy
		}

		# correlation matrix
		correlation[1:n_facets, 1:n_facets] <- t(U_star[1:n_facets, 1:n_facets]) %*% U_star[1:n_facets, 1:n_facets]

		# FACETS: priors for relationship of site-level mean facet to environment
		beta_facet_1[1] ~ dnorm(0, sd = beta_facet_prior_dnorm_sd_1)
		for (i in 2:n_terms_facet_1) {
			beta_facet_1[i] ~ dnorm(0, sd = beta_facet_prior_dnorm_sd) # broad prior
		}

		beta_facet_2[1] ~ dnorm(0, sd = beta_facet_prior_dnorm_sd_1)
		for (i in 2:n_terms_facet_2) {
			beta_facet_2[i] ~ dnorm(0, sd = beta_facet_prior_dnorm_sd) # broad prior
		}

		# prior for sd of values WITHIN a site on lognormal (~ half-Cauchy), ==> vague
		log(sigma_facet_1_within_sites) ~ dnorm(0, sd = sigma_facet_within_sites_log_prior_sd)
		log(sigma_facet_2_within_sites) ~ dnorm(0, sd = sigma_facet_within_sites_log_prior_sd)

		# FACETS: parameters of facet distribution are latent and functions of environment
		# individual plant biomasses are samples from the site-level distribution defined by the site-level distribution (next chunk after this one)
		for (i in 1:n_pheno_sites) {

			# INTEGRATION
			phis_site[i, 1:n_facets] <- c(log_site_facet_1_mu[i], log_site_facet_2_mu[i])
			Phi_site[i, 1:n_facets] ~ dmnorm(phis_site[i, 1:n_facets], cholesky = U[1:n_facets, 1:n_facets], prec_param = 0)

			# relationship of FACETs 1 and 2 to the environment
			log_site_facet_1_mu[i] <- inprod(beta_facet_1[1:n_terms_facet_1], x_by_site_facet_1[i, 1:n_terms_facet_1])
			log_site_facet_2_mu[i] <- inprod(beta_facet_2[1:n_terms_facet_2], x_by_site_facet_2[i, 1:n_terms_facet_2])
			
			# zero-inflation
			logit(psi[i]) <- inprod(beta_psi[1:n_terms_psi], x_by_site_psi[i, 1:n_terms_psi])
			# z_site[i] ~ dbern(psi[i])

			# site-level mean facet value
			mu_facet_1_site[i] <- exp(Phi_site[i, 1])
			mu_facet_2_site[i] <- exp(Phi_site[i, 2])

		}

		log_lik_facet_1 <- sum(log_lik_facet_1_y[1:n_plants])
		log_lik_facet_2 <- sum(log_lik_facet_2_y[1:n_plants])
		log_lik <- log_lik_facet_1 + log_lik_facet_2

	})

	### response code for FACET 1
	#############################
	# Switch out zero-inflated gamma or lognormal, depending on which one we're doing.
	if (resp_distrib_facet_1 == 'hGamma') {
	
		response_distrib_code_facet_1 <- nimbleCode({

			# FACET: parameters of facet distribution are latent and functions of environment
			# individual plant facet values are samples from the site-level distribution defined by the site-level distribution (next chunk after this one)
			for (i in 1:n_pheno_sites) {

				# moment matching to get gamma() parameters
				shape_facet_1[i] <- mu_facet_1_site[i]^2 / sigma_facet_1_within_sites^2
				rate_facet_1[i] <- mu_facet_1_site[i] / sigma_facet_1_within_sites^2

			}

			# FACET: likelihood of individual plants
			for (i in 1:n_plants) {

				# likelihood
				# y_facet_1[i] ~ dZIG(shape = shape_facet_1[site_index_facet[i]], rate = rate_facet_1[site_index_facet[i]], z = z_site[site_index_facet[i]])
				y_facet_1[i] ~ dZIG(shape = shape_facet_1[site_index_facet[i]], rate = rate_facet_1[site_index_facet[i]], psi = psi_site[site_index_facet[i]])

				# simulated values for unconditional DHARMa residuals
				# y_facet_1_sim[i] ~ dZIG(shape = shape_facet_1[site_index_facet[i]], rate = rate_facet_1[site_index_facet[i]], z = z_site[site_index_facet[i]])
				y_facet_1_sim[i] ~ dZIG(shape = shape_facet_1[site_index_facet[i]], rate = rate_facet_1[site_index_facet[i]], psi = psi_site[site_index_facet[i]])

				# log_lik_facet_1_y[i] <- dZIG(y_facet_1[i], shape = shape_facet_1[site_index_facet[i]], rate = rate_facet_1[site_index_facet[i]], z = z_site[site_index_facet[i]])
				log_lik_facet_1_y[i] <- dZIG(y_facet_1[i], shape = shape_facet_1[site_index_facet[i]], rate = rate_facet_1[site_index_facet[i]], psi = psi_site[site_index_facet[i]])
		
			}

		})
	
	} else if (resp_distrib_facet_1 == 'hurdleLN') {
	
		response_distrib_code_facet_1 <- nimbleCode({

			# FACET: likelihood of individual plants
			for (i in 1:n_plants) {

				# likelihood
				# y_facet_1[i] ~ dZILN(meanlog = Phi_site[site_index_facet[i], 1], sdlog = sigma_facet_1_within_sites, z = z_site[site_index_facet[i]])
				y_facet_1[i] ~ dZILN(meanlog = Phi_site[site_index_facet[i], 1], sdlog = sigma_facet_1_within_sites, psi = psi_site[site_index_facet[i]])

				# simulated values for unconditional DHARMa residuals
				# y_facet_1_sim[i] ~ dZILN(meanlog = Phi_site[site_index_facet[i], 1], sdlog = sigma_facet_1_within_sites, z = z_site[site_index_facet[i]])
				y_facet_1_sim[i] ~ dZILN(meanlog = Phi_site[site_index_facet[i], 1], sdlog = sigma_facet_1_within_sites, psi = psi_site[site_index_facet[i]])

				# log_lik_facet_1_y[i] <- dZILN(y_facet_1[i], meanlog = Phi_site[site_index_facet[i], 1], sdlog = sigma_facet_1_within_sites, z = z_site[site_index_facet[i]])
				log_lik_facet_1_y[i] <- dZILN(y_facet_1[i], meanlog = Phi_site[site_index_facet[i], 1], sdlog = sigma_facet_1_within_sites, psi = psi_site[site_index_facet[i]])

			}

		})
		
	}


	### response code for FACET #2
	##############################
	# Switch out zero-inflated gamma or lognormal, depending on which one we're doing.
	if (resp_distrib_facet_2 == 'hGamma') {
	
		response_distrib_code_facet_2 <- nimbleCode({

			# FACET: parameters of facet distribution are latent and functions of environment
			# individual plant facet values are samples from the site-level distribution defined by the site-level distribution (next chunk after this one)
			for (i in 1:n_pheno_sites) {

				# moment matching to get gamma() parameters
				shape_facet_2[i] <- mu_facet_2_site[i]^2 / sigma_facet_2_within_sites^2
				rate_facet_2[i] <- mu_facet_2_site[i] / sigma_facet_2_within_sites^2

			}

			# FACET: likelihood of individual plants
			for (i in 1:n_plants) {

				# likelihood
				# y_facet_2[i] ~ dZIG(shape = shape_facet_2[site_index_facet[i]], rate = rate_facet_2[site_index_facet[i]], z = z_site[site_index_facet[i]])
				y_facet_2[i] ~ dZIG(shape = shape_facet_2[site_index_facet[i]], rate = rate_facet_2[site_index_facet[i]], psi = psi_site[site_index_facet[i]])

				# simulated values for unconditional DHARMa residuals
				# y_facet_2_sim[i] ~ dZIG(shape = shape_facet_2[site_index_facet[i]], rate = rate_facet_2[site_index_facet[i]], z = z_site[site_index_facet[i]])
				y_facet_2_sim[i] ~ dZIG(shape = shape_facet_2[site_index_facet[i]], rate = rate_facet_2[site_index_facet[i]], psi = psi_site[site_index_facet[i]])

				# log_lik_facet_2_y[i] <- dZIG(y_facet_2[i], shape = shape_facet_2[site_index_facet[i]], rate = rate_facet_2[site_index_facet[i]], z = z_site[site_index_facet[i]])
				log_lik_facet_2_y[i] <- dZIG(y_facet_2[i], shape = shape_facet_2[site_index_facet[i]], rate = rate_facet_2[site_index_facet[i]], psi = psi_site[site_index_facet[i]])
		
			}

		})
	
	} else if (resp_distrib_facet_2 == 'hurdleLN') {
	
		response_distrib_code_facet_2 <- nimbleCode({

			# FACET: likelihood of individual plants
			for (i in 1:n_plants) {

				# likelihood
				# y_facet_2[i] ~ dZILN(meanlog = Phi_site[site_index_facet[i], 2], sdlog = sigma_facet_2_within_sites, z = z_site[site_index_facet[i]])
				y_facet_2[i] ~ dZILN(meanlog = Phi_site[site_index_facet[i], 2], sdlog = sigma_facet_2_within_sites, psi = psi_site[site_index_facet[i]])

				# simulated values for unconditional DHARMa residuals
				# y_facet_2_sim[i] ~ dZILN(meanlog = Phi_site[site_index_facet[i], 2], sdlog = sigma_facet_2_within_sites, z = z_site[site_index_facet[i]])
				y_facet_2_sim[i] ~ dZILN(meanlog = Phi_site[site_index_facet[i], 2], sdlog = sigma_facet_2_within_sites, psi = psi_site[site_index_facet[i]])

				# log_lik_facet_2_y[i] <- dZILN(y_facet_2[i], meanlog = Phi_site[site_index_facet[i], 2], sdlog = sigma_facet_2_within_sites, z = z_site[site_index_facet[i]])
				log_lik_facet_2_y[i] <- dZILN(y_facet_2[i], meanlog = Phi_site[site_index_facet[i], 2], sdlog = sigma_facet_2_within_sites, psi = psi_site[site_index_facet[i]])

			}

		})
		
	}

	model_code <- glueNimbleCode(
		model_code,
		response_distrib_code_facet_1,
		response_distrib_code_facet_2,
		model_code_beta_psi_priors
	)
	
	print(model_code)

	say('nimbleModel():', level = 2)
	model <- nimbleModel(
		code = model_code, # our model
		constants = constants, # constants
		data = data, # data
		inits = inits, # initialization values
		check = TRUE, # any errors?
		calculate = FALSE,
		buildDerivs = FALSE # need for Hamiltonian Monte Carlo
	)

	say('$initializeInfo() and $calculate():', level = 2)
	model$initializeInfo()
	calc <- model$calculate()
	say('check_nodes()', level = 2)
	check_nodes(model)
	say('model$calculate(): ', calc)
	if (is.na(calc) | is.infinite(calc)) stop('Likelihood is incalculable.')

	say('configureMCMC():', level = 2)

	# monitors for coefficients that have no indexing
	monitors_coeffs_not_indexed <- c(
		'sigma_facet_1_within_sites', 'sigma_facet_2_within_sites', 'eta'
	)

	# coefficients that have bracketed indexing
	monitors_coeffs_single_index <- c(
		'beta_facet_1', 'beta_facet_2',
		'sigmas',
		'beta_psi'
	)

	monitors_coeffs_double_index <- c(
		'U_star'
	)

	monitors_derived_not_indexed <- c(
		'log_lik', 'log_lik_facet_1', 'log_lik_facet_2'
	)

	monitors_derived_single_index <- c(
		'mu_facet_1_site', 'mu_facet_2_site'
	)
	monitors_derived_double_index <- c(
		'correlation'
	)

	monitors_dharma <- c(
		'y_facet_1_sim', 'y_facet_2_sim'
	)

	monitors_resp_curves <- c()

	monitors <- c(monitors_coeffs_not_indexed, monitors_coeffs_single_index, monitors_coeffs_double_index, monitors_derived_not_indexed, monitors_derived_single_index, monitors_derived_double_index, monitors_dharma, monitors_resp_curves)

	conf <- configureMCMC(
		model,
		monitors = monitors,
		print = TRUE,
		enableWAIC = TRUE
	)

	vars <- 'beta_facet_1'
	conf$removeSamplers(vars)
	conf$addSampler(target = vars, type = 'AF_slice')
	say('AF slice sampler added to ', paste(vars, collapse = ' & '), '.')

	vars <- 'beta_facet_2'
	conf$removeSamplers(vars)
	conf$addSampler(target = vars, type = 'AF_slice')
	say('AF slice sampler added to ', paste(vars, collapse = ' & '), '.')

	vars <- 'beta_psi'
	conf$removeSamplers(vars)
	conf$addSampler(target = vars, type = 'AF_slice')
	say('AF slice sampler added to ', paste(vars, collapse = ' & '), '.')

	### compile/build/run model/save MCMC
	build <- buildMCMC(conf)

	compiled <- compileNimble(model, build, showCompilerOutput = FALSE)

	say('runMCMC() ', date())
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
	say('runMCMC() ', date())

	saveRDS(chains, paste0(out_dir, '/chains.rds'))

	say('PRIORS', level = 1)

	say('constants_shared_facet', level = 2)
	print(constants_shared_facet)

	say('constants_shared_psi', level = 2)
	print(constants_shared_psi)

	say('session info', level = 1)
	print(sessionInfo())

	say(date(), pre = 1)
	sink()

######################################################
### post-modeling analysis for FACET-only models ###
######################################################

	### post-modeling analysis of FACET
	descrip1 <- paste0(facet_1, ' ~ ', resp_distrib_facet_1, '(', ifelse(resp_distrib_facet_1 == 'hGamma', 'exp(', ''), 'MVN(env))', ifelse(resp_distrib_facet_1 == 'hGamma', ')', ''))
	descrip2 <- paste0(facet_2, ' ~ ', resp_distrib_facet_2, '(', ifelse(resp_distrib_facet_2 == 'hGamma', 'exp(', ''), 'MVN(env))', ifelse(resp_distrib_facet_2 == 'hGamma', ')', ''))
	descrip <- paste0(descrip1, ' + ', descrip2)
	facets <- paste0(facet_1, ' + ', facet_2)
	workflow_postmodeling_generic(facet = facets, formulae = formulae, descrip = descrip, out_dir = out_dir)

	# workflow_postmodeling_nonbiomass_single_facet(facet = facet, chains = chains, descrip = descrip, formula_facet = formula_facet, formula_psi = formula_psi, resp_distrib = resp_distrib, transform = transform, crossvalidate = crossvalidate, out_dir = out_dir)

say(date())
say('FINIS!', deco = '+', level = 1)

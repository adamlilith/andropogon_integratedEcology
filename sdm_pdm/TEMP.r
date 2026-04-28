# source('C:/Kaji/R/andropogon_integratedEcology/sdm_pdm/TEMP.r')

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

		# # response curves (general)
		# n_response_curve_values = n_response_curve_values, # number of values in response curve array

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
	# file <- paste0('./outputs_loretta/integrated_sdm_pdm/models_occurrence/[occs_zip[psi~', psi_filename, ']~normal~', occs_filename, '_[bias~1]]/chains.rds')
file <- paste0('C:/Kaji/Research/Andropogon/Andropogon/outputs_loretta/integrated_sdm_pdm/models_occurrence/[occs~hurdle_bio1^2_log(bio12)^2_bio15^2]_[bias~1]_niter_160000/chains.rds')
say('Using stopgap chains file for OCCURRENCES!!! Need to change for real runs!!!')
	chains_inits <- readRDS(file)
	
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
	file <- paste0('./outputs_loretta/integrated_sdm_pdm/models_biomass/[biomass~hurdleln_log(bio12)]/chains.rds')
say('Using stopgap chains file for BIOMASS!!! Need to change for real runs!!!')
	chains_inits <- readRDS(file)
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
			log(lambda[i]) <- exp(Phi_county[i, 1])

			# probability of zero inflation
			logit(psi_county[i]) <- inprod(beta_psi[1:n_terms_psi], x_by_county_psi[i, 1:n_terms_psi])

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

			### BIOMASS
			biomass_site[i] <- exp(Phi_site[i, 2])

		}

		# BIOMASS: likelihood of individual plants
		for (i in 1:n_biomass) {

			# likelihood
			y_biomass[i] ~ dHLN(meanlog = Phi_site[site_index_biomass[i], 2], sdlog = sigma_biomass_within_sites, psi = psi_site[site_index_biomass[i]])

			# simulated values for unconditional DHARMa residuals
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
			# y_facet_XYZ[i] ~ dHLN(meanlog = Phi_site[site_index_facet[i], index_XYZ], sdlog = sigma_facet_within_sites_XYZ, z = z_site[site_index_facet[i]])
			y_facet_XYZ[i] ~ dHLN(meanlog = Phi_site[site_index_facet[i], index_XYZ], sdlog = sigma_facet_within_sites_XYZ, psi = psi_site[site_index_facet[i]])

			# simulated values for unconditional DHARMa residuals... add 2 
			# y_facet_sim_XYZ[i] ~ dHLN(meanlog = Phi_site[site_index_facet[i], index_XYZ], sdlog = sigma_facet_within_sites_XYZ, z = z_site[site_index_facet[i]])
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

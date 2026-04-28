## MODELING ANDROPOGON GERARDI DISTRIBUTION, PHENOTYPE, PHYSIOLOGY, GENOTYPE, and ASSOCIATED MICROBIAL COMMUNITIES
### Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2023-12
###
### This model estimates the occurrence (present/absent), relative abundance, and site-level mean values of one or more traits of Andropogon gerardi. Traits include aboveground individual biomass, height, canopy diameter, and others. All of these "facets" integrate through a shared, non-centered multivariate normal distribution that accounts for correlations among facets. Mean values (or transforms of the means) are drawn from the MVN, then used as parameters for facet-specific distributions (e.g., zero-inflated Poisson for abundance, zero-inflated lognormal for biomass, etc.). It assumes relative abundance follows a hurdle (zero-inflated) Poisson distribution, which is a log-linear function of environmental covariates and an offset (number of Poaceae specimens that are not A. gerardi). To obviate issues with taking the log of zero, 1 is added to the offset. However, this can create an upward bias in coefficients since some AG can be found at sites with no other species found (implying zero search effort, so the abundance of AG must be due to highly suitable climate). To correct for this, a binary dummy variable is used for cases where AG is >0 and number of other Poaceae is 0 (dummy variable = 1 in these cases). The interpretation is then that a site with no other Poaceae has an abundance of AG that is equivalent to exp(bias_coeff) non-AG plants. This correction can be paired with a tightly-regularized prior for the intercept. The probability of (inflated) zero is a function of environmental covariates. teh occurrence/abundance, biomass, and trait models share the same hurdle component. Ergo, if the species is predicted to be absent in a location, all of its facet values will be forced to 0, as well.
###
### source('C:/Kaji/R/andropogon_integratedEcology/sdm_pdm/sdm_pdm_08d_model_fully_integrated_occ~hurdle_offset_biomass~hurdle_traits~hurdle.r')
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

	chain <- 1 # also used as seed

	trial <- TRUE # TRUE for testing
	# trial <- FALSE # TRUE for testing

	# calib <- TRUE # use just counties with non-NA Poaceae for calibration region
	calib <- FALSE # use all of North America for calibration region

	# do cross-validation?
	crossvalidate <- TRUE
	# crossvalidate <- FALSE

	if (!trial) {

		### MCMC settings
		# ~1 hr to do 8000 iterations with 2 non-biomass facets
		niter <- 160000
		# nchains <- 1
		nchains <- 4

	} else {
		
		### MCMC settings FOR TESTING
		niter <- 2000
		nchains <- 2

	}
	nburnin <- niter / 2
	thin <- (niter - nburnin) / 1000

	### formula for how aspects of species responds to environment

	formula_occs <- ~ 1 + bio1 + bio12_log10p1 + bio15 + I(bio1^2) + I(bio12_log10p1^2) + I(bio15^2) # response of occurrence to climate and soil
	occs_filename <- 'bio1^2_log(bio12)^2_bio15^2'
	file_occs <- paste0('C:/Kaji/Research/Andropogon/Andropogon/outputs_loretta/integrated_sdm_pdm/models_occurrence/[occs~hurdle_', occs_filename, ']_[bias~poaceae_offset]/chains.rds')

	formula_psi <- formula_occs
	psi_filename <- occs_filename

	formula_biomass <- ~ 1 + bio12_log10p1 # response of biomass to environment
	biomass_filename <- 'log(bio12)'
	resp_distrib_biomass <- 'hurdleLN'
	transform_biomass <- if (resp_distrib_biomass == 'hGamma') { 'exponential' } else if (resp_distrib_biomass == 'hurdleLN') { 'identity' }
	file_biomass <- paste0('./outputs_loretta/integrated_sdm_pdm/models_biomass/[biomass~', tolower(resp_distrib_biomass), '_', biomass_filename, ']/chains.rds')

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
		)
		# cn_ratio = list(
		# 	formula = ~ 1 + bio12,
		# 	filename = 'bio12',
		# 	resp_distrib = 'hurdleLN',
		# 	transform = 'identity'
		# ),
		# height = list(
		# 	formula = ~ 1 + bio1 + bio12 + bio1:bio12,
		# 	filename = 'bio1_x_bio12',
		# 	resp_distrib = 'hurdleLN',
		# 	transform = 'identity'
		# )#,
		# # # internal_co2 = list(
		# 	# formula = ~ 1 + bio12 + nitrogen + bio12:nitrogen,
		# 	# filename = 'bio12_x_nitrogen',
		# 	# resp_distrib = 'hurdleLN',
		# 	# transform = 'identity'
		# # ),
		# # leaf_thickness = list(
		# # 	formula = ~ 1 + bio12,
		# # 	filename = 'bio12',
		# # 	resp_distrib = 'hurdleLN',
		# # 	transform = 'identity'
		# # ),
		# # n_concentration = list(
		# # 	formula = ~ 1 + bio12,
		# # 	filename = 'bio12',
		# # 	resp_distrib = 'hurdleLN',
		# # 	transform = 'identity'
		# # ),
		# photosynthetic_rate = list(
		# 	formula = ~ 1 + bio1 + bio12 + bio1:bio12,
		# 	filename = 'bio1_x_bio12',
		# 	resp_distrib = 'hurdleLN',
		# 	transform = 'identity'
		# ),
		# spad = list(
		# 	formula = ~ 1 + bio12 + insolation_2000_growing_season_kWh_per_m2,
		# 	filename = 'bio12_srad',
		# 	resp_distrib = 'hurdleLN',
		# 	transform = 'identity'
		# ),
		# stomatal_conductance = list(
		# 	formula = ~ 1 + bio1 + bio12_log10p1 + I(bio1^2),
		# 	filename = 'bio1^2_log(bio12)',
		# 	resp_distrib = 'hurdleLN',
		# 	transform = 'identity'
		# ),
		# transpiration_rate = list(
		# 	formula = ~ 1 + bio1 + bio12 + I(bio12^2),
		# 	filename = 'bio1_bio12^2',
		# 	resp_distrib = 'hurdleLN',
		# 	transform = 'identity'
		# )
	)

	nonbiomass_facets <- nonbiomass_facets[sort(names(nonbiomass_facets))]

	nonbiomass_facet_names <- names(nonbiomass_facets)
	nonbiomass_facet_names_short <- gsub(nonbiomass_facet_names, pattern = '_', replacement = '')
	nonbiomass_facet_names_short <- paste(substr(nonbiomass_facet_names_short, 1, 3), collapse = '_')

	for (f in seq_along(nonbiomass_facets)) {

		facet <- names(nonbiomass_facets)[f]
		form <- nonbiomass_facets[[facet]]$formula
		terms <- attr(terms(form), 'term.labels')
		nonbiomass_facets[[facet]]$n_terms <- length(terms) + 1

	}

	### output folder and bias formula
	out_dir <- paste0('./outputs_loretta/integrated_sdm_pdm/models_integrated/', ifelse(trial, 'TRIAL_', ''), 'non_centered_MVN([occs~hurdlepoisson_offset]_[biomass~hurdleln]_[', nonbiomass_facet_names_short, '])_chain_', chain)

#############
### model ###
#############

	if (!trial) if (file.exists(out_dir)) stop('Output folder already exists.')
	dirCreate(out_dir)

	sink(paste0(out_dir, '/runtime_log.txt'), split = TRUE)
	say('MODELING OCCURRENCE, BIOMASS, and NON-BIOMASS FACETS')
	say('Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org |  ', date(), post = 1)

	### data collation
	##################
	say('data collation', post = 1)

	say('This model estimates the occurrence (present/absent), relative abundance, and site-level mean values of two traits of Andropogon gerardi. Traits include aboveground individual biomass, height, canopy diameter, and others. All of these "facets" integrate through a shared, non-centered multivariate normal distribution that accounts for correlations among facets. Mean values (or transforms of the means) are drawn from the MVN, then used as parameters for facet-specific distributions (e.g., zero-inflated Poisson for abundance, zero-inflated lognormal for biomass, etc.). It assumes relative abundance follows a hurdle (zero-inflated) Poisson distribution, which is a log-linear function of environmental covariates and an offset (number of Poaceae specimens that are not A. gerardi). To obviate issues with taking the log of zero, 1 is added to the offset. However, this can create an upward bias in coefficients since some AG can be found at sites with no other species found (implying zero search effort, so the abundance of AG must be due to highly suitable climate). To correct for this, a binary dummy variable is used for cases where AG is >0 and number of other Poaceae is 0 (dummy variable = 1 in these cases). The interpretation is then that a site with no other Poaceae has an abundance of AG that is equivalent to exp(bias_coeff) non-AG plants. This correction can be paired with a tightly-regularized prior for the intercept. The probability of (inflated) zero is a function of environmental covariates. teh occurrence/abundance, biomass, and trait models share the same hurdle component. Ergo, if the species is predicted to be absent in a location, all of its facet values will be forced to 0, as well. This script models TWO non-biomass facets.', breaks = 60, post = 1)

	say('MCMC settings:', level = 2)
	say('trial ........................ ', trial)
	say('niter ........................ ', niter)
	say('nburnin ...................... ', nburnin)
	say('thin ......................... ', thin)
	say('nchains ...................... ', nchains)
	say('formula_occs ................. ', paste(as.character(formula_occs), collapse = ' '))
	say('formula_psi .................. ', paste(as.character(formula_psi), collapse = ' '))
	say('bias_offset .................. poaceae')
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
		formula_psi = formula_psi,
		formula_biomass = formula_biomass,
		nonbiomass_facets = nonbiomass_facets
	)
	saveRDS(formulae, paste0(out_dir, '/formulae.rds'))

	# table indicating which index belongs to which facet (abundance = 1, biomass = 2, traits = 3, 4, ...)
	facet_table <- data.table(
		index = 1:(2 + length(nonbiomass_facets)),
		facet = c('occurrence', 'biomass', names(nonbiomass_facets))
	)

	fwrite(facet_table, paste0(out_dir, '/!facet_codes.csv'))

	########################
	### data preparation ###
	########################

	say('preparing data for occurrences ', date(), level = 2)

	# data for OCCURRENCES at counties
	data_occs_counties <- prepare_occurrence_data(formula_occs = formula_occs, formula_bias = ~ 1, n_response_curve_values = n_response_curve_values, psa_quant = psa_quant, calib = calib)

	# data for OCCURRENCES using site-level environment
	data_occs_sites <- prepare_biomass_data(formula_biomass = formula_occs, n_response_curve_values = n_response_curve_values, calib = calib)

	say('preparing data for psi ', date(), level = 2)

	# data for ZERO-INFLATION at counties
	data_psi_counties <- prepare_occurrence_data(formula_occs = formula_psi, formula_bias = ~ 1, n_response_curve_values = n_response_curve_values, psa_quant = psa_quant, calib = calib)

	# data for psi at sites
	data_psi_sites <- prepare_biomass_data(formula_biomass = formula_psi, n_response_curve_values = n_response_curve_values, calib = calib)

	say('preparing data for biomass ', date(), level = 2)
	# data for BIOMASS at sites
	data_biomass_sites <- prepare_biomass_data(formula_biomass = formula_biomass, n_response_curve_values = n_response_curve_values, calib = calib)

	# data for BIOMASS using county-level environment
	data_biomass_counties <- prepare_occurrence_data(formula_occs = formula_biomass, formula_bias = ~ 1, n_response_curve_values = n_response_curve_values, psa_quant = psa_quant, calib = calib)

	# data for NON-BIOMASS FACETS at sites and counties
	data_nonbiomass_sites <- data_nonbiomass_counties <- list()
	for (f in seq_along(nonbiomass_facets)) {

		facet <- names(nonbiomass_facets)[f]
		say('preparing data for ', facet, ' ', date(), level = 2)
		
		formula_facet <- nonbiomass_facets[[facet]]$formula

		data_nonbiomass_sites[[f]] <- prepare_nonbiomass_data(facet = facet, formula_facet = formula_facet, n_response_curve_values = n_response_curve_values, calib = calib)

		data_nonbiomass_counties[[f]] <- prepare_occurrence_data(formula_occs = formula_facet, formula_bias = ~ 1, n_response_curve_values = n_response_curve_values, psa_quant = psa_quant, calib = calib)

	}

	names(data_nonbiomass_sites) <- names(data_nonbiomass_counties) <- names(nonbiomass_facets)

	#########################
	### inputs for nimble ###
	#########################

	say('inputs ', date(), level = 2)
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

	raw_bias_offset <- ceiling(data_occs_counties$ag_vect_sq$n_poaceae) - data_occs_counties$ag_vect_sq$n_andropogon_gerardi
	# log_bias_offset <- log1p(bias_offset)
	bias_offset <- log(raw_bias_offset + exp(-1))
	zero_bias_correction <- as.numeric(raw_bias_offset == 0)

	constants <- list(
		
		### integration
		n_facets = n_facets,

		### occurrences
		n_counties_occs_calib = data_occs_counties$n_counties_occs_calib, # number of counties in calibration region
		n_terms_occs = data_occs_counties$n_terms, # number of terms in formula for occurrence model (including intercept)
		
		bias_offset = bias_offset,
		zero_bias_correction = zero_bias_correction,

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
	alpha_occs_inits <- 0
	beta_occs_mu_inits <- rep(0, constants$n_terms_occs)
	beta_psi_inits <- rep(0, constants$n_terms_psi)
	
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
		beta_occs = beta_occs_mu_inits, # occurrence ~ environment coefficients (including intercept)
		alpha_occs = alpha_occs_inits, # intercept, area, # of Poaceae
		# y_n_ag_sim = data_occs_counties$y_n_ag, # simulated values of observed number of AG (for DHARMa residuals)
	
		### biomass
		beta_biomass = beta_biomass_inits,
		sigma_biomass_within_sites_log = log_sigma_biomass_within_sites_inits,
		# y_biomass_sim = data_biomass_sites$y_biomass, # simulated values for biomass (for DHARMa residuals)
	
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

		# y_facet_sim <- data_nonbiomass_sites[[facet]]$y_facet # simulated values for facet (for DHARMa residuals)

		inits_facet <- list(
			beta_facet = beta_facet_inits,
			log_sigma_facet_within_sites = log(sigma_facet_within_sites_init)#,
			# y_facet_sim = y_facet_sim
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

	stopifnot(
		dim(inits$U_star)[1] == constants$n_facets,
		dim(inits$U_star)[2] == constants$n_facets
	)

	stopifnot(
		ncol(constants$x_by_county_occs) == constants$n_terms_occs
	)

	stopifnot(
		length(inits$log_sigmas) == constants$n_facets
	)

	### define model
	say('nimbleCode():', level = 2)

	## for faster debugging, un-oblige data collation and load from previous run

	if (trial) {
		saveRDS(data, 'C:/!scratch/ag_fully_integrated_2_facets_data.rds')
		saveRDS(constants, 'C:/!scratch/ag_fully_integrated_2_facets_constants.rds')
		saveRDS(inits, 'C:/!scratch/ag_fully_integrated_2_facets_inits.rds')
	}

	if (FALSE) {
		data <- readRDS('C:/!scratch/ag_fully_integrated_2_facets_data.rds')
		constants <- readRDS('C:/!scratch/ag_fully_integrated_2_facets_constants.rds')
		inits <- readRDS('C:/!scratch/ag_fully_integrated_2_facets_inits.rds')
	}

	model_code_base <- nimbleCode({

		# INTEGRATION: LJK prior for standard deviations and correlations between latent occurrence and biomass processes
		eta ~ dgamma(2, 1)
		U_star[1:n_facets, 1:n_facets] ~ dlkj_corr_cholesky(eta = eta, p = n_facets)
	  
		# INTEGRATION: standard deviations of latent occurrence and biomass
		log(sigmas[1]) ~ dnorm(0, sd = lambda_sigma_prior_sd) # lognormal prior
		log(sigmas[2]) ~ dnorm(0, sd = sigma_biomass_among_sites_log_prior_sd) # lognormal prior

		for (i in 3:n_facets) {
			log(sigmas[i]) ~ dnorm(0, sd = sigma_facet_among_sites_log_prior_sd) # lognormal prior
		}

		# INTEGRATION: county-level
		# NB We have measurements of biomass at the site level, but not county. We thus assume that the estimates of biomass using county-level covariates are indicative of virtual sites that have the same environments as counties. Biomass is still estimated at the site level. Same for the non-biomass trait facets.
		
		# relationship between abundance, biomass, and facets at COUNTY level with environment
		phi_occs_county[1:n_counties_occs_calib] <- (x_by_county_occs[1:n_counties_occs_calib, 1:n_terms_occs] %*% beta_occs[1:n_terms_occs])[ , 1]
		phi_biomass_county[1:n_counties_occs_calib] <- (x_by_county_biomass[1:n_counties_occs_calib, 1:n_terms_biomass] %*% beta_biomass[1:n_terms_biomass])[ , 1]
		# phi_facet_county_XYZ[1:n_counties_occs_calib] <- ... # included for each non-biomass facet in generic code blocks below

		# # matrix of standard normals... including in loop be low row-by-row to enhance mixing
		# z_county[1:n_counties_occs_calib, 1:n_facets] ~ dIIDStandardNorm(n_row = n_counties_occs_calib, n_col = n_facets)

		# non-centered MVN method to avoid sampling directly from MVN
		# stands in for: Phi_county[i, 1:n_facets] ~ dmnorm(phis_county[i, 1:n_facets], cholesky = U[1:n_facets, 1:n_facets], prec_param = 0)
		Phi_county[1:n_counties_occs_calib, 1:n_facets] <- 
			phis_county[1:n_counties_occs_calib, 1:n_facets] + 
			(z_county[1:n_counties_occs_calib, 1:n_facets] %*% U_star[1:n_facets, 1:n_facets]) %*% diag(sigmas[1:n_facets])

		for (i in 1:n_counties_occs_calib) {

			phis_county[i, 1] <- phi_occs_county[i]
			phis_county[i, 2] <- phi_biomass_county[i]
			# phis_county[i, >=3] <- phi_facet_county_XYZ[i] # defined in generic code blocks below for non-biomass facets

			# row of matrix of standard normals
			# z_county[i, 1:n_facets] ~ dmnorm(mean = zeroes[1:n_facets], cholesky = ones_diag[1:n_facets, 1:n_facets], prec_param = 1)
			z_county[i, 1:n_facets] ~ dIIDStandardNorm_row(n_col = n_facets)

			# use offset method to avoid sampling directly from MVN
			# stands in for: Phi_county[i, 1:n_facets] ~ dmnorm(phis_county[i, 1:n_facets], cholesky = U[1:n_facets, 1:n_facets], prec_param = 0)
			# latent non-centered MNV at county level
			# for (j in 1:n_facets) {

			# 	# latent non-centered MNV at county level
			# 	Phi_county[i, j] <- phis_county[i, j] + sigmas[j] * inprod(U_star[1:n_facets, j], z_county[i, 1:n_facets])
			
			# }

		}

		# INTEGRATION: site-level calculations
		# Abundance is naturally estimated at the county level, so we assume that the site-level environmental covariates act like faux counties for abundance, but are indicative of the sites for phenotypic variables.

		# relationship between abundance, biomass, and facets at SITE level with environment
		phi_occs_site[1:n_pheno_sites] <- (x_by_site_occs[1:n_pheno_sites, 1:n_terms_occs] %*% beta_occs[1:n_terms_occs])[ , 1]
		phi_biomass_site[1:n_pheno_sites] <- (x_by_site_biomass[1:n_pheno_sites, 1:n_terms_biomass] %*% beta_biomass[1:n_terms_biomass])[ , 1]
		# phi_facet_site_XYZ[1:n_pheno_sites] <- ... # included for each non-biomass facet in generic code blocks below

		# matrix of standard normals
		z_site[1:n_pheno_sites, 1:n_facets] ~ dIIDStandardNorm(n_row = n_pheno_sites, n_col = n_facets)

		for (i in 1:n_pheno_sites) {

			phis_site[i, 1] <- phi_occs_site[i]
			phis_site[i, 2] <- phi_biomass_site[i]
			# phis_site[i, >=3] <- phi_facet_site_XYZ[i] # defined in generic code blocks below for non-biomass facets

			# use offset method to avoid sampling from non-centered distribution MVN
			# stands in for: Phi_site[i, 1:n_facets] ~ dmnorm(phis_site[i, 1:n_facets], cholesky = U[1:n_facets, 1:n_facets], prec_param = 0)
			# z_site[i, 1:n_facets] ~ dmnorm(mean = zeroes[1:n_facets], cholesky = ones_diag[1:n_facets, 1:n_facets], prec_param = 1)
			for (j in 1:n_facets) {

				# latent non-centered MNV at site level
				Phi_site[i, j] <- phis_site[i, j] + sigmas[j] * inprod(U_star[1:n_facets, j], z_site[i, 1:n_facets])
				
			}

		}

		### OCCURRENCE
		log(lambda_star[1:n_counties_occs_calib]) <-
			Phi_county[1:n_counties_occs_calib, 1] +
			bias_offset[1:n_counties_occs_calib] +
			alpha_occs * zero_bias_correction[1:n_counties_occs_calib]
			
		# zero component
		psi_county_star[1:n_counties_occs_calib] <-
			(x_by_county_psi[1:n_counties_occs_calib, 1:n_terms_psi] %*% beta_psi[1:n_terms_psi])[ , 1]
		logit(psi_county[1:n_counties_occs_calib]) <- psi_county_star[1:n_counties_occs_calib]

		# OCCURRENCE: likelihood
		for (i in 1:n_counties_occs_calib) {
			
			### observed number of AG and sampling bias
			y_n_ag[i] ~ dHurdlePoisson(lambda = lambda_star[i], psi = psi_county[i])

			# # simulate observations for DHARMa residuals
			# y_n_ag_sim[i] ~ dHurdlePoisson(lambda = lambda_star[i], psi = psi_county[i])

		}

		### BIOMASS
		# prior for sd of individual plant biomass on lognormal, ==> vague
		sigma_biomass_within_sites_log ~ dnorm(0, sd = sigma_biomass_within_sites_log_prior_sd)
		sigma_biomass_within_sites <- exp(sigma_biomass_within_sites_log)

		### PROBABILITY OF ZERO INFLATION at site level
		psi_site_star[1:n_pheno_sites] <-
			(x_by_site_psi[1:n_pheno_sites, 1:n_terms_psi] %*% beta_psi[1:n_terms_psi])[ , 1]
		logit(psi_site[1:n_pheno_sites]) <- psi_site_star[1:n_pheno_sites]

		### BIOMASS: likelihood of individual plants
		for (i in 1:n_biomass) {

			# likelihood
			y_biomass[i] ~ dHLN(meanlog = Phi_site[site_index_biomass[i], 2], sdlog = sigma_biomass_within_sites, psi = psi_site[site_index_biomass[i]])

			# # simulated values for unconditional DHARMa residuals
			# y_biomass_sim[i] ~ dHLN(meanlog = Phi_site[site_index_biomass[i], 2], sdlog = sigma_biomass_within_sites, psi = psi_site[site_index_biomass[i]])

		}

	})

	### model code for NON-BIOMASS FACETS
	#####################################

	# Written generically so we can model code each facet in a loop, but the code is the same for each facet except for the index of the facet in the Phi_site matrix and the number of terms in the formula for that facet. We add 2 to the index of the facet in the Phi_site matrix to skip the occurrence and biomass terms.

	code_facets_raw	 <- '{

		# NON-BIOMASS FACET XYZ: prior on coefficients
		beta_facet_XYZ[1] ~ dnorm(0, sd = beta_facet_prior_dnorm_sd_1)
		for (i in 2:n_terms_facet_XYZ) {
			beta_facet_XYZ[i] ~ dnorm(0, sd = beta_facet_prior_dnorm_sd) # broad prior
		}

		# NON-BIOMASS FACET XYZ: relationship between abundance, biomass, and facets at COUNTY level with environment
		phis_county[1:n_counties_occs_calib, XYZ_plus_2] <-
			(x_by_county_facet_XYZ[1:n_counties_occs_calib, 1:n_terms_facet_XYZ] %*% beta_facet_XYZ[1:n_terms_facet_XYZ])[ , 1]

		# NON-BIOMASS FACET XYZ: relationship between abundance, biomass, and facets at SITE level with environment
		phis_site[1:n_pheno_sites, XYZ_plus_2] <-
			(x_by_site_facet_XYZ[1:n_pheno_sites, 1:n_terms_facet_XYZ] %*% beta_facet_XYZ[1:n_terms_facet_XYZ])[ , 1]

		# NON-BIOMASS FACET XYZ: prior for sd of individual plant facet value on lognormal, ==> vague
		log(sigma_facet_within_sites_XYZ) ~ dnorm(0, sd = sigma_facet_within_sites_log_prior_sd)

		# NON-BIOMASS FACET XYZ: likelihood of individual plant non-biomass trait
		for (i in 1:n_non_biomass) {

			### FACET XYZ
			y_facet_XYZ[i] ~ dHLN(meanlog = Phi_site[site_index_facet[i], XYZ_plus_2], sdlog = sigma_facet_within_sites_XYZ, psi = psi_site[site_index_facet[i]])

			# # simulated values for unconditional DHARMa residuals
			# y_facet_sim_XYZ[i] ~ dHLN(meanlog = Phi_site[site_index_facet[i], XYZ_plus_2], sdlog = sigma_facet_within_sites_XYZ, psi = psi_site[site_index_facet[i]])

		}

	}'

	model_code_non_biomass_facets <- list()
	for (i in seq_along(nonbiomass_facets)) {

		code_facets_raw_this_facet <- gsub(pattern = 'XYZ_plus_2', replacement = i + 2, x = code_facets_raw)
		code_facets_raw_this_facet <- gsub(pattern = 'XYZ', replacement = i, x = code_facets_raw_this_facet)
		model_code_non_biomass_facets[[i]] <- code_facets_raw_this_facet

	}
	model_code_non_biomass_facets <- lapply(model_code_non_biomass_facets, function(x) parse(text = x)[[1]])

	model_code <- glueNimbleCode(
		model_code_base,
		model_code_beta_occs_priors,
		model_code_beta_psi_priors,
		model_code_alpha_occs_offset_correction_priors,
		model_code_beta_biomass_priors
	)

	model_code <- Reduce(
		f = function(code_a, code_b) glueNimbleCode(code_a, code_b),
		x = model_code_non_biomass_facets,
		init = model_code
	)

	print(model_code)

	say('nimbleModel()', level = 2)
	model <- nimbleModel(
		code = model_code,
		constants = constants,
		data = data,
		inits = inits,
		check = TRUE,
		calculate = FALSE,
		buildDerivs = TRUE # need for Hamiltonian Monte Carlo
		# buildDerivs = FALSE
	)

	say('initializeInfo() and $calculate() ', date(), level = 2)
	model$initializeInfo()
	calc <- model$calculate()
	check_nodes(model)
	say('model$calculate(): ', calc)
	if (is.na(calc) || is.infinite(calc)) stop('Impossible likelihood.')

	say('configureMCMC() ', date(), level = 2)

	monitors_coeffs_not_indexed <- c('alpha_occs', 'eta', 'sigma_biomass_within_sites')
	for (f in seq_along(nonbiomass_facets)) monitors_coeffs_not_indexed <- c(monitors_coeffs_not_indexed, paste0('sigma_facet_within_sites_', f))
	monitors_coeffs_single_index <- c('beta_occs', 'beta_biomass', 'beta_psi', 'sigmas')
	for (f in seq_along(nonbiomass_facets)) monitors_coeffs_single_index <- c(monitors_coeffs_single_index, paste0('beta_facet_', f))
	monitors_coeffs_double_index <- c('U_star')

	monitors_derived_not_indexed <- c()
	monitors_derived_single_index <- c()
	monitors_derived_double_index <- c()

	monitors_dharma <- c(
		# 'y_n_ag_sim', 'y_biomass_sim', 'y_facet_sim_1', 'y_facet_sim_2'
	)

	monitors_debug <- c(
		# 'lambda', 'psi', 'Phi_county', 'Phi_site'
	)

	monitors <- c(monitors_coeffs_not_indexed, monitors_coeffs_single_index, monitors_coeffs_double_index, monitors_derived_not_indexed, monitors_derived_single_index, monitors_derived_double_index, monitors_dharma, monitors_debug)

	conf <- configureMCMC(
		model,
		monitors = monitors,
		print = FALSE,
		enableWAIC = FALSE
	)

	# HMC samplers for continuous parameters, block sampling for correlated parameters
	hmc_targets <- c('sigmas', 'eta', 'alpha_occs', 'beta_occs', 'beta_psi', 'beta_biomass', 'sigma_biomass_within_sites')
	for (f in seq_along(nonbiomass_facets)) {
		hmc_targets <- c(hmc_targets, paste0('beta_facet_', f))
	}

	conf$removeSamplers(hmc_targets)

	# NUTS sampler — one joint block per group of correlated params
	addHMC(conf,
		target = hmc_targets,
		type = 'NUTS', # auto-tunes; better than 'HMC'
		control = list(
			warmupMode = 'default', # dual-averaging step-size adaptation
			M = 'diagonal', # mass matrix: 'diagonal' or 'dense'
			maxTreeDepth = 10 # safety cap on trajectory length
		)
	)

	# Scalar continuous params can share a second NUTS block or use 'NUTS_classic'
	addHMC(conf, target = scalar_continuous, type = 'NUTS')

	# vars <- c('beta_occs', 'beta_psi')
	# for (var in vars) conf$removeSamplers(var)
	# conf$addSampler(target = vars, type = 'AF_slice')
	# say('Added AF_slice sampler to ', paste(vars, collapse = ' & '), '.')

	# var <- 'beta_biomass'
	# conf$removeSamplers(var)
	# conf$addSampler(target = var, type = 'AF_slice')
	# say('AF_slice sampler added to ', paste(var, collapse = ' & '), '.')

	# for (f in seq_along(nonbiomass_facets)) {

	# 	facet <- names(nonbiomass_facets)[f]
	# 	var <- paste0('beta_facet_', f)

	# 	conf$removeSamplers(var)
	# 	conf$addSampler(target = var, type = 'AF_slice')
	# 	say('AF_slice sampler added to ', paste(var, collapse = ' & '), '.')

	# }

	print(conf)

	### compile/build/run model/save MCMC
	say('buildMCMC() ', date(), level = 2)
	build <- buildMCMC(conf)

	say('compileNimble() ', date(), level = 2)
	compiled <- compileNimble(model, build, showCompilerOutput = TRUE)

	say('runMCMC() ', date(), level = 2)
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
		formula_bias = ~ 1,
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
			formula_bias = ~ 1,
			formula_psi = formula_psi,
			formula_biomass = formula_biomass,
			nonbiomass_facets = nonbiomass_facets,
			out_dir = out_dir
		)
	
	}

say(date())
say('FINIS!', deco = '+', level = 1)

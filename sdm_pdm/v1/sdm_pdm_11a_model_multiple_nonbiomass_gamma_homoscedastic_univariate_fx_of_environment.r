### MODELING ANDROPOGON GERARDI DISTRIBUTION, PHENOTYPE, PHYSIOLOGY, GENOTYPE, and ASSOCIATED MICROBIAL COMMUNITIES
### Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2023-12
###
### This model estimates the distribution of TWO ORE MORE non-biomass traits site-level mean values. The distribution of trait values among ramets at a site is assumed to follow a gamma distribution with individual ramet values drawn from The gamma's parameters are provided by moment matching to a multivariate normal distribution, where the means of each trait's site-level values are functions of ONE environmental covariates (possibly with higher-order terms). Covariances in the multivariate normal can be non-zero. The model is "explicitly" homoscedastic because the standard deviation of each trait's distribution at a site is not is assumed to be constant across sites.
###
### NB Code in this script often uses "traits" (vs "trait") in variable names, but this is for ease in portability. This workflow only models one trait.
###
### source('C:/Kaji/R/andropogon_integratedEcology/sdm_pdm/sdm_pdm_11a_model_multiple_nonbiomass_gamma_homoscedastic_univariate_fx_of_environment.r')

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

	### MCMC settings
	#################

		# ### MCMC settings FOR NON-BIOMASS TRAITS
		# niter <- 20000
		# nburnin <- 10000
		# thin <- 20
		# nchains <- 4
		# waic <- TRUE

		# ### MCMC settings FOR TESTING
		niter <- 1100
		nburnin <- 100
		thin <- 2
		nchains <- 2
		waic <- TRUE

	### formulae for how aspects of species responds to environment
	formula_occs <- ~ 1 # not used except for prepare_occurrences()
	formula_occs_bias <- ~ 1 # not used except for prepare_occurrences()
	
	### traits
	traits <- c(
		'height',
		'canopy_diameter'
	)

	### formulae
	# linear forms and univariate quadratic forms
	formulae_traits <- list(

		height = ~ 1 + bio12,
		canopy_diameter = ~ 1 + bio12 + I(bio12^2)
	
	)

	### main modeling loop
	######################
	# loop over each trait
	# loop over each valid formula

	out_dir <- paste0('./outputs_loretta/integrated_sdm_pdm/models_multiple_traits/', paste(traits, collapse = '_'), ifelse(trial, '_TRIAL', ''), '/')
	if (!trial & file.exists(out_dir)) stop('Output folder already exists.')
	dirCreate(out_dir)
	
	formulae <- formulae_traits

	### start log
	#############
	sink(paste0(out_dir, '/runtime_log.txt'), split = TRUE)
	say('MODELING ', paste(traits, collapse = ' '))
	say('Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | ', date(), post = 1)

	say('data collation', post = 1)

	say('This model estimates the distribution of TWO ORE MORE non-biomass traits site-level mean values. The distribution of trait values among ramets at a site is assumed to follow a gamma distribution with individual ramet values drawn from The gamma\'s parameters are provided by moment matching to a multivariate normal distribution, where the means of each trait\'s site-level values are functions of ONE environmental covariates (possibly with higher-order terms). Covariances in the multivariate normal can be non-zero. The model is "explicitly" homoscedastic because the standard deviation of each trait\'s distribution at a site is not is assumed to be constant across sites.', breaks = 60, post = 1)

	say('Settings:', level = 2)
	say('traits ....................... ', paste(traits, collapse = ' '))
	say('niter ........................ ', niter)
	say('nburnin ...................... ', nburnin)
	say('thin ......................... ', thin)
	say('nchains ...................... ', nchains)
	say('trial ........................ ', trial, post = 2)

	say('formulae_traits:')
	print(formulae_traits)

	say('out_dir:')
	say(out_dir, post = 2)

	########################s
	### data preparation ###
	########################

	data_multitraits <- list()
	for (i in seq_along(traits)) {

		trait <- traits[i]

		data_multitraits[[length(data_multitraits) + 1]] <-
			prepare_nonbiomass_traits(trait = trait, formula = formulae_traits[[trait]], n_response_curve_values = n_response_curve_values, calib = FALSE)

		names(data_multitraits)[length(data_multitraits)] <- trait

	}

	data_occs <- prepare_occurrences(formula_occs = formula_occs, formula_occs_bias = formula_occs_bias, n_response_curve_values = n_response_curve_values, psa_quant = psa_quant, calib = FALSE)

	#########################
	### inputs for nimble ###
	#########################

	say('Inputs:', level = 2)
	data <- list(
		y_height = data_multitraits$height$y_trait,
		y_canopy_diameter = data_multitraits$canopy_diameter$y_trait
	)

	constants <- list(
		
		n_responses = 2, # number of traits

		### traits
		n_traits = length(traits), # number of traits to model
		n_trait_values = data_multitraits[[1]]$n_trait_values, # number of traits observations
		site_index_traits = data_multitraits[[1]]$site_index_traits, # index of sampled site for each row in traits data

		# counties
		n_counties = data_occs$n_counties, # number of counties in the dataset

		# sites sampled for phenotyping
		n_pheno_sites = data_multitrait[[1]]$n_pheno_sites, # number of phenotype sample sites

		# response curves (general)
		n_response_curve_values = n_response_curve_values # number of values in response curve array

	)

	for (trait in traits) {
	
		this_constants <- list(
		
			### traits
			x_by_site_traits = data_multitraits[[trait]]$x_by_site_traits, # MM with covariates for traits (scaled)
			n_terms_traits = data_multitraits[[trait]]$n_terms_traits, # number of terms in formula for traits model (including intercept)

			resp_curves_x = data_multitraits[[trait]]$resp_curves_x_traits, # response curve array for traits

			counties_x_trait_sq = data_multitraits[[trait]]$counties_x_traits_sq,
			counties_x_trait_ssp245_2041_2070 = data_multitraits[[trait]]$counties_x_traits_ssp245_2041_2070,
			counties_x_trait_ssp245_2071_2100 = data_multitraits[[trait]]$counties_x_traits_ssp245_2071_2100,
			counties_x_trait_ssp370_2041_2070 = data_multitraits[[trait]]$counties_x_traits_ssp370_2041_2070,
			counties_x_trait_ssp370_2071_2100 = data_multitraits[[trait]]$counties_x_traits_ssp370_2071_2100

		)
		names(this_constants) <- paste0(names(this_constants), '_', trait)
		constants <- c(constants, this_constants)
	
	}

	inits <- list()

	for (trait in traits) {

		y_mean <- mean(data_multitraits[[trait]]$y_traits[ , 1])
	
		this_inits <- list(

			y_traits_sim = data_multitraits[[trait]]$y_traits[ , 1], # simulated values for trait (for DHARMa residuals)

			log_site_traits_mu = rep(log(y_mean), data_multitraits[[trait]]$n_trait_values),

			traits_mu_county_sq = rep(y_mean, data_multitraits[[trait]]$n_counties),
			traits_mu_county_ssp245_2041_2070 = rep(y_mean, data_multitraits[[trait]]$n_counties),
			traits_mu_county_ssp245_2071_2100 = rep(y_mean, data_multitraits[[trait]]$n_counties),
			traits_mu_county_ssp370_2041_2070 = rep(y_mean, data_multitraits[[trait]]$n_counties),
			traits_mu_county_ssp370_2071_2100 = rep(y_mean, data_multitraits[[trait]]$n_counties),

			cross_site_traits_sigma_log = 1,
			phi_traits_sigma_log = 1,
			beta_traits_mu = rep(0, data_multitraits[[trait]]$n_terms_traits)

		)
		names(this_inits) <- paste0(names(this_inits), '_', trait)
		inits <- c(inits, this_inits)

	}

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
	
		# NON-BIOMASS TRAITS: priors for relationship of mean and variance to environment
		for (i in 1:n_terms_traits_height) {
			beta_traits_height[i] ~ dnorm(0, sd = 10) # broad prior
		}

		for (i in 1:n_terms_traits_canopy_diameter) {
			beta_traits_canopy_diameter[i] ~ dnorm(0, sd = 10) # broad prior
		}

		# prior for sd of site-level trait values on lognormal (~ half-Cauchy), ==> vague
		cross_site_traits_sigma_log_height ~ dnorm(0, sd = 2.5)
		cross_site_traits_sigma_log_canopy_diameter ~ dnorm(0, sd = 2.5)
		
		cross_site_traits_sigma_height <- exp(cross_site_traits_sigma_log_height)
		cross_site_traits_sigma_canopy_diameter <- exp(cross_site_traits_sigma_log_canopy_diameter)

		# # prior for sd of normal distribution for site-level mean trait values
		# phi_traits_sigma_log_height ~ dnorm(0, sd = 2.5)
		# phi_traits_sigma_log_canopy_diameter ~ dnorm(0, sd = 2.5)

		# phi_traits_sigma_height <- exp(phi_traits_sigma_log_height)
		# phi_traits_sigma_canopy_diameter <- exp(phi_traits_sigma_log_canopy_diameter)

		### INTEGRATION

		# LKJ prior on correlation matrix between non-biomass traits
		L_traits[1:n_responses, 1:n_responses] ~ dlkj_corr_cholesky(eta = eta, p = n_responses)

		# MVN precision matrix
		Omega_traits[1:n_responses, 1:n_responses] <- L_traits[1:n_responses, 1:n_responses] %*% t(L_traits[1:n_responses, 1:n_responses])

		# prior for hyperparameter of Cholesky
		eta_traits ~ dgamma(2, 1)

		# NON-BIOMASS TRAITS: parameters of trait distribution are latent and functions of environment
		# individual plant trait values are samples from the site-level distribution (next chunk after this one)
		for (i in 1:n_pheno_sites) {

			# relationship of biomass to the environment
			phi_traits_mu_height[i] <- inprod(beta_traits_height[1:n_terms_traits_height], x_by_site_traits_height[i, 1:n_terms_traits_height])
			phi_traits_mu_canopy_diameter[i] <- inprod(beta_traits_canopy_diameter[1:n_terms_traits_canopy_diameter], x_by_site_traits_canopy_diameter[i, 1:n_terms_traits_canopy_diameter])

			# indexing of traits is paramount!
			# 1: height
			# 2: canopy_diameter
			log_site_traits_mu[i, 1:n_responses] ~ dmnorm(mean = mvnorm_mu[i, 1:n_responses], cholesky = Omega[1:n_responses, 1:n_responses], prec_param = 0)
			mvnorm_mu[i, 1:n_responses] <- c(phi_traits_mu_height[i], phi_traits_mu_canopy_diameter[i])

			# log_site_traits_mu_height[i] ~ dnorm(phi_traits_mu_height[i], sd = phi_traits_sigma_height)
			# log_site_traits_mu_canopy_diameter[i] ~ dnorm(phi_traits_mu_canopy_diameter[i], sd = phi_traits_sigma_canopy_diameter)

			# site-level mean biomass
			site_traits_mu_height[i] <- exp(log_site_traits_mu[i, 1])
			site_traits_mu_canopy_diameter[i] <- exp(log_site_traits_mu[i, 2])

			# moment matching to get dgamma() parameters
			shape_traits_height[i] <- site_traits_mu_height[i]^2 / cross_site_traits_sigma_height^2
			shape_traits_canopy_diameter[i] <- site_traits_mu_canopy_diameter[i]^2 / cross_site_traits_sigma_canopy_diameter^2

			rate_traits_height[i] <- site_traits_mu_height[i] / cross_site_traits_sigma_height^2
			rate_traits_canopy_diameter[i] <- site_traits_mu_canopy_diameter[i] / cross_site_traits_sigma_canopy_diameter^2

		}

		# NON-BIOMASS TRAITS: likelihood of individual plants
		for (i in 1:n_trait_values) {

			# add/subtract a common value to the rate parameter specific to each plant to reflect within-plant correlation between traits on the same plant... assumes correlations between traits are all of the same sign
			log_rate_traits_plant_delta[i] ~ dnorm(0, sd = 10)
			rate_traits_plant_delta[i] <- exp(log_rate_traits_plant_delta[i])

			rate_traits_height_plant[i] <- rate_traits_height[site_index_traits[i]] + rate_traits_plant_delta[i]
			rate_traits_canopy_diameter_plant[i] <- rate_traits_canopy_diameter[site_index_traits[i]] + rate_traits_plant_delta[i]

			rate_traits_height_plant[i] ~ dconstraint(rate_traits_height_plant[i] > 0)
			rate_traits_canopy_diameter_plant[i] ~ dconstraint(rate_traits_canopy_diameter_plant[i] > 0)

			# likelihood
			# NB we assume "y" is in matrix form... easier for porting to new models with multiple traits and for use in auxiliary functions
			y_traits_height[i, 1] ~ dgamma(shape = shape_traits_height[site_index_traits[i]], rate = rate_traits_height_plant[i])
			y_traits_canopy_diameter[i, 1] ~ dgamma(shape = shape_traits_canopy_diameter[site_index_traits[i]], rate = rate_traits_canopy_diameter_plant[i])

			# simulated values for unconditional DHARMa residuals
			y_traits_sim_height[i] ~ dgamma(shape = shape_traits_height[site_index_traits[i]], rate = rate_traits_height_plant[i])
			y_traits_sim_canopy_diameter[i] ~ dgamma(shape = shape_traits_canopy_diameter[site_index_traits[i]], rate = rate_traits_canopy_diameter_plant[i])
	
			log_lik_y_traits_height[i] <- dgamma(y_traits_height[i, 1], shape = shape_traits_height[site_index_traits[i]], rate = rate_traits_height_plant[i], log = 1)
			log_lik_y_traits_canopy_diameter[i] <- dgamma(y_traits_canopy_diameter[i, 1], shape = shape_traits_canopy_diameter[site_index_traits[i]], rate = rate_traits_canopy_diameter_plant[i], log = 1)

		}

		log_lik_height <- sum(log_lik_y_traits_height[1:n_trait_values])
		log_lik_canopy_diameter <- sum(log_lik_y_traits_canopy_diameter[1:n_trait_values])
		log_lik <- log_lik_height + log_lik_canopy_diameter

		# # NON-BIOMASS TRAITS: posterior samplers for predictions of to counties in status quo and future
		# for (i in 1:n_counties) {

		# 	# non-biomass traits: status quo
		# 	log_traits_county_mu_sq[i] <-
		# 		inprod(beta_traits_mu[1:n_terms_traits], counties_x_traits_sq[i, 1:n_terms_traits])
		# 	exp_traits_mu_county_sq[i] <- exp(log_traits_county_mu_sq[i])

		# 	shape_traits_county_mu_sq[i] <- exp_traits_mu_county_sq[i]^2 / cross_site_traits_sigma^2
		# 	rate_traits_county_mu_sq[i] <- exp_traits_mu_county_sq[i] / cross_site_traits_sigma^2
		# 	traits_mu_county_sq[i] ~ dgamma(shape = shape_traits_county_mu_sq[i], rate = rate_traits_county_mu_sq[i])

		# 	# non-biomass traits: future ssp245_2041_2070
		# 	log_traits_county_mu_ssp245_2041_2070[i] <-
		# 		inprod(beta_traits_mu[1:n_terms_traits], counties_x_traits_ssp245_2041_2070[i, 1:n_terms_traits])
		# 	exp_traits_mu_county_ssp245_2041_2070[i] <- exp(log_traits_county_mu_ssp245_2041_2070[i])
			
		# 	shape_traits_county_mu_ssp245_2041_2070[i] <- exp_traits_mu_county_ssp245_2041_2070[i]^2 / cross_site_traits_sigma^2
		# 	rate_traits_county_mu_ssp245_2041_2070[i] <- exp_traits_mu_county_ssp245_2041_2070[i] / cross_site_traits_sigma^2
		# 	traits_mu_county_ssp245_2041_2070[i] ~ dgamma(shape = shape_traits_county_mu_ssp245_2041_2070[i], rate = rate_traits_county_mu_ssp245_2041_2070[i])

		# 	# non-biomass traits: future ssp245_2071_2100
		# 	log_traits_county_mu_ssp245_2071_2100[i] <-
		# 		inprod(beta_traits_mu[1:n_terms_traits], counties_x_traits_ssp245_2071_2100[i, 1:n_terms_traits])
		# 	exp_traits_mu_county_ssp245_2071_2100[i] <- exp(log_traits_county_mu_ssp245_2071_2100[i])
			
		# 	shape_traits_county_mu_ssp245_2071_2100[i] <- exp_traits_mu_county_ssp245_2071_2100[i]^2 / cross_site_traits_sigma^2
		# 	rate_traits_county_mu_ssp245_2071_2100[i] <- exp_traits_mu_county_ssp245_2071_2100[i] / cross_site_traits_sigma^2
		# 	traits_mu_county_ssp245_2071_2100[i] ~ dgamma(shape = shape_traits_county_mu_ssp245_2071_2100[i], rate = rate_traits_county_mu_ssp245_2071_2100[i])

		# 	# non-biomass traits: future ssp370_2041_2070
		# 	log_traits_county_mu_ssp370_2041_2070[i] <-
		# 		inprod(beta_traits_mu[1:n_terms_traits], counties_x_traits_ssp370_2041_2070[i, 1:n_terms_traits])
		# 	exp_traits_mu_county_ssp370_2041_2070[i] <- exp(log_traits_county_mu_ssp370_2041_2070[i])

		# 	shape_traits_county_mu_ssp370_2041_2070[i] <- exp_traits_mu_county_ssp370_2041_2070[i]^2 / cross_site_traits_sigma^2
		# 	rate_traits_county_mu_ssp370_2041_2070[i] <- exp_traits_mu_county_ssp370_2041_2070[i] / cross_site_traits_sigma^2
		# 	traits_mu_county_ssp370_2041_2070[i] ~ dgamma(shape = shape_traits_county_mu_ssp370_2041_2070[i], rate = rate_traits_county_mu_ssp370_2041_2070[i])

		# 	# non-biomass traits: future ssp370_2071_2100
		# 	log_traits_county_mu_ssp370_2071_2100[i] <-
		# 		inprod(beta_traits_mu[1:n_terms_traits], counties_x_traits_ssp370_2071_2100[i, 1:n_terms_traits])
		# 	exp_traits_mu_county_ssp370_2071_2100[i] <- exp(log_traits_county_mu_ssp370_2071_2100[i])

		# 	shape_traits_county_mu_ssp370_2071_2100[i] <- exp_traits_mu_county_ssp370_2071_2100[i]^2 / cross_site_traits_sigma^2
		# 	rate_traits_county_mu_ssp370_2071_2100[i] <- exp_traits_mu_county_ssp370_2071_2100[i] / cross_site_traits_sigma^2
		# 	traits_mu_county_ssp370_2071_2100[i] ~ dgamma(shape = shape_traits_county_mu_ssp370_2071_2100[i], rate = rate_traits_county_mu_ssp370_2071_2100[i])

		# }

		# # NON-BIOMASS TRAITS: posterior predictive sampler for response curves: site-level mean
		# # NB we assume just one predictor, so the response curve "x" is a matrix, not an array
		# for (i in 1:n_response_curve_values) {
				
		# 	log_response_curves_traits_mu[i] <-
		# 		inprod(beta_traits_mu[1:n_terms_traits], resp_curves_x_traits[i, 1:n_terms_traits])
			
		# 	response_curves_traits_mu[i] <- exp(log_response_curves_traits_mu[i])
			
		# }

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
	print(model$initializeInfo())

	calc <- model$calculate()
	say('model$calculate(): ', calc)
	if (is.na(calc) | is.infinite(calc)) stop('Likelihood is incalculable.')

	say('configureMCMC():', level = 2)

	# monitors for coefficients that have no indexing
	monitors_coeffs_not_indexed <- c(
		'cross_site_traits_sigma_height', 'cross_site_traits_sigma_canopy_diameter'
		# 'phi_traits_sigma'
	)

	# coefficients that have bracketed indexing
	monitors_coeffs_single_index <- c(
		'beta_traits_height',
		'beta_traits_canopy_diameter'
	)
	monitors_coeffs_double_index <- c(
	)

	monitors_derived_not_indexed <- c(
		'log_lik_height',
		'log_lik_canopy_diameter',
		'log_lik'
	)
	monitors_derived_single_index <- c(
		'site_traits_mu_height',
		'site_traits_mu_canopy_diameter'
	)
	monitors_derived_double_index <- c(
	)

	monitors_geog <- c(
		'traits_mu_county_sq', 'traits_mu_county_ssp245_2041_2070', 	
		'traits_mu_county_ssp245_2071_2100', 'traits_mu_county_ssp370_2041_2070', 'traits_mu_county_ssp370_2071_2100'
	)

	monitors_dharma <- c(
		'y_traits_sim_height',
		'y_traits_sim_canopy_diameter'
	)

	monitors_resp_curves <- c(
		'response_curves_traits_mu_height',
		'response_curves_traits_mu_canopy_diameter'
	)

	monitors <- c(monitors_coeffs_not_indexed, monitors_coeffs_single_index, monitors_coeffs_double_index, monitors_derived_not_indexed, monitors_derived_single_index, monitors_derived_double_index, monitors_geog, monitors_dharma, monitors_resp_curves)

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

	conf$removeSamplers('beta_traits_mu')
	conf$addSampler(target = c('beta_traits_mu'), type = 'AF_slice')
	say('AF slice sampler added to beta_traits_mu.')

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

	### collate all predictions into one SpatVector
	###############################################

		pred_vect_nam <- data_occs$ag_vect_sq

		for (var in monitors_geog) {
		
			preds <- hammer_extract(chains, param = var, j = TRUE, stat = 'mean')
			pred_vect_nam[[var]] <- preds
			this_var <- sub(var, pattern = 'traits_', replacement = paste0(trait, '_'))
			names(pred_vect_nam)[ncol(pred_vect_nam)] <- this_var

		}

		writeVector(pred_vect_nam, paste0(out_dir, '/prediction_vector.gpkg'), overwrite = TRUE)

	### post-modeling analysis of BIOMASS
	descrip <- paste0(trait, ': gamma~normal homoscedastic')
	workflow_postmodeling_generic(facet = trait, formulae = formulae, descrip = descrip, out_dir = out_dir)
	workflow_postmodeling_nonbiomass_single_trait(trait = trait, formula_traits = formula_traits, descrip = descrip, homoscedastic = TRUE, pred_vect_nam = pred_vect_nam, out_dir = out_dir)

say(date())
say('FINIS!', deco = '+', level = 1)

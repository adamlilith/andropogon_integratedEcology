### MODELING ANDROPOGON GERARDI DISTRIBUTION, PHENOTYPE, PHYSIOLOGY, GENOTYPE, and ASSOCIATED MICROBIAL COMMUNITIES
### Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2023-12
###
### This model estimates the distribution of site-level mean values of non-biomass traits ("facets") of Andropogon gerardi ramets (i.e., height, blade width, SPAD, etc.). The distribution of values among ramets at a site follows a hurdle (zero-inflated) gamma or hurdle lognormal. The site-level mean of the facet is a function of soil/climate. The probability of a non-zero value is a function of the same covariate(s) and has the same functional form as the submodel of the mean.
###
### source('C:/Kaji/R/andropogon_integratedEcology/sdm_pdm/sdm_pdm_06a_model_nonbiomass_single_facet_zero_inflated.r')

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
	trial <- FALSE # FALSE for running for real

	# do cross-validation?
	crossvalidate <- TRUE
	# crossvalidate <- FALSE

	# # ALL non-biomass facets to model
	# facets <- c('n_concentration', 'cn_ratio', 'height', 'blade_width', 'leaf_thickness', 'spad', 'canopy_diameter', 'photosynthetic_rate', 'stomatal_conductance', 'internal_co2', 'transpiration_rate', 'water_potential')

	# # facets that have positive-only values
	facets <- c(
		'blade_width',
		'canopy_diameter',
		'cn_ratio',
		'height',
		'internal_co2',
		'leaf_thickness',
		'n_concentration',
		'photosynthetic_rate',
		'spad',
		'stomatal_conductance',
		'transpiration_rate'
	)
	# facets <- c('cn_ratio')
	
	# # facets that can have positive and negative values
	# facets <- c('water_potential')
	
	# facets <- 'n_concentration'
	# facets <- 'cn_ratio'
	# facets <- 'height'
	# facets <- 'blade_width'
	# facets <- 'leaf_thickness'
	# facets <- 'spad'
	# facets <- 'canopy_diameter'
	# facets <- 'photosynthetic_rate'
	# facets <- 'stomatal_conductance'
	# facets <- 'internal_co2'
	# facets <- 'transpiration_rate'
	# facets <- 'water_potential' # normal distribution

	# distribution of response variable
	resp_distribs <- c('hGamma', 'hurdleLN')

	### MCMC settings
	if (!trial) {

		niter <- 160000
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

	for (facet in facets) {

		for (resp_distrib in resp_distribs) {

			names(resp_distrib) <- facet
			transform <- if (resp_distrib == 'hGamma') { 'exponential' } else if (resp_distrib == 'hurdleLN') { 'identity' }

			### formula for how facets of species responds to environment
			formulae_facet <- list(

				# precip logged
				'log(bio12)' = ~ 1 + bio12_log10p1,
				'log(bio12)^2' = ~ 1 + bio12_log10p1 + I(bio12_log10p1^2),
				'bio1^2_log(bio12)' = ~ 1 + bio1 + bio12_log10p1 + I(bio1^2),
				'bio1^2_x_log(bio12)' = ~ 1 + bio1 + bio12_log10p1 + I(bio1^2) + bio1:bio12_log10p1,
				'bio1_log(bio12)' = ~ 1 + bio1 + bio12_log10p1,
				'bio1_log(bio12)^2' = ~ 1 + bio1 + bio12_log10p1 + I(bio12_log10p1^2),
				'bio1_x_log(bio12)' = ~ 1 + bio1 + bio12_log10p1 + bio1:bio12_log10p1,
				'log(bio12)_sand' = ~ 1 + bio12_log10p1 + sand,
				'log(bio12)_x_sand' = ~ 1 + bio12_log10p1 + sand + bio12_log10p1:sand,
				'log(bio12)_nitrogen' = ~ 1 + bio12_log10p1 + nitrogen,
				'log(bio12)_x_nitrogen' = ~ 1 + bio12_log10p1 + nitrogen + bio12_log10p1:nitrogen,
				'log(bio12)_sand_nitrogen' = ~ 1 + bio12_log10p1 + sand + nitrogen,

				# precip not logged
				'bio12' = ~ 1 + bio12,
				'bio12^2' = ~ 1 + bio12 + I(bio12^2),
				'bio12_sand' = ~ 1 + bio12 + sand,
				'bio12_x_sand' = ~ 1 + bio12 + sand + bio12:sand,
				'bio12_nitrogen' = ~ 1 + bio12 + nitrogen,
				'bio12_x_nitrogen' = ~ 1 + bio12 + nitrogen + bio12:nitrogen,
				'bio12_sand_nitrogen' = ~ 1 + bio12 + sand + nitrogen,
				'aridity' = ~ 1 + aridity,
				'aridity^2' = ~ 1 + aridity + I(aridity^2),
				'aridity_sand' = ~ 1 + aridity + sand,
				'aridity_x_sand' = ~ 1 + aridity + sand + aridity:sand,
				'bio1' = ~ 1 + bio1,
				'bio1^2' = ~ 1 + bio1 + I(bio1^2),
				'bio1^2_bio12' = ~ 1 + bio1 + bio12 + I(bio1^2),
				'bio1^2_x_bio12' = ~ 1 + bio1 + bio12 + I(bio1^2) + bio1:bio12,
				'bio1_bio12' = ~ 1 + bio1 + bio12,
				'bio1_bio12^2' = ~ 1 + bio1 + bio12 + I(bio12^2),
				'bio1_x_bio12' = ~ 1 + bio1 + bio12 + bio1:bio12,
				'nitrogen' = ~ 1 + nitrogen,
				'ph' = ~ 1 + ph,
				'ph^2' = ~ 1 + ph + I(ph^2),
				'sand' = ~ 1 + sand,
				'sand^2' = ~ 1 + sand + I(sand^2)
			)

			if (facet %in% c('internal_co2', 'photosynthetic_rate', 'spad', 'stomatal_conductance', 'transpiration_rate')) {

				formulae_facet <- c(
					formulae_facet,
					list(
						'srad' = ~ 1 + insolation_2000_growing_season_kWh_per_m2,
						'bio12_srad' = ~ 1 + bio12 + insolation_2000_growing_season_kWh_per_m2,
						'log(bio12)_srad' = ~ 1 + bio12_log10p1 + insolation_2000_growing_season_kWh_per_m2,
						'bio1_srad' = ~ 1 + bio1 + insolation_2000_growing_season_kWh_per_m2,
						'bio12_x_srad' = ~ 1 + bio12 + insolation_2000_growing_season_kWh_per_m2,
						'log(bio12)_x_srad' = ~ 1 + bio12_log10p1 + insolation_2000_growing_season_kWh_per_m2 + bio12_log10p1:insolation_2000_growing_season_kWh_per_m2,
						'aridity_srad' = ~ 1 + aridity + insolation_2000_growing_season_kWh_per_m2,
						'aridity_x_srad' = ~ 1 + aridity + insolation_2000_growing_season_kWh_per_m2 + aridity:insolation_2000_growing_season_kWh_per_m2,
						'bio1_bio12_srad' = ~ 1 + bio1 + bio12 + insolation_2000_growing_season_kWh_per_m2,
						'bio1^2_bio12_srad' = ~ 1 + bio1 + bio12 + insolation_2000_growing_season_kWh_per_m2 + I(bio1^2),
						'bio1^2_x_bio12_srad' = ~ 1 + bio1 + bio12 + insolation_2000_growing_season_kWh_per_m2 + I(bio1^2) + bio1:bio12
					)
				)
			}
		
			for (count_formula in seq_along(formulae_facet)) {

				formula_facet <- formulae_facet[[count_formula]]
				formula_psi <- formula_facet

				filename_facet <- names(formulae_facet)[count_formula]

				out_dir <- paste0('./outputs_loretta/integrated_sdm_pdm/models_', facet, '/', ifelse(trial, 'TRIAL_', ''), '[', facet, '~', tolower(resp_distrib), '(', filename_facet, ')]')

				ok <- trial | !file.exists(out_dir)
				if (!ok) {
					say(out_dir)
					warning('Output folder already exists.')
				} else {
					
					dirCreate(out_dir)

					sink(paste0(out_dir, '/runtime_log.txt'), split = TRUE)
					say('MODELING ', toupper(facet))
					say('Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | ', date(), post = 2)

					say('This model estimates the distribution of site-level mean values of non-biomass traits ("facets") of Andropogon gerardi ramets (i.e., height, blade width, SPAD, etc.). The distribution of values among ramets at a site follows a hurdle (zero-inflated) gamma or hurdle lognormal. The site-level mean of the facet is a function of soil/climate. The probability of a non-zero value is a function of the same covariate(s) and has the same functional form as the submodel of the mean.', breaks = 60, post = 1)

					say('MCMC settings:', level = 2)
					say('facet ........................ ', facet)
					say('trial ........................ ', trial)
					say('niter ........................ ', niter)
					say('nburnin ...................... ', nburnin)
					say('thin ......................... ', thin)
					say('nchains ...................... ', nchains)
					say('formula_facet ................ ', paste(as.character(formula_facet), collapse = ' '))
					say('formula_psi .................. ', paste(as.character(formula_psi), collapse = ' '))
					say('resp_distrib ................. ', unname(resp_distrib))
					say('transform .................... ', transform)
					say('zero_inflated ................ ', zero_inflated)
					say('calib ........................ ', calib, post = 2)

					say('out_dir:')
					say(out_dir)

					formulae <- list(
						formula_facet = formula_facet,
						formula_psi = formula_psi
					)
					saveRDS(formulae, paste0(out_dir, '/formulae.rds'))

					########################s
					### data preparation ###
					########################

					data_facet <- prepare_nonbiomass_data(facet = facet, formula_facet = formula_facet, n_response_curve_values = n_response_curve_values, calib = calib)

					#########################
					### inputs for nimble ###
					#########################

					say('Inputs:', level = 2)
					data <- list(
						y_facet = data_facet$y_facet # trait of individual plants
					)

					constants <- list(
						
						# sampled sites
						n_pheno_sites = data_facet$n_pheno_sites, # number of phenotype sample sites

						### facet
						x_by_site_facet = data_facet$x_by_site, # MM with covariates for facet (scaled)
						n_plants = data_facet$n_plants, # number of facet observations
						site_index_facet = data_facet$site_index_facet, # index of sampled site for each row in facet data
						n_terms_facet = data_facet$n_terms, # number of terms in formula for facet model (including intercept)
						n_terms_psi = data_facet$n_terms # number of terms in formula for facet model (including intercept)

					)

					constants <- c(constants, constants_shared_facet)
					if (!is.null(formula_psi)) constants <- c(constants, constants_shared_psi)

					# use GLM to get initial values for betas
					if (resp_distrib == 'hGamma') {
						fam <- Gamma(link = 'log')
						y <- data_facet$y_facet
					} else if (resp_distrib == 'hurdleLN' ) {
						fam <- gaussian(link = 'identity')
						y <- log(data_facet$y_facet)
					}
					

					indices <- rep(1:26, each = 6)
					x <- constants$x_by_site
					x <- x[indices, ]
					prelim_model <- glm.fit(y = y, x = x, family = fam)
					beta_inits <- prelim_model$coefficients

					# manual selection based on what yields non-infinite likelihood using just initializations
					log_sigma_facet_within_sites_init <- if (facet %in% c('height', 'n_concentration')) {
						1
					} else if (facet == 'internal_co2') {
						3
					} else if (facet == 'canopy_diameter') {
						10
					} else {
						0.1
					}

					inits <- list(

						log_sigma_facet_within_sites = log_sigma_facet_within_sites_init,
						y_facet_sim = data_facet$y_facet, # simulated values for facet (for DHARMa residuals)
						beta_facet = beta_inits,
						beta_psi = rep(0, constants$n_terms_psi)

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
					
						# FACET: priors for relationship of site-level mean facet to environment
						beta_facet[1] ~ dnorm(0, sd = beta_facet_prior_dnorm_sd_1)
						for (i in 2:n_terms_facet) {
							beta_facet[i] ~ dnorm(0, sd = beta_facet_prior_dnorm_sd) # broad prior
						}

						# # prior for sd of mean of facet at a site on lognormal (~ half-Cauchy), ==> vague
						# log(sigma_facet_among_sites) ~ dnorm(0, sd = sigma_facet_among_sites_log_prior_sd)

						# prior for sd of individual plant facet value on lognormal (~ half-Cauchy), ==> vague
						log(sigma_facet_within_sites) ~ dnorm(0, sd = sigma_facet_within_sites_log_prior_sd)

						# FACET: parameters of facet distribution are latent and functions of environment
						# individual plant biomasses are samples from the site-level distribution defined by the site-level distribution (next chunk after this one)
						for (i in 1:n_pheno_sites) {

							# # relationship of facet to the environment
							# log_mu_site_facet[i] ~ dnorm(log_site_facet_mu_mean[i], sd = sigma_facet_among_sites)
							# log_site_facet_mu_mean[i] <- inprod(beta_facet[1:n_terms_facet], x_by_site_facet[i, 1:n_terms_facet])
							
							# relationship of facet to the environment
							log_mu_site_facet[i] <- inprod(beta_facet[1:n_terms_facet], x_by_site_facet[i, 1:n_terms_facet])
							
							# zero-inflation
							logit(psi_site[i]) <- inprod(beta_psi[1:n_terms_facet], x_by_site_facet[i, 1:n_terms_facet])
							# z_site[i] ~ dbern(psi[i])

							# site-level mean facet value
							mu_facet_site[i] <- exp(log_mu_site_facet[i])

						}

						# log_lik <- sum(log_lik_facet[1:n_plants])

					})

					### response code
					#################
					# Switch out zero-inflated gamma or lognormal, depending on which one we're doing.

					if (resp_distrib == 'hGamma') {
					
						response_distrib_code <- nimbleCode({

							# FACET: parameters of facet distribution are latent and functions of environment
							# individual plant facet values are samples from the site-level distribution defined by the site-level distribution (next chunk after this one)
							for (i in 1:n_pheno_sites) {

								# moment matching to get gamma() parameters
								shape_facet[i] <- mu_facet_site[i]^2 / sigma_facet_within_sites^2
								rate_facet[i] <- mu_facet_site[i] / sigma_facet_within_sites^2

							}

							# FACET: likelihood of individual plants
							for (i in 1:n_plants) {

								# likelihood
								# y_facet[i] ~ dHurdleGamma(shape = shape_facet[site_index_facet[i]], rate = rate_facet[site_index_facet[i]], z = z_site[site_index_facet[i]])
								y_facet[i] ~ dHurdleGamma(shape = shape_facet[site_index_facet[i]], rate = rate_facet[site_index_facet[i]], psi = psi_site[site_index_facet[i]])

								# simulated values for unconditional DHARMa residuals
								# y_facet_sim[i] ~ dHurdleGamma(shape = shape_facet[site_index_facet[i]], rate = rate_facet[site_index_facet[i]], z = z_site[site_index_facet[i]])
								y_facet_sim[i] ~ dHurdleGamma(shape = shape_facet[site_index_facet[i]], rate = rate_facet[site_index_facet[i]], psi = psi_site[site_index_facet[i]])

								# log_lik_facet[i] <- dHurdleGamma(y_facet[i], shape = shape_facet[site_index_facet[i]], rate = rate_facet[site_index_facet[i]], z = z_site[site_index_facet[i]])
								# log_lik_facet[i] <- dHurdleGamma(y_facet[i], shape = shape_facet[site_index_facet[i]], rate = rate_facet[site_index_facet[i]], psi = psi_site[site_index_facet[i]])
						
							}

						})
					
					} else if (resp_distrib == 'hurdleLN') {
					
						response_distrib_code <- nimbleCode({

							# FACET: likelihood of individual plants
							for (i in 1:n_plants) {

								# likelihood
								# y_facet[i] ~ dHLN(meanlog = log_mu_site_facet[site_index_facet[i]], sdlog = sigma_facet_within_sites, z = z_site[site_index_facet[i]])
								y_facet[i] ~ dHLN(meanlog = log_mu_site_facet[site_index_facet[i]], sdlog = sigma_facet_within_sites, psi = psi_site[site_index_facet[i]])

								# simulated values for unconditional DHARMa residuals
								# y_facet_sim[i] ~ dHLN(meanlog = log_mu_site_facet[site_index_facet[i]], sdlog = sigma_facet_within_sites, z = z_site[site_index_facet[i]])
								y_facet_sim[i] ~ dHLN(meanlog = log_mu_site_facet[site_index_facet[i]], sdlog = sigma_facet_within_sites, psi = psi_site[site_index_facet[i]])

								# log_lik_facet[i] <- dHLN(y_facet[i], meanlog = log_mu_site_facet[site_index_facet[i]], sdlog = sigma_facet_within_sites, z = z_site[site_index_facet[i]])
								# log_lik_facet[i] <- dHLN(y_facet[i], meanlog = log_mu_site_facet[site_index_facet[i]], sdlog = sigma_facet_within_sites, psi = psi_site[site_index_facet[i]])

							}

						})
						
					}

					model_code <- glueNimbleCode(
						model_code,
						response_distrib_code,
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
						'sigma_facet_within_sites'#, 'sigma_facet_among_sites'
					)

					# coefficients that have bracketed indexing
					monitors_coeffs_single_index <- c(
						'beta_facet', 'beta_psi'
					)

					monitors_coeffs_double_index <- c(
					)

					monitors_derived_not_indexed <- c(
						# 'log_lik'
					)

					monitors_derived_single_index <- c(
						'mu_facet_site'
					)
					monitors_derived_double_index <- c(
					)

					monitors_geog_conus <- c(
					)

					monitors_dharma <- c(
						'y_facet_sim'
					)

					monitors_resp_curves <- c(
						# 'response_curves_facet_mu', 'response_curves_psi'
					)

					monitors_debug <- c(
						'psi_site', 'log_mu_site_facet'
					)
					# if (resp_distrib == 'hurdleGamma') monitors_debug <- c(monitors_debug, 'shape_facet', 'rate_facet')

					monitors <- c(monitors_coeffs_not_indexed, monitors_coeffs_single_index, monitors_coeffs_double_index, monitors_derived_not_indexed, monitors_derived_single_index, monitors_derived_double_index, monitors_dharma, monitors_resp_curves, monitors_debug)

					conf <- configureMCMC(
						model,
						monitors = monitors,
						print = TRUE,
						enableWAIC = TRUE
					)

					vars <- 'beta_facet'
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
					descrip <- paste0(facet, ' ~ ', resp_distrib, '(', ifelse(resp_distrib == 'hGamma', 'exp(', ''), 'normal(env))', ifelse(resp_distrib == 'hGamma', ')', ''))
					workflow_postmodeling_generic(facet = facet, formulae = formulae, descrip = descrip, out_dir = out_dir)

					workflow_postmodeling_nonbiomass_single_facet(facet = facet, chains = chains, descrip = descrip, formula_facet = formula_facet, formula_psi = formula_psi, resp_distrib = resp_distrib, transform = transform, crossvalidate = crossvalidate, out_dir = out_dir)

				} # output folder already exists?

			} # next formula

		} # next response type (ZIG or 'hurdleLN')

	} # next facet

say(date())
say('FINIS!', deco = '+', level = 1)

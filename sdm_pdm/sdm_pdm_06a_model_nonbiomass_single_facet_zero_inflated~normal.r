### MODELING ANDROPOGON GERARDI DISTRIBUTION, PHENOTYPE, PHYSIOLOGY, GENOTYPE, and ASSOCIATED MICROBIAL COMMUNITIES
### Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2023-12
###
### This model estimates the distribution of site-level mean values of non-biomass traits ("facets") of Andropogon gerardi ramets (i.e., height, blade width, SPAD, etc.). The distribution of values among ramets at a site follows a zero-inflated gamma or zero-inflated lognormal. The site-level mean of the facet is a function of soil/climate. The probability of a zero value is a function of the same covariate(s) and has the same functional form as the submodel of the mean.
###
### source('C:/Kaji/R/andropogon_integratedEcology/sdm_pdm/sdm_pdm_06a_model_nonbiomass_single_facet_zero_inflated~normal.r')

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
	facets <- c('blade_width', 'canopy_diameter', 'cn_ratio', 'height', 'internal_co2', 'leaf_thickness', 'n_concentration', 'photosynthetic_rate', 'spad', 'stomatal_conductance', 'transpiration_rate')
	
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
	# ZIG = zero-inflated gamma
	# ZILN = zero-inflated lognormal
	# resp_distribs <- c('ZIG', 'ZILN')
	resp_distribs <- c('ZILN', 'ZIG')
	# resp_distribs <- 'ZIG'

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

	for (facet in facets) {

		for (resp_distrib in resp_distribs) {

			names(resp_distrib) <- facet
			transform <- if (resp_distrib == 'ZIG') { 'exponential' } else if (resp_distrib == 'ZILN') { 'identity' }

			for (log_precip in c(TRUE, FALSE)) {

				if (log_precip) {
				
					### formula for how facets of species responds to environment
					formulae_facet <- list(
						'bio12' = ~ 1 + bio12,
						'bio12^2' = ~ 1 + bio12 + I(bio12^2),
						'bio12_sand' = ~ 1 + bio12 + sand,
						'bio12_x_sand' = ~ 1 + bio12 + sand + bio12:sand,
						'bio12_nitrogen' = ~ 1 + bio12 + site_nitrogen,
						'bio12_x_nitrogen' = ~ 1 + bio12 + site_nitrogen + bio12:site_nitrogen,
						'bio12_sand_nitrogen' = ~ 1 + bio12 + sand + site_nitrogen,
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
						'ph' = ~ 1 + ph,
						'ph^2' = ~ 1 + ph + I(ph^2),
						'sand' = ~ 1 + sand,
						'sand^2' = ~ 1 + sand + I(sand^2)
					)
				
				} else if (!log_precip) {
				
					### formula for how facets of species responds to environment
					formulae_facet <- list(
						'bio1^2_bio12' = ~ 1 + bio1 + bio12 + I(bio1^2),
						'bio1^2_x_bio12' = ~ 1 + bio1 + bio12 + I(bio1^2) + bio1:bio12,
						'bio1_bio12' = ~ 1 + bio1 + bio12,
						'bio1_bio12^2' = ~ 1 + bio1 + bio12 + I(bio12^2),
						'bio1_x_bio12' = ~ 1 + bio1 + bio12 + bio1:bio12,
						'bio12' = ~ 1 + bio12,
						'bio12^2' = ~ 1 + bio12 + I(bio12^2),
						'bio12_sand' = ~ 1 + bio12 + sand,
						'bio12_x_sand' = ~ 1 + bio12 + sand + bio12:sand,
						'nitrogen' = ~ 1 + site_nitrogen,
						'bio12_nitrogen' = ~ 1 + bio12 + site_nitrogen,
						'bio12_x_nitrogen' = ~ 1 + bio12 + site_nitrogen + bio12:site_nitrogen,
						'bio12_sand_nitrogen' = ~ 1 + bio12 + sand + site_nitrogen
					)
		
				}

				for (count_formula in seq_along(formulae_facet)) {

					formula_facet <- formulae_facet[[count_formula]]
					formula_psi <- formula_facet

					preds_filename <- names(formulae_facet)[count_formula]

					out_dir <- paste0('./outputs_loretta/integrated_sdm_pdm/models_', facet, '/', ifelse(trial, 'TRIAL_', ''), '[', facet, '_', tolower(resp_distrib), '~normal~', preds_filename, ']', ifelse(log_precip, '_log_precip', ''), '/')

					ok <- trial | !file.exists(out_dir)
					if (!ok) {
						say(out_dir)
						warning('Output folder already exists.')
					} else {
						
						dirCreate(out_dir)

						sink(paste0(out_dir, '/runtime_log.txt'), split = TRUE)
						say('MODELING ', toupper(facet))
						say('Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | ', date(), post = 2)

						say('This model estimates the distribution of site-level mean values of non-biomass traits ("facets") of Andropogon gerardi ramets (i.e., height, blade width, SPAD, etc.). The distribution of values among ramets at a site follows a zero-inflated gamma or zero-inflated lognormal. The site-level mean of the facet is a function of soil/climate. The probability of a zero value is a function of the same covariate(s) and has the same functional form as the submodel of the mean.', breaks = 60, post = 1)

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

						data_facet <- prepare_nonbiomass_data(facet = facet, formula_facet = formula_facet, log_precip = log_precip, n_response_curve_values = n_response_curve_values, calib = calib)

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
							x_by_site_facet = data_facet$x_by_site_facet, # MM with covariates for facet (scaled)
							n_plants = data_facet$n_plants, # number of facet observations
							site_index_facet = data_facet$site_index_facet, # index of sampled site for each row in facet data
							n_terms_facet = data_facet$n_terms_facet, # number of terms in formula for facet model (including intercept)
							n_terms_psi = data_facet$n_terms_facet, # number of terms in formula for facet model (including intercept)

							n_covariates_facet = data_facet$n_covariates_facet,
							resp_curves_x_facet = data_facet$resp_curves_x_facet, # response curve array for facet

							# response curves (general)
							n_response_curve_values = n_response_curve_values # number of values in response curve array

						)

						constants <- c(constants, constants_shared_facet)
						if (!is.null(formula_psi)) constants <- c(constants, constants_shared_psi)

						# use GLM to get initial values for betas
						if (resp_distrib == 'ZIG') {
							fam <- Gamma(link = 'log')
							y <- data_facet$y_facet
						} else if (resp_distrib == 'ZILN' ) {
							fam <- gaussian(link = 'identity')
							y <- log(data_facet$y_facet)
						}

						indices <- rep(1:26, each = 6)
						x <- constants$x_by_site_facet
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

							log_site_facet_mu = log(unlist(data_facet$site_vect[[paste0(facet, '_mean')]])),
							log_sigma_facet_within_sites = log_sigma_facet_within_sites_init,
							log_sigma_facet_among_sites = 1,
							y_facet_sim = data_facet$y_facet, # simulated values for facet (for DHARMa residuals)
							beta_facet = beta_inits,
							beta_psi = beta_inits,
							z_site = rep(1, data_facet$n_pheno_sites)

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
							for (i in 1:n_terms_facet) {
								beta_facet[i] ~ dnorm(0, sd = beta_facet_prior_dnorm_sd) # broad prior
							}

							# FACET: priors for probability of presence
							beta_psi[1] ~ dnorm(0, sd = beta_psi_prior_dnorm_sd_1)
							for (i in 2:n_terms_psi) {
								beta_psi[i] ~ ddexp(0, rate = beta_psi_prior_dnorm_sd)
							}

							# prior for sd of mean of facet at a site on lognormal (~ half-Cauchy), ==> vague
							log(sigma_facet_among_sites) ~ dnorm(0, sd = sigma_facet_among_sites_log_prior_sd)

							# prior for sd of individual plant facet value on lognormal (~ half-Cauchy), ==> vague
							log(sigma_facet_within_sites) ~ dnorm(0, sd = sigma_facet_within_sites_log_prior_sd)

							# FACET: parameters of facet distribution are latent and functions of environment
							# individual plant biomasses are samples from the site-level distribution defined by the site-level distribution (next chunk after this one)
							for (i in 1:n_pheno_sites) {

								# relationship of facet to the environment
								log_site_facet_mu[i] ~ dnorm(log_site_facet_mu_mean[i], sd = sigma_facet_among_sites)
								log_site_facet_mu_mean[i] <- inprod(beta_facet[1:n_terms_facet], x_by_site_facet[i, 1:n_terms_facet])
								
								# zero-inflation
								logit(psi[i]) <- inprod(beta_psi[1:n_terms_facet], x_by_site_facet[i, 1:n_terms_facet])
								z_site[i] ~ dbern(psi[i])

								# site-level mean facet value
								mu_facet_site[i] <- exp(log_site_facet_mu[i])

							}

							log_lik <- sum(log_lik_facet[1:n_plants])

						})

						### response code
						#################
						# Switch out zero-inflated gamma or lognormal, depending on which one we're doing.

						if (resp_distrib == 'ZIG') {
						
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
									y_facet[i] ~ dZIG(shape = shape_facet[site_index_facet[i]], rate = rate_facet[site_index_facet[i]], z = z_site[site_index_facet[i]])

									# simulated values for unconditional DHARMa residuals
									y_facet_sim[i] ~ dZIG(shape = shape_facet[site_index_facet[i]], rate = rate_facet[site_index_facet[i]], z = z_site[site_index_facet[i]])

									log_lik_facet[i] <- dZIG(y_facet[i], shape = shape_facet[site_index_facet[i]], rate = rate_facet[site_index_facet[i]], z = z_site[site_index_facet[i]])
							
								}

							})
						
						} else if (resp_distrib == 'ZILN') {
						
							response_distrib_code <- nimbleCode({

								# FACET: likelihood of individual plants
								for (i in 1:n_plants) {

									# likelihood
									y_facet[i] ~ dZILN(meanlog = log_site_facet_mu[site_index_facet[i]], sdlog = sigma_facet_within_sites, z = z_site[site_index_facet[i]])

									# simulated values for unconditional DHARMa residuals
									y_facet_sim[i] ~ dZILN(meanlog = log_site_facet_mu[site_index_facet[i]], sdlog = sigma_facet_within_sites, z = z_site[site_index_facet[i]])

									log_lik_facet[i] <- dZILN(y_facet[i], meanlog = log_site_facet_mu[site_index_facet[i]], sdlog = sigma_facet_within_sites, z = z_site[site_index_facet[i]])

								}

							})
							
						}

						### response curve model code
						#############################
						# Which one we use depends on how many covariates are in the formula for the facet.
						# If just one, the response curve is a matrix. If >1, an array.
						if (data_facet$n_covariates_facet == 1) {

							### univariate
							response_curve_code <- nimbleCode({

								# FACET: posterior predictive sampler for response curves: site-level mean
								# NB we assume just one predictor for facet, so the response curve "x" is a matrix, not an array
								for (i in 1:n_response_curve_values) {
										
									log_response_curves_facet_mu[i] <-
										inprod(beta_facet[1:n_terms_facet], resp_curves_x_facet[i, 1:n_terms_facet])
									response_curves_facet_mu[i] <- exp(log_response_curves_facet_mu[i])

									logit(response_curves_psi[i]) <-
										inprod(beta_psi[1:n_terms_facet], resp_curves_x_facet[i, 1:n_terms_facet])
									
								}

							})
						
						} else {
						
							### multivariate
							response_curve_code <- nimbleCode({

								# FACET: posterior predictive sampler for response curves: site-level mean
								# NB we assume just one predictor for facet, so the response curve "x" is a matrix, not an array
								for (i in 1:n_covariates_facet) {
									
									for (j in 1:n_response_curve_values) {
											
										log_response_curves_facet_mu[j, i] <-
											inprod(beta_facet[1:n_terms_facet], resp_curves_x_facet[j, 1:n_terms_facet, i])
										response_curves_facet_mu[j, i] <- exp(log_response_curves_facet_mu[j, i])

										logit(response_curves_psi[j, i]) <-
											inprod(beta_psi[1:n_terms_facet], resp_curves_x_facet[j, 1:n_terms_facet, i])
										
									}
								}

							})
						
						}

						model_code <- glueNimbleCode(model_code, response_distrib_code, response_curve_code)

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
							'sigma_facet_within_sites', 'sigma_facet_among_sites'
						)

						# coefficients that have bracketed indexing
						monitors_coeffs_single_index <- c(
							'beta_facet', 'beta_psi'
						)

						monitors_coeffs_double_index <- c(
						)

						monitors_derived_not_indexed <- c(
							'log_lik'
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
							'response_curves_facet_mu', 'response_curves_psi'
						)

						monitors <- c(monitors_coeffs_not_indexed, monitors_coeffs_single_index, monitors_coeffs_double_index, monitors_derived_not_indexed, monitors_derived_single_index, monitors_derived_double_index, monitors_dharma, monitors_resp_curves)

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

						say('session info', level = 2)
						print(sessionInfo())

						say(date(), pre = 1)
						sink()

					######################################################
					### post-modeling analysis for FACET-only models ###
					######################################################

						### post-modeling analysis of FACET
						descrip <- paste0(facet, ' ~ ', resp_distrib, '(', ifelse(resp_distrib == 'ZIG', 'exp(', ''), 'normal(env))', ifelse(resp_distrib == 'ZIG', ')', ''))
						workflow_postmodeling_generic(facet = facet, formulae = formulae, descrip = descrip, out_dir = out_dir)

						workflow_postmodeling_nonbiomass_single_facet(facet = facet, chains = chains, descrip = descrip, formula_facet = formula_facet, formula_psi = formula_psi, resp_distrib = resp_distrib, transform = transform, log_precip = log_precip, crossvalidate = crossvalidate, out_dir = out_dir)

					} # output folder already exists?

				} # next formula

			} # log_precip = T/F

		} # next response type (ZIG or ZILN)

	} # next facet

say(date())
say('FINIS!', deco = '+', level = 1)

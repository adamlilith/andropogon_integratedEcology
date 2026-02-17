### MODELING ANDROPOGON GERARDI DISTRIBUTION, PHENOTYPE, PHYSIOLOGY, GENOTYPE, and ASSOCIATED MICROBIAL COMMUNITIES
### Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2023-12
###
### This model estimates the distribution of site-level mean Andropogon gerardi ramet biomass. The distribution of biomasses among ramets at a site follows a zero-inflated gamma or zero-inflated lognormal distribution with individual ramet biomasses drawn from this distribution. The site-level mean biomass is a function of soil/climate. The probability of presence (non-zero) is a function of the same covariate(s) and has the same functional form as the submodel of the mean.
###
### source('C:/Kaji/R/andropogon_integratedEcology/sdm_pdm/sdm_pdm_04b_model_biomass_zero_inflated~normal.r')

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

	# do cross-validation?
	# crossvalidate <- FALSE
	crossvalidate <- TRUE

	if (!trial) {

		### MCMC settings
		# niter <- 1600000
		# nburnin <- niter / 2
		# thin <- 800
		# nchains <- 4
		# waic <- TRUE

		niter <- 400000
		nchains <- 4

	} else {
		
		### MCMC settings FOR TESTING
		niter <- 2000
		nchains <- 2

	}
	nburnin <- niter / 2
	thin <- (niter - nburnin) / 1000

	zero_inflated <- TRUE # SPECIFIC TO THIS SCRIPT--SHOULD NOT BE CHANGED

	### formula for how aspects of species responds to environment

	resp_distribs <- c('ZIG', 'ZILN')
	for (resp_distrib in resp_distribs) {

		transform <- if (resp_distrib == 'ZIG') { 'exponential' } else if (resp_distrib == 'ZILN' ) { 'identity' }

		for (log_precip in c(FALSE, TRUE)) {

			if (!log_precip) {

				all_formulae <- list(
					'bio12' = ~ 1 + bio12,
					'bio1' = ~ 1 + bio1,
					'bio1^2' = ~ 1 + bio1 + I(bio1^2),
					'bio1_bio12' = ~ 1 + bio1 + bio12,
					'bio1_x_bio12' = ~ 1 + bio1 + bio12 + bio1:bio12,
					'bio1^2_bio12' = ~ 1 + bio1 + bio12 + I(bio1^2),
					'bio1^2_x_bio12' = ~ 1 + bio1 + bio12 + I(bio1^2) + bio1:bio12,
					'aridity' = ~ 1 + aridity,
					'ph' = ~ 1 + ph,
					'ph^2' = ~ 1 + ph + I(ph^2),
					'sand' = ~ 1 + sand,
					'bio12_sand' = ~ 1 + bio12 + sand,
					'bio12_x_sand' = ~ 1 + bio12 + sand + bio12:sand,
					'aridity_sand' = ~ 1 + aridity + sand,
					'aridity_x_sand' = ~ 1 + aridity + sand + aridity:sand,
					'nitrogen' = ~ 1 + site_nitrogen,
					'bio1_nitrogen' = ~ 1 + bio1 + site_nitrogen,
					'bio1_x_nitrogen' = ~ 1 + bio1 + site_nitrogen + bio1:site_nitrogen,
					'bio12_nitrogen' = ~ 1 + bio12 + site_nitrogen,
					'bio12_x_nitrogen' = ~ 1 + bio12 + site_nitrogen + bio12:site_nitrogen
				)

			} else if (log_precip) {

				all_formulae <- list(
					'bio12' = ~ 1 + bio12,
					'bio1_bio12' = ~ 1 + bio1 + bio12,
					'bio1_x_bio12' = ~ 1 + bio1 + bio12 + bio1:bio12,
					'bio1^2_bio12' = ~ 1 + bio1 + bio12 + I(bio1^2),
					'bio1^2_x_bio12' = ~ 1 + bio1 + bio12 + I(bio1^2) + bio1:bio12,
					'bio12_sand' = ~ 1 + bio12 + sand,
					'bio12_x_sand' = ~ 1 + bio12 + sand + bio12:sand,
					'bio12_nitrogen' = ~ 1 + bio12 + site_nitrogen,
					'bio12_x_nitrogen' = ~ 1 + bio12 + site_nitrogen + bio12:site_nitrogen
				)

			}

			for (count_formula in seq_along(all_formulae)) {

				formula_biomass <- all_formulae[[count_formula]]
				preds_filename <- names(all_formulae)[count_formula]

				out_dir <- paste0('./outputs_loretta/integrated_sdm_pdm/models_biomass/', ifelse(trial, 'TRIAL_', ''), '[biomass_', tolower(resp_distrib), '~normal~', preds_filename, ']', ifelse(log_precip, '_log_precip', ''), '/')

				say(out_dir)

				if (trial | (!trial & !file.exists(out_dir))) {

					# do not change
					formula_psi <- formula_biomass

					if (file.exists(out_dir)) warning('Output folder already exists.')
					dirCreate(out_dir)

					sink(paste0(out_dir, '/runtime_log.txt'), split = TRUE)
					say('MODELING BIOMASS')
					say('Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | ', date(), post = 1)

					### data collation
					##################
					say('data collation', post = 1)

					say('This model estimates the distribution of site-level mean Andropogon gerardi ramet biomass. The distribution of biomasses among ramets at a site follows a zero-inflated gamma or zero-inflated lognormal distribution with individual ramet biomasses drawn from this distribution. The site-level mean biomass is a function of soil/climate. The probability of presence (non-zero) is a function of the same covariate(s) and has the same functional form as the submodel of the mean.', breaks = 60, post = 1)

					say('MCMC settings:', level = 2)
					say('trial ........................ ', trial)
					say('niter ........................ ', niter)
					say('nburnin ...................... ', nburnin)
					say('thin ......................... ', thin)
					say('nchains ...................... ', nchains)
					say('resp_distrib ................. ', resp_distrib)
					say('transform .................... ', transform)
					say('formula_biomass .............. ', paste(as.character(formula_biomass), collapse = ' '))
					say('formula_psi .................. ', paste(as.character(formula_psi), collapse = ' '))
					say('zero_inflated ................ ', zero_inflated)
					say('calib ........................ ', calib, post = 2)

					say('out_dir:')
					say(out_dir)

					formulae <- list(
						formula_biomass = formula_biomass,
						formula_psi = formula_psi
					)

					########################s
					### data preparation ###
					########################

					data_biomass <- prepare_biomass_data(formula_biomass = formula_biomass, log_precip = log_precip, n_response_curve_values = n_response_curve_values, calib = calib)

					#########################
					### inputs for nimble ###
					#########################

					say('Inputs:', level = 2)
					data <- list(
						y_biomass = data_biomass$y_biomass		# biomass of individual plants
					)

					constants <- list(
						
						# sampled sites
						n_pheno_sites = data_biomass$n_pheno_sites, # number of phenotype sample sites

						### biomass
						x_by_site_biomass = data_biomass$x_by_site_biomass, # MM with covariates for biomass (scaled)
						n_biomass = data_biomass$n_biomass, # number of biomass observations
						site_index_biomass = data_biomass$site_index_biomass, # index of sampled site for each row in biomass data
						n_terms_biomass = data_biomass$n_terms_biomass, # number of terms in formula for biomass model (including intercept)
						n_terms_psi = data_biomass$n_terms_biomass, # number of terms in formula for psi model (including intercept)

						n_covariates_biomass = data_biomass$n_covariates_biomass,
						resp_curves_x_biomass = data_biomass$resp_curves_x_biomass, # response curve array for biomass

						# response curves (general)
						n_response_curve_values = n_response_curve_values # number of values in response curve array

					)

					constants <- c(constants, constants_shared_biomass, constants_shared_psi)

					beta_biomass_inits <- rep(0, constants$n_terms_biomass)
					beta_psi_inits <- -1 * beta_biomass_inits

					inits <- list(

						log_site_biomass_mu = log(data_biomass$site_vect_biomass$biomass_mean),
						log_sigma_biomass_within_sites = 0.1,
						log_sigma_biomass_among_sites = 1,
						y_biomass_sim = data_biomass$y_biomass, # simulated values for biomass (for DHARMa residuals)
						beta_biomass = beta_biomass_inits,
						beta_psi = beta_psi_inits,
						z_site = rep(1, data_biomass$n_pheno_sites)

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
					
						# BIOMASS: priors for relationship of site-level mean biomass to environment
						for (i in 1:n_terms_biomass) {
							beta_biomass[i] ~ dnorm(0, sd = beta_biomass_prior_dnorm_sd) # broad prior
						}

						# BIOMASS: priors for probability of presence
						beta_psi[1] ~ dnorm(0, sd = beta_psi_prior_dnorm_sd_1)
						for (i in 2:n_terms_psi) {
							beta_psi[i] ~ ddexp(0, rate = beta_psi_prior_dnorm_sd)
						}

						# prior for sd of mean of biomass at a site on lognormal (~ half-Cauchy), ==> vague
						log(sigma_biomass_among_sites) ~ dnorm(0, sd = sigma_biomass_among_sites_log_prior_sd)

						# prior for sd of individual plant biomass on lognormal (~ half-Cauchy), ==> vague
						log(sigma_biomass_within_sites) ~ dnorm(0, sd = sigma_biomass_within_sites_log_prior_sd)

						# BIOMASS: parameters of biomass distribution are latent and functions of environment
						# individual plant biomasses are samples from the site-level distribution defined by the site-level distribution (next chunk after this one)
						for (i in 1:n_pheno_sites) {

							# relationship of biomass to the environment
							log_site_biomass_mu[i] ~ dnorm(log_site_biomass_latent_mu[i], sd = sigma_biomass_among_sites)
							log_site_biomass_latent_mu[i] <- inprod(beta_biomass[1:n_terms_biomass], x_by_site_biomass[i, 1:n_terms_biomass])
							
							# zero-inflation
							logit(psi[i]) <- inprod(beta_psi[1:n_terms_biomass], x_by_site_biomass[i, 1:n_terms_biomass])
							z_site[i] ~ dbern(psi[i])

							# site-level mean biomass
							mu_biomass_site[i] <- exp(log_site_biomass_mu[i])

						}

						log_lik <- sum(log_lik_biomass[1:n_biomass])

					})

					### response likelihood
					#######################
					if (resp_distrib == 'ZIG') { ### zero-inflated gamma
					
						resp_distrib_code <- nimbleCode({

							for (i in 1:n_pheno_sites) {

								# moment matching to get gamma() parameters
								shape_biomass[i] <- mu_biomass_site[i]^2 / sigma_biomass_within_sites^2
								rate_biomass[i] <- mu_biomass_site[i] / sigma_biomass_within_sites^2

							}

							# BIOMASS: likelihood of individual plants
							for (i in 1:n_biomass) {

								# likelihood
								y_biomass[i] ~ dZIG(shape = shape_biomass[site_index_biomass[i]], rate = rate_biomass[site_index_biomass[i]], z = z_site[site_index_biomass[i]])

								# simulated values for unconditional DHARMa residuals
								y_biomass_sim[i] ~ dZIG(shape = shape_biomass[site_index_biomass[i]], rate = rate_biomass[site_index_biomass[i]], z = z_site[site_index_biomass[i]])

								log_lik_biomass[i] <- dZIG(y_biomass[i], shape = shape_biomass[site_index_biomass[i]], rate = rate_biomass[site_index_biomass[i]], z = z_site[site_index_biomass[i]])
						
							}

						})
					
					} else if (resp_distrib == 'ZILN') { ### zero-inflated lognormal

						resp_distrib_code <- nimbleCode({

							# BIOMASS: likelihood of individual plants
							for (i in 1:n_biomass) {

								# likelihood
								y_biomass[i] ~ dZILN(meanlog = log_site_biomass_mu[site_index_biomass[i]], sdlog = sigma_biomass_within_sites, z = z_site[site_index_biomass[i]])

								# simulated values for unconditional DHARMa residuals
								y_biomass_sim[i] ~ dZILN(meanlog = log_site_biomass_mu[site_index_biomass[i]], sdlog = sigma_biomass_within_sites, z = z_site[site_index_biomass[i]])

								# log likelihood
								log_lik_biomass[i] <- dZILN(y_biomass[i], meanlog = log_site_biomass_mu[site_index_biomass[i]], sdlog = sigma_biomass_within_sites, z = z_site[site_index_biomass[i]])

							}

						})
					
					}
				
					### response curves
					###################
					if (data_biomass$n_covariates_biomass == 1) {

						### univariate
						response_curve_code <- nimbleCode({

							# BIOMASS: posterior predictive sampler for response curves: site-level mean
							# NB we assume just one predictor for biomass, so the response curve "x" is a matrix, not an array
							for (i in 1:n_response_curve_values) {
									
								log_response_curves_biomass_mu[i] <-
									inprod(beta_biomass[1:n_terms_biomass], resp_curves_x_biomass[i, 1:n_terms_biomass])
								response_curves_biomass_mu[i] <- exp(log_response_curves_biomass_mu[i])

								logit(response_curves_psi[i]) <-
									inprod(beta_psi[1:n_terms_biomass], resp_curves_x_biomass[i, 1:n_terms_biomass])
								
							}

						})
					
					} else {
					
						### multivariate
						response_curve_code <- nimbleCode({

							# BIOMASS: posterior predictive sampler for response curves: site-level mean
							# NB we assume just one predictor for biomass, so the response curve "x" is a matrix, not an array
							for (i in 1:n_covariates_biomass) {
								
								for (j in 1:n_response_curve_values) {
										
									log_response_curves_biomass_mu[j, i] <-
										inprod(beta_biomass[1:n_terms_biomass], resp_curves_x_biomass[j, 1:n_terms_biomass, i])
									response_curves_biomass_mu[j, i] <- exp(log_response_curves_biomass_mu[j, i])

									logit(response_curves_psi[j, i]) <-
										inprod(beta_psi[1:n_terms_biomass], resp_curves_x_biomass[j, 1:n_terms_biomass, i])
									
								}
							}

						})
					
					}

					model_code <- glueNimbleCode(model_code, resp_distrib_code, response_curve_code)

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
					)

					say('$initializeInfo() and $calculate():', level = 2)
					model$initializeInfo()
					calc <- model$calculate()
					say('model$calculate(): ', calc)
					if (is.na(calc) | is.infinite(calc)) stop('Likelihood is incalculable.')

					say('configureMCMC():', level = 2)

					# monitors for coefficients that have no indexing
					monitors_coeffs_not_indexed <- c(
						'sigma_biomass_within_sites', 'sigma_biomass_among_sites'
					)

					# coefficients that have bracketed indexing
					monitors_coeffs_single_index <- c(
						'beta_biomass', 'beta_psi'
					)

					monitors_coeffs_double_index <- c(
					)

					monitors_derived_not_indexed <- c(
						'log_lik'
					)

					monitors_derived_single_index <- c(
						# 'mu_biomass_site'
					)
					monitors_derived_double_index <- c(
					)

					monitors_geog_conus <- c(
					)

					monitors_dharma <- c(
						'y_biomass_sim'
					)

					monitors_resp_curves <- c(
						'response_curves_biomass_mu', 'response_curves_psi'
					)

					monitors <- c(monitors_coeffs_not_indexed, monitors_coeffs_single_index, monitors_coeffs_double_index, monitors_derived_not_indexed, monitors_derived_single_index, monitors_derived_double_index, monitors_dharma, monitors_resp_curves)

					conf <- configureMCMC(
						model,
						monitors = monitors,
						print = TRUE,
						enableWAIC = TRUE
					)

					vars <- 'beta_biomass'
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
				### post-modeling analysis for BIOMASS-only models ###
				######################################################

					### post-modeling analysis of BIOMASS
					descrip <- paste0('biomass ~ ', resp_distrib, '(', ifelse(resp_distrib == 'ZIG', 'exp(', ''), 'normal(env))', ifelse(resp_distrib == 'ZIG', ')', ''))

					workflow_postmodeling_generic(facet = 'biomass', formulae = formulae, descrip = descrip, out_dir = out_dir)

					workflow_postmodeling_biomass(chains = chains, descrip = descrip, formula_biomass = formula_biomass, formula_biomass_sigma = formula_biomass_sigma, formula_psi = formula_psi, resp_distrib = resp_distrib, transform = transform, log_precip = log_precip, crossvalidate = crossvalidate, out_dir = out_dir)

				} # output folder does not exist or we're doing a TRIAL model

			} # next formula

		} # log_precip

	} # next response distribution

say(date())
say('FINIS!', deco = '+', level = 1)

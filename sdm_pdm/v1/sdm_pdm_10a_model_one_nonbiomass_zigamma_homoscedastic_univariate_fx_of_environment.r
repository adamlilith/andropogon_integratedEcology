### MODELING ANDROPOGON GERARDI DISTRIBUTION, PHENOTYPE, PHYSIOLOGY, GENOTYPE, and ASSOCIATED MICROBIAL COMMUNITIES
### Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2023-12
###
### This model estimates the distribution of a SINGLE non-biomass trait's site-level mean value. The distribution of trait values among ramets at a site is assumed to follow a zero-inflated gamma distribution with individual ramet values drawn from this distribution. The site-level mean and the probability of a zero value are both functions of environmental covariates, and have the same functional form (but not coefficients). The model is "explicitly" homoscedastic because the standard deviation is not assumed to vary by site. The mean of the gamma is a draw from the a normal distribution with a mean given by the function of the environment. The function allows for ONE climatic and/or edaphic predictor (possibly with higher-order terms).
###
### NB Code in this script often uses "traits" (vs "trait") in variable names, but this is for ease in portability. This workflow only models one trait.
###
### source('C:/Kaji/R/andropogon_integratedEcology/sdm_pdm/sdm_pdm_10a_model_one_nonbiomass_zigamma_homoscedastic_univariate_fx_of_environment.r')

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

	### MCMC settings
	#################

		### MCMC settings FOR NON-BIOMASS TRAITS
		niter <- 20000
		nburnin <- 10000
		thin <- 20
		nchains <- 4

		# ### MCMC settings FOR TESTING
		# niter <- 1100
		# nburnin <- 100
		# thin <- 2
		# nchains <- 2

	### formulae for how aspects of species responds to environment
	formula_occs <- ~ 1 # not used except for prepare_occurrences()
	formula_occs_bias <- ~ 1 # not used except for prepare_occurrences()
	
	### traits
	traits <- c(
		'height',
		'blade_width',
		'leaf_thickness',
		'spad',
		'canopy_diameter',
		'photosynthetic_rate',
		'stomatal_conductance',
		'internal_co2',
		'transpiration_rate',
		'n_concentration',
		'cn_ratio'
	)

traits <- traits[5]

	### formulae
	# linear forms and univariate quadratic forms
	formulae_traits <- list(
	
		~ 1 + bio12,
		~ 1 + bio12 + I(bio12^2),
		~ 1 + bio18,
		~ 1 + bio18 + I(bio18^2),
		~ 1 + aridity,
		~ 1 + aridity + I(aridity^2),
		~ 1 + bio5,
		~ 1 + bio5 + I(bio5^2),
		~ 1 + bio10,
		~ 1 + bio10 + I(bio10^2),
		~ 1 + bio15,
		~ 1 + bio15 + I(bio15^2)
	
	)

formulae_traits <- formulae_traits[2]

	### main modeling loop
	######################
	# loop over each trait
	# loop over each valid formula

	for (trait in traits) {

		for (count_form in seq_along(formulae_traits)) {

			formula_traits <- formulae_traits[[count_form]]
			formula_traits_char <- paste(as.character(formula_traits), collapse = ' ')

			say(trait, ' ', formula_traits_char, level = 1)

			# create folder name
			terms <- gsub('~ 1 \\+ ', formula_traits_char, replacement = '')
			terms <- strsplit(terms, split = ' ')[[1]]
			terms <- terms[terms != '+']

			quads <- which(grepl(terms, pattern = '\\^2\\)'))
			if (length(quads) > 0) terms <- terms[-quads]

			if (length(quads) > 0) terms <- paste0(terms, '^2')
			terms <- sort(terms)
			terms <- paste(terms, collapse = ' ')

			out_dir <- paste0('./outputs_loretta/integrated_sdm_pdm/models_', trait, '/', trait, '_zigamma~normal_homoscedastic_', terms, ifelse(trial, '_TRIAL', ''), '/')
			if (!trial & file.exists(out_dir)) stop('Output folder already exists.')
			dirCreate(out_dir)
			
			### start log
			#############
			sink(paste0(out_dir, '/runtime_log.txt'), split = TRUE)
			say('MODELING ', trait, ' ', formula_traits_char)
			say('Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | ', date(), post = 1)

			say('data collation', post = 1)

			say('This model estimates the distribution of a SINGLE non-biomass trait\'s site-level mean value. The distribution of trait values among ramets at a site is assumed to follow a zero-inflated gamma distribution with individual ramet values drawn from this distribution. The site-level mean and the probability of a zero value are both functions of environmental covariates, and have the same functional form (but not coefficients). The model is "explicitly" homoscedastic because the standard deviation is not assumed to vary by site. The mean of the gamma is a draw from the a normal distribution with a mean given by the function of the environment. The function allows for ONE climatic and/or edaphic predictor (possibly with higher-order terms).', breaks = 60, post = 1)

			say('Settings:', level = 2)
			say('trait ........................ ', trait)
			say('niter ........................ ', niter)
			say('nburnin ...................... ', nburnin)
			say('thin ......................... ', thin)
			say('nchains ...................... ', nchains)
			say('formula_traits ............... ', paste(as.character(formula_traits), collapse = ' '))
			say('trial ........................ ', trial, post = 2)

			say('out_dir:')
			say(out_dir, post = 2)

			########################s
			### data preparation ###
			########################

			data_traits <- prepare_nonbiomass_traits(trait = trait, formula = formula_traits, n_response_curve_values = n_response_curve_values, calib = FALSE)
			# data_biomass <- prepare_biomass(formula = formula_biomass, n_response_curve_values = n_response_curve_values, calib = FALSE)

			data_occs <- prepare_occurrences(formula_occs = formula_occs, formula_occs_bias = formula_occs_bias, n_response_curve_values = n_response_curve_values, psa_quant = psa_quant, calib = FALSE)
			
			formulae <- list(
				formula_traits = formula_traits
			)

			#########################
			### inputs for nimble ###
			#########################

			say('Inputs:', level = 2)
			data <- list(
				y_traits = data_traits$y_trait		# 1-column matrix of trait values of individual plants
			)

			constants <- list(
				
				# sampled sites
				n_pheno_sites = data_traits$n_pheno_sites, # number of phenotype sample sites

				### traits
				x_by_site_traits = data_traits$x_by_site_traits, # MM with covariates for traits (scaled)
				# n_traits = data_traits$n_traits, # number of traits observations, not needed for ONE-trait model
				n_trait_values = data_traits$n_trait_values, # number of traits observations
				site_index_traits = data_traits$site_index_traits, # index of sampled site for each row in traits data
				n_terms_traits = data_traits$n_terms_traits, # number of terms in formula for traits model (including intercept)

				# n_covariates_traits = data_traits$n_covariates_traits, # do not need if just one predictor for traits
				resp_curves_x_traits = data_traits$resp_curves_x_traits, # response curve array for traits

				counties_x_traits_sq = data_traits$counties_x_traits_sq,
				counties_x_traits_ssp245_2041_2070 = data_traits$counties_x_traits_ssp245_2041_2070,
				counties_x_traits_ssp245_2071_2100 = data_traits$counties_x_traits_ssp245_2071_2100,
				counties_x_traits_ssp370_2041_2070 = data_traits$counties_x_traits_ssp370_2041_2070,
				counties_x_traits_ssp370_2071_2100 = data_traits$counties_x_traits_ssp370_2071_2100,

				# counties
				n_counties = data_occs$n_counties, # number of counties in the dataset

				# response curves (general)
				n_response_curve_values = n_response_curve_values # number of values in response curve array

			)

			y_mean <- mean(data_traits$y_traits)
			inits <- list(

				y_traits_sim = data_traits$y_traits, # simulated values for trait (for DHARMa residuals)

				log_site_traits_mu = rep(log(y_mean), data_traits$n_trait_values),

				traits_mu_county_sq = rep(y_mean, data_traits$n_counties),
				traits_mu_county_ssp245_2041_2070 = rep(y_mean, data_traits$n_counties),
				traits_mu_county_ssp245_2071_2100 = rep(y_mean, data_traits$n_counties),
				traits_mu_county_ssp370_2041_2070 = rep(y_mean, data_traits$n_counties),
				traits_mu_county_ssp370_2071_2100 = rep(y_mean, data_traits$n_counties),

				cross_site_traits_sigma_log = 3, # carefully selected so dgamma(..., log = 1) != 0
				phi_traits_sigma_log = 1,
				beta_traits_mu = rep(0, data_traits$n_terms_traits),
				beta_traits_pzero = rep(-2, data_traits$n_terms_traits)

			)

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
			
				# NON-BIOMASS TRAITS
				# prior for relationship of mean to environment
				for (i in 1:n_terms_traits) {
					beta_traits_mu[i] ~ dnorm(0, sd = 10) # broad prior
				}

				# prior for probability of trait having zero value
				for (i in 1:n_terms_traits) {
					beta_traits_pzero[i] ~ dnorm(0, sd = 100) # broad prior
				}

				# prior for sd of site-level trait values on lognormal (~ half-Cauchy), ==> vague
				cross_site_traits_sigma_log ~ dnorm(0, sd = 2.5)
				cross_site_traits_sigma <- exp(cross_site_traits_sigma_log)

				# prior for sd of normal distribution for site-level mean trait values
				phi_traits_sigma_log ~ dnorm(0, sd = 2.5)
				phi_traits_sigma <- exp(phi_traits_sigma_log)

				# NON-BIOMASS TRAITS: parameters of trait distribution are latent and functions of environment
				# individual plant trait values are samples from the site-level distribution (next chunk after this one)
				for (i in 1:n_pheno_sites) {

					# relationship of biomass to the environment
					phi_traits_mu[i] <- inprod(beta_traits_mu[1:n_terms_traits], x_by_site_traits[i, 1:n_terms_traits])
					log_site_traits_mu[i] ~ dnorm(phi_traits_mu[i], sd = phi_traits_sigma)

					# site-level mean biomass
					site_traits_mu[i] <- exp(log_site_traits_mu[i])

					# moment matching to get dgamma() parameters
					shape_traits[i] <- site_traits_mu[i]^2 / cross_site_traits_sigma^2
					rate_traits[i] <- site_traits_mu[i] / cross_site_traits_sigma^2

					logit(pzero_traits[i]) <- inprod(beta_traits_pzero[1:n_terms_traits], x_by_site_traits[i, 1:n_terms_traits])

				}

				# NON-BIOMASS TRAITS: likelihood of individual plants
				for (i in 1:n_trait_values) {

					# likelihood
					y_traits[i] ~ dzigamma(shape = shape_traits[site_index_traits[i]], rate = rate_traits[site_index_traits[i]], pzero = pzero_traits[site_index_traits[i]])

					# simulated values for unconditional DHARMa residuals
					y_traits_sim[i] ~ dzigamma(shape = shape_traits[site_index_traits[i]], rate = rate_traits[site_index_traits[i]], pzero = pzero_traits[site_index_traits[i]])
			
					log_lik_y_traits[i] <- dzigamma(y_traits[i], shape = shape_traits[site_index_traits[i]], rate = rate_traits[site_index_traits[i]], pzero = pzero_traits[site_index_traits[i]], log = 1)

				}

				log_lik <- sum(log_lik_y_traits[1:n_trait_values])

				# NON-BIOMASS TRAITS: posterior samplers for predictions of to counties in status quo and future
				for (i in 1:n_counties) {

					# non-biomass traits: status quo
					log_traits_county_mu_sq[i] <-
						inprod(beta_traits_mu[1:n_terms_traits], counties_x_traits_sq[i, 1:n_terms_traits])
					exp_traits_mu_county_sq[i] <- exp(log_traits_county_mu_sq[i])

					logit(pzero_traits_county_sq[i]) <- inprod(beta_traits_pzero[1:n_terms_traits], counties_x_traits_sq[i, 1:n_terms_traits])

					shape_traits_county_mu_sq[i] <- exp_traits_mu_county_sq[i]^2 / cross_site_traits_sigma^2
					rate_traits_county_mu_sq[i] <- exp_traits_mu_county_sq[i] / cross_site_traits_sigma^2
					traits_mu_county_sq[i] ~ dzigamma(shape = shape_traits_county_mu_sq[i], rate = rate_traits_county_mu_sq[i], pzero = pzero_traits_county_sq[i])

					# non-biomass traits: future ssp245_2041_2070
					log_traits_county_mu_ssp245_2041_2070[i] <-
						inprod(beta_traits_mu[1:n_terms_traits], counties_x_traits_ssp245_2041_2070[i, 1:n_terms_traits])
					exp_traits_mu_county_ssp245_2041_2070[i] <- exp(log_traits_county_mu_ssp245_2041_2070[i])
					
					logit(pzero_traits_county_ssp245_2041_2070[i]) <- inprod(beta_traits_pzero[1:n_terms_traits], counties_x_traits_ssp245_2041_2070[i, 1:n_terms_traits])

					shape_traits_county_mu_ssp245_2041_2070[i] <- exp_traits_mu_county_ssp245_2041_2070[i]^2 / cross_site_traits_sigma^2
					rate_traits_county_mu_ssp245_2041_2070[i] <- exp_traits_mu_county_ssp245_2041_2070[i] / cross_site_traits_sigma^2
					traits_mu_county_ssp245_2041_2070[i] ~ dzigamma(shape = shape_traits_county_mu_ssp245_2041_2070[i], rate = rate_traits_county_mu_ssp245_2041_2070[i], pzero = pzero_traits_county_ssp245_2041_2070[i])

					# non-biomass traits: future ssp245_2071_2100
					log_traits_county_mu_ssp245_2071_2100[i] <-
						inprod(beta_traits_mu[1:n_terms_traits], counties_x_traits_ssp245_2071_2100[i, 1:n_terms_traits])
					exp_traits_mu_county_ssp245_2071_2100[i] <- exp(log_traits_county_mu_ssp245_2071_2100[i])
					
					logit(pzero_traits_county_ssp245_2071_2100[i]) <- inprod(beta_traits_pzero[1:n_terms_traits], counties_x_traits_ssp245_2071_2100[i, 1:n_terms_traits])

					shape_traits_county_mu_ssp245_2071_2100[i] <- exp_traits_mu_county_ssp245_2071_2100[i]^2 / cross_site_traits_sigma^2
					rate_traits_county_mu_ssp245_2071_2100[i] <- exp_traits_mu_county_ssp245_2071_2100[i] / cross_site_traits_sigma^2
					traits_mu_county_ssp245_2071_2100[i] ~ dzigamma(shape = shape_traits_county_mu_ssp245_2071_2100[i], rate = rate_traits_county_mu_ssp245_2071_2100[i], pzero = pzero_traits_county_ssp245_2071_2100[i])

					# non-biomass traits: future ssp370_2041_2070
					log_traits_county_mu_ssp370_2041_2070[i] <-
						inprod(beta_traits_mu[1:n_terms_traits], counties_x_traits_ssp370_2041_2070[i, 1:n_terms_traits])
					exp_traits_mu_county_ssp370_2041_2070[i] <- exp(log_traits_county_mu_ssp370_2041_2070[i])

					logit(pzero_traits_county_ssp370_2041_2070[i]) <- inprod(beta_traits_pzero[1:n_terms_traits], counties_x_traits_ssp370_2041_2070[i, 1:n_terms_traits])

					shape_traits_county_mu_ssp370_2041_2070[i] <- exp_traits_mu_county_ssp370_2041_2070[i]^2 / cross_site_traits_sigma^2
					rate_traits_county_mu_ssp370_2041_2070[i] <- exp_traits_mu_county_ssp370_2041_2070[i] / cross_site_traits_sigma^2
					traits_mu_county_ssp370_2041_2070[i] ~ dzigamma(shape = shape_traits_county_mu_ssp370_2041_2070[i], rate = rate_traits_county_mu_ssp370_2041_2070[i], pzero = pzero_traits_county_ssp370_2041_2070[i])

					# non-biomass traits: future ssp370_2071_2100
					log_traits_county_mu_ssp370_2071_2100[i] <-
						inprod(beta_traits_mu[1:n_terms_traits], counties_x_traits_ssp370_2071_2100[i, 1:n_terms_traits])
					exp_traits_mu_county_ssp370_2071_2100[i] <- exp(log_traits_county_mu_ssp370_2071_2100[i])

					logit(pzero_traits_county_ssp370_2071_2100[i]) <- inprod(beta_traits_pzero[1:n_terms_traits], counties_x_traits_ssp370_2071_2100[i, 1:n_terms_traits])

					shape_traits_county_mu_ssp370_2071_2100[i] <- exp_traits_mu_county_ssp370_2071_2100[i]^2 / cross_site_traits_sigma^2
					rate_traits_county_mu_ssp370_2071_2100[i] <- exp_traits_mu_county_ssp370_2071_2100[i] / cross_site_traits_sigma^2
					traits_mu_county_ssp370_2071_2100[i] ~ dzigamma(shape = shape_traits_county_mu_ssp370_2071_2100[i], rate = rate_traits_county_mu_ssp370_2071_2100[i], pzero = pzero_traits_county_ssp370_2071_2100[i])

				}

				# NON-BIOMASS TRAITS: posterior predictive sampler for response curves: site-level mean
				# NB we assume just one predictor, so the response curve "x" is a matrix, not an array
				for (i in 1:n_response_curve_values) {
						
					log_response_curves_traits_mu[i] <-
						inprod(beta_traits_mu[1:n_terms_traits], resp_curves_x_traits[i, 1:n_terms_traits])
					
					response_curves_traits_mu[i] <- exp(log_response_curves_traits_mu[i])
					
				}

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
				'cross_site_traits_sigma', 'phi_traits_sigma'
			)

			# coefficients that have bracketed indexing
			monitors_coeffs_single_index <- c(
				'beta_traits_mu', 'beta_traits_pzero'
			)
			monitors_coeffs_double_index <- c(
			)

			monitors_derived_not_indexed <- c(
				'log_lik'
			)
			monitors_derived_single_index <- c(
				'site_traits_mu', 'pzero_traits'
			)
			monitors_derived_double_index <- c(
			)

			monitors_geog <- c(
				'traits_mu_county_sq', 'traits_mu_county_ssp245_2041_2070', 	
				'traits_mu_county_ssp245_2071_2100', 'traits_mu_county_ssp370_2041_2070', 'traits_mu_county_ssp370_2071_2100',
				'pzero_traits_county_sq', 'pzero_traits_county_ssp245_2041_2070', 	
				'pzero_traits_county_ssp245_2071_2100', 'pzero_traits_county_ssp370_2041_2070', 'pzero_traits_county_ssp370_2071_2100'
			)

			monitors_dharma <- c(
				'y_traits_sim'
			)

			monitors_resp_curves <- c(
				'response_curves_traits_mu'
			)

			monitors <- c(monitors_coeffs_not_indexed, monitors_coeffs_single_index, monitors_coeffs_double_index, monitors_derived_not_indexed, monitors_derived_single_index, monitors_derived_double_index, monitors_geog, monitors_dharma, monitors_resp_curves)

			conf <- configureMCMC(
				model,
				monitors = monitors,
				print = TRUE,
				enableWAIC = TRUE
			)

			# # add no U-turn sampler (Hamiltonian Monte Carlo)
			# conf$addSampler(target = monitors_coeffs, type = 'NUTS')
			# say('NUTS sampler added to all continuous parameters.')

			# conf$removeSamplers('beta_occs_vs_biomass')
			# conf$addSampler(target = 'beta_occs_vs_biomass', type = 'AF_slice')
			# say('AF slice sampler added to beta_occs_vs_biomass.')

			vars <- c('beta_traits_mu', 'beta_traits_pzero')
			conf$removeSamplers(vars)
			conf$addSampler(target = vars, type = 'AF_slice')
			say('AF slice sampler added to ', paste(vars, collapse = ' & '))

			unsampleds <- conf$getUnsampledNodes()
			if (length(unsampleds) > 0) {
				say('Unsampled nodes:')
				say(unsampleds)
				stop('Unsampled nodes!')
			}

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

			### collate all predictions into one SpatVector
			###############################################

				pred_vect_nam <- data_occs$ag_vect_sq

				for (var in monitors_geog) {
				
					preds <- hammer_extract(chains, param = var, j = TRUE, stat = 'mean')
					pred_vect_nam[[var]] <- preds
					if (grepl(var, pattern = 'traits_mu_county_')) {
						this_var <- sub(var, pattern = 'traits_mu_county_', replacement = paste0(trait, '_mu_county_'))
					} else if (grepl(var, pattern = 'pzero_traits_county_')) {
						this_var <- sub(var, pattern = 'pzero_traits_county_', replacement = paste0(trait, '_pzero_county_'))
					}
					names(pred_vect_nam)[ncol(pred_vect_nam)] <- this_var

				}

				writeVector(pred_vect_nam, paste0(out_dir, '/prediction_vector.gpkg'), overwrite = TRUE)

			### post-modeling analysis of BIOMASS
			descrip <- paste0(trait, ': zigamma homoscedastic')
			workflow_postmodeling_generic(facet = trait, formulae = formulae, descrip = descrip, out_dir = out_dir)
			workflow_postmodeling_nonbiomass_single_trait(trait = trait, formula_traits = formula_traits, descrip = descrip, homoscedastic = TRUE, pred_vect_nam = pred_vect_nam, out_dir = out_dir)

	} # next formula

} # next trait

say(date())
say('FINIS!', deco = '+', level = 1)

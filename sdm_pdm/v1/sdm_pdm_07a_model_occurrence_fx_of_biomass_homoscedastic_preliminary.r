### MODELING ANDROPOGON GERARDI DISTRIBUTION, PHENOTYPE, PHYSIOLOGY, GENOTYPE, and ASSOCIATED MICROBIAL COMMUNITIES
### Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2023-12
###
### This model estimates the distribution of site-level mean Andropogon gerardi ramet biomass. The distribution of biomasses among ramets at a site follows a gamma distribution with individual ramet biomasses drawn from this distribution. The site-level mean biomass is a function of soil/climate. The mean is assumed to respond to one climatic predictor (possibly with higher-order terms).
###
### source('C:/Kaji/R/andropogon_integratedEcology/sdm_pdm/sdm_pdm_07_model_occurrence_fx_of_biomass_homoscedastic.r')

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

	# calib <- TRUE # use just counties with non-NA Poaceae for calibration region
	calib <- FALSE # use all of North America for calibration region

	out_dir <- paste0('./outputs_loretta/integrated_sdm_pdm/models_occurrence_fx_of_biomass/occurrence_as_fx_of_biomass_homoscedastic', ifelse(trial, '_TRIAL', ''), '/')
	if (!trial & file.exists(out_dir)) stop('Output folder already exists.')

	# formula for how aspects of species responds to environment
	formula_biomass_mu <- ~ 1 + bio12 # response of biomass to environment
	formula_occs <- ~ 1
	formula_occs_vs_biomass <- ~ 1 + biomass + I(biomass^2) # response of occurrence to biomass

	### MCMC settings
	niter <- 240000
	nburnin <- 40000
	thin <- 200
	nchains <- 4
	waic <- TRUE

	# ### MCMC settings
	# niter <- 120000
	# nburnin <- 20000
	# thin <- 100
	# nchains <- 4
	# waic <- TRUE

	# # ### MCMC settings FOR TESTING
	# niter <- 2200
	# nburnin <- 200
	# thin <- 2
	# nchains <- 2
	# waic <- TRUE

#############
### model ###
#############

	dirCreate(out_dir)
	sink(paste0(out_dir, '/runtime_log.txt'), split = TRUE)
	say('MODELING OCCURRENCE AS A FUNCTION OF BIOMASS')
	say('Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | ', date(), post = 1)

	formula <- list(
		formula_biomass_mu = formula_biomass_mu,
		formula_occs_vs_biomass = formula_occs_vs_biomass
	)
	saveRDS(formula, paste0(out_dir, '/formula.rds'))

	### data collation
	##################
	say('data collation', post = 1)

	say('This model estimates the distribution of Andropogon gerardi abundance as a function of biomass. Biomass is modeled as the site-level mean ramet biomass. The distribution of biomasses among ramets at a site follows a gamma distribution with individual ramet biomasses drawn from this distribution. The site-level mean biomass is a function of soil/climate. The mean is assumed to respond to one climatic predictor (possibly with higher-order terms).', breaks = 60, post = 1)

	say('MCMC settings:', level = 2)
	say('niter ........................ ', niter)
	say('nburnin ...................... ', nburnin)
	say('thin ......................... ', thin)
	say('nchains ...................... ', nchains)
	say('formula_biomass_mu ........... ', paste(as.character(formula_biomass_mu), collapse = ' '))
	say('formula_occs_vs_biomass ...... ', paste(as.character(formula_occs_vs_biomass), collapse = ' '))
	say('trial ........................ ', trial, post = 2)

	########################s
	### data preparation ###
	########################

	data_biomass_mu <- prepare_biomass(formula = formula_biomass_mu, n_response_curve_values = n_response_curve_values, calib = calib)
	data_occs <- prepare_occurrences(formula = formula_occs, n_response_curve_values = n_response_curve_values, psa_quant = psa_quant, calib = calib)
	data_response_curves_occs_vs_biomass <- create_response_curve_array_occ_vs_biomass(formula = formula_occs_vs_biomass, data_biomass = data_biomass_mu, n_response_curve_values = n_response_curve_values)

	terms_occs_vs_biomass <- attr(terms(formula_occs_vs_biomass), 'term.labels')
	n_terms_occs_vs_biomass <- length(terms_occs_vs_biomass) + 1 # + 1 for intercept

	#########################
	### inputs for nimble ###
	#########################

	say('Inputs:', level = 2)
	data <- list(
		y_biomass = data_biomass_mu$y_biomass,	# biomass of individual plants
		y_n_ag = data_occs$y_n_ag				# number of AG observations in each county
	)

	constants <- list(
		
		### occurrences
		n_counties_occs_calib = data_occs$n_counties_occs_calib, # number of counties in calibration region
		area_km2_log10 = data_occs$area_km2_log10, # log10(area_km2) for each county
		n_poaceae_log10p1 = data_occs$n_poaceae_log10p1, # log10(n_poaceae + 1) for each county

		index_biomass_county_in_calib_county = data_biomass_mu$index_biomass_county_in_calib_county, # index of biomass county in calibration region

		response_curve_x_occs_vs_biomass = data_response_curves_occs_vs_biomass$response_curves_x_occs_vs_biomass, # response curve array for 

		### integration of occurrence ~ fx(biomass)
		site_biomass_mean = data_biomass_mu$site_biomass_mean, # mean biomass for scaling
		site_biomass_sd = data_biomass_mu$site_biomass_sd, # sd of biomass for scaling
		n_terms_occs_vs_biomass = n_terms_occs_vs_biomass, # number of terms in abundance ~ f(biomass)

		### general phenotypic traits
		n_pheno_sites = data_biomass_mu$n_pheno_sites, # number of phenotype sample sites

		### biomass
		x_by_site_biomass = data_biomass_mu$x_by_site_biomass, # MM with covariates for biomass (scaled)
		n_biomass = data_biomass_mu$n_biomass, # number of biomass observations
		site_index_biomass = data_biomass_mu$site_index_biomass, # index of sampled site for each row in biomass data
		n_terms_biomass_mu = data_biomass_mu$n_terms_biomass, # number of terms in formula for biomass model (including intercept)

		resp_curves_x_biomass = data_biomass_mu$resp_curves_x_biomass, # response curve array for biomass

		counties_x_biomass_sq = data_biomass_mu$counties_x_biomass_sq,
		counties_x_biomass_ssp245_2041_2070 = data_biomass_mu$counties_x_biomass_ssp245_2041_2070,
		counties_x_biomass_ssp245_2071_2100 = data_biomass_mu$counties_x_biomass_ssp245_2071_2100,
		counties_x_biomass_ssp370_2041_2070 = data_biomass_mu$counties_x_biomass_ssp370_2041_2070,
		counties_x_biomass_ssp370_2071_2100 = data_biomass_mu$counties_x_biomass_ssp370_2071_2100,

		### general
		n_counties = data_occs$n_counties, # number of counties in the dataset
		n_response_curve_values = n_response_curve_values # number of values in response curve array

	)

	# initial values for occurrences
	N_inits_calib <- data_occs$y_n_ag * 100 + 10
	N_inits_all_counties <- data_occs$ag_vect_sq$n_andropogon_gerardi * 100 + 10

	# initial values for biomass ~ fx(env)
	beta_biomass_mu_inits <- rep(0, constants$n_terms_biomass_mu)

	inits <- list(

		y_n_ag_sim = data_occs$y_n_ag, # simulated values of observed number of AG (for DHARMa residuals)
		alpha_occs = c(-1, 1, 1), # intercept, area, # of Poaceae
		beta_occs_vs_biomass = c(1, -5, -5), # occurrence ~ biomass coefficients

		N = N_inits_calib, # number of latent AG in calibration counties
		N_ag_county_sq = N_inits_all_counties, # number of latent AG in all counties
		N_ag_county_ssp245_2041_2070 = N_inits_all_counties,
		N_ag_county_ssp245_2071_2100 = N_inits_all_counties,
		N_ag_county_ssp370_2041_2070 = N_inits_all_counties,
		N_ag_county_ssp370_2071_2100 = N_inits_all_counties,

		y_biomass_sim = data_biomass_mu$y_biomass, # simulated values for biomass (for DHARMa residuals)

		site_biomass_plant_sigma_log = 3.8, # from preliminary runs
		beta_biomass_mu = beta_biomass_mu_inits

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
	# biomass: gamma distribution with mean a functions of ONE environmental predictor... if more, then need an n-dimensional response array
	model_code <- nimbleCode({
	
		### BIOMASS
		###########

			# BIOMASS: priors for relationship of mean and variance to environment
			for (i in 1:n_terms_biomass_mu) {
				beta_biomass_mu[i] ~ dnorm(0, sd = 10) # broad prior
			}

			# prior for sd of site-level biomass on lognormal (~ half-Cauchy), ==> vague
			site_biomass_plant_sigma_log ~ dnorm(0, sd = 2.5)
			site_biomass_plant_sigma_log <- exp(site_biomass_plant_sigma_log)


			# BIOMASS: parameters of biomass distribution are latent and functions of environment
			# individual plant biomasses are samples from the site-level distribution defined by the site-level distribution (next chunk after this one)
			for (i in 1:n_pheno_sites) {

				# log of biomass and shape parameter of the gamma distribution
				log_site_biomass_mu[i] <- inprod(beta_biomass_mu[1:n_terms_biomass_mu], x_by_site_biomass[i, 1:n_terms_biomass_mu])
				
				# site-level mean and s.d. of biomass
				# s.d. is just estimated for comparison to observed s.d
				mu_biomass_site[i] <- exp(log_site_biomass_mu[i])

				# moment matching to get dgamma() parameters
				shape_biomass[i] <- mu_biomass_site[i]^2 / site_biomass_plant_sigma_log^2
				rate_biomass[i] <- mu_biomass_site[i] / site_biomass_plant_sigma_log^2

			}

			# BIOMASS: likelihood of individual plants
			for (i in 1:n_biomass) {

				# likelihood
				y_biomass[i] ~ dgamma(shape = shape_biomass[site_index_biomass[i]], rate = rate_biomass[site_index_biomass[i]])

				# simulated values for unconditional DHARMa residuals
				y_biomass_sim[i] ~ dgamma(shape = shape_biomass[site_index_biomass[i]], rate = rate_biomass[site_index_biomass[i]])
		
			}

			# BIOMASS: posterior samplers for predictions of to counties in status quo and future
			for (i in 1:n_counties) {

				# biomass: status quo
				log_biomass_county_mu_sq[i] <-
					inprod(beta_biomass_mu[1:n_terms_biomass_mu], counties_x_biomass_sq[i, 1:n_terms_biomass_mu])
				mu_biomass_county_sq[i] <- exp(log_biomass_county_mu_sq[i])

				# biomass: future
				log_biomass_county_mu_ssp245_2041_2070[i] <-
					inprod(beta_biomass_mu[1:n_terms_biomass_mu], counties_x_biomass_ssp245_2041_2070[i, 1:n_terms_biomass_mu])
				mu_biomass_county_ssp245_2041_2070[i] <- exp(log_biomass_county_mu_ssp245_2041_2070[i])
				
				log_biomass_county_mu_ssp245_2071_2100[i] <-
					inprod(beta_biomass_mu[1:n_terms_biomass_mu], counties_x_biomass_ssp245_2071_2100[i, 1:n_terms_biomass_mu])
				mu_biomass_county_ssp245_2071_2100[i] <- exp(log_biomass_county_mu_ssp245_2071_2100[i])
				
				log_biomass_county_mu_ssp370_2041_2070[i] <-
					inprod(beta_biomass_mu[1:n_terms_biomass_mu], counties_x_biomass_ssp370_2041_2070[i, 1:n_terms_biomass_mu])
				mu_biomass_county_ssp370_2041_2070[i] <- exp(log_biomass_county_mu_ssp370_2041_2070[i])

				log_biomass_county_mu_ssp370_2071_2100[i] <-
					inprod(beta_biomass_mu[1:n_terms_biomass_mu], counties_x_biomass_ssp370_2071_2100[i, 1:n_terms_biomass_mu])
				mu_biomass_county_ssp370_2071_2100[i] <- exp(log_biomass_county_mu_ssp370_2071_2100[i])

			}

			# BIOMASS: posterior predictive sampler for response curves: site-level mean
			# NB we assume just one predictor for biomass, so the response curve "x" is a matrix, not an array
			for (i in 1:n_response_curve_values) {
					
				log_response_curves_biomass_mu[i] <-
					inprod(beta_biomass_mu[1:n_terms_biomass_mu], resp_curves_x_biomass[i, 1:n_terms_biomass_mu])
				
				response_curves_biomass_mu[i] <- exp(log_response_curves_biomass_mu[i])
				
			}

		### OCCURRENCE
		##############

			# OCCURRENCE: priors for relationship to biomass
			for (i in 1:n_terms_occs_vs_biomass) {
				beta_occs_vs_biomass[i] ~ dnorm(0, sd = 10)
			}

			# OCCURRENCE: priors for sampling bias
			for (i in 1:3) {
				alpha_occs[i] ~ dnorm(0, sd = 10) # 
			}

			# OCCURRENCE: likelihood
			for (i in 1:n_counties_occs_calib) {
				
				### actual abundance (latent--unobserved)
				N[i] ~ dpois(lambda_sq[i])

				### observed number of AG
				y_n_ag[i] ~ dbin(prob = p[i], size = N[i])
				y_n_ag_sim[i] ~ dbin(prob = p[i], size = N[i]) # simulate values for DHARMa residuals

				# sampling bias
				logit(p[i]) <- alpha_occs[1] + alpha_occs[2] * area_km2_log10[i] + alpha_occs[3] * n_poaceae_log10p1[i]

				# scaled, site-level mean biomass for this county
				biomass_county_mu_sq_calib_scaled[i] <-
					(mu_biomass_county_sq[index_biomass_county_in_calib_county[i]] - site_biomass_mean) / site_biomass_sd

				# relationship between expected (latent) abundance and biomass
				log(lambda_sq[i]) <-
					beta_occs_vs_biomass[1] +
					beta_occs_vs_biomass[2] * biomass_county_mu_sq_calib_scaled[i] +
					beta_occs_vs_biomass[3] * biomass_county_mu_sq_calib_scaled[i]^2

			}

			# OCCURRENCE: posterior samplers for geographic predictions
			for (i in 1:n_counties) {

				N_ag_county_sq[i] ~ dpois(county_lambda_sq[i])
				log(county_lambda_sq[i]) <-
					beta_occs_vs_biomass[1] +
					beta_occs_vs_biomass[2] * biomass_county_mu_sq_scaled[i] +
					beta_occs_vs_biomass[3] * biomass_county_mu_sq_scaled[i]^2
				biomass_county_mu_sq_scaled[i] <-
					(mu_biomass_county_sq[i] - site_biomass_mean) / site_biomass_sd

				N_ag_county_ssp245_2041_2070[i] ~ dpois(county_lambda_ssp245_2041_2070[i])
				log(county_lambda_ssp245_2041_2070[i]) <-
					beta_occs_vs_biomass[1] +
					beta_occs_vs_biomass[2] * biomass_county_mu_ssp245_2041_2070_scaled[i] +
					beta_occs_vs_biomass[3] * biomass_county_mu_ssp245_2041_2070_scaled[i]^2
				biomass_county_mu_ssp245_2041_2070_scaled[i] <-
					(mu_biomass_county_ssp245_2041_2070[i] - site_biomass_mean) / site_biomass_sd

				N_ag_county_ssp245_2071_2100[i] ~ dpois(county_lambda_ssp245_2071_2100[i])
				log(county_lambda_ssp245_2071_2100[i]) <-
					beta_occs_vs_biomass[1] +
					beta_occs_vs_biomass[2] * biomass_county_mu_ssp245_2071_2100_scaled[i] +
					beta_occs_vs_biomass[3] * biomass_county_mu_ssp245_2071_2100_scaled[i]^2
				biomass_county_mu_ssp245_2071_2100_scaled[i] <-
					(mu_biomass_county_ssp245_2071_2100[i] - site_biomass_mean) / site_biomass_sd

				N_ag_county_ssp370_2041_2070[i] ~ dpois(county_lambda_ssp370_2041_2070[i])
				log(county_lambda_ssp370_2041_2070[i]) <-
					beta_occs_vs_biomass[1] +
					beta_occs_vs_biomass[2] * biomass_county_mu_ssp370_2041_2070_scaled[i] +
					beta_occs_vs_biomass[3] * biomass_county_mu_ssp370_2041_2070_scaled[i]^2
				biomass_county_mu_ssp370_2041_2070_scaled[i] <-
					(mu_biomass_county_ssp370_2041_2070[i] - site_biomass_mean) / site_biomass_sd

				N_ag_county_ssp370_2071_2100[i] ~ dpois(county_lambda_ssp370_2071_2100[i])
				log(county_lambda_ssp370_2071_2100[i]) <-
					beta_occs_vs_biomass[1] +
					beta_occs_vs_biomass[2] * biomass_county_mu_ssp370_2071_2100_scaled[i] +
					beta_occs_vs_biomass[3] * biomass_county_mu_ssp370_2071_2100_scaled[i]^2
				biomass_county_mu_ssp370_2071_2100_scaled[i] <-
					(mu_biomass_county_ssp370_2071_2100[i] - site_biomass_mean) / site_biomass_sd

			}

			# OCCURRENCE: posterior predictive sampler for BIOMASS response curves
			for (i in 1:n_response_curve_values) {
				
				log(response_curve_occs_vs_biomass[i]) <-
					inprod(beta_occs_vs_biomass[1:n_terms_occs_vs_biomass], response_curve_x_occs_vs_biomass[i, 1:n_terms_occs_vs_biomass])

				# exp_response_curve_occs_vs_biomass[i] <- exp(response_curve_occs_vs_biomass[i])

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
	model$initializeInfo()
	calc <- model$calculate()
	say('model$calculate(): ', calc)
	if (is.na(calc)) stop('Likelihood is incalculable.')

	say('configureMCMC():', level = 2)

	# monitors for coefficients that have no indexing
	monitors_coeffs_not_indexed <- c(
		'site_biomass_plant_sigma_log'
	)

	# coefficients that have bracketed indexing
	monitors_coeffs_indexed <- c(
		'beta_biomass_mu',
		'beta_occs_vs_biomass', 'alpha_occs'
	)

	monitors_derived_not_indexed <- c(
	)

	monitors_derived_indexed <- c(
		'mu_biomass_site'
	)

	monitors_geog <- c(
		'mu_biomass_county_sq', 'mu_biomass_county_ssp245_2041_2070', 	
		'mu_biomass_county_ssp245_2071_2100', 'mu_biomass_county_ssp370_2041_2070', 'mu_biomass_county_ssp370_2071_2100',

		'N_ag_county_sq', 'N_ag_county_ssp245_2041_2070', 'N_ag_county_ssp245_2071_2100', 'N_ag_county_ssp370_2041_2070', 
		'N_ag_county_ssp370_2071_2100'

	)

	monitors_dharma <- c(
		'y_biomass_sim',
		'y_n_ag_sim', 'lambda_sq'
	)

	monitors_resp_curves <- c(
		'response_curves_biomass_mu',
		'response_curve_occs_vs_biomass'
	)

	monitors <- c(monitors_coeffs_not_indexed, monitors_coeffs_indexed, monitors_derived_not_indexed, monitors_derived_indexed, monitors_geog, monitors_dharma, monitors_resp_curves)

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

	conf$removeSamplers('beta_occs_vs_biomass')
	conf$addSampler(target = 'beta_occs_vs_biomass', type = 'AF_slice')
	say('AF slice sampler added to beta_occs_vs_biomass.')

	conf$removeSamplers('beta_biomass_mu')
	conf$addSampler(target = c('beta_biomass_mu'), type = 'AF_slice')
	say('AF slice sampler added to beta_biomass_mu.')

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

#################################################
### post-modeling analysis for BIOMASS models ###
#################################################

	### collate all predictions into one SpatVector
	###############################################

		pred_vect_nam <- data_occs$ag_vect_sq

		for (var in monitors_geog) {
		
			preds <- hammer_extract(chains, param = var, j = TRUE, stat = 'mean')
			pred_vect_nam[[var]] <- preds
			names(pred_vect_nam)[ncol(pred_vect_nam)] <- var

		}

		writeVector(pred_vect_nam, paste0(out_dir, '/prediction_vector.gpkg'), overwrite = TRUE)

	### post-modeling analysis of BIOMASS
	workflow_postmodeling_biomass(homoscedastic = TRUE, pred_vect_nam = pred_vect_nam)

say('###############################')
say('### trace and density plots ###')
say('###############################')

	mcmc <- hammer_subset(chains, 'beta_occs_vs_biomass', j = TRUE)
	ggs_mcmc <- ggs(mcmc$samples)

	# OCCURRENCE ~ f(BIOMASS): graphing trace and density plots for all beta
	pars <- 'beta_occs_vs_biomass'
	trace <- ggs_traceplot(ggs_mcmc, family = pars)
	density <- ggs_density(ggs_mcmc, family = pars)

	combo <- trace + density

	file <- paste0(out_dir, '/beta_occs_vs_biomass_density_trace.png')
	ggsave(combo, file = file, width = 12, height = 8, dpi = 450, bg = 'white')

say('#############################')
say('### effective sample size ###')
say('#############################')

	mcmc <- chains$samples
	cols <- paste0('beta_occs_vs_biomass[', 1:n_terms_occs_vs_biomass, ']')
	if (homoscedastic) {
		cols <- c(cols, 'site_biomass_plant_sigma_log')
	} else {
		cols <- c(cols, paste0('beta_biomass_sigma[', 1:n_terms_biomass_sigma, ']'))
	}
	for (i in 1:nchains) {
		mcmc[[i]] <- mcmc[[i]][ , cols]
	}

	sink(paste0(out_dir, '/effective_sample_size_occurrence_fx_biomass.txt'), split = TRUE)
		say('EFFECTIVE SAMPLE SIZE', post = 2)
		print(effectiveSize(mcmc))
		say('')
	sink()

say('##############################################')
say('### response curves: OCCURRENCE vs BIOMASS ###')
say('##############################################')

	responses <- graph_response_curves_occurrence_vs_biomass(
		out_dir = out_dir,
		chains = chains,
		data_response_curves_occs_vs_biomass = data_response_curves_occs_vs_biomass,
		data_biomass = NULL,
		quant_threshold = 0.5
	)

say('########################')
say('### cross validation ###')
say('########################')

	if (!trial) {

		# CV for biomass and occurrence
		folds_fx <- folds_for_biomass_and_occurrences # change according to response data type we're using

		cv <- runCrossValidate(
			MCMCconfiguration = conf,
			k = k_folds, # universal setting
			foldFunction = folds_fx,
			lossFunction = 'MSE',
			MCMCcontrol = list(niter = niter, nburnin = nburnin),
			returnSamples = FALSE,
			nCores = 1,
			nBootReps = 200,
			silent = FALSE
		)

		sink(paste0(out_dir, '/cross_validation.txt'), split = TRUE)
		say('GEO-FOLD CROSS VALIDATION')
		say(date(), post = 2)
		say('Number of geo-folds: ', k_folds, post = 2)
		say('MSE for combined biomass and occurrence: ', k_folds, post = 2)
		print(cv)

		# CV for BIOMASS and occurrence
		folds_fx <- folds_for_biomass_and_occurrences # change according to response data type we're using

		cv <- runCrossValidate(
			MCMCconfiguration = conf,
			k = k_folds, # universal setting
			foldFunction = folds_fx,
			lossFunction = 'MSE',
			MCMCcontrol = list(niter = niter, nburnin = nburnin),
			returnSamples = FALSE,
			nCores = 1,
			nBootReps = 200,
			silent = FALSE
		)

		sink(paste0(out_dir, '/cross_validation.txt'), split = TRUE)
		say('GEO-FOLD CROSS VALIDATION')
		say(date(), post = 2)
		say('Number of geo-folds: ', k_folds, post = 2)
		say('MSE for combined biomass and occurrence: ', k_folds, post = 2)
		print(cv)



		sink()

	}

say(date())
say('FINIS!', deco = '+', level = 1)

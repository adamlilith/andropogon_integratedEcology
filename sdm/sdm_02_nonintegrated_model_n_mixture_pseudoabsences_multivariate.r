### MODELING ANDROPOGON GERARDI DISTRIBUTION, PHENOTYPE, PHYSIOLOGY, GENOTYPE, and ASSOCIATED MICROBIAL COMMUNITIES
### Erica Newman | Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2023-12
###
### This script constructs a non-integrated model for AG geographic distribution.
###
### source('C:/Ecology/R/andropogon_integratedEcology/sdm/sdm_02_nonintegrated_model_n_mixture_pseudoabsences_multivariate.r')
### source('C:/Subarashi/R/andropogon_integratedEcology/sdm/sdm_02_nonintegrated_model_n_mixture_pseudoabsences_multivariate.r')
###
### CONTENTS ###
### setup ###
### non-integrated SDM ###
### model diagnostics ###
### map of current distribution ###
### maps of future distribution ###
### response curves ###

#############
### setup ###
#############

	rm(list = ls())

	# drive <- 'C:/Ecology/'
	drive <- 'C:/Subarashi/'

	setwd(paste0(drive, '/Research/Andropogon/Andropogon'))

	library(bayesplot) # graphing
	library(coda) # Bayesian diagnostics
	library(cowplot) # combining ggplots
	library(ggplot2) # plotting
	library(ggspatial) # plotting spatial things
	library(nimble) # Bayes
	# library(nimbleHMC) # Hamiltonian Monte Carlo samplers
	library(omnibus)
	library(scales) # for plotting transparency
	library(terra) # spatial objects

	### user-defined values
	#######################

		climate_predictor_names <- c('bio1', 'bio12', 'bio15') # set for Jack's multivariate analyses, chosen based on simple models and collinearity

		# soil_predictor_names <- c('ph', 'sand', 'soc', 'silt')
		soil_predictor_names <- NULL

		predictor_names <- c(climate_predictor_names, soil_predictor_names)

		psa_quant <- 0.99 # to define pseudoabsences, use this quantile of Poaceae occurrences across counties with Poaceae occurrences 
		# psa_quant <- 0 # to define pseudoabsences, use this quantile of Poaceae occurrences across counties with Poaceae occurrences 

		### MCMC settings
		niter <- 280000
		nburnin <- 40000
		thin <- 200
		nchains <- 4
		trial <- FALSE

		# ### MCMC settings: SHORT RUNS
		# niter <- 140000
		# nburnin <- 40000
		# thin <- 100
		# nchains <- 2
		# trial <- TRUE

		# ### MCMC settings: UNIVARIATE
		# niter <- 240000
		# nburnin <- 40000
		# thin <- 200
		# nchains <- 4
		# trial <- FALSE

		# # # for testing
		# niter <- 90
		# nburnin <- 10
		# thin <- 1
		# nchains <- 2
		# trial <- TRUE

	# out_dir <- paste0('./outputs_loretta/sdm_[TEMP]/')
	# out_dir <- paste0('./outputs_loretta/sdm_[nmixture]_[pseudoabsences_', psa_quant, ']_[', paste(predictor_names, collapse = '_'), '_quad_ia]_[priors_ddnorm]', ifelse(trial, '_TRIAL', ''), '/')
	# out_dir <- paste0('./outputs_loretta/sdm_[nmixture]_[pseudoabsences_', psa_quant, ']_[bio1^2_bio12^2_bio15^2_bio1xbio12_bio1xbio12]_[priors_ddnorm]', ifelse(trial, '_TRIAL', ''), '/')
	# out_dir <- paste0('./outputs_loretta/sdm_[nmixture]_[pseudoabsences_', psa_quant, ']_[bio1^2_bio12^2_bio15^2_bio1xbio12xbio15]_[priors_ddnorm]', ifelse(trial, '_TRIAL', ''), '/')
	# out_dir <- paste0('./outputs_loretta/sdm_[nmixture]_[pseudoabsences_', psa_quant, ']_[bio1^2_bio12^2_bio15^2_bio1xbio12]_[priors_ddnorm]', ifelse(trial, '_TRIAL', ''), '/')
	# out_dir <- paste0('./outputs_loretta/sdm_[nmixture]_[pseudoabsences_', psa_quant, ']_[bio1^2_bio12^2_bio15^2_bio1xbio12xbio15]_[priors_ddnorm]', ifelse(trial, '_TRIAL', ''), '/')
	# out_dir <- paste0('./outputs_loretta/sdm_[nmixture]_[pseudoabsences_', psa_quant, ']_[bio1^2_bio12^2_bio15^2_bio1xbio12_bio1xbio15_bio12xbio15]_[priors_ddnorm]', ifelse(trial, '_TRIAL', ''), '/')
	out_dir <- paste0('./outputs_loretta/sdm_[nmixture]_[pseudoabsences_', psa_quant, ']_[bio1^2_bio12^2_bio15^2_bio12xbio15]_[priors_ddnorm]', ifelse(trial, '_TRIAL', ''), '/')
	dirCreate(out_dir)

	sink(paste0(out_dir, '/runtime_log.txt'), split = TRUE)
	say()

say('##########################')
say('### non-integrated SDM ###')
say('##########################')
say(date(), post = 1)

	say('This model is for the spatial distribution of AG. Currently, it assumes distribution is driven only by climate and soil. Occurrences are at the county level, so county area is used as an offset.', breaks = 60, post = 1)

	say('MCMC settings:', level = 2)
	say('niter ......', niter)
	say('nburnin ... ', nburnin)
	say('thin .......', thin)
	say('nchains ... ', nchains)
	say('psa_quant . ', psa_quant, post = 1)

	### model predictors and terms
	say('Predictors and model formula:', level = 2)

	say('We are using climate predictors: ', paste(climate_predictor_names, collapse = ' '), post = 1)
	say('We are using soil predictors: ', paste(soil_predictor_names, collapse = ' '), post = 1)

	# combos <- combn(climate_predictor_names, 2, simplify = FALSE)

	# form <- c(1, climate_predictor_names, paste0('I(', climate_predictor_names, '^2)'))
	# for (i in seq_along(combos)) form <- c(form, paste(combos[[i]], collapse = ':'))
	# form <- paste(form, collapse = ' + ')
	# form <- paste0('~ ', form)


	# if (!is.null(soil_predictor_names)) form <- paste0(form, ' + ', paste0(soil_predictor_names, collapse = ' + '))
	# form <- '~ 1 + bio1 + bio12 + bio15 + I(bio1^2) + I(bio12^2) + I(bio15^2) + bio1:bio12 + bio1:bio15 + bio12:bio15 + bio1:bio12:bio15'
	# form <- '~ 1 + bio1 + bio12 + bio15 + I(bio1^2) + I(bio12^2) + I(bio15^2) + bio1:bio12 + bio1:bio15 + bio12:bio15'
	form <- '~ 1 + bio1 + bio12 + bio15 + I(bio1^2) + I(bio12^2) + I(bio15^2) + bio12:bio15'
	# create model formula
	say('Model formula:')
	say(form, post = 1, breaks = 80)
	form <- as.formula(form)

	### load AG data
	ag_vect_sq <- vect('./data_from_adam_and_loretta/andropogon_gerardi_occurrences_with_environment_1961_2020.gpkg')
	ag_vect_ssp245_2041_2070 <- vect('./data_from_adam_and_loretta/andropogon_gerardi_occurrences_with_environment_ensemble_8GCMs_ssp245_2041_2070.gpkg')
	ag_vect_ssp245_2071_2100 <- vect('./data_from_adam_and_loretta/andropogon_gerardi_occurrences_with_environment_ensemble_8GCMs_ssp245_2071_2100.gpkg')
	ag_vect_ssp370_2041_2070 <- vect('./data_from_adam_and_loretta/andropogon_gerardi_occurrences_with_environment_ensemble_8GCMs_ssp370_2041_2070.gpkg')
	ag_vect_ssp370_2071_2100 <- vect('./data_from_adam_and_loretta/andropogon_gerardi_occurrences_with_environment_ensemble_8GCMs_ssp370_2071_2100.gpkg')

	ag_vect_ssp245_2041_2070 <- ag_vect_ssp245_2041_2070[!is.infinite(ag_vect_ssp245_2041_2070$aridity) & !is.na(ag_vect_ssp245_2041_2070$aridity)]
	ag_vect_ssp245_2071_2100 <- ag_vect_ssp245_2071_2100[!is.infinite(ag_vect_ssp245_2071_2100$aridity) & !is.na(ag_vect_ssp245_2071_2100$aridity)]
	ag_vect_ssp370_2041_2070 <- ag_vect_ssp370_2041_2070[!is.infinite(ag_vect_ssp370_2041_2070$aridity) & !is.na(ag_vect_ssp370_2041_2070$aridity)]
	ag_vect_ssp370_2071_2100 <- ag_vect_ssp370_2071_2100[!is.infinite(ag_vect_ssp370_2071_2100$aridity) & !is.na(ag_vect_ssp370_2071_2100$aridity)]

	# ag_vect_soil <- vect('./data_from_adam_and_loretta/andropogon_gerardi_occurrences_with_environment_soil.gpkg')
	# soil_names <- names(ag_vect_soil)
	# soil_in_vect <- soil_names[soil_names %in% predictor_names]
	
	# soil_columns <- as.data.frame(ag_vect_soil)
	# soil_columns <- soil_columns[ , soil_in_vect]
	
	# ag_vect_sq <- cbind(ag_vect_sq, soil_columns)
	# ag_vect_ssp245_2041_2070 <- cbind(ag_vect_ssp245_2041_2070, soil_columns)
	# ag_vect_ssp245_2071_2100 <- cbind(ag_vect_ssp245_2071_2100, soil_columns)
	# ag_vect_ssp370_2041_2070 <- cbind(ag_vect_ssp370_2041_2070, soil_columns)
	# ag_vect_ssp370_2071_2100 <- cbind(ag_vect_ssp370_2071_2100, soil_columns)

	fields <- c('area_km2', 'any_ag_quality_1_to_3', 'num_poaceae_records', predictor_names)
	ag_vect_sq <- ag_vect_sq[ , fields]

	# for counties with NA Poaceae and AG, assign a maximal number of Poaceae and 0 AG
	n_pseudoabs <- quantile(ag_vect_sq$num_poaceae_records, psa_quant, na.rm = TRUE)
	ag_vect_sq$num_poaceae_records[is.na(ag_vect_sq$num_poaceae_records)] <- n_pseudoabs
	ag_vect_sq$any_ag_quality_1_to_3[is.na(ag_vect_sq$any_ag_quality_1_to_3)] <- 0

	### collate data
	ag_sq <- as.data.frame(ag_vect_sq)

	### county area... used as an offset to make underlying model fit an IPP
	area_km2 <- ag_sq$area_km2
	log_area_km2 <- log(area_km2)
	log_area_km2_scaled <- scale(log_area_km2)
	log_area_center <- attributes(log_area_km2_scaled)$`scaled:center`
	log_area_scale <- attributes(log_area_km2_scaled)$`scaled:scale`
	log_area_km2_scaled <- log_area_km2_scaled[ , 1]

	### number of Poaceae records... used to model sampling bias
	log_num_poaceae_records <- log1p(ag_sq$num_poaceae_records) # log(x_sq + 1)
	log_num_poaceae_records_scaled <- scale(log_num_poaceae_records)
	log_num_poaceae_records_scaled <- log_num_poaceae_records_scaled[ , 1]

	### subset, scale, and manipulate predictors into model frame
	# status quo: counties with training data
	x_raw <- as.data.frame(ag_sq[ , c('area_km2', predictor_names)])
	log_area_km2 <- log(x_raw$area_km2)
	log_area_km2 <- data.frame(log_area_km2 = log_area_km2)
	x_raw$area_km2 <- NULL
	x_raw <- cbind(log_area_km2, x_raw)
	x_raw_scaled <- scale(x_raw)
	x_centers <- attr(x_raw_scaled, 'scaled:center')
	x_scales <- attr(x_raw_scaled, 'scaled:scale')
	x_raw_scaled <- as.data.frame(x_raw_scaled)
	x_sq <- model.matrix(form, as.data.frame(x_raw_scaled))

	# ssp245_2041_2070
	x_fut <- ag_vect_ssp245_2041_2070
	x_fut <- as.data.frame(x_fut[ , c('area_km2', predictor_names)])
	log_area_km2 <- log(x_fut$area_km2)
	log_area_km2 <- data.frame(log_area_km2 = log_area_km2)
	x_fut$area_km2 <- NULL
	x_fut <- cbind(log_area_km2, x_fut)
	x_fut <- scale(x_fut, center = x_centers, scale = x_scales)
	x_fut <- as.data.frame(x_fut)
	x_ssp245_2041_2070 <- model.matrix(form, as.data.frame(x_fut))

	# ssp245_2071_2100
	x_fut <- ag_vect_ssp245_2071_2100
	x_fut <- as.data.frame(x_fut[ , c('area_km2', predictor_names)])
	log_area_km2 <- log(x_fut$area_km2)
	log_area_km2 <- data.frame(log_area_km2 = log_area_km2)
	x_fut$area_km2 <- NULL
	x_fut <- cbind(log_area_km2, x_fut)
	x_fut <- scale(x_fut, center = x_centers, scale = x_scales)
	x_fut <- as.data.frame(x_fut)
	x_ssp245_2071_2100 <- model.matrix(form, as.data.frame(x_fut))

	# ssp370_2041_2070
	x_fut <- ag_vect_ssp370_2041_2070
	x_fut <- as.data.frame(x_fut[ , c('area_km2', predictor_names)])
	log_area_km2 <- log(x_fut$area_km2)
	log_area_km2 <- data.frame(log_area_km2 = log_area_km2)
	x_fut$area_km2 <- NULL
	x_fut <- cbind(log_area_km2, x_fut)
	x_fut <- scale(x_fut, center = x_centers, scale = x_scales)
	x_fut <- as.data.frame(x_fut)
	x_ssp370_2041_2070 <- model.matrix(form, as.data.frame(x_fut))

	# ssp245_2071_2100
	x_fut <- ag_vect_ssp370_2071_2100
	x_fut <- as.data.frame(x_fut[ , c('area_km2', predictor_names)])
	log_area_km2 <- log(x_fut$area_km2)
	log_area_km2 <- data.frame(log_area_km2 = log_area_km2)
	x_fut$area_km2 <- NULL
	x_fut <- cbind(log_area_km2, x_fut)
	x_fut <- scale(x_fut, center = x_centers, scale = x_scales)
	x_fut <- as.data.frame(x_fut)
	x_ssp370_2071_2100 <- model.matrix(form, as.data.frame(x_fut))

	### response array
	##################

	# This is an array for plotting the response curves of the SDM
	# columns: variables, all but one held at mean value across counties with AG
	# rows: focal variable increments, others held constant
	# depths: identity of focal variable changes

	# number of rows in array--represents number of points along which we'll be able to make the response curve
	n_predictors <- length(predictor_names)
	n_array_rows <- 200

	env_array <- array(
		NA,
		dim = c(n_array_rows, n_predictors, n_predictors),
		dimnames = list(
			1:n_array_rows,
			predictor_names,
			predictor_names
		)
	)

	# calculate means of each variable across counties with AG presences
	ag_pres <- ag_vect_sq[ag_vect_sq$any_ag_quality_1_to_3 > 0]
	ag_pres <- as.data.frame(ag_pres)[ , predictor_names, drop = FALSE]
	means <- colMeans(ag_pres)

	for (i in seq_along(predictor_names)) {
		for (j in seq_along(predictor_names)) {
			env_array[ , i, j] <- means[i]
		}
	}

	mins <- apply(ag_pres, 2, min)
	maxs <- apply(ag_pres, 2, max)

	for (i in seq_along(predictor_names)) {
		env_array[ , i, i] <- seq(mins[i], maxs[i], length.out = n_array_rows)
	}

	# scale
	env_array_scaled <- env_array
	for (i in seq_along(predictor_names)) {
		env_array_scaled[ , , i] <- sweep(env_array_scaled[ , , i], 2, x_centers[predictor_names], '-')
		env_array_scaled[ , , i] <- sweep(env_array_scaled[ , , i], 2, x_scales[predictor_names], '/')
	}

	# make into model matrix
	mm_array <- list()
	for (i in seq_along(predictor_names)) {
		mm_array[[i]] <- model.matrix(form, as.data.frame(env_array_scaled[ , , i, drop = TRUE]))
	}
	nc <- ncol(mm_array[[1]])
	response_curve_x <- array(
		as.numeric(unlist(mm_array)),
		dim = c(n_array_rows, nc, n_predictors),
		dimnames = list(1:n_array_rows, colnames(mm_array[[1]]), predictor_names)
	)

	n_response_curve_rows <- nrow(response_curve_x)

	### inputs for nimble
	say('Inputs:', level = 2)
	data <- list(
		y = ag_sq$any_ag_quality_1_to_3 # number of AG records in each county
	)

	n_counties <- nrow(ag_sq)
	n_counties_future <- nrow(x_ssp370_2071_2100)

	n_sdm_terms <- ncol(x_sq)

	constants <- list(
		n_counties = n_counties,
		n_counties_future = n_counties_future,

		n_sdm_terms = n_sdm_terms,

		x_sq = x_sq,

		n_response_curve_rows = n_response_curve_rows,
		n_predictors = n_predictors,
		response_curve_x = response_curve_x,

		log_area_km2_scaled = log_area_km2_scaled,
		log_num_poaceae_records_scaled = log_num_poaceae_records_scaled,

		x_ssp245_2041_2070 = x_ssp245_2041_2070,
		x_ssp245_2071_2100 = x_ssp245_2071_2100,

		x_ssp370_2041_2070 = x_ssp370_2041_2070,
		x_ssp370_2071_2100 = x_ssp370_2071_2100

	)

	lambda_fut_inits <- rep(mean(ag_sq$any_ag_quality_1_to_3), n_counties_future)
	N_inits <- 10 * ag_sq$any_ag_quality_1_to_3
	lambda_sq_inits <- 1 + ag_sq$any_ag_quality_1_to_3

	beta_inits <- rep(0, n_sdm_terms)
	beta_inits[grepl(colnames(x_sq), pattern = '\\^2')] <- -2
	beta_inits[grepl(colnames(x_sq), pattern = '\\:')] <- 0

	inits <- list(

		beta = beta_inits,
		alpha0_sampling = 0,
		alpha_area = 1,
		alpha_poaceae = 1,
		
		N = N_inits,

		lambda_sq = lambda_sq_inits,

		lambda_ssp245_2041_2070 = lambda_fut_inits,
		lambda_ssp245_2071_2100 = lambda_fut_inits,
		lambda_ssp370_2041_2070 = lambda_fut_inits,
		lambda_ssp370_2071_2100 = lambda_fut_inits

	)

	say('Data:')
	print(str(data))

	say('Constants:', pre = 1)
	print(str(constants))

	say('Initializations:', pre = 1)
	print(str(inits))

	### define model
	say('nimbleCode():', level = 2)
	say('We assume an N-mixture model (latent, real abundance ~ Poisson, and observations ~ binomial draws from latent abundance.', post = 2)
	model_code <- nimbleCode({
	  
		# priors for relationship to environment
		beta[1] ~ dnorm(0, sd = 10) # intercept... not regularized
		for (j in 2:n_sdm_terms) {
			# beta[j] ~ ddexp(0, rate = 1) # regularization toward 0
			beta[j] ~ dnorm(0, sd = 10) # broad prior
		}

		# priors for sampling bias
		alpha0_sampling ~ dnorm(0, sd = 10)
		alpha_area ~ dnorm(0, sd = 10)
		alpha_poaceae ~ dnorm(0, sd = 10)

		# likelihood
		for (i in 1:n_counties) {    # this specifies estimates be made for each county
			
			### actual abundance (latent--unobserved)
			N[i] ~ dpois(lambda_sq[i])

			# relationship between latent abundance and environment
			log(lambda_sq[i]) <- inprod(beta[1:n_sdm_terms], x_sq[i, 1:n_sdm_terms])

			### observed number of AG
			y[i] ~ dbin(prob = p[i], size = N[i])

			# sampling bias
			logit(p[i]) <- alpha0_sampling + alpha_area * log_area_km2_scaled[i] + alpha_poaceae * log_num_poaceae_records_scaled[i]

		}

		# posterior samplers for future predictions
		for (i in 1:n_counties_future) {

			log(lambda_ssp245_2041_2070[i]) <- inprod(beta[1:n_sdm_terms], x_ssp245_2041_2070[i, 1:n_sdm_terms])
			log(lambda_ssp245_2071_2100[i]) <- inprod(beta[1:n_sdm_terms], x_ssp245_2071_2100[i, 1:n_sdm_terms])
			
			log(lambda_ssp370_2041_2070[i]) <- inprod(beta[1:n_sdm_terms], x_ssp370_2041_2070[i, 1:n_sdm_terms])
			log(lambda_ssp370_2071_2100[i]) <- inprod(beta[1:n_sdm_terms], x_ssp370_2071_2100[i, 1:n_sdm_terms])

		}

		# posterior predictive sampler for response curves
		for (j in 1:n_predictors) {
			for (i in 1:n_response_curve_rows) {
				
				log(response_curves_sdm_lambda[i, j]) <-
					inprod(beta[1:n_sdm_terms], response_curve_x[i, 1:n_sdm_terms, j])

				exp_response_curves_sdm_lambda[i, j] <- exp(response_curves_sdm_lambda[i, j])

			}
		}

	})

	print(model_code)

	say('nimbleModel():', level = 2)
	model_species <- nimbleModel(
		code = model_code, # our model
		constants = constants, # constants
		data = data, # data
		inits = inits, # initialization values
		check = TRUE, # any errors?
		calculate = FALSE
		# buildDerivs = TRUE # need for Hamiltonian Monte Carlo
	)

	say('$initializeInfo() and $calculate():', level = 2)
	model_species$initializeInfo()
	model_species$calculate()

	say('configureMCMC():', level = 2)
	monitors <- c(
		'beta', 'alpha0_sampling', 'alpha_area', 'alpha_poaceae',
		'exp_response_curves_sdm_lambda',
		'lambda_sq',	
		'lambda_ssp245_2041_2070', 'lambda_ssp245_2071_2100', 'lambda_ssp370_2041_2070', 'lambda_ssp370_2071_2100'
	)
	
	conf <- configureMCMC(
		model_species,
		monitors = monitors,
		print = TRUE,
		enableWAIC = FALSE
	)

	# # add no U-turn sampler (Hamiltonian Monte Carlo)
	# conf$addSampler(target = c('beta', 'alpha_area', 'alpha_poaceae'), type = 'NUTS')

	### compile/build/run model/save MCMC
	build <- buildMCMC(conf)

	compiled <- compileNimble(model_species, build, showCompilerOutput = FALSE)

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

	saveRDS(chains, paste0(out_dir, '/sdm_nmixture_chains.rds'))

	say('session info', level = 2)
	print(sessionInfo())

	say(date(), pre = 1)
	sink()

say('#########################')
say('### model diagnostics ###')
say('#########################')

	chains <- readRDS(paste0(out_dir, '/sdm_nmixture_chains.rds'))
	mcmc <- chains$samples

	for (i in 1:nchains) {
		cols <- c(paste0('beta[', 1:n_sdm_terms, ']'), 'alpha0_sampling', 'alpha_poaceae', 'alpha_area')
		mcmc[[i]] <- mcmc[[i]][ , cols]
	}

	# graphing trace plots for all betas
	pars <- paste0('beta[', 1:n_sdm_terms, ']')
	file <- paste0(out_dir, '/beta_trace.png')
	ggsave(mcmc_trace(mcmc, pars = pars), file = file, width = 10, height = 8, dpi = 450, bg = 'white')

	# graphing trace plots for all "extra" betas
	pars <- c('alpha0_sampling', 'alpha_poaceae', 'alpha_area')
	file <- paste0(out_dir, '/alpha_trace.png')
	ggsave(mcmc_trace(mcmc, pars = pars), file = file, width = 10, height = 4, dpi = 450, bg = 'white')

	# graphing density plots for all betas
	pars <- paste0('beta[', 1:n_sdm_terms, ']')
	file <- paste0(out_dir, '/beta_density.png')
	ggsave(mcmc_dens_overlay(mcmc, pars = pars), file = file, width = 12, height = 8, dpi = 450, bg = 'white')

	# graphing trace plots for all "extra" betas
	pars <- c('alpha0_sampling', 'alpha_poaceae', 'alpha_area')
	file <- paste0(out_dir, '/alpha_density.png')
	ggsave(mcmc_dens_overlay(mcmc, pars = pars), file = file, width = 10, height = 4, dpi = 450, bg = 'white')

	# Gelman-Rubin statistic
	rhats <- gelman.diag(mcmc, autoburnin = FALSE, multivariate = TRUE)

	sink(paste0(out_dir, '/convergence.txt'), split = TRUE)
	say('GELMAN-RUBIN STATISTICS')
	say(date(), post = 2)
	print(rhats)
	sink()

say('###################################', pre = 1)
say('### map of current distribution ###')
say('###################################')

	chains <- readRDS(paste0(out_dir, '/sdm_nmixture_chains.rds'))

	# subset chain summary to just the lambdas
	summary <- chains$summary$all.chains

	# just counties with data
	which_lambda <- grepl(rownames(summary), pattern = 'lambda_sq')
	lambda <- summary[which_lambda, ]

	ag_vect_sq$lambda_mean <- lambda[ , 'Mean']
	ag_vect_sq$lambda_0.05ci <- lambda[ , '95%CI_low']
	ag_vect_sq$lambda_0.95ci <- lambda[ , '95%CI_upp']
	ag_vect_sq$lambda_ci <- ag_vect_sq$lambda_0.95ci - ag_vect_sq$lambda_0.05ci

	lambda_sq_quants <- quantile(ag_vect_sq$lambda_mean, c(0.25, 0.5, 0.75, 0.90, 0.95))
	quant_labels <- c('[0 - 0.25)', '[0.25 - 0.50)', '[0.50 - 0.75)', '[0.75 - 0.90)', '[0.90 - 0.95)', '[0.95 - 1]')
	
	ag_vect_sq$quant_col <- NA
	ag_vect_sq$quant_col[ag_vect_sq$lambda_mean < lambda_sq_quants[1]] <- quant_labels[1]
	ag_vect_sq$quant_col[ag_vect_sq$lambda_mean >= lambda_sq_quants[1] & ag_vect_sq$lambda_mean < lambda_sq_quants[2]] <- quant_labels[2]
	ag_vect_sq$quant_col[ag_vect_sq$lambda_mean >= lambda_sq_quants[2] & ag_vect_sq$lambda_mean < lambda_sq_quants[3]] <- quant_labels[3]
	ag_vect_sq$quant_col[ag_vect_sq$lambda_mean >= lambda_sq_quants[3] & ag_vect_sq$lambda_mean < lambda_sq_quants[4]] <- quant_labels[4]
	ag_vect_sq$quant_col[ag_vect_sq$lambda_mean >= lambda_sq_quants[4] & ag_vect_sq$lambda_mean < lambda_sq_quants[5]] <- quant_labels[5]
	ag_vect_sq$quant_col[ag_vect_sq$lambda_mean >= lambda_sq_quants[5]] <- quant_labels[6]

	nam <- vect(paste0(drive, '/Research Data/GADM/Version 4.1/High Res North America Level 1 sans Great Lakes SpatVector WGS84.gpkg'))
	nam <- project(nam, ag_vect_sq)

	ag_vect_sq_pres <- ag_vect_sq[ag_vect_sq$any_ag_quality_1_to_3 > 0]
	extent <- ext(ag_vect_sq_pres)
	extent <- as.vector(extent)
	x_range <- (extent[2] - extent[1])
	y_range <- (extent[4] - extent[3])
	extent[1] <- extent[1] + 0.15 * x_range
	extent[3] <- extent[3] + 0.125 * y_range
	extent[4] <- extent[4] - 0.2 * y_range

	cents_with_ag <- ag_vect_sq[ag_vect_sq$any_ag_quality_1_to_3 > 0]
	cents_with_ag <- centroids(cents_with_ag)

	fill_scale <- c(
		'[0 - 0.25)' = 'gray85',
		'[0.25 - 0.5)' = alpha('forestgreen', 0.20),
		'[0.50 - 0.75)' = alpha('forestgreen', 0.40),
		'[0.75 - 0.90)' = alpha('forestgreen', 0.65),
		'[0.90 - 0.95)' = alpha('forestgreen', 0.75),
		'[0.95 - 1]' = 'forestgreen'
	)

	map <- ggplot() +
		layer_spatial(ag_vect_sq, aes(fill = quant_col), color = NA) +
		layer_spatial(nam, color = 'gray30', fill = NA, linewidth = 0.1) +
		scale_fill_manual(
			name = 'Quantile\n of λ',
			values = fill_scale
		) +
		layer_spatial(cents_with_ag, pch = 16, color = alpha('gray20', 0.5), size = 0.35) +
		guides(fill = guide_legend(reverse = TRUE)) +
		xlim(extent[1], extent[2]) + ylim(extent[3], extent[4]) +
		ggtitle(expression('Present-day distribution of ' * italic('Andropogon gerardi')), subtitle = '1961-2020 | N-mixture model') +
		theme(
			plot.title = element_text(size = 16),
			plot.subtitle = element_text(size = 14)
		)

	ggsave(plot = map, filename = paste0(out_dir, '/sdm_nmixture_lambda_status_quo.png'), width = 12, height = 10, dpi = 600)

	writeVector(ag_vect_sq, paste0(out_dir, '/sdm_nmixture_1961_2020_climate_focus.gpkg'), overwrite = TRUE)

say('###################################')
say('### maps of future distribution ###')
say('###################################')

	futs <- c(
		'ssp245_2041_2070',
		'ssp245_2071_2100',
		'ssp370_2041_2070',
		'ssp370_2071_2100'
	)

	nam <- vect(paste0(drive, '/Research Data/GADM/Version 4.1/High Res North America Level 1 sans Great Lakes SpatVector WGS84.gpkg'))
	nam <- project(nam, ag_vect_sq)

	chains <- readRDS(paste0(out_dir, '/sdm_nmixture_chains.rds'))

	# subset chain summary to just the lambdas
	summary <- chains$summary$all.chains

	for (fut in futs) {

		which_lambda <- grepl(rownames(summary), pattern = paste0('lambda_', fut))
		lambda <- summary[which_lambda, ]

		this_ag_vect <- get(paste0('ag_vect_', fut))

		# get just counties with data
		this_ag_vect$lambda_mean <- lambda[ , 'Mean']
		this_ag_vect$lambda_0.05ci <- lambda[ , '95%CI_low']
		this_ag_vect$lambda_0.95ci <- lambda[ , '95%CI_upp']
		this_ag_vect$lambda_ci <- this_ag_vect$lambda_0.95ci - this_ag_vect$lambda_0.05ci

		quant_labels <- c('[0 - 0.25)', '[0.25 - 0.50)', '[0.50 - 0.75)', '[0.75 - 0.90)', '[0.90 - 0.95)', '[0.95 - 1]') # should be same as above
		# NB lambda_sq_quants is from above!!!
		this_ag_vect$quant_col <- NA
		this_ag_vect$quant_col[this_ag_vect$lambda_mean < lambda_sq_quants[1]] <- quant_labels[1]
		this_ag_vect$quant_col[this_ag_vect$lambda_mean >= lambda_sq_quants[1] & this_ag_vect$lambda_mean < lambda_sq_quants[2]] <- quant_labels[2]
		this_ag_vect$quant_col[this_ag_vect$lambda_mean >= lambda_sq_quants[2] & this_ag_vect$lambda_mean < lambda_sq_quants[3]] <- quant_labels[3]
		this_ag_vect$quant_col[this_ag_vect$lambda_mean >= lambda_sq_quants[3] & this_ag_vect$lambda_mean < lambda_sq_quants[4]] <- quant_labels[4]
		this_ag_vect$quant_col[this_ag_vect$lambda_mean >= lambda_sq_quants[4] & this_ag_vect$lambda_mean < lambda_sq_quants[5]] <- quant_labels[5]
		this_ag_vect$quant_col[this_ag_vect$lambda_mean >= lambda_sq_quants[5]] <- quant_labels[6]

		pretty_title <- paste0('SSP ', substr(fut, 4, 6), ' ', substr(fut, 8, 11), '-', substr(fut, 13, 16), ' | N-mixture model')

		fill_scale <- c(
			'[0 - 0.25)' = 'gray85',
			'[0.25 - 0.5)' = alpha('forestgreen', 0.20),
			'[0.50 - 0.75)' = alpha('forestgreen', 0.40),
			'[0.75 - 0.90)' = alpha('forestgreen', 0.65),
			'[0.90 - 0.95)' = alpha('forestgreen', 0.75),
			'[0.95 - 1]' = 'forestgreen'
		)

		
		map <- ggplot() +
			layer_spatial(this_ag_vect, aes(fill = quant_col), color = NA) +
			layer_spatial(nam, color = 'gray50', fill = NA, linewidth = 0.1) +
			scale_fill_manual(
				name = 'Quantile\n of λ',
				values = fill_scale
			) +
			# layer_spatial(cents_with_ag, pch = 16, col = 'black', size = 0.2) +
			guides(fill = guide_legend(reverse = TRUE)) +
			xlim(extent[1], extent[2]) + ylim(extent[3], extent[4]) +
			ggtitle(expression('Future equilibrial distribution of ' * italic('Andropogon gerardi')), subtitle = pretty_title) +
			theme(
				plot.title = element_text(size = 16),
				plot.subtitle = element_text(size = 14)
			)
				

		ggsave(plot = map, filename = paste0(out_dir, '/sdm_nmixture_lambda_', fut, '.png'), width = 12, height = 10, dpi = 600)

		writeVector(this_ag_vect, paste0(out_dir, '/sdm_nmixture_', fut, '.gpkg'), overwrite = TRUE)

	} # next future

say('###########################')
say('### parameter estimates ###')
say('###########################')

	chains <- readRDS(paste0(out_dir, '/sdm_nmixture_chains.rds'))

	# subset chain summary to just the lambdas
	summary <- chains$summary$all.chains

	model_term <- attr(terms(form), 'term.labels')
	model_term_nice <- model_term
	model_term_nice <- gsub(model_term_nice, pattern = 'I\\(', replacement = '')
	model_term_nice <- gsub(model_term_nice, pattern = '\\)', replacement = '')
	model_term_nice <- c('Intercept', model_term_nice)

	pars <- paste0('beta[', 1:(1 + length(model_term)), ']')

	coeff_chains <- chains$samples
	for (i in 1:nchains) {
		coeff_chains[[i]] <- coeff_chains[[i]][ , pars]
		colnames(coeff_chains[[i]]) <- model_term_nice
	}

	caterpillars <- mcmc_intervals(coeff_chains, pars = model_term_nice)

	# plot posterior of coefficient estimates
	caterpillars <- caterpillars +
		xlab('Estimated value') +
		theme(plot.background = element_rect(fill = 'white'))

	ggsave(caterpillars, filename = paste0(out_dir, '/sdm_nmixture_coefficients.png'), width = 10, height = 12, dpi = 300)

say('#######################')
say('### response curves ###')
say('#######################')

	chains <- readRDS(paste0(out_dir, '/sdm_nmixture_chains.rds'))

	# make a plot for how AG SDM and each ancestral population responds to each predictor
	responses <- list()
	for (i in 1:n_predictors) {

		pred <- predictor_names[i]

		# unscaled predictor value
		x <- env_array[ , pred, pred]
		
		# create data frame with SDM response
		pars <- paste0('exp_response_curves_sdm_lambda[', 1:n_response_curve_rows, ', ', i, ']')
		sdm_response_mean <- chains$summary$all.chains[pars, 'Mean']
		sdm_response_lower <- chains$summary$all.chains[pars, '95%CI_low']
		sdm_response_upper <- chains$summary$all.chains[pars, '95%CI_upp']

		# data frames to hold predictions in long format
		df_mean <- data.frame(
			x = x,
			response = log(sdm_response_mean)
		)

		df_lower <- data.frame(
			x = x,
			response = log(sdm_response_lower)
		)

		df_upper <- data.frame(
			x = x,
			response = log(sdm_response_upper)
		)

		this_df_ci <- data.frame(
			x = c(x, rev(x)),
			y = c(df_upper$response, rev(df_lower$response))
		)

		lambda_max <- max(-Inf, df_upper$response[!is.infinite(df_upper$response)])

		if (pred == 'aridity') {
			nice_title <- 'Aridity ((temperature + 10) / (precipitation / 1000))'
			nice_axis <- 'Aridity'
		} else if (pred == 'bio1') {
			nice_title <- 'Mean annual temperature (BIO01)'
			nice_axis <- 'Mean annual temperature (°C)'
		} else if (pred == 'bio5') {
			nice_title <- 'Temperature of the hottest month (BIO05)'
			nice_axis <- 'Temperature of the hottest month (°C)'
		} else if (pred == 'bio6') {
			nice_title <- 'Temperature of the coldest month (BIO06)'
			nice_axis <- 'Temperature of the coldest month (°C)'
		} else if (pred == 'bio7') {
			nice_title <- 'Temperature annual range (BIO07)'
			nice_axis <- 'Temperature annual range (°C)'
		} else if (pred == 'bio12') {
			nice_title <- 'Total annual precipitation (BIO12)'
			nice_axis <- 'Total annual precipitation (mm)'
		} else if (pred == 'bio15') {
			nice_title <- 'Precipitation seasonality (BIO15)'
			nice_axis <- 'Precipitation seasonality'
		} else if (pred == 'bio18') {
			nice_title <- 'Precipitation of warmest quarter (BIO18)'
			nice_axis <- 'Precipitation of warmest quarter (mm)'
		} else if (pred == 'ph') {
			nice_title <- 'Soil pH'
			nice_axis <- 'pH'
		} else if (pred == 'sand') {
			nice_title <- 'Soil proportion sand'
			nice_axis <- 'Proportion sand'
		} else if (pred == 'silt') {
			nice_title <- 'Soil proportion silt'
			nice_axis <- 'Proportion silt'
		} else if (pred == 'soc') {
			nice_title <- 'Soil organic matter'
			nice_axis <- 'Soil organic matter'
		}

		response <- ggplot() +
			geom_polygon(
				data = this_df_ci,
				mapping = aes(x = x, y = y),
				color = NA,
				fill = alpha('blue', 0.1)
			) +
			geom_line(
				data = df_mean,
				mapping = aes(x = x, y = response)
			) +
			xlab(nice_axis) +
			ylab('Lambda') +
			ggtitle(nice_title) +
			theme(
				plot.title = element_text(size = 10)
			)

		responses[[pred]] <- response

	} # next predictor

	for (i in 1:n_predictors) {
		responses[[i]] <- responses[[i]] + ylim(0, lambda_max)
	}

	if (length(predictor_names) == 1) {
		nrow <- 1
		width <- 10
		height <- 8
	} else if (length(predictor_names) <= 4) {
		nrow <- 1
		width <- 14
		height <- 4
	} else {
		nrow <- 2
		width <- 16
		height <- 10	
	}

	responses <- plot_grid(plotlist = responses, nrow = nrow)

	ggsave(plot = responses, filename = paste0(out_dir, '/response_curves.png'), width = width, height = height, dpi = 600)

say(date())
say('FINIS!', deco = '+', level = 1)

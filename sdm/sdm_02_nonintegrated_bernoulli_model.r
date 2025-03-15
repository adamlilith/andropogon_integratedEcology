### MODELING ANDROPOGON GERARDI DISTRIBUTION, PHENOTYPE, PHYSIOLOGY, GENOTYPE, and ASSOCIATED MICROBIAL COMMUNITIES
### Erica Newman | Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2023-12
###
### This script constructs a non-integrated model for AG geographic distribution.
###
### source('C:/Ecology/R/andropogon_integratedEcology/sdm/sdm_02_nonintegrated_bernoulli_model.r')
### source('E:/Adam/R/andropogon_integratedEcology/sdm/sdm_02_nonintegrated_bernoulli_model.r')
###
### CONTENTS ###
### setup ###
### non-integrated SDM ###
### model diagnostics ###
### map of current distribution ###
### maps of future distribution ###

#############
### setup ###
#############

	rm(list = ls())

	drive <- 'C:/Ecology/'
	# drive <- 'E:/Adam/'

	setwd(paste0(drive, '/Research/Andropogon/Andropogon'))

	library(bayesplot) # graphing
	library(coda) # Bayesian diagnostics
	library(ggplot2) # plotting
	library(ggspatial) # plotting spatial things
	library(nimble) # Bayes
	# library(nimbleHMC) # Hamiltonian Monte Carlo samplers
	library(omnibus)
	library(scales) # for plotting transparency
	library(terra) # spatial objects

	out_dir <- paste0('./outputs_loretta/sdm_[bernoulli]/')
	dirCreate(out_dir)

	sink(paste0(out_dir, '/runtime_log.txt'), split = TRUE)
	say()

say('##########################')
say('### non-integrated SDM ###')
say('##########################')
say(date(), post = 1)

	say('This model is for the spatial distribution of AG. Currently, it assumes distribution is driven only by climate and soil. Occurrences are at the county level, so county area is used as an offset.', breaks = 60, post = 1)

	### MCMC settings
	niter <- 130000
	nburnin <- 10000
	thin <- 120
	nchains <- 4

	# # for testing
	# niter <- 140
	# nburnin <- 40
	# thin <- 1
	# nchains <- 2

	say('MCMC settings:', level = 2)
	say('niter ......', niter)
	say('nburnin ... ', nburnin)
	say('thin .......', thin)
	say('nchains ... ', nchains, post = 1)

	### model predictors and terms
	say('Predictors and model formula:', level = 2)

	# predictor_names <- c('aridity', 'bio7', 'bio2', 'sand', 'bio12', 'bio15', 'cec', 'soc', 'silt')
	predictor_names <- c('aridity', 'bio7', 'ph', 'sand', 'bio12', 'bio15', 'cec', 'soc', 'silt')
	say('We are using predictors: ', paste(predictor_names, collapse = ' '), post = 1)

	preselected_model_terms <- read.csv('./outputs_loretta/sig_coeffs_elastic_net_2024_06_04.csv')

	# get vector of linear predictors... we need to ensure all predictors appear at least as linear terms to respect marginality
	terms <- preselected_model_terms$term
	terms <- terms[!(terms %in% c('(Intercept)', 'log_area_km2'))]
	linear_terms <- gsub(terms, pattern = 'I\\(', replacement = '')
	linear_terms <- gsub(linear_terms, pattern = '\\^2)', replacement = '')
	linear_terms <- strsplit(linear_terms, ':')
	linear_terms <- unlist(linear_terms)
	linear_terms <- unique(linear_terms)
	linear_terms <- sort(linear_terms)

	# get non-linear terms
	terms <- terms[!(terms %in% linear_terms)]

	# create model formula
	form <- paste0(' ~ 1 + ', paste(c(linear_terms, terms), collapse = ' + '))
	say('Model formula:')
	say(form, post = 1, breaks = 80)
	form <- as.formula(form)

	### load AG data
	ag_vect_sq <- vect('./data_from_adam_and_loretta/andropogon_gerardi_occurrences_with_environment_1961_2020.gpkg')
	ag_vect_ssp245_2041_2070 <- vect('./data_from_adam_and_loretta/andropogon_gerardi_occurrences_with_environment_ensemble_8GCMs_ssp245_2041_2070.gpkg')
	ag_vect_ssp245_2071_2100 <- vect('./data_from_adam_and_loretta/andropogon_gerardi_occurrences_with_environment_ensemble_8GCMs_ssp245_2071_2100.gpkg')
	ag_vect_ssp370_2041_2070 <- vect('./data_from_adam_and_loretta/andropogon_gerardi_occurrences_with_environment_ensemble_8GCMs_ssp370_2041_2070.gpkg')
	ag_vect_ssp370_2071_2100 <- vect('./data_from_adam_and_loretta/andropogon_gerardi_occurrences_with_environment_ensemble_8GCMs_ssp370_2071_2100.gpkg')

	ag_vect_soil <- vect('./data_from_adam_and_loretta/andropogon_gerardi_occurrences_with_environment_soil.gpkg')
	soil_names <- names(ag_vect_soil)
	soil_in_vect <- soil_names[soil_names %in% predictor_names]
	
	soil_columns <- as.data.frame(ag_vect_soil)
	soil_columns <- soil_columns[ , soil_in_vect]
	
	ag_vect_sq <- cbind(ag_vect_sq, soil_columns)
	ag_vect_ssp245_2041_2070 <- cbind(ag_vect_ssp245_2041_2070, soil_columns)
	ag_vect_ssp245_2071_2100 <- cbind(ag_vect_ssp245_2071_2100, soil_columns)
	ag_vect_ssp370_2041_2070 <- cbind(ag_vect_ssp370_2041_2070, soil_columns)
	ag_vect_ssp370_2071_2100 <- cbind(ag_vect_ssp370_2071_2100, soil_columns)

	fields <- c('area_km2', 'any_ag_quality1to3', 'num_poaceae_records', predictor_names)
	ag_vect_sq <- ag_vect_sq[ , fields]
	ag_vect_sq$num_ag_quality1to3 <- ag_vect_sq$any_ag_quality1to3 # becomes a sampling predictor
	ag_vect_sq$any_ag_quality1to3 <- as.integer(ag_vect_sq$any_ag_quality1to3 > 0) # becomes the response

	# separate counties with/without training data
	ag <- as.data.frame(ag_vect_sq)
	completes <- complete.cases(as.data.frame(ag_vect_sq))
	ag_focus <- ag[completes, ]
	ag_vect_focus <- ag_vect_sq[completes, ]
	
	# complement of counties with training data for full-continent predictions
	ag_vect_focus_complement <- ag_vect_sq[!completes, ]
	ag_focus_complement <- as.data.frame(ag_vect_focus_complement)

	### collate data

	### county area... used as an offset to make underlying model fit an IPP
	area_km2 <- ag_focus$area_km2
	log_area_km2 <- log(area_km2)
	log_area_km2_scaled <- scale(log_area_km2)
	log_area_center <- attributes(log_area_km2_scaled)$`scaled:center`
	log_area_scale <- attributes(log_area_km2_scaled)$`scaled:scale`
	log_area_km2_scaled <- log_area_km2_scaled[ , 1]

	### number of Poaceae records... used to model sampling bias
	log_num_poaceae_records <- log1p(ag_focus$num_poaceae_records) # log(x + 1)
	log_num_poaceae_records_scaled <- scale(log_num_poaceae_records)
	log_num_poaceae_records_scaled <- log_num_poaceae_records_scaled[ , 1]

	### number of AG records... used to model sampling bias
	log_num_ag_records <- log1p(ag_focus$num_ag_quality1to3) # log(x + 1)
	log_num_ag_records_scaled <- scale(log_num_ag_records)
	log_num_ag_records_scaled <- log_num_ag_records_scaled[ , 1]

	### subset, scale, and manipulate predictors into model frame
	# status quo: counties with training data
	x_raw <- as.data.frame(ag_focus[ , predictor_names])
	# x_raw <- as.data.frame(ag_focus[ , c('area_km2', predictor_names)])
	# log_area_km2 <- log(x_raw$area_km2)
	# log_area_km2 <- data.frame(log_area_km2 = log_area_km2)
	# x_raw$area_km2 <- NULL
	# x_raw <- cbind(log_area_km2, x_raw)
	x_raw_scaled <- scale(x_raw)
	x_centers <- attr(x_raw_scaled, 'scaled:center')
	x_scales <- attr(x_raw_scaled, 'scaled:scale')
	x_raw_scaled <- as.data.frame(x_raw_scaled)
	x_sq <- model.matrix(form, as.data.frame(x_raw_scaled))

	# status quo: counties without training data
	# x_sq_complement <- as.data.frame(ag_focus_complement[ , c('area_km2', predictor_names)])
	x_sq_complement <- as.data.frame(ag_focus_complement[ , predictor_names])
	# log_area_km2 <- log(x_sq_complement$area_km2)
	# log_area_km2 <- data.frame(log_area_km2 = log_area_km2)
	# x_sq_complement$area_km2 <- NULL
	# x_sq_complement <- cbind(log_area_km2, x_sq_complement)
	x_sq_complement <- scale(x_sq_complement, center = x_centers, scale = x_scales)
	x_sq_complement <- as.data.frame(x_sq_complement)
	x_sq_complement <- model.matrix(form, as.data.frame(x_sq_complement))

	# ssp245_2041_2070
	x_fut <- ag_vect_ssp245_2041_2070
	x_fut <- as.data.frame(x_fut[ , predictor_names])
	# x_fut <- as.data.frame(x_fut[ , c('area_km2', predictor_names)])
	# log_area_km2 <- log(x_fut$area_km2)
	# log_area_km2 <- data.frame(log_area_km2 = log_area_km2)
	# x_fut$area_km2 <- NULL
	# x_fut <- cbind(log_area_km2, x_fut)
	x_fut <- scale(x_fut, center = x_centers, scale = x_scales)
	x_fut <- as.data.frame(x_fut)
	x_ssp245_2041_2070 <- model.matrix(form, as.data.frame(x_fut))

	# ssp245_2071_2100
	x_fut <- ag_vect_ssp245_2071_2100
	x_fut <- as.data.frame(x_fut[ , predictor_names])
	# x_fut <- as.data.frame(x_fut[ , c('area_km2', predictor_names)])
	# log_area_km2 <- log(x_fut$area_km2)
	# log_area_km2 <- data.frame(log_area_km2 = log_area_km2)
	# x_fut$area_km2 <- NULL
	# x_fut <- cbind(log_area_km2, x_fut)
	x_fut <- scale(x_fut, center = x_centers, scale = x_scales)
	x_fut <- as.data.frame(x_fut)
	x_ssp245_2071_2100 <- model.matrix(form, as.data.frame(x_fut))

	# ssp370_2041_2070
	x_fut <- ag_vect_ssp370_2041_2070
	x_fut <- as.data.frame(x_fut[ , predictor_names])
	# x_fut <- as.data.frame(x_fut[ , c('area_km2', predictor_names)])
	# log_area_km2 <- log(x_fut$area_km2)
	# log_area_km2 <- data.frame(log_area_km2 = log_area_km2)
	# x_fut$area_km2 <- NULL
	# x_fut <- cbind(log_area_km2, x_fut)
	x_fut <- scale(x_fut, center = x_centers, scale = x_scales)
	x_fut <- as.data.frame(x_fut)
	x_ssp370_2041_2070 <- model.matrix(form, as.data.frame(x_fut))

	# ssp245_2071_2100
	x_fut <- ag_vect_ssp370_2071_2100
	x_fut <- as.data.frame(x_fut[ , predictor_names])
	# x_fut <- as.data.frame(x_fut[ , c('area_km2', predictor_names)])
	# log_area_km2 <- log(x_fut$area_km2)
	# log_area_km2 <- data.frame(log_area_km2 = log_area_km2)
	# x_fut$area_km2 <- NULL
	# x_fut <- cbind(log_area_km2, x_fut)
	x_fut <- scale(x_fut, center = x_centers, scale = x_scales)
	x_fut <- as.data.frame(x_fut)
	x_ssp370_2071_2100 <- model.matrix(form, as.data.frame(x_fut))

	### inputs for nimble
	say('Inputs:', level = 2)
	data <- list(
		y = ag_focus$any_ag_quality1to3 # any AG in a county?
	)

	n_counties <- nrow(ag_focus)
	n_counties_complement <- nrow(x_sq_complement)
	n_counties_future <- nrow(x_ssp370_2071_2100)

	n_sdm_terms <- ncol(x_sq)

	constants <- list(
	  n_counties = n_counties,
	  n_counties_complement = n_counties_complement,
	  n_counties_future = n_counties_future,

	  n_sdm_terms = n_sdm_terms,

	  x_sq = x_sq,
	  x_sq_complement = x_sq_complement,

	  log_area_km2_scaled = log_area_km2_scaled,
	  log_num_poaceae_records_scaled = log_num_poaceae_records_scaled,
	  log_num_ag_records_scaled = log_num_ag_records_scaled,

	  x_ssp245_2041_2070 = x_ssp245_2041_2070,
	  x_ssp245_2071_2100 = x_ssp245_2071_2100,
	  
	  x_ssp370_2041_2070 = x_ssp370_2041_2070,
	  x_ssp370_2071_2100 = x_ssp370_2071_2100

	)

	inits <- list()
	for (i in seq_len(nchains)) {

		# lambda_sq_inits <- 1 + rpois(n_counties, ag_focus$any_ag_quality1to3)
		# lambda_sq_complement_inits <- 1 + rpois(n_counties_complement, mean(ag_focus$any_ag_quality1to3))
		# lambda_fut_inits <- rpois(n_counties_future, mean(ag_focus$any_ag_quality1to3))
		
		lambda_sq_inits <- 1 + rpois(n_counties, ag_focus$any_ag_quality1to3)
		lambda_sq_complement_inits <- 1 + rpois(n_counties_complement, mean(ag_focus$any_ag_quality1to3))
		lambda_fut_inits <- rpois(n_counties_future, mean(ag_focus$any_ag_quality1to3))
		
		y_inits <- sample(ag_focus$any_ag_quality1to3, n_counties)
		psi_inits <- runif(n_counties, 0.05, 0.5)

		inits[[i]] <- list(

			# beta = runif(n_sdm_terms, -0.5, 0.5),
			# # alpha0_sampling = runif(1, -0.5, 0.5),
			# alpha_area = runif(1, 0, 1),
			# alpha_poaceae = runif(1, -0.5, 0.5),
			# alpha_ag = runif(1, -0, 1),
			# alpha_ag = runif(1, -0, 1),
			
			beta = rep(0, n_sdm_terms),
			# alpha0_sampling = runif(1, -0.5, 0.5),
			alpha_area = 0,
			alpha_poaceae = 0,
			alpha_ag = 0,
			alpha_ag = 0,
			
			y = y_inits,
			psi = psi_inits,

			lambda_sq = lambda_sq_inits,
			lambda_sq_complement = lambda_sq_complement_inits,

			lambda_ssp245_2041_2070 = lambda_fut_inits,
			lambda_ssp245_2071_2100 = lambda_fut_inits,
			lambda_ssp370_2041_2070 = lambda_fut_inits,
			lambda_ssp370_2071_2100 = lambda_fut_inits
		
		)

	}

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
			beta[j] ~ ddexp(0, rate = 1) # regularization toward 0
			# beta[j] ~ dnorm(0, sd = 10) # broad prior
		}

		# priors for sampling bias
		# alpha0_sampling ~ dnorm(0, sd = 10)
		alpha_area ~ dnorm(0, sd = 10)
		alpha_poaceae ~ dnorm(0, sd = 10)
		alpha_ag ~ dnorm(0, sd = 10)

		# likelihood
		for (i in 1:n_counties) {    # this specifies estimates be made for each county
			
			### observed AG or not
			y[i] ~ dbern(psi[i])

			cloglog(psi[i]) <- log(lambda_sq[i]) + log(p[i])

			# relationship between latent abundance and environment
			log(lambda_sq[i]) <- inprod(beta[1:n_sdm_terms], x_sq[i, 1:n_sdm_terms])

			# sampling bias
			logit(p[i]) <- # alpha0_sampling +
				alpha_area * log_area_km2_scaled[i] +
				alpha_poaceae * log_num_poaceae_records_scaled[i] +
				alpha_ag * log_num_ag_records_scaled[i]

		}

		# posterior samplers for complement of status quo predictions
		for (i in 1:n_counties_complement) {

			log(lambda_sq_complement[i]) <- inprod(beta[1:n_sdm_terms], x_sq_complement[i, 1:n_sdm_terms])

		}
		
		# posterior samplers for future predictions
		for (i in 1:n_counties_future) {

			log(lambda_ssp245_2041_2070[i]) <- inprod(beta[1:n_sdm_terms], x_ssp245_2041_2070[i, 1:n_sdm_terms])
			log(lambda_ssp245_2071_2100[i]) <- inprod(beta[1:n_sdm_terms], x_ssp245_2071_2100[i, 1:n_sdm_terms])
			
			log(lambda_ssp370_2041_2070[i]) <- inprod(beta[1:n_sdm_terms], x_ssp370_2041_2070[i, 1:n_sdm_terms])
			log(lambda_ssp370_2071_2100[i]) <- inprod(beta[1:n_sdm_terms], x_ssp370_2071_2100[i, 1:n_sdm_terms])

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
	# monitors <- c('beta', 'alpha0_sampling', 'alpha_area', 'alpha_poaceae', 'alpha_ag', 'lambda_sq', 'lambda_sq_complement', 'lambda_ssp245_2041_2070', 'lambda_ssp245_2071_2100', 'lambda_ssp370_2041_2070', 'lambda_ssp370_2071_2100')
	monitors <- c('beta', 'alpha_area', 'alpha_poaceae', 'alpha_ag', 'lambda_sq', 'lambda_sq_complement', 'lambda_ssp245_2041_2070', 'lambda_ssp245_2071_2100', 'lambda_ssp370_2041_2070', 'lambda_ssp370_2071_2100')
	
	conf <- configureMCMC(
	  model_species,
	  monitors = monitors,
	  print = TRUE,
	  enableWAIC = FALSE
	)

	# # add no U-turn sampler (Hamiltonian Monte Carlo)
	# conf$addSampler(target = c('beta', 'alpha0_sampling', 'alpha_area', 'alpha_poaceae', 'alpha_ag'), type = 'NUTS')

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
	  WAIC = FALSE,
	  perChainWAIC = FALSE
	)

	saveRDS(chains, paste0(out_dir, '/sdm_bernoulli_chains.rds'))

	say('session info', level = 2)
	print(sessionInfo())

	say(date(), pre = 1)
	sink()

say('#########################')
say('### model diagnostics ###')
say('#########################')

	chains <- readRDS(paste0(out_dir, '/sdm_bernoulli_chains.rds'))
	mcmc <- chains$samples

	for (i in 1:nchains) {
		# cols <- c(paste0('beta[', 1:n_sdm_terms, ']'), 'alpha0_sampling', 'alpha_area', 'alpha_poaceae', 'alpha_ag')
		cols <- c(paste0('beta[', 1:n_sdm_terms, ']'), 'alpha_area', 'alpha_poaceae', 'alpha_ag')
		mcmc[[i]] <- mcmc[[i]][ , cols]
	}

	# graphing trace plots for all betas
	pars <- paste0('beta[', 1:n_sdm_terms, ']')
	file <- paste0(out_dir, '/sdm_bernoulli_beta_trace.png')
	ggsave(mcmc_trace(mcmc, pars = pars), file = file, width = 10, height = 8, dpi = 450, bg = 'white')

	# graphing trace plots for all sampling coefficients
	# pars <- c('alpha0_sampling', 'alpha_area', 'alpha_poaceae', 'alpha_ag')
	pars <- c('alpha_area', 'alpha_poaceae', 'alpha_ag')
	file <- paste0(out_dir, '/sdm_bernoulli_alpha_trace.png')
	ggsave(mcmc_trace(mcmc, pars = pars), file = file, width = 10, height = 4, dpi = 450, bg = 'white')

	# graphing density plots for all betas
	pars <- paste0('beta[', 1:n_sdm_terms, ']')
	file <- paste0(out_dir, '/sdm_bernoulli_beta_density.png')
	ggsave(mcmc_dens_overlay(mcmc, pars = pars), file = file, width = 12, height = 8, dpi = 450, bg = 'white')

	# graphing trace plots for all sampling coefficients
	# pars <- c('alpha0_sampling', 'alpha_area', 'alpha_poaceae', 'alpha_ag')
	pars <- c('alpha_area', 'alpha_poaceae', 'alpha_ag')
	file <- paste0(out_dir, '/sdm_bernoulli_alpha_density.png')
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

	chains <- readRDS(paste0(out_dir, '/sdm_bernoulli_chains.rds'))

	# subset chain summary to just the lambdas
	summary <- chains$summary$all.chains

	# just counties with data
	which_lambda <- grepl(rownames(summary), pattern = 'lambda_sq')
	lambda <- summary[which_lambda, ]

	ag_vect_focus$lambda_mean <- lambda[ , 'Mean']
	ag_vect_focus$lambda_0.05ci <- lambda[ , '95%CI_low']
	ag_vect_focus$lambda_0.95ci <- lambda[ , '95%CI_upp']
	ag_vect_focus$lambda_ci <- ag_vect_focus$lambda_0.95ci - ag_vect_focus$lambda_0.05ci

	lambda_sq_quants <- quantile(ag_vect_focus$lambda_mean, c(0.25, 0.5, 0.75, 0.95))
	quant_labels <- c('[0 - 0.25)', '[0.25 - 0.5)', '[0.5 - 0.75)', '[0.75 - 0.95)', '[0.95-1]')
	
	ag_vect_focus$quant_col <- NA
	ag_vect_focus$quant_col[ag_vect_focus$lambda_mean < lambda_sq_quants[1]] <- quant_labels[1]
	ag_vect_focus$quant_col[ag_vect_focus$lambda_mean >= lambda_sq_quants[1] & ag_vect_focus$lambda_mean < lambda_sq_quants[2]] <- quant_labels[2]
	ag_vect_focus$quant_col[ag_vect_focus$lambda_mean >= lambda_sq_quants[2] & ag_vect_focus$lambda_mean < lambda_sq_quants[3]] <- quant_labels[3]
	ag_vect_focus$quant_col[ag_vect_focus$lambda_mean >= lambda_sq_quants[3] & ag_vect_focus$lambda_mean < lambda_sq_quants[4]] <- quant_labels[4]
	ag_vect_focus$quant_col[ag_vect_focus$lambda_mean >= lambda_sq_quants[4]] <- quant_labels[5]

	# complement of counties with data
	which_lambda <- grepl(rownames(summary), pattern = 'lambda_sq_complement')
	lambda <- summary[which_lambda, ]

	ag_vect_focus_complement$lambda_mean <- lambda[ , 'Mean']
	ag_vect_focus_complement$lambda_0.05ci <- lambda[ , '95%CI_low']
	ag_vect_focus_complement$lambda_0.95ci <- lambda[ , '95%CI_upp']
	ag_vect_focus_complement$lambda_ci <- ag_vect_focus_complement$lambda_0.95ci - ag_vect_focus_complement$lambda_0.05ci

	ag_vect_focus_complement$quant_col <- NA
	ag_vect_focus_complement$quant_col[ag_vect_focus_complement$lambda_mean < lambda_sq_quants[1]] <- quant_labels[1]
	ag_vect_focus_complement$quant_col[ag_vect_focus_complement$lambda_mean >= lambda_sq_quants[1] & ag_vect_focus_complement$lambda_mean < lambda_sq_quants[2]] <- quant_labels[2]
	ag_vect_focus_complement$quant_col[ag_vect_focus_complement$lambda_mean >= lambda_sq_quants[2] & ag_vect_focus_complement$lambda_mean < lambda_sq_quants[3]] <- quant_labels[3]
	ag_vect_focus_complement$quant_col[ag_vect_focus_complement$lambda_mean >= lambda_sq_quants[3] & ag_vect_focus_complement$lambda_mean < lambda_sq_quants[4]] <- quant_labels[4]
	ag_vect_focus_complement$quant_col[ag_vect_focus_complement$lambda_mean >= lambda_sq_quants[4]] <- quant_labels[5]

	nam <- vect(paste0(drive, '/Research Data/GADM/Version 4.1/High Res North America Level 1 sans Great Lakes SpatVector WGS84.gpkg'))
	nam <- project(nam, ag_vect_sq)

	extent <- ext(ag_vect_focus)
	extent <- as.vector(extent)
	x_range <- (extent[2] - extent[1])
	y_range <- (extent[4] - extent[3])
	extent[1] <- extent[1] + 0.15 * x_range
	extent[3] <- extent[3] + 0.125 * y_range
	extent[4] <- extent[4] - 0.2 * y_range

	cents_with_ag <- ag_vect_focus[ag_vect_focus$any_ag_quality1to3 > 0]
	cents_with_ag <- centroids(cents_with_ag)

	fill_scale <- c('[0 - 0.25)' = 'gray90', '[0.25 - 0.5)' = alpha('forestgreen', 0.15), '[0.5 - 0.75)' = alpha('forestgreen', 0.3), '[0.75 - 0.95)' = alpha('forestgreen', 0.6), '[0.95-1]' = 'forestgreen')

	map <- ggplot() +
		layer_spatial(ag_vect_focus, aes(fill = quant_col), color = NA) +
		layer_spatial(ag_vect_focus_complement, aes(fill = quant_col), color = NA) +
		layer_spatial(nam, color = 'gray30', fill = NA, linewidth = 0.1) +
		scale_fill_manual(
			name = 'Quantile\n of λ',
			values = fill_scale
		) +
		layer_spatial(cents_with_ag, pch = 16, color = alpha('gray20', 0.5), size = 0.35) +
		guides(fill = guide_legend(reverse = TRUE)) +
		xlim(extent[1], extent[2]) + ylim(extent[3], extent[4]) +
		ggtitle(expression('Present-day distribution of ' * italic('Andropogon gerardi')), subtitle = '1961-2020 | Bernoulli model') +
		theme(
			plot.title = element_text(size = 16),
			plot.subtitle = element_text(size = 14)
		)

	ggsave(plot = map, filename = paste0(out_dir, '/sdm_bernoulli_lambda_status_quo.png'), width = 12, height = 10, dpi = 600)

	writeVector(ag_vect_focus, paste0(out_dir, '/sdm_bernoulli_1961_2020_climate_focus.gpkg'), overwrite = TRUE)
	writeVector(ag_vect_focus_complement, paste0(out_dir, '/sdm_bernoulli_1961_2020_climate_focus_complement.gpkg'), overwrite = TRUE)

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

	chains <- readRDS(paste0(out_dir, '/sdm_bernoulli_chains.rds'))

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

		quant_labels <- c('[0 - 0.25)', '[0.25 - 0.5)', '[0.5 - 0.75)', '[0.75 - 0.95)', '[0.95-1]') # should be same as above
		# NB lambda_sq_quants is from above!!!
		this_ag_vect$quant_col <- NA
		this_ag_vect$quant_col[this_ag_vect$lambda_mean < lambda_sq_quants[1]] <- quant_labels[1]
		this_ag_vect$quant_col[this_ag_vect$lambda_mean >= lambda_sq_quants[1] & this_ag_vect$lambda_mean < lambda_sq_quants[2]] <- quant_labels[2]
		this_ag_vect$quant_col[this_ag_vect$lambda_mean >= lambda_sq_quants[2] & this_ag_vect$lambda_mean < lambda_sq_quants[3]] <- quant_labels[3]
		this_ag_vect$quant_col[this_ag_vect$lambda_mean >= lambda_sq_quants[3] & this_ag_vect$lambda_mean < lambda_sq_quants[4]] <- quant_labels[4]
		this_ag_vect$quant_col[this_ag_vect$lambda_mean >= lambda_sq_quants[4]] <- quant_labels[5]

		pretty_title <- paste0('SSP ', substr(fut, 4, 6), ' ', substr(fut, 8, 11), '-', substr(fut, 13, 16), ' | Bernoulli model')

		fill_scale <- c('[0 - 0.25)' = 'gray90', '[0.25 - 0.5)' = alpha('forestgreen', 0.15), '[0.5 - 0.75)' = alpha('forestgreen', 0.3), '[0.75 - 0.95)' = alpha('forestgreen', 0.6), '[0.95-1]' = 'forestgreen')
		
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
				

		ggsave(plot = map, filename = paste0(out_dir, '/sdm_bernoulli_lambda_', fut, '.png'), width = 12, height = 10, dpi = 600)

		writeVector(this_ag_vect, paste0(out_dir, '/sdm_bernoulli_', fut, '.gpkg'), overwrite = TRUE)

	} # next future

say('FINIS!', deco = '+', level = 1)

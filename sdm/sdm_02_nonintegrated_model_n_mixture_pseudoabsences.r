### MODELING ANDROPOGON GERARDI DISTRIBUTION, PHENOTYPE, PHYSIOLOGY, GENOTYPE, and ASSOCIATED MICROBIAL COMMUNITIES
### Erica Newman | Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2023-12
###
### This script constructs a non-integrated model for AG geographic distribution.
###
### source('C:/Ecology/R/andropogon_integratedEcology/sdm/sdm_02_nonintegrated_model_n_mixture_pseudoabsences.r')
### source('E:/Adam/R/andropogon_integratedEcology/sdm/sdm_02_nonintegrated_model_n_mixture_pseudoabsences.r')
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

	### user-defined values
	#######################

		psa_quant <- 0.99 # to define pseudoabsences, use this quantile of Poaceae occurrences across counties with Poaceae occurrences 

		### MCMC settings
		niter <- 130000
		nburnin <- 10000
		thin <- 120
		nchains <- 4

		# # for testing
		# niter <- 90
		# nburnin <- 10
		# thin <- 1
		# nchains <- 2

	out_dir <- paste0('./outputs_loretta/sdm_[nmixture]_[pseudoabsences_', psa_quant, ']_[climate]_[priors_ddexp]/')
	# out_dir <- paste0('./outputs_loretta/sdm_[nmixture]_[pseudoabsences_', psa_quant, ']_[climate_soil]_[priors_dnorm]/')
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

	# predictor_names <- c('aridity', 'bio7', 'bio2', 'sand', 'bio12', 'bio15', 'cec', 'soc', 'silt')
	# predictor_names <- c('aridity', 'bio7', 'ph', 'sand', 'bio12', 'bio15', 'cec', 'soc', 'silt')
	predictor_names <- c('aridity', 'bio7', 'bio12', 'bio15')
	say('We are using predictors: ', paste(predictor_names, collapse = ' '), post = 1)

	# preselected_model_terms <- read.csv('./outputs_loretta/sig_coeffs_elastic_net_2024_06_04.csv')

	# # get vector of linear predictors... we need to ensure all predictors appear at least as linear terms to respect marginality
	# terms <- preselected_model_terms$term
	# terms <- terms[!(terms %in% c('(Intercept)', 'log_area_km2'))]
	# linear_terms <- gsub(terms, pattern = 'I\\(', replacement = '')
	# linear_terms <- gsub(linear_terms, pattern = '\\^2)', replacement = '')
	# linear_terms <- strsplit(linear_terms, ':')
	# linear_terms <- unlist(linear_terms)
	# linear_terms <- unique(linear_terms)
	# linear_terms <- sort(linear_terms)


	# # get non-linear terms
	# terms <- terms[!(terms %in% linear_terms)]

	combos <- combn(predictor_names, 2, simplify = FALSE)

	form <- c(1, predictor_names, paste0('I(', predictor_names, '^2)'))
	for (i in seq_along(combos)) form <- c(form, paste(combos[[i]], collapse = ':'))
	form <- paste(form, collapse = ' + ')
	form <- paste0('~ ', form)

	# create model formula
	# form <- paste0(' ~ 1 + ', paste(c(linear_terms, terms), collapse = ' + '))
	say('Model formula:')
	say(form, post = 1, breaks = 80)
	form <- as.formula(form)

	### load AG data
	ag_vect_sq <- vect('./data_from_adam_and_loretta/andropogon_gerardi_occurrences_with_environment_1961_2020.gpkg')
	ag_vect_ssp245_2041_2070 <- vect('./data_from_adam_and_loretta/andropogon_gerardi_occurrences_with_environment_ensemble_8GCMs_ssp245_2041_2070.gpkg')
	ag_vect_ssp245_2071_2100 <- vect('./data_from_adam_and_loretta/andropogon_gerardi_occurrences_with_environment_ensemble_8GCMs_ssp245_2071_2100.gpkg')
	ag_vect_ssp370_2041_2070 <- vect('./data_from_adam_and_loretta/andropogon_gerardi_occurrences_with_environment_ensemble_8GCMs_ssp370_2041_2070.gpkg')
	ag_vect_ssp370_2071_2100 <- vect('./data_from_adam_and_loretta/andropogon_gerardi_occurrences_with_environment_ensemble_8GCMs_ssp370_2071_2100.gpkg')

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

	fields <- c('area_km2', 'any_ag_quality1to3', 'num_poaceae_records', predictor_names)
	ag_vect_sq <- ag_vect_sq[ , fields]

	# for counties with NA Poaceae and AG, assign a maximal number of Poaceae and 0 AG
	n_pseudoabs <- quantile(ag_vect_sq$num_poaceae_records, psa_quant, na.rm = TRUE)
	ag_vect_sq$num_poaceae_records[is.na(ag_vect_sq$num_poaceae_records)] <- n_pseudoabs
	ag_vect_sq$any_ag_quality1to3[is.na(ag_vect_sq$any_ag_quality1to3)] <- 0


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

	### inputs for nimble
	say('Inputs:', level = 2)
	data <- list(
		y = ag_sq$any_ag_quality1to3 # number of AG records in each county
	)

	n_counties <- nrow(ag_sq)
	n_counties_future <- nrow(x_ssp370_2071_2100)

	n_sdm_terms <- ncol(x_sq)

	constants <- list(
	  n_counties = n_counties,
	  n_counties_future = n_counties_future,

	  n_sdm_terms = n_sdm_terms,

	  x_sq = x_sq,

	  log_area_km2_scaled = log_area_km2_scaled,
	  log_num_poaceae_records_scaled = log_num_poaceae_records_scaled,

	  x_ssp245_2041_2070 = x_ssp245_2041_2070,
	  x_ssp245_2071_2100 = x_ssp245_2071_2100,
	  
	  x_ssp370_2041_2070 = x_ssp370_2041_2070,
	  x_ssp370_2071_2100 = x_ssp370_2071_2100

	)

	lambda_fut_inits <- rep(mean(ag_sq$any_ag_quality1to3), n_counties_future)
	inits <- list()
	for (i in seq_len(nchains)) {

		N_inits <- 10 * ag_sq$any_ag_quality1to3
		lambda_sq_inits <- 1 + ag_sq$any_ag_quality1to3

		inits[[i]] <- list(
			
			beta = runif(n_sdm_terms, -0.5, 0.5),
			alpha0_sampling = runif(1, -0.5, 0.5),
			alpha_area = runif(1, -0.5, 0.5),
			alpha_poaceae = runif(1, -0.5, 0.5),
			
			N = N_inits,

			lambda_sq = lambda_sq_inits,

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
	monitors <- c('beta', 'alpha0_sampling', 'alpha_area', 'alpha_poaceae', 'lambda_sq',		'lambda_ssp245_2041_2070', 'lambda_ssp245_2071_2100', 'lambda_ssp370_2041_2070', 'lambda_ssp370_2071_2100')
	
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
	  WAIC = FALSE,
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
	file <- paste0(out_dir, '/sdm_nmixture_beta_trace.png')
	ggsave(mcmc_trace(mcmc, pars = pars), file = file, width = 10, height = 8, dpi = 450, bg = 'white')

	# graphing trace plots for all "extra" betas
	pars <- c('alpha0_sampling', 'alpha_poaceae', 'alpha_area')
	file <- paste0(out_dir, '/sdm_nmixture_alpha_trace.png')
	ggsave(mcmc_trace(mcmc, pars = pars), file = file, width = 10, height = 4, dpi = 450, bg = 'white')

	# graphing density plots for all betas
	pars <- paste0('beta[', 1:n_sdm_terms, ']')
	file <- paste0(out_dir, '/sdm_nmixture_beta_density.png')
	ggsave(mcmc_dens_overlay(mcmc, pars = pars), file = file, width = 12, height = 8, dpi = 450, bg = 'white')

	# graphing trace plots for all "extra" betas
	pars <- c('alpha0_sampling', 'alpha_poaceae', 'alpha_area')
	file <- paste0(out_dir, '/sdm_nmixture_alpha_density.png')
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

	lambda_sq_quants <- quantile(ag_vect_sq$lambda_mean, c(0.25, 0.5, 0.75, 0.95))
	quant_labels <- c('[0 - 0.25)', '[0.25 - 0.5)', '[0.5 - 0.75)', '[0.75 - 0.95)', '[0.95-1]')
	
	ag_vect_sq$quant_col <- NA
	ag_vect_sq$quant_col[ag_vect_sq$lambda_mean < lambda_sq_quants[1]] <- quant_labels[1]
	ag_vect_sq$quant_col[ag_vect_sq$lambda_mean >= lambda_sq_quants[1] & ag_vect_sq$lambda_mean < lambda_sq_quants[2]] <- quant_labels[2]
	ag_vect_sq$quant_col[ag_vect_sq$lambda_mean >= lambda_sq_quants[2] & ag_vect_sq$lambda_mean < lambda_sq_quants[3]] <- quant_labels[3]
	ag_vect_sq$quant_col[ag_vect_sq$lambda_mean >= lambda_sq_quants[3] & ag_vect_sq$lambda_mean < lambda_sq_quants[4]] <- quant_labels[4]
	ag_vect_sq$quant_col[ag_vect_sq$lambda_mean >= lambda_sq_quants[4]] <- quant_labels[5]

	nam <- vect(paste0(drive, '/Research Data/GADM/Version 4.1/High Res North America Level 1 sans Great Lakes SpatVector WGS84.gpkg'))
	nam <- project(nam, ag_vect_sq)

	ag_vect_sq_pres <- ag_vect_sq[ag_vect_sq$any_ag_quality1to3 > 0]
	extent <- ext(ag_vect_sq_pres)
	extent <- as.vector(extent)
	x_range <- (extent[2] - extent[1])
	y_range <- (extent[4] - extent[3])
	extent[1] <- extent[1] + 0.15 * x_range
	extent[3] <- extent[3] + 0.125 * y_range
	extent[4] <- extent[4] - 0.2 * y_range

	cents_with_ag <- ag_vect_sq[ag_vect_sq$any_ag_quality1to3 > 0]
	cents_with_ag <- centroids(cents_with_ag)

	fill_scale <- c('[0 - 0.25)' = 'gray85', '[0.25 - 0.5)' = alpha('forestgreen', 0.25), '[0.5 - 0.75)' = alpha('forestgreen', 0.45), '[0.75 - 0.95)' = alpha('forestgreen', 0.7), '[0.95-1]' = 'forestgreen')

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

		quant_labels <- c('[0 - 0.25)', '[0.25 - 0.5)', '[0.5 - 0.75)', '[0.75 - 0.95)', '[0.95-1]') # should be same as above
		# NB lambda_sq_quants is from above!!!
		this_ag_vect$quant_col <- NA
		this_ag_vect$quant_col[this_ag_vect$lambda_mean < lambda_sq_quants[1]] <- quant_labels[1]
		this_ag_vect$quant_col[this_ag_vect$lambda_mean >= lambda_sq_quants[1] & this_ag_vect$lambda_mean < lambda_sq_quants[2]] <- quant_labels[2]
		this_ag_vect$quant_col[this_ag_vect$lambda_mean >= lambda_sq_quants[2] & this_ag_vect$lambda_mean < lambda_sq_quants[3]] <- quant_labels[3]
		this_ag_vect$quant_col[this_ag_vect$lambda_mean >= lambda_sq_quants[3] & this_ag_vect$lambda_mean < lambda_sq_quants[4]] <- quant_labels[4]
		this_ag_vect$quant_col[this_ag_vect$lambda_mean >= lambda_sq_quants[4]] <- quant_labels[5]

		pretty_title <- paste0('SSP ', substr(fut, 4, 6), ' ', substr(fut, 8, 11), '-', substr(fut, 13, 16), ' | N-mixture model')

		fill_scale <- c('[0 - 0.25)' = 'gray85', '[0.25 - 0.5)' = alpha('forestgreen', 0.25), '[0.5 - 0.75)' = alpha('forestgreen', 0.45), '[0.75 - 0.95)' = alpha('forestgreen', 0.7), '[0.95-1]' = 'forestgreen')
		
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


say('FINIS!', deco = '+', level = 1)

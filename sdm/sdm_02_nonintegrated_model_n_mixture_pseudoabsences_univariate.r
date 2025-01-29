### MODELING ANDROPOGON GERARDI DISTRIBUTION, PHENOTYPE, PHYSIOLOGY, GENOTYPE, and ASSOCIATED MICROBIAL COMMUNITIES
### Erica Newman | Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2023-12
###
### This script constructs a non-integrated model for AG geographic distribution.
###
### source('C:/Ecology/R/andropogon_integratedEcology/sdm/sdm_02_nonintegrated_model_n_mixture_pseudoabsences_univariate.r')
### source('C:/Subarashi/R/andropogon_integratedEcology/sdm/sdm_02_nonintegrated_model_n_mixture_pseudoabsences_univariate.r')
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

	drive <- 'C:/Ecology/'
	# drive <- 'C:/Subarashi/'

	setwd(paste0(drive, '/Research/Andropogon/Andropogon'))

	library(bayesplot) # graphing
	library(coda) # Bayesian diagnostics
	library(cowplot) # combining ggplots
	library(enmSdmX) # GIS and SDMing
	library(ggplot2) # plotting
	library(ggspatial) # plotting spatial things
	library(nimble) # Bayes
	# library(nimbleHMC) # Hamiltonian Monte Carlo samplers
	library(omnibus) # helper functions
	library(readxl) # open Excel files
	library(scales) # for plotting transparency
	library(terra) # spatial objects

	### user-defined values
	#######################

		# to define pseudoabsences, use this quantile of Poaceae occurrences across counties with Poaceae occurrences 
		# psa_quant <- 0.99
		psa_quant <- 0 

		### MCMC settings
		
		# for all variables except aridity
		predictor_names <- c('bio1', 'bio5', 'bio7', 'bio12', 'bio15')
		niter <- 250000
		nburnin <- 50000
		thin <- 200
		nchains <- 4

		# # for aridity
		# predictor_names <- 'aridity'
		# niter <- 1100000
		# nburnin <- 100000
		# thin <- 1000
		# nchains <- 4

		# # for testing
		# niter <- 90
		# nburnin <- 10
		# thin <- 1
		# nchains <- 2
		
		# to generate range across which to calculate response, expand the range by this proportion of observed range on the lower and upper sides
		predictor_range_expand_factor <- 0.1

		# number of values along which to estimate response curves
		n_response_curve_rows <- 800

	# out_dir <- paste0('./outputs_loretta/sdm_univariate_[TEMP]/')
	# out_dir <- paste0('./outputs_loretta/sdm_univariate_[nmixture]_[pseudoabsences_', psa_quant, ']_[climate]_[priors_ddexp]/')
	out_dir <- paste0('./outputs_loretta/sdm_univariate_[nmixture]_[pseudoabsences_', psa_quant, ']_[climate]_[priors_dnorm]/')
	dirCreate(out_dir)

# say('######################')
# say('### univariate SDM ###')
# say('######################')

	# ### load AG data
	# ag_vect_sq <- vect('./data_from_adam_and_loretta/andropogon_gerardi_occurrences_with_environment_1961_2020.gpkg')

	# fields <- c('area_km2', 'any_ag_quality_1_to_3', 'num_poaceae_records', predictor_names)
	# ag_vect_sq <- ag_vect_sq[ , fields]

	# # for counties with NA Poaceae and AG, assign a maximal number of Poaceae and 0 AG
	# n_pseudoabs <- quantile(ag_vect_sq$num_poaceae_records, psa_quant, na.rm = TRUE)
	# ag_vect_sq$num_poaceae_records[is.na(ag_vect_sq$num_poaceae_records)] <- n_pseudoabs
	# ag_vect_sq$any_ag_quality_1_to_3[is.na(ag_vect_sq$any_ag_quality_1_to_3)] <- 0

	# ### collate data
	# ag_sq <- as.data.frame(ag_vect_sq)

	# ### county area... used as an offset to make underlying model fit an IPP
	# area_km2 <- ag_sq$area_km2
	# log_area_km2 <- log(area_km2)
	# log_area_km2_scaled <- scale(log_area_km2)
	# log_area_center <- attributes(log_area_km2_scaled)$`scaled:center`
	# log_area_scale <- attributes(log_area_km2_scaled)$`scaled:scale`
	# log_area_km2_scaled <- log_area_km2_scaled[ , 1]

	# ### number of Poaceae records... used to model sampling bias
	# log_num_poaceae_records <- log1p(ag_sq$num_poaceae_records) # log(x_sq + 1)
	# log_num_poaceae_records_scaled <- scale(log_num_poaceae_records)
	# log_num_poaceae_records_scaled <- log_num_poaceae_records_scaled[ , 1]

	# # # ### environment at sampled sites
	# # # site_data <- read_xlsx('./data_from_loretta/!plant_sitelevel_data_24OCT2024 - Google Sheets [aggregated by Erica] before adding metadata 2024-12-28.xlsx', sheet = 'plantmaster_bysite_19NOV2024')
	
	# # # site_vect <- vect(site_data, geom = c('LONGITUDE', 'LATITUDE'), crs = getCRS('WGS84'))
	# # # env_rasts <- rast(paste0(drive, '/Research Data/ClimateNA/v 7.3 AdaptWest/1961-2020/bioclim_variables_1961_2020.tif'))
	
	# # # site_vect <- project(site_vect, env_rasts)
	# # # env_at_samples <- extract(env_rasts, site_vect)
		
# 	for (predictor_name in predictor_names) {

# 		say(predictor_name, level = 1)

# 		sink(paste0(out_dir, '/', predictor_name, '_runtime_log.txt'), split = TRUE)

# 		say(predictor_name)
# 		say(date(), post = 1)

# 		say('This project creates a series of univariate climatic niche models for AG. Occurrences are at the county level, so county area is used as an offset. Number of Poaceae samples, interpreted as an index of sampling effort, is also used as a offset.', breaks = 60, post = 1)

# 		say('MCMC settings:', level = 2)
# 		say('niter ....... ', niter)
# 		say('nburnin ..... ', nburnin)
# 		say('thin ........ ', thin)
# 		say('nchains ..... ', nchains)
# 		say('psa_quant ... ', psa_quant, post = 1)

# 		# create model formula
# 		# form <- paste0(' ~ 1 + ', paste(c(linear_terms, terms), collapse = ' + '))
# 		say('Model formula:')
# 		form <- paste0(' ~ ', predictor_name, ' + I(', predictor_name, '^2)')
# 		say(form, post = 1, breaks = 80)
# 		form <- as.formula(form)

# 		### subset, scale, and manipulate predictors into model frame
# 		# status quo: counties with training data
# 		x_raw <- as.data.frame(ag_sq[ , c('area_km2', predictor_name)])
# 		log_area_km2 <- log(x_raw$area_km2)
# 		log_area_km2 <- data.frame(log_area_km2 = log_area_km2)
# 		x_raw$area_km2 <- NULL
# 		x_raw <- cbind(log_area_km2, x_raw)
# 		x_raw_scaled <- scale(x_raw)
# 		x_centers <- attr(x_raw_scaled, 'scaled:center')
# 		x_scales <- attr(x_raw_scaled, 'scaled:scale')
# 		x_raw_scaled <- as.data.frame(x_raw_scaled)
# 		x_sq <- model.matrix(form, as.data.frame(x_raw_scaled))

# 		### response array
# 		##################

# 		# This is a matrix for plotting the response curves of the SDM
# 		# columns: variables, each increasing from lowest to highest value across counties with morphological/physiological samples
# 		# rows: focal variable increments, others held constant

# 		pred_values_at_sites <- x_raw[ , predictor_name]
		
# 		pred_range_at_sites <- range(pred_values_at_sites)
# 		pred_range_at_sites <- diff(pred_range_at_sites)
# 		pred_lower_at_sites <- min(pred_values_at_sites) - predictor_range_expand_factor * pred_range_at_sites
# 		pred_upper_at_sites <- max(pred_values_at_sites) + predictor_range_expand_factor * pred_range_at_sites

# 		if (predictor_name %in% c('bio7', 'bio12', 'bio15')) {
# 			pred_lower_at_sites <- max(0, pred_lower_at_sites)
# 		}

# 		# number of rows in array--represents number of points along which we'll be able to make the response curve
# 		n_predictors <- length(predictor_names)

# 		sequence <- seq(pred_lower_at_sites, pred_upper_at_sites, length.out = n_response_curve_rows)
# 		sequence_scaled <- (sequence - x_centers[predictor_name]) / x_scales[predictor_name]
# 		env_array <- matrix(sequence, nrow = n_response_curve_rows, ncol = 1)
# 		env_array_scaled <- matrix(sequence_scaled, nrow = n_response_curve_rows, ncol = 1)
# 		colnames(env_array) <- colnames(env_array_scaled) <- predictor_name

# 		# make into model matrix
# 		response_curve_x <- model.matrix(form, as.data.frame(env_array_scaled))
		
# 	### inputs
# 	##########

# 	say('inputs:', level = 2)
	
# 		data <- list(
# 			y = ag_sq$any_ag_quality_1_to_3 # number of AG records in each county
# 		)

# 		n_counties <- nrow(ag_sq)
# 		n_sdm_terms <- ncol(x_sq)

# 		constants <- list(
# 			n_counties = n_counties,

# 			n_sdm_terms = n_sdm_terms,

# 			x_sq = x_sq,

# 			n_response_curve_rows = n_response_curve_rows,
# 			response_curve_x = response_curve_x,

# 			log_area_km2_scaled = log_area_km2_scaled,
# 			log_num_poaceae_records_scaled = log_num_poaceae_records_scaled

# 		)

# 		N_inits <- 10 * ag_sq$any_ag_quality_1_to_3
# 		lambda_sq_inits <- 1 + ag_sq$any_ag_quality_1_to_3

# 		inits <- list(
			
# 			beta = rep(0, n_sdm_terms),
# 			alpha0_sampling = 0,
# 			alpha_area = 0,
# 			alpha_poaceae = 0,
			
# 			N = N_inits,

# 			lambda_sq = lambda_sq_inits
			
# 		)

# 		say('Data:')
# 		print(str(data))

# 		say('Constants:', pre = 1)
# 		print(str(constants))

# 		say('Initializations:', pre = 1)
# 		print(str(inits))

# 		### define model
# 		say('nimbleCode():', level = 2)
# 		say('We assume an N-mixture model (latent, real abundance ~ Poisson, and observations ~ binomial draws from latent abundance.', post = 2)
# 		model_code <- nimbleCode({
		  
# 			# priors for relationship to environment
# 			beta[1] ~ dnorm(0, sd = 10) # intercept... not regularized
# 			for (j in 2:n_sdm_terms) {
# 				# beta[j] ~ ddexp(0, rate = 1) # regularization toward 0
# 				beta[j] ~ dnorm(0, sd = 10) # broad prior
# 			}

# 			# priors for sampling bias
# 			alpha0_sampling ~ dnorm(0, sd = 10)
# 			alpha_area ~ dnorm(0, sd = 10)
# 			alpha_poaceae ~ dnorm(0, sd = 10)

# 			# likelihood
# 			for (i in 1:n_counties) {    # this specifies estimates be made for each county
				
# 				### actual abundance (latent--unobserved)
# 				N[i] ~ dpois(lambda_sq[i])

# 				# relationship between latent abundance and environment
# 				log(lambda_sq[i]) <- inprod(beta[1:n_sdm_terms], x_sq[i, 1:n_sdm_terms])

# 				### observed number of AG
# 				y[i] ~ dbin(prob = p[i], size = N[i])

# 				# sampling bias
# 				logit(p[i]) <- alpha0_sampling + alpha_area * log_area_km2_scaled[i] + alpha_poaceae * log_num_poaceae_records_scaled[i]

# 			}

# 			# posterior predictive sampler for response curves
# 			for (i in 1:n_response_curve_rows) {
				
# 				log(response_curves_sdm_lambda[i]) <-
# 					inprod(beta[1:n_sdm_terms], response_curve_x[i, 1:n_sdm_terms])

# 			}
			
# 		})

# 		print(model_code)

# 		say('nimbleModel():', level = 2)
# 		model_species <- nimbleModel(
# 			code = model_code, # our model
# 			constants = constants, # constants
# 			data = data, # data
# 			inits = inits, # initialization values
# 			check = TRUE, # any errors?
# 			calculate = FALSE
# 			# buildDerivs = TRUE # need for Hamiltonian Monte Carlo
# 		)

# 		say('$initializeInfo() and $calculate():', level = 2)
# 		model_species$initializeInfo()
# 		model_species$calculate()

# 		say('configureMCMC():', level = 2)
# 		monitors <- c(
# 			'beta', 'alpha0_sampling', 'alpha_area', 'alpha_poaceae',
# 			'response_curves_sdm_lambda'
# 		)

# 		conf <- configureMCMC(
# 			model_species,
# 			monitors = monitors,
# 			print = TRUE,
# 			enableWAIC = FALSE
# 		)

# 		# # add no U-turn sampler (Hamiltonian Monte Carlo)
# 		# conf$addSampler(target = c('beta', 'alpha_area', 'alpha_poaceae'), type = 'NUTS')

# 		### compile/build/run model/save MCMC
# 		build <- buildMCMC(conf)

# 		compiled <- compileNimble(model_species, build, showCompilerOutput = FALSE)

# 		chains <- runMCMC(
# 			compiled$build,
# 			niter = niter,
# 			nburnin = nburnin,
# 			thin = thin,
# 			nchains = nchains,
# 			inits = inits,
# 			progressBar = TRUE,
# 			samplesAsCodaMCMC = TRUE,
# 			summary = TRUE,
# 			WAIC = FALSE,
# 			perChainWAIC = FALSE
# 		)

# 		saveRDS(chains, paste0(out_dir, '/', predictor_name, '_sdm_nmixture_chains.rds'))

# 	### model diagnostics
# 	#####################

# 	say('model diagnostics', level = 2)

# 		# chains <- readRDS(paste0(out_dir, '/', predictor_name, '_sdm_nmixture_chains.rds'))
# 		mcmc <- chains$samples

# 		for (i in 1:nchains) {
# 			cols <- c(paste0('beta[', 1:n_sdm_terms, ']'), 'alpha0_sampling', 'alpha_poaceae', 'alpha_area')
# 			mcmc[[i]] <- mcmc[[i]][ , cols]
# 		}

# 		# graphing trace plots for all betas
# 		pars <- paste0('beta[', 1:n_sdm_terms, ']')
# 		file <- paste0(out_dir, '/', predictor_name, '_sdm_nmixture_beta_trace.png')
# 		ggsave(mcmc_trace(mcmc, pars = pars), file = file, width = 14, height = 6, dpi = 450, bg = 'white')

# 		# graphing trace plots for all alphas
# 		pars <- c('alpha0_sampling', 'alpha_poaceae', 'alpha_area')
# 		file <- paste0(out_dir, '/', predictor_name, '_sdm_nmixture_alpha_trace.png')
# 		ggsave(mcmc_trace(mcmc, pars = pars), file = file, width = 14, height = 6, dpi = 450, bg = 'white')

# 		# graphing density plots for all betas
# 		pars <- paste0('beta[', 1:n_sdm_terms, ']')
# 		file <- paste0(out_dir, '/', predictor_name, '_sdm_nmixture_beta_density.png')
# 		ggsave(mcmc_dens_overlay(mcmc, pars = pars), file = file, width = 14, height = 6, dpi = 450, bg = 'white')

# 		# graphing density plots for all alphas
# 		pars <- c('alpha0_sampling', 'alpha_poaceae', 'alpha_area')
# 		file <- paste0(out_dir, '/', predictor_name, '_sdm_nmixture_alpha_density.png')
# 		ggsave(mcmc_dens_overlay(mcmc, pars = pars), file = file, width = 14, height = 6, dpi = 450, bg = 'white')

# 		# Gelman-Rubin statistic
# 		rhats <- gelman.diag(mcmc, autoburnin = FALSE, multivariate = TRUE)

# 		say('GELMAN-RUBIN STATISTICS')
# 		print(rhats)

# 	### parameter estimates
# 	#######################

# 	say('parameter estimates', level = 2)

# 		# chains <- readRDS(paste0(out_dir, '/', predictor_name, '_sdm_nmixture_chains.rds'))

# 		# subset chain summary to just the lambdas
# 		summary <- chains$summary$all.chains

# 		model_term <- attr(terms(form), 'term.labels')
# 		model_term_nice <- model_term
# 		model_term_nice <- gsub(model_term_nice, pattern = 'I\\(', replacement = '')
# 		model_term_nice <- gsub(model_term_nice, pattern = '\\)', replacement = '')
# 		model_term_nice <- c('Intercept', model_term_nice)

# 		pars <- paste0('beta[', 1:(1 + length(model_term)), ']')

# 		coeff_chains <- chains$samples
# 		for (i in 1:nchains) {
# 			coeff_chains[[i]] <- coeff_chains[[i]][ , pars]
# 			colnames(coeff_chains[[i]]) <- model_term_nice
# 		}

# 		caterpillars <- mcmc_intervals(coeff_chains, pars = model_term_nice)

# 		# plot posterior of coefficient estimates
# 		caterpillars <- caterpillars +
# 			xlab('Estimated value') +
# 			theme(plot.background = element_rect(fill = 'white'))

# 		ggsave(caterpillars, filename = paste0(out_dir, '/', predictor_name, '_sdm_nmixture_coefficients.png'), width = 10, height = 12, dpi = 300)

# 	### session information
# 	#######################

# 	say('session info', level = 2)
	
# 		print(sessionInfo())

# 		say(date(), pre = 1)
# 		sink()

# 	} # next predictor

say('##########################################################')
say('### plot response curves with mean and extreme markers ###')
say('##########################################################')

	predictor_names <- c('aridity', 'bio1', 'bio5', 'bio7', 'bio12', 'bio15')

	ag_vect_sq <- vect('./data_from_adam_and_loretta/andropogon_gerardi_occurrences_with_environment_1961_2020.gpkg')

	ag_vect_sq <- ag_vect_sq[ , c('any_ag_quality_1_to_3', predictor_names)]
	ag_sq <- as.data.frame(ag_vect_sq)

	# value of predictors at sampled sites
	env_at_sites <- read.csv('./data_from_loretta/!admixture_data_collated/admixture_with_env.csv')
	env_at_sites <- env_at_sites[!duplicated(env_at_sites$population), ]
	env_at_sites <- env_at_sites[ , c('longitude', 'latitude')]

	env_at_sites <- terra::vect(env_at_sites, geom = c('longitude', 'latitude'), crs = '+proj=longlat +datum=WGS84')
	bcs <- rast(paste0(drive, '/Research Data/ClimateNA/v 7.3 AdaptWest/1961-2020/bioclim_variables_1961_2020.tif'))
	env_at_sites <- project(env_at_sites, bcs)

	env_at_sites <- extract(bcs, env_at_sites, ID = FALSE)
	env_at_sites$aridity <- (env_at_sites$bio1 + 10) / (env_at_sites$bio12 / 1000)

	responses <- list() # one element per plot
	collated <- data.frame() # DF to hold all response curve data for exporting
	for (predictor_name in predictor_names) {
	
		say(predictor_name)

		chains <- readRDS(paste0(out_dir, '/', predictor_name, '_sdm_nmixture_chains.rds'))
	
		# get raw predictor values
		### subset, scale, and manipulate predictors into model frame
		# status quo: counties with training data
		pred_values_at_sites <- ag_sq[ , predictor_name]
		
		pred_range_at_sites <- range(pred_values_at_sites)
		pred_range_at_sites <- diff(pred_range_at_sites)
		pred_lower_at_sites <- min(pred_values_at_sites) - predictor_range_expand_factor * pred_range_at_sites
		pred_upper_at_sites <- max(pred_values_at_sites) + predictor_range_expand_factor * pred_range_at_sites

		if (predictor_name %in% c('bio7', 'bio12', 'bio15')) {
			pred_lower_at_sites <- max(0, pred_lower_at_sites)
		}

		# number of rows in array--represents number of points along which we'll be able to make the response curve
		sequence <- seq(pred_lower_at_sites, pred_upper_at_sites, length.out = n_response_curve_rows)

		# for response along this variable, get mean predicted value and 95% CI's
		mean <- chains$summary$all.chains[paste0('response_curves_sdm_lambda[', 1:n_response_curve_rows, ']'), 'Mean']
		lower <- chains$summary$all.chains[paste0('response_curves_sdm_lambda[', 1:n_response_curve_rows, ']'), '95%CI_low']
		upper <- chains$summary$all.chains[paste0('response_curves_sdm_lambda[', 1:n_response_curve_rows, ']'), '95%CI_upp']

		# collate into data frame
		response <- data.frame(
			predictor = predictor_name,
			x = sequence,
			lambda = mean,
			type = 'mean'
		)

		n <- nrow(response)
		response <- rbind(
			response,
			data.frame(
				predictor = predictor_name,
				x = c(sequence, rev(sequence)),
				lambda = c(upper, rev(lower)),
				type = c(rep('upper CI', n), rep('lower CI', n))
			)
		)

		collated <- rbind(collated, response)

		# limit x-range
		max_value <- max(response$lambda[response$type == 'mean'])
		min_value <- min(response$lambda[response$type == 'mean'])
		cutoff_value <- max_value * 0.01
		response <- response[response$lambda >= cutoff_value, ]

		# limit y-range
		max_value <- max(response$lambda[response$type == 'upper CI'])
		ylim <- c(0, 1.05 * max_value)

		if (predictor_name == 'aridity') {
			nice_title <- 'Aridity ((temperature + 10) / (precipitation / 1000))'
			nice_axis <- 'Aridity'
		} else if (predictor_name == 'bio1') {
			nice_title <- 'Mean annual temperature (BIO01)'
			nice_axis <- 'Mean annual temperature (°C)'
		} else if (predictor_name == 'bio5') {
			nice_title <- 'Temperature of the hottest month (BIO05)'
			nice_axis <- 'Temperature of the hottest month (°C)'
		} else if (predictor_name == 'bio7') {
			nice_title <- 'Temperature annual range (BIO07)'
			nice_axis <- 'Temperature annual range (°C)'
		} else if (predictor_name == 'bio12') {
			nice_title <- 'Total annual precipitation (BIO12)'
			nice_axis <- 'Total annual precipitation (mm)'
		} else if (predictor_name == 'bio15') {
			nice_title <- 'Precipitation seasonality (BIO15)'
			nice_axis <- 'Precipitation seasonality'
		} else if (predictor_name == 'ph') {
			nice_title <- 'Soil pH'
			nice_axis <- 'pH'
		} else if (predictor_name == 'sand') {
			nice_title <- 'Soil proportion sand'
			nice_axis <- 'Proportion sand'
		} else if (predictor_name == 'silt') {
			nice_title <- 'Soil proportion silt'
			nice_axis <- 'Proportion silt'
		} else if (predictor_name == 'soc') {
			nice_title <- 'Soil organic matter'
			nice_axis <- 'Soil organic matter'
		}

		# beta1 <- chains$summary$all.chains['beta[1]', 'Mean']
		# beta2 <- chains$summary$all.chains['beta[2]', 'Mean']
		# beta3 <- chains$summary$all.chains['beta[3]', 'Mean']

		# maximum
		max_value_at <- which.max(response$lambda[response$type == 'mean'])
		max_at_x <- response$x[max_value_at]

		# SD

		weighted_mean_x <- sum(response$lambda[response$type == 'mean'] * response$x[response$type == 'mean']) / sum(response$lambda[response$type == 'mean'])

		w <- response$lambda[response$type == 'mean']
		x <- response$x[response$type == 'mean']
		M <- length(x)
		sd <- sqrt( sum(w * (x - weighted_mean_x)^2) / (((M - 1) / M) * sum(w)))

		# sd <- sqrt(sum((response$x[response$type == 'mean'] - weighted_mean_x)^2) / (sum(response$type == 'mean') - 1))

		# rug plot of occurrences
		ag_rug <- data.frame(x = ag_sq[ag_sq$any_ag_quality_1_to_3 > 0, predictor_name])

		# rug plot of sampled sites
		sampled_sites_rug <- data.frame(x = env_at_sites[ , predictor_name])

		responses[[length(responses) + 1]] <- ggplot(
				response[response$type == 'mean', ],
				aes(x = x, y = lambda),
				color = 'forestgreen'
			) +
			geom_line() +
			geom_polygon(
				data = response[response$type %in% c('upper CI', 'lower CI'), ],
				mapping = aes(x = x, y = lambda),
				color = NA,
				fill = 'forestgreen',
				alpha = 0.2
			) +
			geom_rug(data = ag_rug, aes(x = x, y = x), sides = 'b') +
			geom_rug(data = sampled_sites_rug, aes(x = x, y = x), sides = 't', color = 'forestgreen') +
			ylim(ylim[1], ylim[2]) +
			xlab(nice_axis) +
			ylab('Lambda') +
			ggtitle(nice_title) +
			theme(
				plot.title = element_text(size = 11),
				axis.title = element_text(size = 9.5)
			)

			if (predictor_name != 'bio7') {
			
				responses[[length(responses)]] <- responses[[length(responses)]] +
					geom_vline(xintercept = max_at_x, color = 'red', linetype = 'solid') +
					geom_vline(xintercept = max_at_x + sd, color = 'red', linetype = 'dashed') +
					geom_vline(xintercept = max_at_x - sd, color = 'red', linetype = 'dashed')
			
			}

	}

	responses <- plot_grid(plotlist = responses, nrow = 2)

	ggsave(plot = responses, filename = paste0(out_dir, '/univariate_sdm_response_curve_estimates.png'), width = 12, height = 6, dpi = 600)		
	write.csv(collated, paste0(out_dir, '/univariate_sdm_response_curve_estimates.csv'), row.names = FALSE)

say('FINIS!', deco = '+', level = 1)

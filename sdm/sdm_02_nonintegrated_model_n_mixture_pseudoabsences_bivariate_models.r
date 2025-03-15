### MODELING ANDROPOGON GERARDI DISTRIBUTION, PHENOTYPE, PHYSIOLOGY, GENOTYPE, and ASSOCIATED MICROBIAL COMMUNITIES
### Adam B. Smith | Erica Newman | Missouri Botanical Garden | adam.smith@mobot.org | 2023-12
###
### This script constructs a non-integrated model for AG geographic distribution.
###
### source('C:/Ecology/R/andropogon_integratedEcology/sdm/sdm_02_nonintegrated_model_n_mixture_pseudoabsences_bivariate_models.r')
### source('C:/Subarashi/R/andropogon_integratedEcology/sdm/sdm_02_nonintegrated_model_n_mixture_pseudoabsences_bivariate_models.r')
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

		# names of predictors for which to construct bivariate models
		# all_predictor_names <- c('bio1', 'bio5', 'bio6', 'bio7', 'bio12', 'bio15', 'bio18', 'aridity')
		# all_predictor_names <- c('bio18', 'aridity')
		# all_predictor_names <- c('bio12', 'aridity')
		# all_predictor_names <- c('bio6', 'bio12')
		# all_predictor_names <- c('bio5', 'bio18')
		all_predictor_names <- c('bio1', 'bio12')

		# predictors with correlation above this threshold will not be modeled together
		correlation_threshold <- 0.7

		psa_quant <- 0.99 # to define pseudoabsences, use this quantile of Poaceae occurrences across counties with Poaceae occurrences 
		# psa_quant <- 0 # to define pseudoabsences, use this quantile of Poaceae occurrences across counties with Poaceae occurrences 

		### MCMC settings: for stubborn pairs
		niter <- 8 * 240000
		nburnin <- 8 * 40000
		thin <- 8 * 200
		nchains <- 4
		trial <- FALSE

		# ### MCMC settings: default for all variable pairs
		# niter <- 240000
		# nburnin <- 40000
		# thin <- 200
		# nchains <- 4
		# trial <- FALSE

		# # for testing
		# niter <- 90
		# nburnin <- 10
		# thin <- 1
		# nchains <- 2
		# trial <- TRUE

		# out_dir <- paste0('./outputs_loretta/sdm_[nmixture]_[pseudoabsences_', psa_quant, ']_[climate_bivariate_quad_ia]_[priors_dnorm]', ifelse(trial, '_TRIAL', ''), '/')
		out_dir <- paste0('./outputs_loretta/sdm_[nmixture]_[pseudoabsences_', psa_quant, ']_[climate_bivariate_quad]_[priors_dnorm]', ifelse(trial, '_TRIAL', ''), '/')
		dirCreate(out_dir)

		ag_vect_sq <- vect('./data_from_adam_and_loretta/andropogon_gerardi_occurrences_with_environment_1961_2020.gpkg')

#######################################################################
### cycle through each pairwise combination of predictors and model ###
#######################################################################

	for (count_predictor_1 in 1:(length(all_predictor_names) - 1)) {
	# for (count_predictor_1 in 1) {

		predictor_name_1 <- all_predictor_names[count_predictor_1]

		for (count_predictor_2 in (count_predictor_1 + 1):length(all_predictor_names)) {

			predictor_name_2 <- all_predictor_names[count_predictor_2]

			predictor_names <- c(predictor_name_1, predictor_name_2)

			x <- ag_vect_sq[[c(predictor_name_1, predictor_name_2)]]
			correlation <- cor(x, method = 'spearman', use = 'complete.obs')[1, 2]

			if (abs(correlation) <= correlation_threshold) {

				out_tag <- paste0(predictor_name_1, '_', predictor_name_2)

				sink(paste0(out_dir, '/runtime_', out_tag, '.txt'), split = TRUE)
				say('BIVARIATE MODELS')
				say(date(), post = 2)
				say('predictor 1 ............. ', predictor_name_1)
				say('predictor 2 ............. ', predictor_name_2)
				say('niter ................... ', niter)
				say('nburnin ................. ', nburnin)
				say('thin .................... ', thin)
				say('nchains ................. ', nchains)
				say('psa_quant ............... ', psa_quant)
				say('correlation_threshold ... ', correlation_threshold, post = 1)

				# create model formula
				# form <- c(1, predictor_name_1, predictor_name_2, paste0('I(', predictor_names, '^2)'), paste0(predictor_names, collapse = ':'))
				form <- c(1, predictor_name_1, predictor_name_2, paste0('I(', predictor_names, '^2)'))
				form <- paste(form, collapse = ' + ')
				form <- paste0('~ ', form)

				say('Model formula:')
				say(form, post = 1, breaks = 80)
				form <- as.formula(form)

				fields <- c('area_km2', 'any_ag_quality_1_to_3', 'num_poaceae_records', predictor_names)
				this_ag_vect_sq <- ag_vect_sq[ , fields]

				# for counties with NA Poaceae and AG, assign a maximal number of Poaceae and 0 AG
				n_pseudoabs <- quantile(this_ag_vect_sq$num_poaceae_records, psa_quant, na.rm = TRUE)
				this_ag_vect_sq$num_poaceae_records[is.na(this_ag_vect_sq$num_poaceae_records)] <- n_pseudoabs
				this_ag_vect_sq$any_ag_quality_1_to_3[is.na(this_ag_vect_sq$any_ag_quality_1_to_3)] <- 0

				### collate data
				ag_sq <- as.data.frame(this_ag_vect_sq)

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

				### response array
				##################

				# This is an array for plotting the response curves of the SDM
				# columns: variables, all but one held at mean value across counties with AG
				# rows: focal variable increments, others held constant
				# depths: identity of focal variable changes

				# number of rows in array--represents number of points along which we'll be able to make the response curve
				n_predictors <- 2
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
				ag_pres <- this_ag_vect_sq[this_ag_vect_sq$any_ag_quality_1_to_3 > 0]
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

				n_sdm_terms <- ncol(x_sq)

				constants <- list(
					n_counties = n_counties,
					n_sdm_terms = n_sdm_terms,
					x_sq = x_sq,

					n_response_curve_rows = n_response_curve_rows,
					n_predictors = n_predictors,
					response_curve_x = response_curve_x,

					log_area_km2_scaled = log_area_km2_scaled,
					log_num_poaceae_records_scaled = log_num_poaceae_records_scaled

				)

				N_inits <- 10 * ag_sq$any_ag_quality_1_to_3
				lambda_sq_inits <- 1 + ag_sq$any_ag_quality_1_to_3

				beta_inits <- rep(0, n_sdm_terms)
				beta_inits[grepl(colnames(x_sq), pattern = '\\^2')] <- -2
				beta_inits[grepl(colnames(x_sq), pattern = '\\:')] <- 0

				inits <- list(

					beta = beta_inits,
					alpha_0 = -1,
					alpha_area = 2,
					alpha_poaceae = 2,
					
					N = N_inits,

					lambda_sq = lambda_sq_inits

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
					alpha_0 ~ dnorm(0, sd = 10)
					alpha_area ~ dnorm(0, sd = 10)
					alpha_poaceae ~ dnorm(0, sd = 10)

					# likelihood
					for (i in 1:n_counties) {
						
						### actual abundance (latent--unobserved)
						N[i] ~ dpois(lambda_sq[i])

						# relationship between latent abundance and environment
						log(lambda_sq[i]) <- inprod(beta[1:n_sdm_terms], x_sq[i, 1:n_sdm_terms])

						### observed number of AG
						y[i] ~ dbin(prob = p[i], size = N[i])

						# sampling bias
						logit(p[i]) <- alpha_0 + alpha_area * log_area_km2_scaled[i] + alpha_poaceae * log_num_poaceae_records_scaled[i]

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
					'beta', 'alpha_0', 'alpha_area', 'alpha_poaceae',
					'exp_response_curves_sdm_lambda',
					'lambda_sq'
				)
				
				conf <- configureMCMC(
					model_species,
					monitors = monitors,
					print = TRUE,
					enableWAIC = TRUE
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

				saveRDS(chains, paste0(out_dir, '/chains_', out_tag, '.rds'))

			say('#########################')
			say('### model diagnostics ###')
			say('#########################')

				mcmc <- chains$samples

				for (i in 1:nchains) {
					cols <- c(paste0('beta[', 1:n_sdm_terms, ']'), 'alpha_0', 'alpha_poaceae', 'alpha_area')
					mcmc[[i]] <- mcmc[[i]][ , cols]
				}

				# graphing trace plots for all betas
				pars <- paste0('beta[', 1:n_sdm_terms, ']')
				file <- paste0(out_dir, '/trace_beta_', out_tag, '.png')
				graph <- mcmc_trace(mcmc, pars = pars) +
					ggtitle(paste0(toupper(predictor_names), collapse = ' × '))
				ggsave(graph, file = file, width = 10, height = 8, dpi = 450, bg = 'white')

				# graphing trace plots for all "extra" betas
				pars <- c('alpha_0', 'alpha_poaceae', 'alpha_area')
				file <- paste0(out_dir, '/trace_alpha_', out_tag, '.png')
				graph <- mcmc_trace(mcmc, pars = pars) +
					ggtitle(paste0(toupper(predictor_names), collapse = ' × '))
				ggsave(graph, file = file, width = 10, height = 4, dpi = 450, bg = 'white')

				# graphing density plots for all betas
				pars <- paste0('beta[', 1:n_sdm_terms, ']')
				file <- paste0(out_dir, '/density_beta_', out_tag, '.png')
				graph <- mcmc_dens_overlay(mcmc, pars = pars) +
					ggtitle(paste0(toupper(predictor_names), collapse = ' × '))
				ggsave(graph, file = file, width = 12, height = 8, dpi = 450, bg = 'white')

				# graphing trace plots for all "extra" betas
				pars <- c('alpha_0', 'alpha_poaceae', 'alpha_area')
				file <- paste0(out_dir, '/density_alpha_', out_tag, '.png')
				graph <- mcmc_dens_overlay(mcmc, pars = pars) +
					ggtitle(paste0(toupper(predictor_names), collapse = ' × '))
				ggsave(graph, file = file, width = 10, height = 4, dpi = 450, bg = 'white')

				# Gelman-Rubin statistic
				rhats <- tryCatch(
					gelman.diag(mcmc, autoburnin = FALSE, multivariate = TRUE),
					error = function(cond) FALSE
				)

				say('GELMAN-RUBIN STATISTICS', level = 2)
				if (!is.logical(rhats)) {
					print(rhats)
				} else {
					say('Incalculable!')
				}

			say('###################################', pre = 1)
			say('### map of current distribution ###')
			say('###################################')

				# subset chain summary to just the lambdas
				summary <- chains$summary$all.chains

				# just counties with data
				which_lambda <- grepl(rownames(summary), pattern = 'lambda_sq')
				lambda <- summary[which_lambda, ]

				this_ag_vect_sq$lambda_mean <- lambda[ , 'Mean']

				lambda_sq_quants <- quantile(this_ag_vect_sq$lambda_mean, c(0.25, 0.5, 0.75, 0.90, 0.95))
				quant_labels <- c('[0 - 0.25)', '[0.25 - 0.50)', '[0.50 - 0.75)', '[0.75 - 0.90)', '[0.90 - 0.95)', '[0.95 - 1]')
				
				this_ag_vect_sq$quant_col <- NA
				this_ag_vect_sq$quant_col[this_ag_vect_sq$lambda_mean < lambda_sq_quants[1]] <- quant_labels[1]
				this_ag_vect_sq$quant_col[this_ag_vect_sq$lambda_mean >= lambda_sq_quants[1] & this_ag_vect_sq$lambda_mean < lambda_sq_quants[2]] <- quant_labels[2]
				this_ag_vect_sq$quant_col[this_ag_vect_sq$lambda_mean >= lambda_sq_quants[2] & this_ag_vect_sq$lambda_mean < lambda_sq_quants[3]] <- quant_labels[3]
				this_ag_vect_sq$quant_col[this_ag_vect_sq$lambda_mean >= lambda_sq_quants[3] & this_ag_vect_sq$lambda_mean < lambda_sq_quants[4]] <- quant_labels[4]
				this_ag_vect_sq$quant_col[this_ag_vect_sq$lambda_mean >= lambda_sq_quants[4] & this_ag_vect_sq$lambda_mean < lambda_sq_quants[5]] <- quant_labels[5]
				this_ag_vect_sq$quant_col[this_ag_vect_sq$lambda_mean >= lambda_sq_quants[5]] <- quant_labels[6]

				nam <- vect(paste0(drive, '/Research Data/GADM/Version 4.1/High Res North America Level 1 sans Great Lakes SpatVector WGS84.gpkg'))
				nam <- project(nam, this_ag_vect_sq)

				ag_vect_sq_pres <- this_ag_vect_sq[this_ag_vect_sq$any_ag_quality_1_to_3 > 0]
				extent <- ext(ag_vect_sq_pres)
				extent <- as.vector(extent)
				x_range <- (extent[2] - extent[1])
				y_range <- (extent[4] - extent[3])
				extent[1] <- extent[1] + 0.15 * x_range
				extent[3] <- extent[3] + 0.125 * y_range
				extent[4] <- extent[4] - 0.2 * y_range

				cents_with_ag <- this_ag_vect_sq[this_ag_vect_sq$any_ag_quality_1_to_3 > 0]
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
					layer_spatial(this_ag_vect_sq, aes(fill = quant_col), color = NA) +
					layer_spatial(nam, color = 'gray30', fill = NA, linewidth = 0.1) +
					scale_fill_manual(
						name = 'Quantile\n of λ',
						values = fill_scale
					) +
					layer_spatial(cents_with_ag, pch = 16, color = alpha('gray20', 0.5), size = 0.35) +
					guides(fill = guide_legend(reverse = TRUE)) +
					xlim(extent[1], extent[2]) + ylim(extent[3], extent[4]) +
					ggtitle(
						expression('Present-day distribution of ' * italic('Andropogon gerardi')),
						subtitle = paste0('1961-2020 | N-mixture model | ', toupper(predictor_name_1), ' × ', toupper(predictor_name_2))
					) +
					theme(
						plot.title = element_text(size = 16),
						plot.subtitle = element_text(size = 14)
					)

				ggsave(plot = map, filename = paste0(out_dir, '/map_1960_2020_', out_tag, '.png'), width = 12, height = 10, dpi = 600)

				writeVector(this_ag_vect_sq, paste0(out_dir, '/predictions_1961_2020_bivariate_', out_tag, '.gpkg'), overwrite = TRUE)

			say('###########################')
			say('### parameter estimates ###')
			say('###########################')

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
					ggtitle(paste0(toupper(predictor_name_1), ' × ', toupper(predictor_name_2))) +
					theme(plot.background = element_rect(fill = 'white'))

				ggsave(caterpillars, filename = paste0(out_dir, '/coefficients_', out_tag, '.png'), width = 8, height = 8, dpi = 300)

			say('#######################')
			say('### response curves ###')
			say('#######################')

				# make a plot for how AG SDM and each ancestral population responds to each predictor
				responses <- list()
				for (count_pred in 1:2) {

					pred <- predictor_names[count_pred]

					# unscaled predictor value
					x <- env_array[ , pred, pred]
					
					# create data frame with SDM response
					pars <- paste0('exp_response_curves_sdm_lambda[', 1:n_response_curve_rows, ', ', count_pred, ']')
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
					} else if (pred == 'bio18') {
						nice_title <- 'Precipitation of warmest quarter (BIO18)'
						nice_axis <- 'Precipitation of warmest quarter (mm)'
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

					responses[[count_pred]] <- response

				} # next predictor

				for (i in 1:n_predictors) {
					responses[[i]] <- responses[[i]] + ylim(0, lambda_max)
				}

				responses <- plot_grid(plotlist = responses, nrow = 1)

				ggsave(plot = responses, filename = paste0(out_dir, '/response_curves_', out_tag, '.png'), width = 12.8, height = 7.2, dpi = 600)

				sink()

			} # if predictors are not too correlated

		} # next predictor 2

	} # next predictor 1

say(date())
say('FINIS!', deco = '+', level = 1)

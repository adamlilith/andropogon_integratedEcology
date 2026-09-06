### MODELING ANDROPOGON GERARDI DISTRIBUTION, PHENOTYPE, PHYSIOLOGY, GENOTYPE, and ASSOCIATED MICROBIAL COMMUNITIES
### Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2023-12
###
### This script compiles results of multiple facet and SDM models for Andropogon gerardi.
###
### source('C:/Kaji/R/andropogon_integratedEcology/sdm_pdm/sdm_pdm_90_compile_model_results.r')
###
### CONTENTS ###
### setup ###
### collate occurrence/abundance-only models ###
### compile non-occurrence model results ###
### compile integrated model chains ###

#############
### setup ###
#############

	rm(list = ls())

	setwd('C:/Kaji/Research/Andropogon/Andropogon')
	source('C:/Kaji/R/andropogon_integratedEcology/sdm_pdm/sdm_pdm_00_shared_functions_and_variables.r')

# say('################################################################')
# say('### collate diagnostics for occurrence/abundance-only models ###')
# say('################################################################')

# 	# Collate performance statistics on occurrence-only models (rhats, WAIC, etc.).

# 	model_based_dir <- './outputs_loretta/integrated_sdm_pdm/models_occurrence'
# 	model_dirs <- listFiles(model_based_dir)

# 	model_dirs <- model_dirs[!grepl(model_dirs, pattern = 'laplace')]

# 	# read through all "occurrence meta" files to get a comprehensive list of coefficients
# 	collated <- data.table()
# 	for (i in seq_along(model_dirs)) {

# 		model_dir <- model_dirs[i]
# 		generic <- readRDS(paste0(model_dir, '/!meta_generic.rds'))
# 		occ <- readRDS(paste0(model_dir, '/!meta_occs.rds'))

# 		# this_coeffs <- generic$coeffs

# 		if (grepl(model_dir, pattern = 'poaceae_offset_plus_one_and_zero_correction_covariate')) {
# 			offset_simple <- 'bias ~ offset + 1 & zero correction covariate'
# 		} else if (grepl(model_dir, pattern = 'bias~poaceae_offset_plus_exp_neg_1_WITH_zero_correction_covariate')) {
# 			offset_simple <- 'bias ~ offset + exp(-1) & zero correction covariate'
# 		} else if (grepl(model_dir, pattern = 'bias~poaceae_offset_plus_exp_neg_1_SANS_zero_correction_covariate_large_intercept_prior')) {
# 			offset_simple <- 'bias ~ offset + exp(-1) & zero correction covariate & large intercept prior'
# 		} else if (grepl(model_dir, pattern = 'bias~poaceae_offset_plus_exp_neg_1_SANS_zero_correction_covariate')) {
# 			offset_simple <- 'bias ~ offset + exp(-1) & zero correction covariate'
# 		} else if (grepl(model_dir, pattern = 'bias~poaceae_offset_+_exp_neg_1_SANS_zero_correction_covariate')) {
# 			offset_simple <- 'bias ~ offset + exp(-1)'
# 		} else if (grepl(model_dir, pattern = 'bias~poaceae_offset_plus_one_and_zero_correction_covariate')) {
# 			offset_simple <- 'bias ~ offset + 1 & zero correction covariate'
# 		} else if (grepl(model_dir, pattern = 'bias~poaceae_offset_plus_one_and_correction_covariate')) {
# 			offset_simple <- 'bias ~ offset + 1 & correction covariate (duplicate?)'
# 		} else {
# 			offset_simple <- '(modeled through covariates)'
# 		}

# 		dharma <- occ$dharma_resids
# 		if (!inherits(dharma, 'data.table')) {
		
# 			dharma_sac <-
# 				dharma_uniformity <-
# 				dharma_dispersion <-
# 				dharma_outliers <-
# 				dharma_quantiles_overall <-
# 				dharma_quantiles_upper <-
# 				dharma_quantiles_middle <-
# 				dharma_quantiles_lower <- NA
		
# 		} else {
		
# 			dharma_sac <- dharma$significant[dharma$test == 'spatial autocorrelation']
# 			dharma_uniformity <- dharma$significant[dharma$test == 'uniformity']
# 			dharma_dispersion <- dharma$significant[dharma$test == 'dispersion']
# 			dharma_outliers <- dharma$significant[dharma$test == 'outliers']
# 			dharma_quantiles_overall <- dharma$significant[dharma$test == 'quantiles, overall']
# 			dharma_quantiles_upper <- dharma$significant[dharma$test == 'quantiles, upper']
# 			dharma_quantiles_middle <- dharma$significant[dharma$test == 'quantiles, middle']
# 			dharma_quantiles_lower <- dharma$significant[dharma$test == 'quantiles, lower']
		
# 		}

# 		rhat_max_of_mean <- max(generic$coeffs$rhat)
# 		ess_min <- min(generic$coeffs$eff_sample_size)

# 		waic <- if (is.numeric(generic$waic)) { generic$waic } else { generic$waic$WAIC }

# 		this_collated <- data.table(
# 			facet = generic$facet,
# 			descrip = generic$descrip,
# 			model_dir = model_dir,
# 			homoscedastic = occ$homoscedastic,
# 			zero_inflated = occ$zero_inflated,
# 			formula_occs = paste(generic$formulae$formula_occs, collapse = ' '),
# 			formula_bias = paste(generic$formulae$formula_bias, collapse = ' '),
# 			offset_simple = offset_simple,
# 			waic = waic,
# 			# lppd = generic$waic$lppd,
# 			# pwaic = generic$waic$pWAIC,
# 			rhat_max_of_mean = rhat_max_of_mean,
# 			ess_min = ess_min,
# 			dharma_sac = dharma_sac,
# 			dharma_uniformity = dharma_uniformity,
# 			dharma_dispersion = dharma_dispersion,
# 			dharma_outliers = dharma_outliers,
# 			dharma_quantiles_overall = dharma_quantiles_overall,
# 			dharma_quantiles_upper = dharma_quantiles_upper,
# 			dharma_quantiles_middle = dharma_quantiles_middle,
# 			dharma_quantiles_lower = dharma_quantiles_lower
# 		)
		
# 		# this_collated <- this_collated[rep(1, nrow(this_coeffs))]
# 		# this_collated <- cbind(this_collated, this_coeffs)

# 		valid <-
# 			rhat_max_of_mean <= 1.1 &
# 			ess_min >= 1000
# 		this_collated$valid <- valid

# 		collated <- rbind(collated, this_collated)
	
# 	}

# 	collated <- collated[order(waic, decreasing = FALSE)]

# 	fwrite(collated, './outputs_loretta/integrated_sdm_pdm/summary_of_ALL_occurrence_models.csv')

# say('############################################')
# say('### compile non-occurrence model results ###')
# say('############################################')

# 	facets <- c(
# 		'biomass',
# 		'blade_width',
# 		'canopy_diameter',
# 		'cn_ratio',
# 		'height',
# 		'internal_co2',
# 		'leaf_thickness',
# 		'n_concentration',
# 		'photosynthetic_rate',
# 		'spad',
# 		'stomatal_conductance',
# 		'transpiration_rate'
# 	)

# 	# facets <- c(
# 		# 'canopy_diameter'
# 		# 'blade_width',
# 		# 'cn_ratio'
# 	# )

# 	for (facet in facets) {
	
# 		say(facet, level = 2)

# 		results <- data.table()
# 		folders <- listFiles(paste0('./outputs_loretta/integrated_sdm_pdm/models_', facet), pattern = paste0('\\[', facet, '~'))
# 		for (folder in folders) {

# 			say(folder)

# 			### generic meta
# 			################

# 			meta <- readRDS(paste0(folder, '/!meta_generic.rds'))

# 			descrip <- meta$descrip

# 			# formulae
# 			form_facet <- if (facet == 'biomass') {
# 					meta$formulae[paste0('formula_', facet)][[1]]
# 			} else {
# 					meta$formulae['formula_facet'][[1]]
# 			}

# 			form_psi <- if (any(names(meta$formula) == 'formula_psi')) {
# 					meta$formulae['formula_psi'][[1]]
# 			} else {
# 					NA_character_
# 			}

# 			terms <- attr(terms(form_facet), 'term.labels')
# 			n_terms <- length(terms)

# 			# do CIs for higher-order terms not include 0?
# 			if (any(grepl(terms, pattern = '\\^') | grepl(terms, pattern = '\\:'))) {

# 				# do CIs of higher order terms contain 0?
# 				higher_order_terms <- TRUE
# 				higher_order_indices <- which(grepl(terms, pattern = '\\^') | grepl(terms, pattern = '\\:'))
# 				if (facet == 'biomass') {
# 					coeff <- paste0('beta_biomass[', higher_order_indices, ']')
# 				} else {
# 					coeff <- paste0('beta_facet[', higher_order_indices, ']')
# 				}
# 				higher_order_sig <- all(meta$coeffs$significant[meta$coeffs$coeff %in% coeff] == '*')
# 				lower_order_sig <- NA
				
# 			} else {
			
# 				higher_order_terms <- FALSE
# 				higher_order_sig <- NA

# 				# do CIs of lower-order terms contain 0?
# 				coeff <- if (facet == 'biomass') {
# 					paste0('beta_biomass[', 1 + (1:n_terms), ']')
# 				} else {
# 					paste0('beta_facet[', 1 + (1:n_terms), ']')
# 				}
# 				lower_order_sig <- all(meta$coeffs$significant[meta$coeffs$coeff %in% coeff] == '*')

# 			}
				
# 			form_facet <- as.character(form_facet)[2]
# 			if (!identical(form_psi, NA_character_)) {
# 				form_psi <- as.character(form_psi)[2]
# 			} else {
# 				form_psi <- NA_character_
# 			}

# 			### specific meta
# 			#################
# 			meta_trait <- readRDS(paste0(folder, '/!meta_', facet, '.rds'))

# 			# # # if (is.null(meta_trait$descrip)) {
				
# 			# # # 	descrip <- descrip <- paste0('biomass ~ ', meta_trait$resp_distrib, '(env)')
# 			# # # 	meta_trait$descrip <- descrip
# 			# # # 	saveRDS(meta_trait, paste0(folder, '/!meta_', facet, '.rds'))

# 			# # # }

# 			results <- rbind(
# 				results,
# 				data.table(
# 					facet = facet,
# 					descrip = meta_trait$descrip,
# 					resp_distrib = meta_trait$resp_distrib,
# 					transform = meta_trait$transform,
# 					formula_facet = form_facet,
# 					formula_psi = form_psi,
# 					n_terms = n_terms,
# 					lower_order_sig = lower_order_sig,
# 					higher_order_terms = higher_order_terms,
# 					higher_order_sig = higher_order_sig,
# 					rhat_max_of_mean = max(meta$coeffs$rhat),
# 					ess_min = min(meta$coeffs$eff_sample_size),
# 					waic = meta$waic$WAIC,
# 					pwaic = meta$waic$pWAIC,
# 					lppd = meta$waic$lppd,
# 					correl_lower = meta_trait$crossvalidation[k == 'means', 'correl_lower'][[1]],
# 					correl_mean = meta_trait$crossvalidation[k == 'means', 'correl_mean'][[1]],
# 					correl_upper = meta_trait$crossvalidation[k == 'means', 'correl_upper'][[1]],
# 					mape_lower = meta_trait$crossvalidation[k == 'means', 'mape_lower'][[1]],
# 					mape_mean = meta_trait$crossvalidation[k == 'means', 'mape_mean'][[1]],
# 					mape_upper = meta_trait$crossvalidation[k == 'means', 'mape_upper'][[1]],
# 					obs_quants_prop_in_inner_90 = meta_trait$crossvalidation[k == 'means', 'obs_quants_prop_in_inner_90'][[1]],
# 					dharma_sac = meta_trait$dharma_resids$significant[meta_trait$dharma_resids$test == 'spatial autocorrelation'],
# 					dharma_uniformity = meta_trait$dharma_resids$significant[meta_trait$dharma_resids$test == 'uniformity'],
# 					dharma_dispersion = meta_trait$dharma_resids$significant[meta_trait$dharma_resids$test == 'dispersion'],
# 					dharma_outliers = meta_trait$dharma_resids$significant[meta_trait$dharma_resids$test == 'outliers'],
# 					dharma_quantiles_overall = meta_trait$dharma_resids$significant[meta_trait$dharma_resids$test == 'quantiles, overall'],
# 					dharma_quantiles_upper = meta_trait$dharma_resids$significant[meta_trait$dharma_resids$test == 'quantiles, upper'],
# 					dharma_quantiles_middle = meta_trait$dharma_resids$significant[meta_trait$dharma_resids$test == 'quantiles, middle'],
# 					dharma_quantiles_lower = meta_trait$dharma_resids$significant[meta_trait$dharma_resids$test == 'quantiles, lower']
# 				)
# 			)

# 			} # next folder for this facet

# 		results[ , valid := TRUE]
# 		results$valid[results$rhat_max_of_mean > 1.1] <- FALSE
# 		results$valid[results$ess_min < 1000] <- FALSE
# 		results$valid[results$higher_order_terms & !results$higher_order_sig] <- FALSE
# 		results$valid[!results$lower_order_sig] <- FALSE
# 		# results$valid[results$dharma_sac == '*' | results$dharma_uniformity == '*' | results$dharma_outliers == '*' | results$dharma_quantiles_overall == '*'] <- FALSE
# 		results$valid[results$dharma_sac == '*' | results$dharma_uniformity == '*' | results$dharma_outliers == '*'] <- FALSE

# 		results <- results[order(!valid)]
# 		results <- results[order(correl_mean, decreasing = TRUE)] # bio 12 * sand BUT weird FUT and 1930s

# 		results <- results[order(!valid, waic, decreasing = FALSE)] # bio 12 * sand BUT weird FUT and 1930s
# 		results <- results[order(!valid, mape_mean, decreasing = FALSE)] # bio 12 * sand BUT weird FUT and 1930s
# 		fwrite(results, paste0('./outputs_loretta/integrated_sdm_pdm/summary_of_', facet, '_models.csv'))

# 	}

# say('########################################')
# say('### report top models for each facet ###')
# say('########################################')

# 	facets <- c(
# 		'biomass',
# 		'height',
# 		'blade_width',
# 		'cn_ratio',
# 		'leaf_thickness',
# 		'spad',
# 		'canopy_diameter',
# 		'photosynthetic_rate',
# 		'stomatal_conductance',
# 		'internal_co2',
# 		'transpiration_rate',
# 		'n_concentration'
# 	)

# 	top_models <- data.table()
# 	for (facet in facets) {
	
# 		results <- fread(paste0('./outputs_loretta/integrated_sdm_pdm/summary_of_', facet, '_models.csv'))
# 		index_valid <- which(results$valid)
# 		results <- results[index_valid]
		
# 		results <- results[order(correl_mean, decreasing = TRUE)]
# 		top_models <- rbind(top_models, results[1:2])
	
# 		results <- results[order(mape_mean, decreasing = FALSE)]
# 		top_models <- rbind(top_models, results[1:2])
	
# 		results <- results[order(waic, decreasing = FALSE)]
# 		top_models <- rbind(top_models, results[1:2])
	
# 	}

# 	fwrite(top_models, paste0('./outputs_loretta/integrated_sdm_pdm/summary_of_ALL_top_models.csv'))

say('#######################################')
say('### compile integrated model chains ###')
say('#######################################')

	# This chunk collates the separate chains run for an integrated model. It assumes that different R processes have been used to run each chain, and that samples are stored as matrices in "chunks" (e.g., of 1000 iterations each).  Since chains are run in parallel, until they are all complete, some will have more sets of completed chunks than others. The custom predictXYZ() functions assume that each chain has an equal number of iterations, so the output may be a set of chains with different numbers of iterations.

	# # # name of model folders sans "_chain_X" suffixes
	# model_dir <- './outputs_loretta/integrated_sdm_pdm/models_integrated/non_centered_MVN([occs~hurdlepoisson_offset]_[biomass~hurdleln]_[bla_can_hei])'
	# chains_to_do <- 1:4

	# # # name of model folders sans "_chain_X" suffixes
	# model_dir <- './outputs_loretta/integrated_sdm_pdm/models_integrated/non_centered_MVN([occs~hurdlepoisson_offset]_[biomass~hurdleln]_[bla_can_hei])_eta=2'
	# chains_to_do <- 1:4
	# drops <- c('U_star[1, 1]', 'U_star[2, 1]', 'U_star[3, 1]', 'U_star[4, 1]', 'U_star[5, 1]', 'U_star[3, 2]', 'U_star[4, 2]', 'U_star[5, 2]', 'U_star[4, 3]', 'U_star[5, 3]', 'U_star[5, 4]')
	# facets <- c('abundance', 'biomass', 'blade width', 'canopy diameter', 'height')

	model_dir <- './outputs_loretta/integrated_sdm_pdm/models_integrated/non_centered_MVN([occs~hurdlepoisson_offset]_[biomass~hurdleln]_[cnr_pho_spa_sto_tra])_eta=2'
	chains_to_do <- 1:5
	drops <- c('U_star[1, 1]', 'U_star[2, 1]', 'U_star[3, 1]', 'U_star[3, 2]', 'U_star[4, 1]', 'U_star[4, 2]', 'U_star[4, 3]', 'U_star[5, 1]', 'U_star[5, 2]', 'U_star[5, 3]', 'U_star[5, 4]', 'U_star[6, 1]', 'U_star[6, 2]', 'U_star[6, 3]', 'U_star[6, 4]', 'U_star[6, 5]', 'U_star[7, 1]', 'U_star[7, 2]', 'U_star[7, 3]', 'U_star[7, 4]', 'U_star[7, 5]', 'U_star[7, 6]')
	facets <- c('abundance', 'biomass', 'C:N', 'photosynthesis', 'SPAD', 'stomatal conductance', 'transpiration rate')

	say('Collating chains ', paste(chains_to_do, collapse = ', '), ' from \n   ', paste(model_dir, collapse = '\n   '), '...')

	# see how many set files the chain that has made the least progress has made
	smallest_number_set_files <- Inf
	for (chain_to_do in chains_to_do) {

		this_folder <- paste0(model_dir, '_chain_', chain_to_do)
		set_files <- listFiles(this_folder, pattern = 'chains_set_')
		smallest_number_set_files <- min(smallest_number_set_files, length(set_files))

	}

	# load from the second half of the set files (discarding the first half as burnin)
	if (smallest_number_set_files %% 2 == 1) {
		start_set <- (smallest_number_set_files - 1) / 2
	} else {
		start_set <- smallest_number_set_files / 2
	}

	say('Chain with fewest iterations has ', smallest_number_set_files, ' set files. Discarding first ', start_set, ' as burnin.')

	chains <- list(samples = list())
	for (chain_to_do in chains_to_do) {

		say('Collating chain ', chain_to_do, '...')

		this_folder <- paste0(model_dir, '_chain_', chain_to_do)
		set_files <- listFiles(this_folder, pattern = 'chains_set_')

		if (exists('chain')) rm(chain)
		for (set in start_set:smallest_number_set_files) {

			chain_chunk <- readRDS(set_files[set])
			if (exists('chain')) {
				chain <- rbind(chain, chain_chunk)
			} else {
				chain <- chain_chunk
			}

		}

		# thin so we get 1000 samples
		keeps <- round(seq(1, nrow(chain), length.out = 1000))
		chain <- chain[keeps, ]

		chains$samples[[chain_to_do]] <- as.mcmc(chain, start = 1, end = 1000, thin = 1)

	}

	chains$samples <- as.mcmc.list(chains$samples)
	chains <- mc_resummarize(chains)

	dirCreate(model_dir)
	saveRDS(chains, paste0(model_dir, '/chains.rds'))

	### copy other files
	file.copy(paste0(model_dir, '_chain_1/!facet_codes.csv'), paste0(model_dir, '/facet_codes.csv'), overwrite = TRUE)
	file.copy(paste0(model_dir, '_chain_1/.constants.rds'), paste0(model_dir, '/.constants.rds'), overwrite = TRUE)
	file.copy(paste0(model_dir, '_chain_1/.data.rds'), paste0(model_dir, '/.data.rds'), overwrite = TRUE)
	file.copy(paste0(model_dir, '_chain_1/.inits.rds'), paste0(model_dir, '/.inits.rds'), overwrite = TRUE)
	file.copy(paste0(model_dir, '_chain_1/model_code.txt'), paste0(model_dir, '/model_code.txt'), overwrite = TRUE)
	file.copy(paste0(model_dir, '_chain_1/monitors.txt'), paste0(model_dir, '/monitors.txt'), overwrite = TRUE)
	file.copy(paste0(model_dir, '_chain_1/formulae.rds'), paste0(model_dir, '/formulae.rds'), overwrite = TRUE)
	file.copy(paste0(model_dir, '_chain_1/runtime_log.txt'), paste0(model_dir, '/runtime_log.txt'), overwrite = TRUE)

	### calculate facet-facet correlations from U_star and add to chains object
	###########################################################################

	# replace U_star with elements from correlation matrix
	indices <- which(grepl(colnames(chains$samples[[1]]), pattern = 'U_star'))
	chains_corr <- list(samples = vector('list', length(chains_to_do)))
	n_facets <- length(facets)

	nms <- expand.grid(facets, facets)
	nms <- nms[nms[[1]] != nms[[2]], , drop = FALSE] # remove same facet--same-facet names
	nms <- apply(nms, 1, function(x) paste0('cor_', x[1], '_vs_', x[2]))
	nms <- gsub(nms, pattern = ' ', replacement = '_')

	# collate correlation matrices into mcmc-like object for further analysis
	for (chain in chains_to_do) {
		for (i in 1:1000) {

			# correlation
			U_star <- chains$samples[[chain]][i, indices]
			U_star <- matrix(U_star, nrow = n_facets, ncol = n_facets)
			this_corr <- t(U_star) %*% U_star
			diag(this_corr) <- NA

			this_corr <- c(this_corr)
			this_corr <- this_corr[!is.na(this_corr)]
			this_corr <- matrix(this_corr, nrow = 1)
			colnames(this_corr) <- nms
			if (!is.matrix(chains_corr$samples[[chain]])) {
				chains_corr$samples[[chain]] <- this_corr
			} else {
				chains_corr$samples[[chain]] <- rbind(chains_corr$samples[[chain]], this_corr)
			}

		}
	}

	# remove cases of "facet2-vs-facet1" duplicating "facet1-vs-facet2" and "facet1-vs-facet1"
	for (i in chains_to_do) chains_corr$samples[[i]] <- as.mcmc(chains_corr$samples[[i]], start = 1, end = 1000, thin = 1)
	chains_corr$samples <- as.mcmc.list(chains_corr$samples)

	nms <- expand.grid(facets, facets)
	nms <- nms[nms[[1]] == nms[[2]] | duplicated(apply(nms, 1, function(x) paste(sort(x), collapse = '_'))), , drop = FALSE]
	nms <- apply(nms, 1, function(x) paste0('cor_', x[1], '_vs_', x[2]))
	nms <- gsub(nms, pattern = ' ', replacement = '_')
	for (nm in nms) chains_corr <- mc_subset(chains_corr, param = nm, keep = FALSE, resummarize = FALSE)

	# drop U_star because convergence should be calculated on the correlations
	chains <- mc_subset(chains, param = 'U_star', j = TRUE, k = TRUE, keep = FALSE)

	# add correlations
	for (chain in chains_to_do) {
		chains$samples[[chain]] <- cbind(chains$samples[[chain]], chains_corr$samples[[chain]])
	}
	chains <- mc_resummarize(chains)

	say('trace plots', level = 2)
	vars <- rownames(chains$summary$all.chains)
	x <- 1:1000
	for (var in vars) {

		say(var)

		png(paste0(model_dir, '/trace_', var, '.png'), width = 800, height = 1200, res = 120)
		par(mfrow = c(max(chains_to_do), 1), cex.main = 1.6, cex.axis = 1.4, cex.lab = 1.4)

		ylim <- c(Inf, -Inf)
		for (i in 1:length(chains$samples)) {
			y <- as.numeric(chains$samples[[i]][ , var])
			ylim[1] <- min(ylim[1], min(y))
			ylim[2] <- max(ylim[2], max(y))
		}


		for (i in 1:length(chains$samples)) {
			y <- as.numeric(chains$samples[[i]][ , var])
			if (length(y) < n) y <- c(y, rep(NA_real_, n - length(y)))
			plot(x, y, type = 'l', main = paste0('Chain ', i, ' - ', var), xlab = 'Iteration', ylab = 'Value', ylim = ylim)
		}

		dev.off()

	}

	# Gelman-Rubin R-hat and effective sample size
	rhats <- gelman.diag(chains$samples, autoburnin = FALSE, multivariate = FALSE)
	ess <- effectiveSize(chains$samples)

	sink(paste0(model_dir, '/!meta.txt'), split = TRUE)
		
		say('RHATS', post = 2)
		print(rhats)

		say('EFFECTIVE SAMPLE SIZE', pre = 1, post = 2)
		print(as.data.frame(ess))

	sink()

# say('#######################################################################')
# say('### post-modeling diagnostics and predictions for integrated models ###')
# say('#######################################################################')

# 	# Create trace/density plots and calculate relevant diagnostics for integrated models. For other models, these are done within each model script, but integrated models are run in parallel so need collated together first (see code chunk above).

# 	# # name of model folder
# 	# model_dir <- './outputs_loretta/integrated_sdm_pdm/models_integrated/non_centered_MVN([occs~hurdlepoisson_offset]_[biomass~hurdleln]_[bla_can_hei])'

# 	model_dir <- './outputs_loretta/integrated_sdm_pdm/models_integrated/non_centered_MVN([occs~hurdlepoisson_offset]_[biomass~hurdleln]_[cnr_pho_spa_sto_tra])'

# 	# descrip <- 'non_centered_MVN([occs~hurdlepoisson_offset]_[biomass~hurdleln])'
# 	# for (f in seq_along(nonbiomass_facets)) {
# 	# 	descrip <- paste0(descrip, ', ', names(nonbiomass_facets)[f], ' ~ ', tolower(nonbiomass_facets[[f]]$resp_distrib), '(MVN)')
# 	# }
# 	# facets_nice <- paste0('occurrence + biomass + ', paste(names(nonbiomass_facets), collapse = ' + '))

# 	chains <- readRDS(paste0(model_dir, '/chains.rds'))
# 	formulae <- readRDS(paste0(model_dir, '/formulae.rds'))
# 	nonbiomass_facets <- formulae$nonbiomass_facets

# 	workflow_postmodeling_fully_integrated(
# 		formula_occs = formulae$formula_occs,
# 		formula_bias = formulae$formula_bias,
# 		formula_psi = formulae$formula_psi,
# 		formula_biomass = formulae$formula_biomass,
# 		nonbiomass_facets = nonbiomass_facets,
# 		out_dir = model_dir
# 	)


say('DONE', level = 1, deco = '@')

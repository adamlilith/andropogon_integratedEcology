### MODELING ANDROPOGON GERARDI DISTRIBUTION, PHENOTYPE, PHYSIOLOGY, GENOTYPE, and ASSOCIATED MICROBIAL COMMUNITIES
### Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2023-12
###
### This script compiles results of multiple facet and SDM models for Andropogon gerardi.
###
### source('C:/Kaji/R/andropogon_integratedEcology/sdm_pdm/sdm_pdm_90_compile_model_results.r')
###
### CONTENTS ###
### setup ###
### collate occurrence-only models ###
### compile non-occurrence model results ###

#############
### setup ###
#############

	rm(list = ls())

	setwd('C:/Kaji/Research/Andropogon/Andropogon')
	source('C:/Kaji/R/andropogon_integratedEcology/sdm_pdm/sdm_pdm_00_shared_functions_and_variables.r')

say('######################################')
say('### collate occurrence-only models ###')
say('######################################')

	# Collate performance statistics on occurrence-only models (rhats, WAIC, etc.).

	model_based_dir <- './outputs_loretta/integrated_sdm_pdm/models_occurrence'
	model_dirs <- listFiles(model_based_dir)

	model_dirs <- model_dirs[!grepl(model_dirs, pattern = 'laplace')]

	# read through all "occurrence meta" files to get a comprehensive list of coefficients
	collated <- data.table()
	for (i in seq_along(model_dirs)) {

		model_dir <- model_dirs[i]
		generic <- readRDS(paste0(model_dir, '/!meta_generic.rds'))
		occ <- readRDS(paste0(model_dir, '/!meta_occs.rds'))

		this_coeffs <- generic$coeffs

		if (grepl(model_dir, pattern = 'poaceae_offset_plus_one_and_zero_correction_covariate')) {
			offset_simple <- 'bias ~ offset + 1 & zero correction covariate'
		} else if (grepl(model_dir, pattern = 'bias~poaceae_offset_plus_exp_neg_1_WITH_zero_correction_covariate')) {
			offset_simple <- 'bias ~ offset + exp(-1) & zero correction covariate'
		} else if (grepl(model_dir, pattern = 'bias~poaceae_offset_plus_exp_neg_1_SANS_zero_correction_covariate_large_intercept_prior')) {
			offset_simple <- 'bias ~ offset + exp(-1) & zero correction covariate & large intercept prior'
		} else if (grepl(model_dir, pattern = 'bias~poaceae_offset_plus_exp_neg_1_SANS_zero_correction_covariate')) {
			offset_simple <- 'bias ~ offset + exp(-1) & zero correction covariate'
		} else if (grepl(model_dir, pattern = 'bias~poaceae_offset_+_exp_neg_1_SANS_zero_correction_covariate')) {
			offset_simple <- 'bias ~ offset + exp(-1)'
		} else if (grepl(model_dir, pattern = 'bias~poaceae_offset_plus_one_and_zero_correction_covariate')) {
			offset_simple <- 'bias ~ offset + 1 & zero correction covariate'
		} else if (grepl(model_dir, pattern = 'bias~poaceae_offset_plus_one_and_correction_covariate')) {
			offset_simple <- 'bias ~ offset + 1 & correction covariate (duplicate?)'
		} else {
			offset_simple <- '(modeled through covariates)'
		}

		dharma <- occ$dharma_resids
		if (!inherits(dharma, 'data.table')) {
		
			dharma_sac <-
				dharma_uniformity <-
				dharma_dispersion <-
				dharma_outliers <-
				dharma_quantiles_overall <-
				dharma_quantiles_upper <-
				dharma_quantiles_middle <-
				dharma_quantiles_lower <- NA
		
		} else {
		
			dharma_sac <- dharma$significant[dharma$test == 'spatial autocorrelation']
			dharma_uniformity <- dharma$significant[dharma$test == 'uniformity']
			dharma_dispersion <- dharma$significant[dharma$test == 'dispersion']
			dharma_outliers <- dharma$significant[dharma$test == 'outliers']
			dharma_quantiles_overall <- dharma$significant[dharma$test == 'quantiles, overall']
			dharma_quantiles_upper <- dharma$significant[dharma$test == 'quantiles, upper']
			dharma_quantiles_middle <- dharma$significant[dharma$test == 'quantiles, middle']
			dharma_quantiles_lower <- dharma$significant[dharma$test == 'quantiles, lower']
		
		}

		rhat_max_of_mean <- max(generic$coeffs$rhat)
		ess_min <- min(generic$coeffs$eff_sample_size)

		this_collated <- data.table(
			facet = generic$facet,
			descrip = generic$descrip,
			model_dir = model_dir,
			homoscedastic = occ$homoscedastic,
			zero_inflated = occ$zero_inflated,
			formula_occs = paste(generic$formulae$formula_occs, collapse = ' '),
			formula_bias = paste(generic$formulae$formula_bias, collapse = ' '),
			offset_simple = offset_simple,
			waic = generic$waic$WAIC,
			lppd = generic$waic$lppd,
			pwaic = generic$waic$pWAIC,
			rhat_max_of_mean = rhat_max_of_mean,
			ess_min = ess_min,
			dharma_sac = dharma_sac,
			dharma_uniformity = dharma_uniformity,
			dharma_dispersion = dharma_dispersion,
			dharma_outliers = dharma_outliers,
			dharma_quantiles_overall = dharma_quantiles_overall,
			dharma_quantiles_upper = dharma_quantiles_upper,
			dharma_quantiles_middle = dharma_quantiles_middle,
			dharma_quantiles_lower = dharma_quantiles_lower
		)
		
		this_collated <- this_collated[rep(1, nrow(this_coeffs))]
		this_collated <- cbind(this_collated, this_coeffs)

		valid <-
			rhat_max_of_mean <= 1.1 &
			ess_min >= 1000 &
			dharma_sac != '*' &
			dharma_uniformity != '*' &
			dharma_dispersion != '*' &
			dharma_outliers != '*' &
			dharma_quantiles_overall != '*' &
			dharma_quantiles_upper != '*' &
			dharma_quantiles_middle != '*' &
			dharma_quantiles_lower != '*'
		this_collated$valid <- valid

		collated <- rbind(collated, this_collated)
	
	}

	collated <- collated[order(waic, decreasing = FALSE)]

	fwrite(collated, './outputs_loretta/integrated_sdm_pdm/summary_of_ALL_occurrence_models.csv')

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

say('DONE', level = 1, deco = '@')

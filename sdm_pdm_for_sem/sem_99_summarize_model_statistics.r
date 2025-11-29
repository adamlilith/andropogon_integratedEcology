### MODELING ANDROPOGON GERARDI DISTRIBUTION, PHENOTYPE, PHYSIOLOGY, GENOTYPE, and ASSOCIATED MICROBIAL COMMUNITIES
### Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2023-12
###
### This script collates model statistics across model types.
###
### source('C:/Kaji/R/andropogon_integratedEcology/sdm_pdm_for_sem/sem_99_summarize_model_statistics.r')
### 
#############
### setup ###
#############

	rm(list = ls())

	drive <- 'C:/Kaji/'

	setwd(paste0(drive, '/Research/Andropogon/Andropogon'))
	source(paste0(drive, '/R/andropogon_integratedEcology/sdm_pdm_for_sem/sem_00_shared_functions_and_variables.r'))

###########################
### user-defined values ###
###########################


say('######################################')
say('### collate occurrence-only models ###')
say('#######################################')

	# Collate performance statistics on occurrence-only models (rhats, WAIC, etc.).

	model_based_dir <- './outputs_loretta/sdm_pdm_for_sem/models_occurrence'
	model_dirs <- listFiles(model_based_dir)

	model_dirs <- model_dirs[!grepl(model_dirs, pattern = 'laplace')]

	# read through all "occurrence meta" files to get a comprehensive list of coefficients
	coeffs <- data.table()
	for (i in seq_along(model_dirs)) {

		model_dir <- model_dirs[i]
		generic <- readRDS(paste0(model_dir, '/!meta_generic.rds'))
		occ <- readRDS(paste0(model_dir, '/!meta_occs.rds'))

		model_coeffs <- generic$coeffs$coeff
		
		this_coeffs <- data.table(
			coeff = model_coeffs,
			term = NA_character_,
			submodel = NA_character_
		)

		# match coefficients to terms: occurrence
		form <- generic$formulae$formula_occs
		terms <- terms(form)
		terms <- attr(terms, 'term.labels')
		terms <- c('intercept', terms)

		this_coeffs$term[grepl(this_coeffs$coeff, pattern = 'beta_occs_mu')] <- terms
		this_coeffs$submodel[grepl(this_coeffs$coeff, pattern = 'beta_occs_mu')] <- 'occurrence'

		# match coefficients to terms: occurrence bias
		form <- generic$formulae$formula_occs_bias
		terms <- terms(form)
		terms <- attr(terms, 'term.labels')
		terms <- c('intercept', terms)

		this_coeffs$term[grepl(this_coeffs$coeff, pattern = 'alpha_occs')] <- terms
		this_coeffs$submodel[grepl(this_coeffs$coeff, pattern = 'alpha_occs')] <- 'bias'

		# sigma
		if (occ$homoscedastic) {
			this_coeffs$term[this_coeffs$coeff == 'lambda_sigma'] <- 'occurrence s.d.'
			this_coeffs$submodel[this_coeffs$coeff == 'lambda_sigma'] <- 'occurrence s.d.'
		} else {
			form <- generic$formulae$formula_occs
			terms <- terms(form)
			terms <- attr(terms, 'term.labels')
			terms <- c('intercept', terms)

			this_coeffs$term[grepl(this_coeffs$coeff, pattern = 'beta_occs_sigma')] <- terms
			this_coeffs$submodel[grepl(this_coeffs$coeff, pattern = 'beta_occs_sigma')] <- 'occurrence s.d.'
		}

		# zero-inflated
		if (occ$zero_inflated) {
			
			form <- generic$formulae$formula_occs
			terms <- terms(form)
			terms <- attr(terms, 'term.labels')
			terms <- c('intercept', terms)

			this_coeffs$term[grepl(this_coeffs$coeff, pattern = 'beta_pzero')] <- terms
			this_coeffs$submodel[grepl(this_coeffs$coeff, pattern = 'beta_pzero')] <- 'occurrence pzero'
		}

		coeffs <- rbind(coeffs, this_coeffs)

	}

	coeffs <- coeffs[!duplicated(coeffs[ , c('term', 'submodel')])]
	coeffs <- coeffs[order(submodel)]
	coeffs$coeff <- NULL
	coeffs[ , c('lower', 'mean', 'upper') := NA_real_]
	
	### collate statistics for each model
	collated <- data.table()
	for (i in seq_along(model_dirs)) {
	
		model_dir <- model_dirs[i]
		say(model_dir)

		bias_simple <- grepl(model_dir, pattern = 'p~dunif')

		generic <- readRDS(paste0(model_dir, '/!meta_generic.rds'))
		occ <- readRDS(paste0(model_dir, '/!meta_occs.rds'))

		ll <- generic$log_lik
		ll_lower <- ll[ll == min(ll)]
		ll_upper <- ll[ll == max(ll)]
		ll_mean <- ll[ll != max(ll) & ll != min(ll)]

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

		collated <- rbind(
			collated,
			data.table(
				facet = generic$facet,
				descrip = generic$descrip,
				model_dir = model_dir,
				homoscedastic = occ$homoscedastic,
				zero_inflated = occ$zero_inflated,
				formula_occs = paste(generic$formulae$formula_occs, collapse = ' '),
				formula_bias = paste(generic$formulae$formula_occs_bias, collapse = ' '),
				bias_simple = bias_simple,
				waic = generic$waic$WAIC,
				lppd = generic$waic$lppd,
				pwaic = generic$waic$pWAIC,
				log_lik_lower = ll_lower,
				log_lik_mean = ll_mean,
				log_lik_upper = ll_upper,
				rhat_max_of_mean = max(generic$coeffs$rhat),
				ess_min = min(generic$coeffs$eff_sample_size),
				dharma_sac = dharma_sac,
				dharma_uniformity = dharma_uniformity,
				dharma_dispersion = dharma_dispersion,
				dharma_outliers = dharma_outliers,
				dharma_quantiles_overall = dharma_quantiles_overall,
				dharma_quantiles_upper = dharma_quantiles_upper,
				dharma_quantiles_middle = dharma_quantiles_middle,
				dharma_quantiles_lower = dharma_quantiles_lower
			)
		)
	
	}

	fwrite(collated, './outputs_loretta/sdm_pdm_for_sem/summaries_occurrence_models.csv')


say(date())
say('FINIS!', deco = '+', level = 1)

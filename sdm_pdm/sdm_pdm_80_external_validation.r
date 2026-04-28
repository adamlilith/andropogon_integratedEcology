### MODELING ANDROPOGON GERARDI DISTRIBUTION, PHENOTYPE, PHYSIOLOGY, GENOTYPE, and ASSOCIATED MICROBIAL COMMUNITIES
### Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2026-03
###
### This script conducts external validation of Andropogon gerardi models.
###
### source('C:/Kaji/R/andropogon_integratedEcology/sdm_pdm/sdm_pdm_80_external_validation.r')
###
### CONTENTS ###
### setup ###
### compare models against AG height from Figure 3 in McMillan (1964) ###

#############
### setup ###
#############

	rm(list = ls())

	setwd('C:/Kaji/Research/Andropogon/Andropogon')
	source('C:/Kaji/R/andropogon_integratedEcology/sdm_pdm/sdm_pdm_00_shared_functions_and_variables.r')

say('#########################################################################')
say('### compare models against AG height from Figure 3 in McMillan (1964) ###')
say('#########################################################################')

	library(brms)

	# Compare model predictions at sites from where McMillan collected AG seed to AG height data in Figure 3 of McMillan (1964). We will use the slope of a linear mixed model (obs ~ predicted | site) and the Spearman rank correlation between observed and predicted. Observed values were binned into ordinal classes by McMillan, so the regression assumes a Poisson response.

	mcm <- vect('./data_from_mcmillan/mcmillan_1964_fig_3_ag_prism_1952_1961_soilgrids_top_5_cm.gpkg')
	mcm <- calculate_logged_vars(mcm)
	obs <- mcm[ , c('id', paste0('height', 1:3))]
	obs <- as.data.table(obs)

	model_folders <- listFiles('./outputs_loretta/integrated_sdm_pdm/models_height')

	centers_scales <- fread('./outputs_loretta/integrated_sdm_pdm/centers_and_scales_for_covariates.csv')

	summary <- data.table()
	for (i in seq_along(model_folders)) {

		model_folder <- model_folders[i]
		say(model_folder)

		# posterior samples
		chains <- readRDS(paste0(model_folder, '/chains.rds'))

		# metadata
		meta_facet <- readRDS(paste0(model_folder, '/!meta_height.rds'))
		formula_facet <- meta_facet$formulae$formula_facet
		resp_distrib <- meta_facet$resp_distrib
		transform <- meta_facet$transform

		terms <- terms(formula_facet)
		terms <- attr(terms, 'term.labels')

		covariates <- terms
		covariates <- covariates[!grepl(covariates, pattern = '\\^2')]
		covariates <- covariates[!grepl(covariates, pattern = '\\:')]
		covariates <- covariates[!grepl(covariates, pattern = '\\*')]

		mm <- mcm[ , covariates, drop = FALSE]
		mm <- as.data.table(mm)

		for (covariate in covariates) {

			center <- centers_scales$center_sites[centers_scales$variable == covariate]
			scale <- centers_scales$scale_sites[centers_scales$variable == covariate]

			mm[ , covariate] <- (mm[[covariate]] - center) / scale

		}

		x <- model.matrix(formula_facet, data = mm)

		# predict facet model
		preds <- predict_nonbiomass_single_trait(chains, x = x, resp_distrib = resp_distrib, transform = transform)
		predictions <- colMeans(preds)
		
		# compare facet model predictions to observed height classes
		predicted_scaled <- scale(predictions)[ , 1]
		obs_pred <- cbind(obs, data.table(predicted_scaled = predicted_scaled))

		melted <- melt(
			obs_pred,
			id.vars = setdiff(names(obs_pred), c('height1', 'height2', 'height3')), 
			variable.name = 'plant', value.name = 'obs_height'
		)

		melted <- melted[order(id)]
		melted$id <- factor(melted$id)
		melted <- melted[complete.cases(melted)]

		fit <- brm(
			obs_height ~ predicted_scaled + (1 | id),
			data = melted,
			family = cumulative('logit'),
			cores = 2,
			silent = 0
		)

		slope <- fixef(fit)['predicted_scaled', ]

		# corelation between facet model predictions and observed height classes across posterior samples
		niters <- nrow(preds)
		correl <- rep(NA_real_, niters)
		for (niter in 1:niters) {
			
			predicted <- preds[niter, ]
			obs_pred <- cbind(obs, data.table(predicted = predicted))

			melted <- melt(
				obs_pred,
				id.vars = setdiff(names(obs_pred), c('height1', 'height2', 'height3')), 
				variable.name = 'plant', value.name = 'obs_height'
			)
			melted <- melted[order(id)]
			melted$id <- factor(melted$id)
			melted <- melted[complete.cases(melted)]

			correl[niter] <- cor(melted$predicted, melted$obs_height, method = 'spearman')

		}

		formula_facet_char <- paste(as.character(formula_facet), collapse = ' ')

		summary <- rbind(
			summary,
			data.table(
				facet = 'height',
				model_folder = basename(model_folder),
				formula_facet = formula_facet_char, 
				resp_distrib = resp_distrib,
				transform = transform,
				slope_mean = slope[['Estimate']],
				slope_95CI_lower = slope[['Q2.5']],
				slope_95CI_upper = slope[['Q97.5']],
				correl_mean = mean(correl),
				correl_95CI_lower = quantile(correl, 0.025),
				correl_95CI_upper = quantile(correl, 0.975)
			)
		)


	} # next model

	fwrite(summary, './outputs_loretta/integrated_sdm_pdm/models_height/validation_vs_mcmillan_1964_height.csv')

say('DONE', level = 1, deco = '@')

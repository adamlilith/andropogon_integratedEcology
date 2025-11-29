### MODELING ANDROPOGON GERARDI DISTRIBUTION, PHENOTYPE, PHYSIOLOGY, GENOTYPE, and ASSOCIATED MICROBIAL COMMUNITIES
### Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2023-12
###
### This script compiles results of multiple trait and SDM models for Andropogon gerardi.
###
### source('C:/Kaji/R/andropogon_integratedEcology/sdm_pdm_for_sem/sem_xx_compile_model_results.r')
###
### CONTENTS ###
### setup ###
### compile model results ###

#############
### setup ###
#############

	rm(list = ls())

	drive <- 'C:/Kaji/'

	setwd(paste0(drive, '/Research/Andropogon/Andropogon'))
	source(paste0(drive, '/R/andropogon_integratedEcology/sdm_pdm_for_sem/sem_00_shared_functions_and_variables.r'))

say('#############################')
say('### compile model results ###')
say('#############################')

	model_results <- data.table()

	folders <- listFiles(
		'./outputs_loretta/sdm_pdm_for_sem',
		pattern = 'model_biomass'
	)

	for (folder in folders) {

	  say(folder)
	  
	  meta <- readRDS(paste0(folder, '/!meta_biomass.rds'))
	  chains <- readRDS(paste0(folder, '/chains.rds'))
	  
		# formulae
		forms <- meta$formulae
		formulae <- character()
		for (i in seq_along(forms)) {

			form_name <- names(forms)[i]
			formulae[i] <- paste0(as.character(forms[[i]]), collapse = '')

		}
	  names(formulae) <- names(forms)

	  # remember
	  model_results <- rbind(
			model_results,
			data.table(
				formula_biomass = formulae['formula_biomass_mu'],
				zero_inflated = formulae['formula_pzero'] != '',
				waic = chains$WAIC$WAIC,
				pwaic = chains$WAIC$pWAIC,
				crossvalidation_correl = mean(meta$crossvalidation$correl_mean),
				crossvalidation_mae = mean(meta$crossvalidation$mae_mean),
				crossvalidation_mape = mean(meta$crossvalidation$mape_mean),
				crossvalidation_rmse = mean(meta$crossvalidation$rmse_mean),
				crossvalidation_obs_quants = mean(meta$crossvalidation$obs_quants_mean)
			)
		)
	  
	} # next model

	model_results <- model_results[order(waic, decreasing = FALSE)]
	fwrite(model_results, './outputs_loretta/sdm_pdm_for_sem/compiled_model_results_biomass.csv')


	  
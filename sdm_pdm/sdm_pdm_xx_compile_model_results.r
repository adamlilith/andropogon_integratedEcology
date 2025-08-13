### MODELING ANDROPOGON GERARDI DISTRIBUTION, PHENOTYPE, PHYSIOLOGY, GENOTYPE, and ASSOCIATED MICROBIAL COMMUNITIES
### Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2023-12
###
### This script compiles results of multiple trait and SDM models for Andropogon gerardi.
###
### source('C:/Kaji/R/andropogon_integratedEcology/sdm_pdm/sdm_pdm_xx_compile_model_results.r')
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
	source(paste0(drive, '/R/andropogon_integratedEcology/sdm_pdm/sdm_pdm_00_shared_functions_and_variables.r'))

say('#############################')
say('### compile model results ###')
say('#############################')

	traits <- c(
		'height',
		'blade_width',
		'leaf_thickness',
		'spad',
		'canopy_diameter',
		'photosynthetic_rate',
		'stomatal_conductance',
		'internal_co2',
		'transpiration_rate',
		'n_concentration'
	)

	model_results <- data.table()

	for (trait in traits) {
	
		folders <- listFiles(paste0('./outputs_loretta/integrated_sdm_pdm/models_', trait))
		for (folder in folders) {

			# formulae
			forms <- readRDS(paste0(folder, '/formula.rds'))
			formulae <- character()
			for (i in seq_along(forms)) {

				form_name <- names(forms)[i]
				formulae[i] <- paste0(as.character(forms[[i]]), collapse = '')

			}

			# model fit
			raw <- readLines(paste0(folder, '/model_fit_', trait, '.txt'))
			raw <- raw[grepl(raw, pattern = 'Correlation')]
			raw <- strsplit(raw, split = ': ')
			correl <- as.numeric(raw[[1]][2])
			

			# WAIC
			raw <- readLines(paste0(folder, '/waic.txt'))
			waic <- as.numeric(substr(raw[6], 5, nchar(raw[6])))
			lppd <- as.numeric(substr(raw[8], 5, nchar(raw[8])))
			pwaic <- as.numeric(substr(raw[10], 5, nchar(raw[10])))

			# CV
			raw <- fread(paste0(folder, '/crossvalidation_', trait, '.csv'))
			cv_score <- raw$cv_value[raw$k == 'summary']
			cv_se <- raw$cv_value_se[raw$k == 'summary']

		}


	}




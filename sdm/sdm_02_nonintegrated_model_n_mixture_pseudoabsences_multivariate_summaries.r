### MODELING ANDROPOGON GERARDI DISTRIBUTION, PHENOTYPE, PHYSIOLOGY, GENOTYPE, and ASSOCIATED MICROBIAL COMMUNITIES
### Erica Newman | Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2023-12
###
### This script evaluates alternative multivariate SDMs for AG distribution.
###
### source('C:/Ecology/R/andropogon_integratedEcology/sdm/sdm_02_nonintegrated_model_n_mixture_pseudoabsences_multivariate_summaries.r')
### source('C:/Subarashi/R/andropogon_integratedEcology/sdm/sdm_02_nonintegrated_model_n_mixture_pseudoabsences_multivariate_summaries.r')
###
### CONTENTS ###
### setup ###
### compare models ###

#############
### setup ###
#############

	rm(list = ls())

	# drive <- 'C:/Ecology/'
	drive <- 'C:/Subarashi/'

	setwd(paste0(drive, '/Research/Andropogon/Andropogon'))

	library(coda) # Bayesian diagnostics
	library(data.table) # data wrangling
	library(omnibus)

	### user-defined values
	#######################

		climate_predictor_names <- c('bio1', 'bio12', 'bio15') # set for Jack's multivariate analyses, chosen based on simple models and collinearity

		psa_quant <- 0.99 # to define pseudoabsences, use this quantile of Poaceae occurrences across counties with 

say('######################')
say('### compare models ###')
say('######################')

	models_dirs <- c(
		'sdm_[nmixture]_[pseudoabsences_0.99]_[bio1^2_bio12^2_bio15^2_bio1xbio12]_[priors_ddnorm]',
		'sdm_[nmixture]_[pseudoabsences_0.99]_[bio1^2_bio12^2_bio15^2_bio1xbio12_bio1xbio15]_[priors_ddnorm]',
		'sdm_[nmixture]_[pseudoabsences_0.99]_[bio1^2_bio12^2_bio15^2_bio1xbio12_bio1xbio15_bio12xbio15]_[priors_ddnorm]',
		'sdm_[nmixture]_[pseudoabsences_0.99]_[bio1^2_bio12^2_bio15^2_bio1xbio12_bio12xbio15]_[priors_ddnorm]',
		'sdm_[nmixture]_[pseudoabsences_0.99]_[bio1^2_bio12^2_bio15^2_bio1xbio12xbio15]_[priors_ddnorm]',
		'sdm_[nmixture]_[pseudoabsences_0.99]_[bio1^2_bio12^2_bio15^2_bio1xbio15]_[priors_ddnorm]',
		'sdm_[nmixture]_[pseudoabsences_0.99]_[bio1^2_bio12^2_bio15^2_bio1xbio15_bio12xbio15]_[priors_ddnorm]',
		'sdm_[nmixture]_[pseudoabsences_0.99]_[bio1^2_bio12^2_bio15^2_bio12xbio15]_[priors_ddnorm]'
	)

	results <- data.table()
	for (model in models_dirs) {
	
		chains <- readRDS(paste0('./outputs_loretta/', model, '/sdm_nmixture_chains.rds'))
	
	
	}



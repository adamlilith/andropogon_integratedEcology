### MODELING ANDROPOGON GERARDI DISTRIBUTION, PHENOTYPE, PHYSIOLOGY, GENOTYPE, and ASSOCIATED MICROBIAL COMMUNITIES
### Erica Newman | Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2023-12
###
### This script constructs an integrated species distribution model for Andropogon gerardi, where occurrence is a function of environmental variables and biomass, which is in turn also a function of environmental variables.
###
### source('C:/Ecology/R/andropogon_integratedEcology/sdm_pdm/sdm_pdm_03a_simple_univariate_phenotype_models.r')
### source('C:/Subarashi/R/andropogon_integratedEcology/sdm_pdm/sdm_pdm_03a_simple_univariate_phenotype_models.r')
###
### CONTENTS ###
### setup ###
### user-defined values ###

#############
### setup ###
#############

	rm(list = ls())

	drive <- 'C:/Ecology/'
	# drive <- 'C:/Subarashi/'

	setwd(paste0(drive, '/Research/Andropogon/Andropogon'))
	source(paste0(drive, '/R/andropogon_integratedEcology/sdm_pdm/sdm_pdm_00_shared_functions_and_variables.r'))

###########################
### user-defined values ###
###########################

	predictors <- c('bio1', 'bio12', 'bio15', 'bio18', 'aridity')

	pheno_vars <- c('biomass', 'height', 'canopy diameter', 'blade width', 'leaf thickness', 'leaf nitrogen', 'photosynthetic rate', 'chholorphyll absorb', 'transpitaion rate', 'water use efficency', 'stomatal conductance', 'midway water potential', 'C isotope discrim')
	
	# "R name" is the variable as it appears in the spreadsheet
	# name in quites is what I am calling it
	pheno_vars <- c(
		Biomass = 'biomass',
		Height = 'vegetative height',
		CanopyDiam = 'canopy diameter',
		BladeWidth = 'blade width',
		LeafThick = 'leaf thickness',
		N_conc = 'leaf nitrogen',
		PhotoRate = 'photosynthetic rate',
		SPAD = 'chlorophyll absorbance',
		TranspRate = 'transpiration rate',
		# 'water use efficiency',
		StomCond = 'stomatal conductance',
		WatPot = 'midway water potential',
		Delta13C = 'C isotope discrimination'
	)

	pheno_units <- c(
		Biomass = 'g',
		Height = 'cm',
		CanopyDiam = 'cm',
		BladeWidth = 'cm',
		LeafThick = 'cm',
		N_conc = '?',
		PhotoRate = 'umol CO2m-2s-1',
		SPAD = '?',
		TranspRate = 'mmolH2Om-2s-1',
		# 'water use efficiency',
		StomCond = 'molH2Om-2s-1',
		WatPot = '?',
		Delta13C = 'δ¹³C'
	)

	# ### MCMC settings
	# niter <- 120000
	# nburnin <- 20000
	# thin <- 100
	# nchains <- 4
	# waic <- TRUE
	# trial <- FALSE

	# ### MCMC settings FOR TESTING
	niter <- 220
	nburnin <- 20
	thin <- 2
	nchains <- 2
	waic <- FALSE
	trial <- TRUE

	out_dir <- paste0('./outputs_loretta/sdm_pdm_simple_univariate_phenotype_models', ifelse(trial, '_TRIAL', NULL), '/')
	dirCreate(out_dir)

	# number of values in response curve array used to depict responses of SDM and PDM
	n_response_curve_values <- 200

#############################################################################################################
### model each phenotype variable and its standard deviation as a function of one environmental predictor ###
#############################################################################################################

	### load and collate data
	#########################

		site_data_raw <- readRDS('./data_from_loretta/sdm_pdm_00_merged_site_data_with_climate/sites.rds')
		biomass_data_raw <- readRDS('./data_from_loretta/sdm_pdm_00_merged_site_data_with_climate/biomass.rds')
		morpho_phys_data_raw <- readRDS('./data_from_loretta/sdm_pdm_00_merged_site_data_with_climate/morpho_phys.rds')

		# ensure sites in each array appear in the same order
		site_data_raw <- site_data_raw[order(site_id)]
		biomass_data_raw <- biomass_data_raw[order(SITE)]
		morpho_phys_data_raw <- morpho_phys_data_raw[order(SITE)]

		stopifnot(all(site_data_raw$site_id == unique(biomass_data_raw$SITE)))
		stopifnot(all(site_data_raw$site_id == unique(morpho_phys_data_raw$SITE)))

		n_pheno_sites <- length(unique(site_data_raw$site_id))
		
	### univariate models
	#####################

	for (pheno_var in pheno_vars) {

		pheno_var_ <- gsub(pheno_var, pattern = ' ', replacement = '_')

		sink(paste0(out_dir, '/runtime_log_', pheno_var, '.txt'), split = TRUE)
		say('SIMPLE UNIVARIATE MODEL OF ', toupper(pheno_var))
		say('Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | ', date(), post = 1)

		### data collation
		##################
		say('data collation', post = 1)

		say('This model simply relates the given phenotypic variable to one of several environmental predictors, one per model. For each environmental variable, four models are calibrated:', breaks = 60, post = 1)

		say('Linear model:')
		say('   y ~ N(mu, sigma)')
		say('   mu = beta0_model1_mu + beta1_model1_mu x')
		say('   sigma = beta0_model1_sigma + beta1_model1_sigma x', post = 1)

		say('Quadratic mu model:')
		say('   y ~ N(mu, sigma)')
		say('   mu = beta0_model2_mu + beta1_model2_mu x + beta2_model2_mu x^2')
		say('   sigma = beta0_model2_sigma + beta1_model2_sigma x', post = 1)

		say('Quadratic sigma model:')
		say('   y ~ N(mu, sigma)')
		say('   mu = beta0_model3_mu + beta1_model3_mu x')
		say('   sigma = beta0_model3_sigma + beta1_model3_sigma x + beta2_model3_sigma x^2', post = 1)

		say('Quadratic linear and sigma model:')
		say('   y ~ N(mu, sigma)')
		say('   mu = beta0_model4_mu + beta1_model4_mu x + beta1_model4_mu x^2')
		say('   sigma = beta0_model4_sigma + beta1_model4_sigma x + beta2_model4_sigma x^2', post = 1)

		say('MCMC settings:', level = 2)
		say('niter ....... ', niter)
		say('nburnin ..... ', nburnin)
		say('thin ........ ', thin)
		say('nchains ..... ', nchains, post = 1)

		### collate data
		################

		for (predictor in predictors) {

			# center and scale at site level
			x <- site_data_raw[ , .predictors]
			x <- scale(x)
			x_centers <- attr(x, 'scaled:center')
			x_scales <- attr(x, 'scaled:scale')
			x <- as.data.frame(x)
			
			x_mu_model_1 <- model.matrix(form_1, x)
			x_mu_model_2 <- model.matrix(form_2, x)
			x_mu_model_3 <- model.matrix(form_3, x)
			x_mu_model_4 <- model.matrix(form_4, x)
			
			biomass_by_site_x_sigma <- model.matrix(pheno_form_sigma, x)

			# make vector of which sampled site matches each row in the biomass and morphology/physiology data
			n_biomass <- nrow(biomass_data_raw)

			biomass_site_index <- rep(NA, n_biomass)
			for (i in 1:n_biomass) {
				biomass_site_index[i] <- which(site_data_raw$site_id == biomass_data_raw$SITE[i])
			}

		} # next predictor


	} # next phenotypic variable



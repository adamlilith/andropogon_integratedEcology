### MODELING OF MICROBES RELATED TO ANDROPOGON GERARDI: CONSTRUCT MODEL(S)
### Erica Newman & Adam B. Smith | 2024-07
###
### This script calibrates and evaluates a Hierarchical Modeling of Species Communities (HMSC) model for microbes associated with Andropogon gerardi.
###
### source('C:/Ecology/Research/Andropogon/Andropogon/andropogon_integratedEcology/hmsc_vs_microbes/microbes_03_run_hmsc_model.r')
### source('E:/Adam/Research/Andropogon/Andropogon/andropogon_integratedEcology/hmsc_vs_microbes/microbes_03_run_hmsc_model.r')
###
### CONTENTS ###
### setup ###
### settings ###
### model setup ###
### run model ###
### assess convergence ###
### assess model fit ###
### parameter estimates ###
### predictions ###
### create rasters of predictions ###
### make maps of predicted abundance for each taxon ###

#############
### setup ###
#############

rm(list = ls())
set.seed(1)

# drive <- 'C:/Ecology/'
drive <- 'E:/Adam/'

# workDir <- 'C:/Ecology/Research/Andropogon/Andropogon/'
workDir <- 'E:/Adam/Research/Andropogon/Andropogon/'

setwd(workDir)

library(abind) # for arrays
library(BayesLogit) # sometimes HMSC complains if we don't attach this
library(bayesplot) # for Bayesian model diagnostic plots
library(cowplot) # for combining multiple ggplots
library(data.table) # large data frames
library(enmSdmX) # GIS & SDMing
library(geodata) # geographic data
library(ggplot2) # plotting
library(ggspatial) # plotting spatial objects
library(Hmsc) # workhorse
library(ape) # we need this to construct a taxonomic tree
library(ggplot2) # plots
library(omnibus) # utilities
library(viridis) # colors

################
### settings ###
################

	# output_subfolder <- 'hmsc_[sampling_linear]_[climate_linear_quadratic]_[soil_linear_quadratic]_[ag_linear_quadratic]_[sans_traits]_[sans_phylo]_[poisson]'
	output_subfolder <- 'hmsc_[sampling_linear]_[climate_linear_quadratic_ia]_[soil_linear_quadratic]_[ag]_[sans_traits]_[sans_phylo]_[sans_sac]_[poisson]'

	quant_common_species <- 0.5 # analyze species with sum of all abundances across sites >= this quantile
	# quant_common_species <- 0 # analyze species with sum of all abundances across sites >= this quantile

	# response <- 'lognormal poisson' # gives very extreme values
	response <- 'poisson'

	x_formula <- ~ + sampling_ppt_mm +
		summer_tmean_C +
		annual_ppt_mm + I(annual_ppt_mm^2) +
		summer_tmean_C:annual_ppt_mm +
		ppt_cv + I(ppt_cv^2) +
		pH + I(pH^2) + soc + I(soc^2) + soil_n + I(soil_n^2) + clay + I(clay^2) + sand + I(clay^2) +
		ag_lambda
	predictors <- c('sampling_ppt_mm', 'summer_tmean_C', 'annual_ppt_mm', 'ppt_cv', 'pH', 'soc', 'soil_n', 'clay', 'sand', 'ag_lambda')
	# predictors <- c('sampling_ppt_mm', 'summer_tmean_C', 'annual_ppt_mm', 'ppt_cv', 'pH', 'ag_lambda')
	
	# include_traits <- TRUE
	include_traits <- FALSE

	include_phylogeny <- TRUE # include phylogeny
	# include_phylogeny <- FALSE # do not include

	# number of MCMC iterations in the final result (ie, not number of total MCMC iterations!)
	samples <- 1000

	# burn-in
	# transient <- 1000
	transient <- NULL

	# thinning rate
	thin <- 100

	nParallel <- 1 # default: nParallel = nChains, set to 1 to disable parallel execution, values >1 --> no sampling!?!
	nChains <- 4

	nfolds <- 4 # folds for cross-validation

	dirCreate(paste0('./outputs_sonny/', output_subfolder))
	sink(paste0('./outputs_sonny/', output_subfolder, '/hmsc_runtime.txt'), split = TRUE)

		say('HMSC on microbes associated with the Andropogon gerardi rhizosphere')
		say(date(), post = 2)
		say('quant_common_species ........................ ', quant_common_species)
		say('samples ..................................... ', samples)
		say('transient ................................... ', ifelse(is.null(transient), 'NULL', transient))
		say('thin ........................................ ', thin)
		say('nChains ..................................... ', nChains)
		say('nfolds ...................................... ', nfolds)
		say('include phylogeny ........................... ', include_phylogeny)
		say('include traits .............................. ', include_traits)
		say('response .................................... ', response)

		say('x_formula:')
		say(x_formula)
		say('')

	sink()

say('###################')
say('### model setup ###')
say('###################')

	if (!include_traits) {
		trait_formula <- NULL
	} else {
		trait_formula <- ~ in_bulk
	}

	### site-by-species abundances/occurrence data ###
	##################################################
	say('abundance data', level = 3)

	# Column names are improper--need 'R-friendly' columns (no spaces) that correspond to each taxon.
	# We also need the column names to be 'class_rID~~~', with domain and phylum excluded because the
	# inclusion of a phylogeny will assume column names are as such.
	abundances <- read.csv('./data/data_from_sonny/compiled_for_modeling_with_hmsc/2024_07_26_erica/species_rhz_26JUL2024.csv',  as.is = FALSE, check.names = FALSE)
	
	# # add abundances across columns that represent the same taxon
	# taxon_names <- colnames(abundances)
	# taxon_names <- taxon_names[taxon_names %notin% c('my.id', 'id')]
	# rID_positions <- regexpr(taxon_names, pattern = 'rID_')
	# taxon_names <- substr(taxon_names, 1, rID_positions - 1)
	# taxon_names <- trimws(taxon_names)
	# unique_taxa <- unique(taxon_names)

	# collapsed_abundances <- abundances[ , c('my.id', 'id')]
	# for (i in seq_along(unique_taxa)) {
	
		# taxon <- unique_taxa[i]
		# columns_with_taxon <- which(taxon_names == taxon) + 2 # add 2 because of 'my.id' and 'id' columns
		# this_abund <- abundances[ , columns_with_taxon, drop = FALSE]
		# this_abund <- rowSums(this_abund)
		# this_abund <- data.frame(this_abund)
		# colnames(this_abund) <- taxon
		
		# collapsed_abundances <- cbind(collapsed_abundances, this_abund)
	
	# }

	# make nice taxon names
	taxa <- colnames(abundances)
	taxa <- strsplit(taxa, split = ' ')

	new_names <- rep(NA_character_, length(taxa))
	for (i in seq_along(taxa)) {

		if (length(taxa[[i]]) == 1) {
			new_names[i] <- taxa[[i]]
		} else {
			this_name <- taxa[[i]]
			this_name <- gsub(this_name, pattern = 'd__', replacement = '')
			this_name <- gsub(this_name, pattern = 'p__', replacement = '')
			this_name <- gsub(this_name, pattern = 'c__', replacement = '')
			this_name <- gsub(this_name, pattern = '__', replacement = 'unknown')
			this_name <- gsub(this_name, pattern = '-', replacement = '_')
			this_name <- paste(this_name, collapse = '_')
			new_names[i] <- this_name
		}

	}

	new_names <- trimws(new_names)
	names(abundances) <- new_names

	# create response matrix with abundances
	Y <- abundances[ , colnames(abundances) %notin% c('my.id', 'id')]
	Y <- as.matrix(Y)

	### site-by-environment data
	############################
	say('environmental data', level = 3)

	site_by_env <- fread('./data/data_from_sonny/compiled_for_modeling_with_hmsc/environment_rhz_26JUL2024.csv')

	# names(site_by_env)[names(site_by_env) == 'x-coordinate'] <- 'longitude'
	# names(site_by_env)[names(site_by_env) == 'y-coordinate'] <- 'latitude'
	# names(site_by_env)[names(site_by_env) == 'SoilC'] <- 'soc'
	# names(site_by_env)[names(site_by_env) == 'SoilN'] <- 'soil_n'
	# names(site_by_env)[names(site_by_env) == 'CLAY'] <- 'clay'
	# names(site_by_env)[names(site_by_env) == 'SAND'] <- 'sand'

		### add SDM-estimated suitability of Andropogon gerardi as predictor
		####################################################################

		# load spatial vector with AG presences
		# we will need to 1) burn the suitability values into this
		# and then 2) extract the suitability values to sampled sites,
		# and finally, 3) covert this to a raster so we can make predictions

		terms <- terms(x_formula)
		terms <- attr(terms, 'term.labels')
		if ('ag_lambda' %in% terms) {

			# SDM output (posterior estimates from non-integrated SDM)
			sdm_chains <- readRDS('./outputs_loretta/noinintegrated_sdm_with_bio2_without_pH/nonintegrated_sdm_chains.rds')
			sdm_summary <- sdm_chains$summary$all.chains

			which_lambda <- grepl(rownames(sdm_summary), pattern = 'lambda')
			lambda <- sdm_summary[which_lambda, ]

			# spatial vector into which to burn posteriors
			ag_vect <- vect('./data/occurrence_data/andropogon_gerardi_occurrences_with_environment.gpkg')

			ag <- as.data.frame(ag_vect)
			completes <- complete.cases(as.data.frame(ag_vect))
			ag_focus <- ag[completes, ]
			ag_vect_focus <- ag_vect[completes, ]
			
			# remove table... makes things faster
			for (i in ncol(ag_vect_focus):1) ag_vect_focus[ , i] <- NULL

			# attach mean of posteriors to vector
			ag_vect_focus$lambda_mean <- lambda[ , 'Mean']
			
			# extract AG suitability to sampled sites
			site_by_env_vect <- vect(site_by_env, geom = c('longitude', 'latitude'), crs = getCRS('WGS84'))
			site_by_env_vect <- project(site_by_env_vect, ag_vect_focus)

			sample_site_lambdas <- extract(ag_vect_focus, site_by_env_vect)
			site_by_env$ag_lambda <- sample_site_lambdas$lambda_mean
		
		}

	site_by_env$id <- as.factor(site_by_env$id)
	site_by_env$sampling_ppt_mm <- log1p(site_by_env$sampling_ppt_mm) # unskew

	site_by_env <- site_by_env[ , ..predictors]
	site_by_env <- scale(site_by_env)

	env_centers <- attr(site_by_env, 'scaled:center')
	env_scales <- attr(site_by_env, 'scaled:scale')

	site_by_env <- as.data.frame(site_by_env)
	
	### phylogenetic data ###
	#########################
	say('phylogenetic data', level = 3)

	taxonomy_traits <- read.csv('./data/data_from_sonny/compiled_for_modeling_with_hmsc/2024_07_26_erica/traits_rhz_26JUL2024.csv', as.is = FALSE)
	taxonomy_traits$intercept <- 1

	# # 'species' names must match those from Y
	# taxa <- taxonomy_traits[ , c('Domain', 'Phylum', 'Class')]
	# taxa <- apply(taxa, 2, gsub, pattern = 'd__', replacement = '')
	# taxa <- apply(taxa, 2, gsub, pattern = 'p__', replacement = '')
	# taxa <- apply(taxa, 2, gsub, pattern = 'c__', replacement = '')
	# taxa <- apply(taxa, 2, gsub, pattern = '__', replacement = 'unknown')
	# taxa <- apply(taxa, 2, gsub, pattern = '-', replacement = '_')
	# taxa <- apply(taxa, 2, trimws)
	# taxa_together <- apply(taxa, 1, paste, collapse = '_')
	
	# taxonomy_traits$taxon <- taxa_together

	# taxonomy_traits$Domain <- taxa[ , 'Domain']
	# taxonomy_traits$Phylum <- taxa[ , 'Phylum']
	# taxonomy_traits$Class <- taxa[ , 'Class']

	# # remove rows with duplicated taxa
	# dups <- duplicated(taxa_together)
	# taxonomy_traits <- taxonomy_traits[!dups, ]

	# taxonomy_traits$Domain <- factor(taxonomy_traits$Domain)
	# taxonomy_traits$Phylum <- factor(taxonomy_traits$Phylum)
	# taxonomy_traits$Class <- factor(taxonomy_traits$Class)

	# rownames(taxonomy_traits) <- taxonomy_traits$taxon

	# # check that columns of Y have same names as phylogeny
	# stopifnot(colnames(Y) == taxonomy_traits$taxon) # no need for manual inspection

	### select species
	##################
	say('select species', level = 3)

	total_abundances <- colSums(Y)
	threshold_abundance <- quantile(total_abundances, quant_common_species)

	# select most common species
	selected_species <- which(total_abundances >= threshold_abundance)
	
	Y <- Y[ , selected_species]
	taxonomy_traits <- taxonomy_traits[taxonomy_traits$taxon %in% colnames(Y), ]
	traits <- taxonomy_traits[ , c('intercept', 'in_bulk'), drop = FALSE]

	say('Number of selected taxa: ', ncol(Y))

	### study design
	################
	say('study design', level = 3)

	study_design <- read.csv('./data/data_from_sonny/compiled_for_modeling_with_hmsc/2024_07_26_erica/studyDesign_rhz_26JUL2024.csv')

	study_design$site <- factor(study_design$location)
	study_design$id <- factor(study_design$plant)

	# site as random effect
	site_effect <- HmscRandomLevel(units = levels(study_design$site))

	# ID as random effect if we want taxon associations at that level
	id_effect <- HmscRandomLevel(units = levels(study_design$id))

	### taxonomic tree
	##################
	say('taxonomic tree', level = 3)

	if (include_phylogeny) {

		# Since we don't have a true phylogeny, we'll use a taxonomic tree.
		taxonomic_tree <- as.phylo(~ Domain / Phylum / Class, data = taxonomy_traits, collapse = FALSE)
		taxonomic_tree$edge.length <- rep(1, length(taxonomic_tree$edge))

		png(paste0('./outputs_sonny/', output_subfolder, '/taxon_tree.png'), width = 1200, height = 1200, res = 300)
		plot(taxonomic_tree, cex = 0.5, no.margin = TRUE, show.node.label = TRUE, type = 'cladogram', edge.color = 'cornflowerblue')
		dev.off()
		
	}

	### graphs of abundance vs each environmental variable
	######################################################
	say('graphs of abundance vs each environmental variable', level = 3)
	
	# make a set of plots showing abundance vs each predictor
	# one plot per taxon
	# all plots for the same predictor compiled into a multi-panel plot
	# this multi-panel plot is saved as a file, one per predictor
	
	n_panel_rows <- 5 # number of rows of subpanels
	n_panel_cols <- 8 # number of columns of subpanels
	n_panels <- n_panel_rows * n_panel_cols
	n_actual_panels <- min(n_panels, ncol(Y))
	for (pred in predictors) {
	
		abund_vs_predictors <- list()
		for (i in seq_len(n_actual_panels)) {
		
			taxon <- colnames(Y)[i]
			x <- site_by_env[ , pred]
			x <- fields::unscale(x, env_centers[[pred]], env_scales[[pred]])
			data <- data.frame(x = x, abundance = Y[ , i] + 1)
			
			in_bulk_indicator <- taxonomy_traits$in_bulk[taxon == taxonomy_traits$taxon]
			in_bulk_color <- ifelse(in_bulk_indicator, 'lightgoldenrodyellow', 'lightgreen')
			
			abund_vs_predictors[[i]] <- ggplot(data, aes(x = x, y = abundance)) +
				geom_point() +
				geom_smooth(formula = y ~ x, method = 'lm', color = 'red', se = FALSE) +
				geom_smooth(formula = y ~ x + I(x^2), method = 'lm', se = FALSE) +
				xlab(pred) + ylab('abundance') +
				scale_y_log10() +
				ggtitle(taxon) +
				theme(
					panel.background = element_rect(fill = in_bulk_color),
					title = element_text(size = 4),
					axis.title = element_text(size = 5),
					axis.text = element_text(size = 4)
				)
		
		}
		
		abund_vs_predictors <- plot_grid(plotlist = abund_vs_predictors, ncol = n_panel_cols, nrow = n_panel_rows)
		ggsave(abund_vs_predictors, filename = paste0('./outputs_sonny/', output_subfolder, '/abundance_vs_', pred, '.png'), width = 12, height = 8, dpi = 200)
	
	}

	### construct model
	###################

	say('construct model', level = 3)

	# arguments to Hmsc (always included)
	args <- list(
		Y = Y, XFormula = x_formula, XData = site_by_env, XScale = FALSE,
		distr = response,
		studyDesign = study_design,
		ranLevels = list(site = site_effect, id = id_effect)

	)

	# trait arguments
	if (include_traits) {
		args$TrData <- traits
		args$TrFormula <- trait_formula
	}

	# phylogeny arguments
	if (include_phylogeny) {
		args$phyloTree <- taxonomic_tree
	}

	model <- do.call(Hmsc, args)
	sampleMcmc(model, samples = 2) # does this yield errors?

	# # obviates(?) error (from https://www.helsinki.fi/assets/drupal/2022-08/hmscdevelopment.pdf)
	# effect of using this is unknown
	# sampleMcmc(model, samples = 2, updater = list(GammaEta = FALSE)) 

	saveRDS(model, file = paste0('./outputs_sonny/', output_subfolder, '/unfitted_model.rds'))

say('#################')
say('### run model ###')
say('#################')

	### settings
	############

	# if (is.null(nParallel)) nParallel <- nChains # setting nParallel > 1 causes sampleHmsc() to fail to sample!

	# model <- readRDS('./models/unfitted_model.rds')

	### MCMC
	########
	say('MCMC', level = 3)

	fit <- sampleMcmc(
		hM = model, samples = samples, thin = thin,
		transient = ifelse(is.null(transient), ceiling(0.5 * samples * thin), transient),
		# adaptNf = rep(ceiling(0.4 * samples * thin), model$nr),
		# adaptNf = rep(ceiling(0.4 * samples * thin), 100),
		nChains = nChains, nParallel = 1, # NB nParallel >1 leads to failure to sample!
		initPar = 'fixed effects', # use MLE to generate initial guesses... much faster!
		# verbose = TRUE
		verbose = FALSE
	)

	saveRDS(fit, file = paste0('./outputs_sonny/', output_subfolder, '/fit_model.rds'))

say('##########################')
say('### assess convergence ###')
say('##########################')

	posteriors <- convertToCodaObject(fit, spNamesNumbers = c(TRUE, FALSE), covNamesNumbers = c(TRUE, FALSE))
	nr <- fit$nr

	### Gelman-Rubin diagnostics (R-hat)
	# We want all R-hat + upper CI values to be <=1.1

	sink(paste0('./outputs_sonny/', output_subfolder, '/model_convergence.txt'), split = TRUE)
	say('ASSESSING MODEL CONVERGENCE')
	say(date(), post = 1)

	say('Beta:', pre = 1)
	rhat_beta <- gelman.diag(posteriors$Beta, multivariate = FALSE)$psrf
	say('Mean R-hat: ', mean(rhat_beta[ , 1]))
	say('Maximum R-hat: ', max(rhat_beta[ , 1]))
	say('Maximum R-hat + upper CI: ', max(rowSums(rhat_beta)))

	rhat_beta <- as.data.frame(rhat_beta)
	names(rhat_beta)[1] <- 'estimate'
	rhat_plot <- ggplot(rhat_beta, aes(x = estimate)) +
		geom_histogram(binwidth = 0.5, fill = '#69b3a2', color = '#e9ecef') +
		geom_vline(xintercept = 1.1) +
		ggtitle('Betas')
		
	ggsave(rhat_plot, file = paste0('./outputs_sonny/', output_subfolder, '/rhat_betas.png'), width = 12, height = 8, dpi = 150)
	
	say('Gamma:', pre = 1)
	rhat_gamma <- gelman.diag(posteriors$Gamma, multivariate = FALSE)$psrf
	say('Mean R-hat: ', mean(rhat_gamma[ , 1]))
	say('Maximum R-hat: ', max(rhat_gamma[ , 1]))
	say('Maximum R-hat + upper CI: ', max(rowSums(rhat_gamma)))
	
	rhat_gamma <- as.data.frame(rhat_gamma)
	names(rhat_gamma)[1] <- 'estimate'
	rhat_plot <- ggplot(rhat_beta, aes(x = estimate)) +
		geom_histogram(binwidth = 0.5, fill = '#69b3a2', color = '#e9ecef') +
		geom_vline(xintercept = 1.1) +
		ggtitle('Gammas')

	ggsave(rhat_plot, file = paste0('./outputs_sonny/', output_subfolder, '/rhat_gamma.png'), width = 12, height = 8, dpi = 150)

	if (include_phylogeny) {

		say('Rho (phylogeny):', pre = 1)
		rhat_rho <- gelman.diag(posteriors$Rho, multivariate = FALSE)$psrf
		say('Mean R-hat: ', mean(rhat_rho[ , 1]))
		say('Maximum R-hat: ', max(rhat_rho[ , 1]))
		say('Maximum R-hat + upper CI: ', max(rowSums(rhat_rho)))
		
		rhat_gamma <- as.data.frame(rhat_rho)
		names(rhat_rho)[1] <- 'estimate'
		rhat_plot <- ggplot(rhat_rho, aes(x = estimate)) +
			geom_histogram(binwidth = 0.5, fill = '#69b3a2', color = '#e9ecef') +
			geom_vline(xintercept = 1.1) +
			ggtitle('Rhos')
		ggsave(rhat_plot, file = paste0('./outputs_sonny/', output_subfolder, '/rhat_rho.png'), width = 12, height = 8, dpi = 150)

	}

	if (nr > 0) {
		
		say('Omega (traits)', pre = 1)
		rhat_omega <- gelman.diag(posteriors$Omega[[1]], multivariate = FALSE)$psrf
		say('Mean R-hat: ', mean(rhat_omega[ , 1]))
		say('Maximum R-hat: ', max(rhat_omega[ , 1]))
		say('Maximum R-hat + upper CI: ', max(rowSums(rhat_omega)))
		
		rhat_omega <- as.data.frame(rhat_omega)
		names(rhat_omega)[1] <- 'estimate'
		rhat_plot <- ggplot(rhat_omega, aes(x = estimate)) +
			geom_histogram(binwidth = 0.5, fill = '#69b3a2', color = '#e9ecef') +
			geom_vline(xintercept = 1.1) +
			ggtitle('Omegas')
		ggsave(rhat_plot, file = paste0('./outputs_sonny/', output_subfolder, '/rhat_omega.png'), width = 12, height = 8, dpi = 150)

		# # spatial random effects
		# if (nr > 0) {
		# 	for (j in 1:nr) {
		# 		if (model$ranLevels[[j]]$sDim > 0) {
				
		# 			say(names(model$ranLevels)[j])
		# 			rhat <- gelman.diag(posteriors$Alpha[[j]], multivariate = FALSE)$psrf
		# 			say('Alpha ', j, ':')
		# 			say(rhat[, 1])
					
					# rhat_sre <- as.data.frame(rhat)
					# names(rhat_sre)[1] <- 'estimate'
					# rhat_plot <- ggplot(rhat_sre, aes(x = estimate)) +
						# geom_histogram(binwidth = 0.5, fill = '#69b3a2', color = '#e9ecef') +
						# geom_vline(xintercept = 1.1) +
						# ggtitle(paste('Alpha', j))
					# ggsave(rhat_plot, file = paste0('./outputs_sonny/', output_subfolder, '/rhat_alpha_', j, '.png'), width = 12, height = 8, dpi = 150)

		# 		}
		# 	}
		# }

	}
	
	sink()

say('########################')
say('### assess model fit ###')
say('########################')

	# level at which to conduct cross-validation
	# NB we need to use at least a site-level random effect
	cv.level <- 'site'

	# predictions without/with cross-validation
	partition <- createPartition(fit, nfolds = nfolds, column = cv.level)
	preds <- computePredictedValues(fit, verbose = FALSE)
	preds_cv <- computePredictedValues(fit, partition = partition, verbose = FALSE)

	model_fit <- evaluateModelFit(hM = fit, predY = preds)
	model_fit_cv <- evaluateModelFit(hM = fit, predY = preds_cv)
	waic <- computeWAIC(fit)
	
	say('WAIC: ', waic)
	
	metrics <- c('RMSE', 'O.RMSE', 'C.RMSE', 'SR2', 'C.SR2')
	fit_plots <- list()
	for (metric in metrics) {
		
		# focal values
		vals <- model_fit[[metric]]
		vals_cv <- model_fit_cv[[metric]]
		
		df <- data.frame(x = vals, y = vals_cv)
		
		# axis limits
		lower <- min(vals, vals_cv)
		upper <- max(vals, vals_cv)
		lims <- c(lower, upper)
		lims <- pretty(lims)
		lims <- c(lims[1], lims[length(lims)])

		fit_plots[[length(fit_plots) + 1]] <- ggplot(df) +
			geom_point(mapping = aes(x = x, y = y), size = 3) +
			xlim(lims[1], lims[2]) + ylim(lims[1], lims[2]) + coord_fixed() +
			geom_abline(slope = 1, intercept = 0) +
			xlab('Explanatory power') + ylab('Predictive power') +
			ggtitle(metric)
	
	}
	
	all_fit_plots <- plot_grid(plotlist = fit_plots) +
		theme(plot.background = element_rect(fill = 'white'))

	ggsave(all_fit_plots, file = paste0('./outputs_sonny/', output_subfolder, '/', model_fit.png), width = 14, height = 10, dpi = 300)

say('###########################')
say('### parameter estimates ###')
say('###########################')

	# Alpha is spatial scale parameter (pages 100-101)

	posteriors <- convertToCodaObject(fit, spNamesNumbers = c(TRUE, FALSE), covNamesNumbers = c(TRUE, FALSE))
	
	### Betas
	#########
	samples <- posteriors$Beta
	
	# make one caterpillar plot per model term
	model_terms <- attr(terms(model$XFormula), 'term.labels')
	for (model_term in model_terms) {

		say(model_term)

		pars <- paste0('B[', model_term, ', ', taxonomy_traits$taxon, ']')
		# if (grepl(model_term, pattern = '\\^2')) {
			# regx <- '\\^2)'
			# caterpillars <- mcmc_intervals(samples, pars = pars, regex_pars = regx)
		# } else {
			caterpillars <- mcmc_intervals(samples, pars = pars)
		# }

		model_term_nice <- model_term
		model_term_nice <- gsub(model_term_nice, pattern = 'I\\(', replacement = '')
		model_term_nice <- gsub(model_term_nice, pattern = '\\)', replacement = '')

		# plot posterior of coefficient estimates
		caterpillars <- caterpillars +
			# scale_y_discrete(labels = cols) +
			xlab(model_term_nice) +
			ggtitle(paste0('Betas for ', model_term_nice)) +
			theme(plot.background = element_rect(fill = 'white'))

		file_name <- model_term
		file_name <- gsub(file_name, pattern = 'I\\(', replacement = '')
		file_name <- gsub(file_name, pattern = '\\^2\\)', replacement = '_2')

		ggsave(caterpillars, file = paste0('./outputs_sonny/', output_subfolder, '/coefficient_estimates_', file_name, '.png'), width = 10, height = 8)

	}

	### Omega 1
	###########
	
	# populate n x n matrix of mean Omega1 values
	samples <- posteriors$Omega[[1]]
	
	means <- lapply(samples, colMeans)
	means <- do.call(rbind, means)
	means <- colMeans(means)
	means_cols <- names(means)
	means_cols <- strsplit(means_cols, split = ', ')
	means_cols <- do.call(rbind, means_cols)
	taxa1 <- means_cols[ , 1]
	taxa2 <- means_cols[ , 2]
	
	taxa <- taxonomy_traits$taxon
	
	omegas <- expand.grid(taxon1 = taxa, taxon2 = taxa)
	omegas$omega <- NA_real_
	
	for (i in 1:nrow(omegas)) {
		
		taxon1 <- grepl(taxa1, pattern = omegas$taxon1[i])
		taxon2 <- grepl(taxa2, pattern = omegas$taxon2[i])
		this_col <- taxon1 & taxon2
		omegas$omega[i] <- means[this_col]
	
	}	
	
	omega_1_plot <- ggplot(omegas, aes(x = taxon1, y = taxon2, fill = omega)) +
		geom_tile() +
		scale_fill_gradient2(
			low = 'red',
			mid = 'white',
			high = 'forestgreen'
		) +
		ggtitle('Omega 1') +
		coord_fixed() +
		theme(
			axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1),
			axis.title = element_blank()
		)

	### Omega 2
	###########
	
	# populate n x n matrix of mean Omega1 values
	samples <- posteriors$Omega[[2]]
	
	means <- lapply(samples, colMeans)
	means <- do.call(rbind, means)
	means <- colMeans(means)
	means_cols <- names(means)
	means_cols <- strsplit(means_cols, split = ', ')
	means_cols <- do.call(rbind, means_cols)
	taxa1 <- means_cols[ , 1]
	taxa2 <- means_cols[ , 2]
	
	taxa <- taxonomy_traits$taxon
	
	omegas <- expand.grid(taxon1 = taxa, taxon2 = taxa)
	omegas$omega <- NA_real_
	
	for (i in 1:nrow(omegas)) {
		
		taxon1 <- grepl(taxa1, pattern = omegas$taxon1[i])
		taxon2 <- grepl(taxa2, pattern = omegas$taxon2[i])
		this_col <- taxon1 & taxon2
		omegas$omega[i] <- means[this_col]
	
	}	
	
	omega_2_plot <- ggplot(omegas, aes(x = taxon1, y = taxon2, fill = omega)) +
		geom_tile() +
		scale_fill_gradient2(
			low = 'red',
			mid = 'white',
			high = 'forestgreen'
		) +
		ggtitle('Omega 2') +
		coord_fixed() +
		theme(
			axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1),
			axis.title = element_blank()
		)

	omegas_plot <- plot_grid(omega_1_plot, omega_2_plot) +
		theme(plot.background = element_rect(fill = 'white'))

	ggsave(omegas_plot, file = paste0('./outputs_sonny/', output_subfolder, '/omegas.png'), width = 18, height = 9, dpi = 600)

say('###################')
say('### predictions ###')
say('###################')

	### save data frames of predictions for the mean predicted abundance, SD, and (optionally, 10th and 90th quantiles)

	### obtain rasters for predictors
	conus_env <- readRDS('./outputs_sonny/conus_environment_as_per_prism_aggregated_by_16.rds')
	
	# add AG suitability as a predictor if it's in the model
	terms <- terms(x_formula)
	terms <- attr(terms, 'term.labels')
	if ('ag_lambda' %in% terms) {

		conus_env_vect <- vect(conus_env, geom = c('x', 'y'), crs = getCRS('NAD83'))
		ag_vect_focus_nad83 <- project(ag_vect_focus, conus_env_vect)
		ag_lambda <- extract(ag_vect_focus_nad83, conus_env_vect)
		conus_env$ag_lambda <- ag_lambda$lambda_mean
		
		conus_env <- conus_env[complete.cases(conus_env), ]
	
	}

	# create data tables for predictions	
	# predictions_mean <- predictions_sd <- predictions_quantile_0.1 <- predictions_quantile_0.9 <- conus_env[ , 'cell']
	predictions_mean <- predictions_sd <- conus_env[ , 'cell']

	# assume average amount of "sampling" rainfall
	conus_env[ , sampling_ppt_mm := mean(sampling_ppt_mm)]
	conus_env <- conus_env[ , ..predictors]

	# center and scale
	conus_env <- scale(conus_env, center = env_centers, scale = env_scales)
	conus_env <- as.data.table(conus_env)

	### create one column per taxon in each data table for predictions
	taxa <- colnames(Y)
	for (i in seq_along(abundances)) {

		predictions_mean[ , DUMMY := NA_real_]
		predictions_sd[ , DUMMY := NA_real_]
		# predictions_quantile_0.1[ , DUMMY := NA_real_]
		# predictions_quantile_0.9[ , DUMMY := NA_real_]

		taxon <- taxa[i]
		colnames(predictions_mean)[ncol(predictions_mean)] <- taxon
		colnames(predictions_sd)[ncol(predictions_sd)] <- taxon
		# colnames(predictions_quantile_0.1)[ncol(predictions_quantile_0.1)] <- taxon
		# colnames(predictions_quantile_0.9)[ncol(predictions_quantile_0.9)] <- taxon

	}

	if (any(colnames(predictions_mean) == 'DUMMY')) predictions_mean[ , DUMMY := NULL]
	if (any(colnames(predictions_sd) == 'DUMMY')) predictions_sd[ , DUMMY := NULL]
	# if (any(colnames(predictions_quantile_0.1) == 'DUMMY') predictions_quantile_0.1[ , DUMMY := NULL]
	# if (any(colnames(predictions_quantile_0.9) == 'DUMMY') predictions_quantile_0.9[ , DUMMY := NULL]

	### make predictions in chunks to save on memory (HMSC crashes otherwise)
	size <- round(nrow(conus_env) / 20)
	sets <- ceiling(nrow(conus_env) / size)
	for (set in seq_len(sets)) {

		say(set, ' of ', sets, ' ', date())

		indices <- (1 + size * (set - 1)):min(nrow(conus_env), size * set)

		# select environmental data
		x_data <- conus_env[indices, ..predictors]
		x_data <- as.data.frame(x_data)
		n_cells <- nrow(x_data)

		# random effects
		pred_site_effect <- data.frame(site = rep('KS', n_cells)) # assuming all places are Kansas
		pred_id_effect <- data.frame(id = rep(40, n_cells)) # assuming all samples are #40
		s_new <- list(site_effect = pred_site_effect, id_effect = pred_id_effect)

		# predict
		# creates a list with one element per posterior iteration, rows are cells of raster and columns are taxa
		grad <- prepareGradient(fit, XDataNew = x_data, sDataNew = s_new) # prepareGradient() is for new data
		preds <- predict(fit, Gradient = grad, expected = TRUE)

		# calculate mean, SD, quantiles across all MCMC samples
		preds <- abind(preds, along = 3)
		preds_mean <- apply(preds, c(1, 2), mean)
		preds_sd <- apply(preds, c(1, 2), sd)
		# preds_quantile_0.1 <- apply(preds, c(1, 2), quantile, probs = 0.1)
		# preds_quantile_0.9 <- apply(preds, c(1, 2), quantile, probs = 0.9)

		preds_mean <- as.data.table(preds_mean)
		preds_sd <- as.data.table(preds_sd)
		# preds_quantile_0.1 <- as.data.table(preds_quantile_0.1)
		# preds_quantile_0.9 <- as.data.table(preds_quantile_0.9)

		# remember predictions
		predictions_mean[indices, (taxa) := preds_mean]
		predictions_sd[indices, (taxa) := preds_sd]
		# predictions_quantile_0.1[indices, (taxa) := preds_quantile_0.1]
		# predictions_quantile_0.9[indices, (taxa) := preds_quantile_0.9]

	}

	saveRDS(predictions_mean, paste0('./outputs_sonny/', output_subfolder, '/predictions_to_prism_aggregated_by_16_mean.rds'))
	saveRDS(predictions_sd, paste0('./outputs_sonny/', output_subfolder, '/predictions_to_prism_aggregated_by_16_sd.rds'))
	# saveRDS(predictions_quantile_0.1, paste0('./outputs_sonny/', output_subfolder, '/predictions_to_prism_aggregated_by_16_quantile_0.1.rds'))
	# saveRDS(predictions_quantile_0.9, paste0('./outputs_sonny/', output_subfolder, '/predictions_to_prism_aggregated_by_16_quantile_0.9.rds'))

say('#####################################')
say('### create rasters of predictions ###')
say('#####################################')

	template <- rast(paste0(drive, '/Research Data/PRISM/PRISM_us_dem_800m.tif'))
	template <- aggregate(template, 16)
	template[] <- NA_real_

	for (i in seq_along(taxa)) {

		taxon <- taxa[i]
		say(taxon)

		### make maps of MEAN prediction
		# for infinite predictions to the maximum non-infinite value
		predictions_filtered <- predictions_mean[[taxon]]
		if (any(is.infinite(predictions_filtered))) {
			non_infinite_max <- max(predictions_filtered[!is.infinite(predictions_filtered)])
			predictions_filtered[is.infinite(predictions_filtered)] <- non_infinite_max
		}

		predictions_filtered <- log10(predictions_filtered)
		map_this_taxon <- setValueByCell(template, val = predictions_filtered, cell = predictions_mean[['cell']])
		names(map_this_taxon) <- taxon

		# ### make maps of MEAN prediction
		# # for infinite predictions to the maximum non-infinite value
		# predictions_filtered <- predictions_mean[[taxon]]
		# map_this_taxon <- setValueByCell(template, val = predictions_filtered, cell = predictions_mean[['cell']])
		# names(map_this_taxon) <- taxon

		if (i == 1) {
			maps_mean <- map_this_taxon
		} else {
			maps_mean <- c(maps_mean, map_this_taxon)
		}

		### make maps of SDM of prediction
		# for infinite predictions to the maximum non-infinite value
		predictions_filtered <- predictions_sd[[taxon]]
		if (any(is.infinite(predictions_filtered))) {
			non_infinite_max <- max(predictions_filtered[!is.infinite(predictions_filtered)])
			predictions_filtered[is.infinite(predictions_filtered)] <- non_infinite_max
		}

		predictions_filtered <- log10(predictions_filtered)
		map_this_taxon <- setValueByCell(template, val = predictions_filtered, cell = predictions_sd[['cell']])
		names(map_this_taxon) <- taxon

		if (i == 1) {
			maps_sd <- map_this_taxon
		} else {
			maps_sd <- c(maps_sd, map_this_taxon)
		}

	}

	writeRaster(maps_mean, paste0('./outputs_sonny/', output_subfolder, '/prediction_rasters_mean.tif'), overwrite = TRUE)
	writeRaster(maps_sd, paste0('./outputs_sonny/', output_subfolder, '/prediction_rasters_sd.tif'), overwrite = TRUE)

say('#######################################################')
say('### make maps of predicted abundance for each taxon ###')
say('#######################################################')

	### maps of means
	maps_mean <- rast(paste0('./outputs_sonny/', output_subfolder, '/prediction_rasters_mean.tif'))

	# North American countries
	nam <- gadm(c('CAN', 'USA', 'MEX'), level = 1, path = tempdir(), resolution = 2)
	# nam <- vect(paste0(drive, './Research Data/GADM/Version 4.1/High Res North America Level 1 sans Great Lakes SpatVector WGS84.gpkg'))
	nam <- project(nam, getCRS('North America Lambert'))

	### collate sites coordinates and abundances so we can plot observed abundances at sites
	sites_abundances <- fread('./data/data_from_sonny/compiled_for_modeling_with_hmsc/environment_rhz_26JUL2024.csv')
	colnames(sites_abundances)[colnames(sites_abundances) == 'x-coordinate'] <- 'latitude'
	colnames(sites_abundances)[colnames(sites_abundances) == 'y-coordinate'] <- 'longitude'

	# add abundances
	taxa <- colnames(Y)
	for (i in seq_along(taxa)) {

		taxon <- taxa[i]
		this_abundance <- collapsed_abundances[match(collapsed_abundances$my.id, sites_abundances$my_id), taxon]
		sites_abundances[ , (taxon) := this_abundance]

	}

	sites_abundances <- vect(sites_abundances, geom = c('latitude', 'longitude'), crs = getCRS('NAD83'))
	sites_abundances <- project(sites_abundances, getCRS('North America Lambert'))
	
	# plotting extent... crop at -107 west, plus a bit larger
	extent <- ext(maps_mean)
	extent <- as.vector(extent)
	extent[1] <- -108
	extent <- ext(extent)
	extent <- as.polygons(extent, crs = crs(maps_mean))
	extent <- buffer(extent, 300000)
	maps_mean <- crop(maps_mean, extent)
	
	# ...but actually plot a smaller area (looks nicer)
	extent <- buffer(extent, -600000)
	extent_for_plot <- project(extent, getCRS('North America Lambert'))
	extent <- as.vector(extent_for_plot)

	pdf(paste0('./outputs_sonny/', output_subfolder, '/prediction_maps.pdf'), width = 9, height = 6)
	for (i in seq_along(taxa)) {
	
		taxon <- taxa[i]
		say(taxon)

		# get map and project		
		this_map <- maps_mean[[taxon]]
		this_map <- project(this_map, getCRS('North America Lambert'))
		# this_map <- crop(this_map, extent_for_crop)

		# to define colors, get min/max abundance across sites and the map
		this_abundance <- as.data.frame(sites_abundances[ , taxon])[ , 1, drop = TRUE]
		this_abundance <- log10(this_abundance + 1)
		map_min <- minmax(this_map)[1, ]
		map_max <- minmax(this_map)[2, ]
		min_abund <- min(map_min, this_abundance)
		
		# remove extreme abundances
		max_abund <- globalx(this_map, quantile, probs = 0.99, na.rm = TRUE)
		this_map <- app(this_map, fun = function(x, max_abund) ifelse(x > max_abund, 0.999 * max_abund, x), max_abund = max_abund)

		if (is.infinite(min_abund)) {
		
			neg_inf_abund_masked <- app(this_map, fun = function(x) ifelse(is.infinite(x), max_abund, x))
			min_abund <- globalx(neg_inf_abund_masked, min)
			this_map <- app(this_map, fun = function(x, min_abund) ifelse(x < min_abund, min_abund, x), min_abund = min_abund)
			
		}

		abund_seq <- seq(min_abund, max_abund, length.out = 100)
		this_abundance_scaled <- (this_abundance - min_abund) / (max_abund - min_abund)
		this_abundance_scaled <- round(100 * this_abundance_scaled)

		pallette <- viridis(n = 100)
		site_colors <- pallette[this_abundance_scaled]

		plot(nam, col = 'gray30', ext = extent_for_plot, lwd = 0.1, main = taxon)
		plot(this_map, range = c(min_abund, max_abund), add = TRUE)
		plot(nam, border = 'gray20', lwd = 0.1, add = TRUE)

		plot(sites_abundances, pch = 21, bg = site_colors, cex = 1.6, add = TRUE)
	
	}
	dev.off()

say('DONE', deco = '!', level = 1)

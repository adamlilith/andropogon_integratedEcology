### MODELING OF MICROBES RELATED TO ANDROPOGON GERARDI: CONSTRUCT MODEL(S)
### Erica Newman & Adam B. Smith | 2024-07
###
### This script calibrates and evaluates a Hierarchical Modeling of Species Communities (HMSC) model for microbes associated with Andropogon gerardi.
###
### source('C:/Ecology/R/andropogon_integratedEcology/hmsc_vs_microbes/microbes_02_run_hmsc_model.R')
### source('C:/Subarashi/R/andropogon_integratedEcology/hmsc_vs_microbes/microbes_02_run_hmsc_model.R')
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

	### EXPERIMENT
	# 1 response: lognormal poisson | predictors: raw | phylo: true | STATUS: IP
	# 1 response: lognormal poisson | predictors: PCs | phylo: true | STATUS: WAITING

	# 2 response: poisson | predictors: raw | phylo: true | STATUS: IP
	# 2 response: poisson | predictors: PCs | phylo: true | STATUS: WAITING
	
	# 3 response: normal | predictors: raw | phylo: true | STATUS: IP
	# 3 response: normal | predictors: PCs | phylo: true | STATUS: WAITING
	

	rm(list = ls())
	set.seed(1)

	# drive <- 'C:/Ecology/'
	drive <- 'C:/Subarashi/'
	# drive <- 'E:/!Scratch/'

	workDir <- paste0(drive, './Research/Andropogon/Andropogon')

	setwd(workDir)

	# .libPaths('H:/Global Change Program/libraries')

	library(abind) # for arrays
	library(BayesLogit) # sometimes HMSC complains if we don't attach this
	library(bayesplot) # for Bayesian model diagnostic plots
	library(cowplot) # for combining multiple ggplots
	library(data.table) # large data frames
	library(enmSdmX) # GIS & SDMing
	library(geodata) # geographic data
	library(ggspatial) # plotting spatial objects
	library(Hmsc) # workhorse
	library(ape) # we need this to construct a taxonomic tree
	library(ggplot2) # plots
	library(omnibus) # utilities
	library(parallel) # parallelization
	library(snow) # parallelization
	library(tictoc) # timing
	library(viridis) # colors

################
### settings ###
################

	# run abbreviated model... for debugging/development
	abbreviated <- FALSE # run full model!
	# abbreviated <- TRUE

	# run analysis at phylum or class level
	# phylum_or_class <- 'phylum'
	phylum_or_class <- 'class'

	# make graphs of each environmental predictor vs abundance? (saves times... should be FALSE only for testing script)
	make_abundance_plots <- TRUE
	# make_abundance_plots <- FALSE

	# use this model of Andropogon gerardi
	ag_model <- 'nmixture' # N-mixture model
	# ag_model <- 'bernoilli' # Bernoulli (binary) model

	# response <- 'lognormal poisson' # response is integers >= 0
	response <- 'poisson' # response is integers >= 0
	# response <- 'normal' # response is log of abundance, with 0 values forced to NA

	# model taxa that may occur only in rhizobiome samples?
	include_rhizo <- TRUE
	# include_rhizo <- FALSE

	# model taxa that may occur only in bulk samples? 
	include_bulk <- TRUE
	# include_bulk <- FALSE

	# model taxa that occur *only* in both rhizobiome and bulk? (this forces these taxa to be included, even if they do not meet the abundance criterion below)... seems this should nearly always be TRUE (only has 12 taxa... 14 is "unknowns" included)
	# if this is TRUE, then include_rhizo and include_bulk should be FALSE!
	# just_both <- TRUE
	just_both <- FALSE

	# model "unknown" Classes?
	# include_unknown_classes <- TRUE
	include_unknown_classes <- FALSE

	# analyze taxa with sum of all abundances across sites >= this quantile
	# rhizobiome-only and bulk-only taxa will be filtered separately (ie, if include_rhizobiome is TRUE, then taxa that *only* occur in the rhizobiome and have an abundance >= quant_abund_threshold will be modeled; same for include_bulk)
	# quant_abund_threshold <- 0.95 # 0.95 is good for testing
	# quant_abund_threshold <- 0.90
	quant_abund_threshold <- 0.75
	# quant_abund_threshold <- 0.50
	# quant_abund_threshold <- 0 # value of 0 ==> all taxa

	# minimum number of sites (not plants) taxon must be present in to model
	min_sites_present <- 13 # 26 total sites

	# include phylogeny?
	include_phylogeny <- TRUE
	# include_phylogeny <- FALSE

	# model predictors
	use_pc_axes_as_predictors <- TRUE
	# use_pc_axes_as_predictors <- FALSE

	# number of MCMC iterations in the final result (ie, not number of total MCMC iterations!)
	samples <- if (abbreviated) { 100 } else { 2000 }

	# burn-in
	# transient <- 1000
	transient <- NULL

	# thinning rate
	thin <- 100
	# thin <- 50
	# thin <- 1

	n_parallel <- if (abbreviated) { 1 } else { 1 } # default: n_parallel = n_chains, set to 1 to disable parallel execution
	n_chains <- if (abbreviated) { 2 } else { 4 }

	n_folds <- if (abbreviated) { 2 } else { 4 }

	# raw_predictors <- c('aridity', 'bio7', 'bio12', 'bio15', 'ph_field', 'sand_field', 'soc_field_perc', 'silt_field', 'ag_lambda_nmixture', 'sampling_ppt_mm', 'sampling_tmean_c')
	raw_predictors <- c('bio1', 'bio12', 'bio15', 'ph_field', 'sand_field', 'ag_lambda_nmixture', 'sampling_ppt_mm', 'sampling_tmean_c')
	pc_predictors <- c('pc1', 'pc2', 'pc3', 'rhizobiome_or_bulk')

	if (just_both & (include_rhizo | include_bulk)) stop('If just_both if TRUE, then include_rhizo and include_bulk should be FALSE.')

	if ((include_rhizo & include_bulk) | just_both) {

		raw_predictors <- c(raw_predictors, 'rhizobiome_or_bulk')

		if (use_pc_axes_as_predictors) {

			x_formula <- '~
				pc1 + I(pc1^2) +
				pc2 + I(pc2^2) +
				pc3 + I(pc3^2) +
				rhizobiome_or_bulk'

		} else {

			x_formula <- '~
				sampling_ppt_mm +
				sampling_tmean_c +

				bio1 + I(bio1^2) +
				bio12 + I(bio12^2) +
				bio15 + I(bio12^2) +
				
				ph_field + I(ph_field^2) +
				sand_field + I(sand_field^2) +
				
				rhizobiome_or_bulk'

			if (ag_model == 'nmixture') {
				x_formula <- paste0(x_formula, ' + ag_lambda_nmixture')
			} else if (ag_model == 'bernoulli') {
				x_formula <- paste0(x_formula, ' + ag_lambda_bernoulli')
			}

		}

	} else {
	
		if (use_pc_axes_as_predictors) {

			x_formula <- '~
				pc1 + I(pc1^2) +
				pc2 + I(pc2^2) +
				pc3 + I(pc3^2)'

		} else {

			x_formula <- '~
				sampling_ppt_mm +
				sampling_tmean_c +

				bio1 + I(bio1^2) +
				bio12 + I(bio12^2) +
				bio15 + I(bio12^2) +

				ph_field + I(ph_field^2) +
				sand_field + I(sand_field^2)'

			if (ag_model == 'nmixture') {
				x_formula <- paste0(x_formula, ' + ag_lambda_nmixture')
			} else if (ag_model == 'bernoulli') {
				x_formula <- paste0(x_formula, ' + ag_lambda_bernoulli')
			}

		}
		
	}
	x_formula <- gsub(x_formula, pattern = '\\n', replacement = ' ')
	x_formula <- as.formula(x_formula)

	### maker folder to which to save results

	resp_string <- paste(
		ifelse(include_rhizo, 'rhizo', ''),
		ifelse(include_bulk, 'bulk', ''),
		ifelse(just_both, 'both', ''),
		sep = ' '
	)
	resp_string <- trimws(resp_string)
	resp_string <- gsub(resp_string, pattern = ' ', replacement = '_')

	header <- if (abbreviated) { 'TEMP_' } else { '' }

	unknowns <- if (include_unknown_classes) { 'all_classes'} else { 'sans_unknown_classes' }

	out_dir <- paste0(
		'./outputs_sonny/', header, 'hmsc_[', phylum_or_class, ']_[',
		ifelse(use_pc_axes_as_predictors, 'pcs', 'climate'),
		']_[abund_',
		quant_abund_threshold, '_min_sites_', min_sites_present, ']_[',
		resp_string,
		']_[',
		ifelse(include_phylogeny, 'phylo', 'sans_phylo'),
		']_[ag_', ag_model, ']_[',
		sub(response, pattern = '_', replacement = '_'),
		']_[', unknowns, ']'

	)

	dirCreate(out_dir)
	sink(paste0(out_dir, '/hmsc_settings.txt'), split = TRUE)

		say('HMSC on microbes associated with the Andropogon gerardi rhizosphere')
		say(date(), post = 2)
		say('ag_model .................................... ', ag_model)
		say('include_rhizo ............................... ', include_rhizo)
		say('include_bulk ................................ ', include_bulk)
		say('just_both ................................... ', just_both)
		say('include_unknown_classes ..................... ', include_unknown_classes)
		say('quant_abund_threshold ....................... ', quant_abund_threshold)
		say('min_sites_present ........................... ', min_sites_present)
		say('samples ..................................... ', samples)
		say('transient ................................... ', ifelse(is.null(transient), 'NULL', transient))
		say('thin ........................................ ', thin)
		say('n_chains .................................... ', n_chains)
		say('n_folds ..................................... ', n_folds)
		say('include phylogeny ........................... ', include_phylogeny)
		say('response .................................... ', response)
		say('use_pc_axes_as_predictors ................... ', use_pc_axes_as_predictors, post = 2)

		say('raw_predictors:')
		say(paste(raw_predictors, collapse = ' '), post = 2)

		say('x_formula:')
		say(x_formula, post = 2)

	sink()

######################
### data collation ###
######################

	say('data collation')

	abund_combined <- fread(paste0('./data_from_sonny/collated_for_hmsc_collapsed_to_', phylum_or_class, '/abundances_site_by_taxon_combined.csv'))
	env_combined <- fread(paste0('./data_from_sonny/collated_for_hmsc_collapsed_to_', phylum_or_class, '/environment_combined.csv'))
	study_design_combined <- fread(paste0('./data_from_sonny/collated_for_hmsc_collapsed_to_', phylum_or_class, '/study_design_combined.csv'))
	taxa_combined <- fread(paste0('./data_from_sonny/collated_for_hmsc_collapsed_to_', phylum_or_class, '/taxa_combined.csv'))

	# get list of taxa we want to analyze
	if (!include_rhizo & !include_bulk & just_both) {

		keeps <- taxa_combined$in_rhizobiome & taxa_combined$in_bulk
		taxa <- taxa_combined[(keeps)]
		abund <- abund_combined # subset columns below
		study_design <- study_design_combined
		env <- env_combined

	} else if (include_rhizo & include_bulk & !just_both) {

		taxa <- taxa_combined
		abund <- abund_combined
		study_design <- study_design_combined
		env <- env_combined

	} else if (include_rhizo & !include_bulk & !just_both) {

		taxa <- taxa_combined[(in_rhizobiome)]
		abund <- abund_combined[rhizobiome_or_bulk == 'rhizobiome']
		study_design <- study_design_combined[study_design_combined$rhizobiome_or_bulk == 'rhizobiome']
		env <- env_combined[rhizobiome_or_bulk == 'rhizobiome']

	} else if (!include_rhizo & include_bulk & !just_both) {

		taxa <- taxa_combined[(in_bulk)]
		abund <- abund_combined[rhizobiome_or_bulk == 'bulk']
		study_design <- study_design_combined[rhizobiome_or_bulk == 'bulk']
		env <- env_combined[rhizobiome_or_bulk == 'bulk']

	}

	site_xy <- env[  , c('longitude', 'latitude')]

	# remove Eukaryota
	names <- names(abund)
	keeps <- !grepl(names, pattern = 'Eukaryota')
	abund <- abund[ , ..keeps]

	keeps <- taxa$domain != 'Eukaryota'
	taxa <- taxa[keeps, ]

	# remove taxa with unknown phylum
	names <- names(abund)
	keeps <- !grepl(names, pattern = '_unknown_')
	abund <- abund[ , ..keeps]
	
	taxa <- taxa[taxa$phylum != 'unknown']

	# if not retaining "unknown" taxa
	if (!include_unknown_classes) {

		# remove unknown classes from abundance matrix
		names <- names(abund)
		nc <- nchar(names)
		names <- substr(names, nc - nchar('_unknown') + 1, nc)
		keeps <- which(names != '_unknown')
		abund <- abund[ , ..keeps]

		# remove unknown classes from taxon matrix
		if (phylum_or_class == 'class') taxa <- taxa[taxa$class != 'unknown']

	}

	# remove taxa in fewer than minimum number of sites
	abund_combined_by_site <- aggregate(abund, by = list(abund_combined$location), mean)
	abund_combined_by_site$location <- NULL
	names(abund_combined_by_site)[1] <- 'location'
	discards <- c('location', 'index', 'plant', 'sample', 'rhizobiome_or_bulk')
	abund_combined_by_site <- abund_combined_by_site[ , !(colnames(abund_combined_by_site) %in% discards)]
	abund_combined_by_site <- abund_combined_by_site > 0
	taxa_sites_present <- colSums(abund_combined_by_site)

	keep_taxa <- names(taxa_sites_present)[taxa_sites_present > min_sites_present]
	abund <- abund_combined[ , ..keep_taxa]

	taxa <- taxa[taxa$taxon %in% keep_taxa]

	# thin abundance matrix to selected taxa
	cond <- colnames(abund) %in% taxa$taxon
	abund <- abund[ , ..cond]

	if (response == 'normal') abund[abund == 0] <- NA

	# filter by abundance threshold
	if (quant_abund_threshold > 0) {
	
		n <- colSums(abund, na.rm = TRUE)	
		threshold_n <- quantile(n, quant_abund_threshold, na.rm = TRUE)

		keeps <- which(n >= threshold_n)
		keep_taxa <- colnames(abund)[keeps]
		abund <- abund[ , ..keep_taxa]
		taxa <- taxa[taxon %in% keep_taxa]
	
	}

	# coerce character factor levels to integer
	env[ , rhizobiome_or_bulk := factor(rhizobiome_or_bulk)]
	env[ , rhizobiome_or_bulk := as.integer(rhizobiome_or_bulk) - 1]

	say('modeling ', nrow(taxa), ' taxa')

	### transform predictors
	########################

	env <- env[ , ..raw_predictors]
	
	cond <- colnames(env) %notin% 'rhizobiome_or_bulk'
	env_sans_categorical <- env[ , ..cond]

	if (use_pc_axes_as_predictors) {
		
		pca <- prcomp(env_sans_categorical, center = TRUE, scale. = TRUE)
		saveRDS(pca, paste0(out_dir, '/pca.rds'))

		sites_by_env <- predict(pca, env)
		sites_by_env <- as.data.table(sites_by_env)
		colnames(sites_by_env) <- tolower(colnames(sites_by_env))
		if ('rhizobiome_or_bulk' %in% raw_predictors) {
			sites_by_env[ , rhizobiome_or_bulk := factor(env$rhizobiome_or_bulk)]
		}

	} else {
		sites_by_env <- env_sans_categorical
	}
	sites_by_env <- as.data.table(sites_by_env)

	terms <- attr(terms(x_formula), 'term.labels')
	terms <- terms[!grepl(terms, pattern = ':')]
	terms <- terms[!grepl(terms, pattern = '\\^2')]

	sites_by_env[ , rhizobiome_or_bulk := env$rhizobiome_or_bulk]
	sites_by_env <- sites_by_env[ , ..terms]

	### study design
	################

	study_design[ , site := factor(as.integer(factor(study_design$location)))]
	study_design[ , id := factor(as.integer(factor(study_design$plant)))]
	study_design_coerced <- study_design[ , c('site', 'id')]

	# site as random effect
	site_effect <- HmscRandomLevel(units = levels(study_design_coerced$site))

	# plant as random effect if we want taxon associations at that level
	plant_effect <- HmscRandomLevel(units = levels(study_design_coerced$id))

	### taxonomic tree
	##################

	if (include_phylogeny) {

		# Since we don't have a true phylogeny, we'll use a taxonomic tree.
		taxa$domain <- factor(taxa$domain)
		taxa$phylum <- factor(taxa$phylum)
		taxa$taxon <- factor(taxa$taxon)

		taxonomic_tree <- as.phylo(~ domain / phylum / taxon, data = taxa, collapse = FALSE)
		taxonomic_tree$edge.length <- rep(1, length(taxonomic_tree$edge))

		png(paste0(out_dir, '/taxon_tree.png'), width = 1200, height = 1400, res = 300)
		plot(taxonomic_tree, cex = 0.4, no.margin = TRUE, show.node.label = TRUE, type = 'cladogram', edge.color = 'cornflowerblue')
		dev.off()
		
	}

	write.csv(abund, paste0(out_dir, '/abund.csv'), row.names = FALSE)
	write.csv(sites_by_env, paste0(out_dir, '/sites_by_env.csv'), row.names = FALSE)
	write.csv(study_design, paste0(out_dir, '/study_design.csv'), row.names = FALSE)

##########################################################
### graphs of abundance vs each environmental variable ###
##########################################################
		
	if (make_abundance_plots & !use_pc_axes_as_predictors) {

		say('graphs of abundance vs each environmental variable')
		
		# make a set of plots showing abundance vs each predictor
		# one plot per taxon
		# all plots for the same predictor compiled into a multi-panel plot
		# this multi-panel plot is saved as a file, one per predictor
		
		n_panel_rows <- 5 # number of rows of subpanels
		n_panel_cols <- 8 # number of columns of subpanels
		n_panels <- n_panel_rows * n_panel_cols
		n_actual_panels <- min(n_panels, ncol(abund))
		these_predictors <- if ('rhizobiome_or_bulk' %in% raw_predictors) {
			raw_predictors[raw_predictors != 'rhizobiome_or_bulk']
		} else {
			raw_predictors
		}
		
		# decide on taxa to plot... sampling regularly along abundance Gradient
		total_abunds <- colSums(abund, na.rm = TRUE)
		total_abunds <- sort(total_abunds)
		taxa_to_plot <- unique(round(seq(1, ncol(abund), length.out = n_panels)))

		if ('rhizobiome_or_bulk' %in% raw_predictors) rhizo_bulk <- ifelse(env$rhizobiome_or_bulk == 1, 'rhizo', 'bulk')
		for (pred in these_predictors) {
		
			abund_vs_predictors <- list()
			# for (i in seq_len(n_actual_panels)) {
			for (i in seq_along(taxa_to_plot)) {
			
				taxon <- colnames(abund)[taxa_to_plot[i]]
				x <- env[[pred]]
				data <- data.frame(x = x, abundance = abund[[i]] + 1)
				if ('rhizobiome_or_bulk' %in% raw_predictors) data$rhizo_bulk <- rhizo_bulk
				
				in_bulk_indicator <- taxa$in_bulk[taxa$taxon == taxon] & !taxa$in_rhizobiome[taxa$taxon == taxon]
				in_bulk_color <- ifelse(in_bulk_indicator, 'lightgoldenrodyellow', 'lightgreen')

					if ('rhizobiome_or_bulk' %in% raw_predictors) {
						
						plot <- ggplot(data, aes(x = x, y = abundance, fill = rhizo_bulk)) +
							geom_point(pch = 21) +
							scale_fill_manual(
								guide = 'none',
								values = c('rhizo' = 'lightgreen', 'bulk' = 'lightgoldenrodyellow')
							)

					} else {
						plot <- ggplot(data, aes(x = x, y = abundance)) +
							geom_point()
					}

				abund_vs_predictors[[i]] <- plot +
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
			ggsave(abund_vs_predictors, filename = paste0(out_dir, '/abundance_vs_', pred, '.png'), width = 12, height = 8, dpi = 200)
		
		}
		
		### same type of graphs for PC axes
		if (use_pc_axes_as_predictors) {
			
			these_predictors <- colnames(sites_by_env)
			these_predictors <- if ('rhizobiome_or_bulk' %in% raw_predictors) {
				these_predictors[these_predictors != 'rhizobiome_or_bulk']
			} else {
				these_predictors
			}

			for (pred in these_predictors) {
			
				abund_vs_predictors <- list()
				for (i in seq_len(n_actual_panels)) {
				
					taxon <- colnames(abund)[i]
					x <- sites_by_env[[pred]]
					data <- data.frame(x = x, abundance = abund[[i]] + 1)
					if ('rhizobiome_or_bulk' %in% raw_predictors) data[ , rhizo_bulk := rhizo_bulk]
					
					in_bulk_indicator <- taxa$in_bulk[taxa$taxon == taxon]
					in_bulk_color <- ifelse(in_bulk_indicator, 'lightgoldenrodyellow', 'lightgreen')
					
					if ('rhizobiome_or_bulk' %in% raw_predictors) {
						
						plot <- ggplot(data, aes(x = x, y = abundance, fill = rhizo_bulk))
							geom_point(pch = 21) +
							scale_fill_manual(
								guide = 'none',
								values = c('rhizo' = 'lightgreen', 'bulk' = 'lightgoldenrodyellow')
							)

					} else {
						plot <- ggplot(data, aes(x = x, y = abundance)) +
							geom_point()
					}
					
					abund_vs_predictors[[i]] <- plot +
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
				ggsave(abund_vs_predictors, filename = paste0(out_dir, '/abundance_vs_', pred, '.png'), width = 12, height = 8, dpi = 200)
			
			}

		}

	}

#######################
### construct model ###
#######################

	say('construct model')

	# arguments to Hmsc (always included)
	args <- list(
		Y = abund,
		XFormula = x_formula,
		XData = as.data.frame(sites_by_env),
		XScale = TRUE,
		distr = response,
		studyDesign = study_design_coerced,
		ranLevels = list(site = site_effect, id = plant_effect),
		YScale = TRUE # only has effect if response = 'normal'
	)

	# phylogeny arguments
	if (include_phylogeny) args$phyloTree <- taxonomic_tree

	model <- do.call(Hmsc, args)
	sampleMcmc(model, samples = 2) # does this yield errors?

	# # obviates(?) error (from https://www.helsinki.fi/assets/drupal/2022-08/hmscdevelopment.pdf)
	# effect of using this is unknown
	# sampleMcmc(model, samples = 2, updater = list(GammaEta = FALSE)) 

#################
### run model ###
#################

	say('MCMC')

	# if (is.null(n_parallel)) n_parallel <- n_chains # setting n_parallel > 1 causes sampleHmsc() to fail to sample!

	tic()
	fit <- sampleMcmc(
		hM = model, samples = samples, thin = thin,
		transient = ifelse(is.null(transient), ceiling(0.5 * samples * thin), transient),
		# adaptNf = rep(ceiling(0.4 * samples * thin), model$nr),
		# adaptNf = rep(ceiling(0.4 * samples * thin), 100),
		nChains = n_chains,
		# nParallel = n_parallel, # NB n_parallel >1 leads to failure to sample!
		nParallel = 1, # NB n_parallel >1 leads to failure to sample!
		initPar = 'fixed effects', # use MLE to generate initial guesses... much faster!
		# verbose = TRUE
		verbose = TRUE
	)
	toc()

	saveRDS(fit, file = paste0(out_dir, '/fit_model.rds'))

##########################
### assess convergence ###
##########################

	say('assess convergence')

	posteriors <- convertToCodaObject(fit, spNamesNumbers = c(TRUE, FALSE), covNamesNumbers = c(TRUE, FALSE))
	nr <- fit$nr

	### Gelman-Rubin diagnostics (R-hat)
	# We want all R-hat + upper CI values to be <=1.1

	sink(paste0(out_dir, '/model_convergence.txt'), split = TRUE)
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
		
	ggsave(rhat_plot, file = paste0(out_dir, '/rhat_betas.png'), width = 12, height = 8, dpi = 150)
	
	say('Gamma:', pre = 1)
	rhat_gamma <- gelman.diag(posteriors$Gamma, multivariate = FALSE)$psrf
	say('Mean R-hat: ', mean(rhat_gamma[ , 1]))
	say('Maximum R-hat: ', max(rhat_gamma[ , 1]))
	say('Maximum R-hat + upper CI: ', max(rowSums(rhat_gamma)))
	
	rhat_gamma <- as.data.frame(rhat_gamma)
	names(rhat_gamma)[1] <- 'estimate'
	rhat_plot <- ggplot(rhat_gamma, aes(x = estimate)) +
		geom_histogram(binwidth = 0.5, fill = '#69b3a2', color = '#e9ecef') +
		geom_vline(xintercept = 1.1) +
		ggtitle('Gammas')

	ggsave(rhat_plot, file = paste0(out_dir, '/rhat_gamma.png'), width = 12, height = 8, dpi = 150)

	if (include_phylogeny) {

		say('Rho (phylogeny):', pre = 1)
		rhat_rho <- gelman.diag(posteriors$Rho, multivariate = FALSE)$psrf
		say('Mean R-hat: ', mean(rhat_rho[ , 1]))
		say('Maximum R-hat: ', max(rhat_rho[ , 1]))
		say('Maximum R-hat + upper CI: ', max(rowSums(rhat_rho)))
		
		rhat_rho <- as.data.frame(rhat_rho)
		names(rhat_rho)[1] <- 'estimate'
		rhat_plot <- ggplot(rhat_rho, aes(x = estimate)) +
			geom_histogram(binwidth = 0.5, fill = '#69b3a2', color = '#e9ecef') +
			geom_vline(xintercept = 1.1) +
			ggtitle('Rhos')
		ggsave(rhat_plot, file = paste0(out_dir, '/rhat_rho.png'), width = 12, height = 8, dpi = 150)

	}

	if (nr > 0) {
		
		say('Omega', pre = 1)
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
	
		ggsave(rhat_plot, file = paste0(out_dir, '/rhat_omega.png'), width = 12, height = 8, dpi = 150)

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
					# ggsave(rhat_plot, file = paste0('./outputs_sonny/', out_dir, '/rhat_alpha_', j, '.png'), width = 12, height = 8, dpi = 150)

		# 		}
		# 	}
		# }

	}
	
	sink()

########################
### assess model fit ###
########################

	say('assess model fit')

	# level at which to conduct cross-validation
	# NB we need to use at least a site-level random effect
	cv.level <- 'site'

	# predictions without/with cross-validation
	partition <- createPartition(fit, nfolds = n_folds, column = cv.level)
	preds <- computePredictedValues(fit, verbose = FALSE)
	preds_cv <- computePredictedValues(fit, partition = partition, verbose = FALSE)

	model_fit <- evaluateModelFit(hM = fit, predY = preds)
	model_fit_cv <- evaluateModelFit(hM = fit, predY = preds_cv)
	waic <- computeWAIC(fit)

	sink(paste0(out_dir, '/waic.txt'), split = TRUE)
	say(date())
	say('WAIC: ', waic)
	sink()
	
	# response <- 'lognormal poisson' # response is integers >= 0
	# response <- 'poisson' # response is integers >= 0
	response <- 'normal' # response is log of abundance, with 0 values forced to NA

	metrics <- if (response %in% c('lognormal', 'poisson')) {
		c('RMSE', 'O.RMSE', 'C.RMSE', 'SR2', 'C.SR2')
	} else if (response == 'normal') {
		c('RMSE', 'SR2', 'O.AUC', 'O.TjurR2', 'O.RMSE', 'C.SR2', 'C.RMSE')
	}

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

	ggsave(all_fit_plots, file = paste0(out_dir, '/model_fit.png'), width = 14, height = 10, dpi = 300)

###########################
### parameter estimates ###
###########################

	say('parameter estimates')

	# Alpha is spatial scale parameter (pages 100-101)

	posteriors <- convertToCodaObject(fit, spNamesNumbers = c(TRUE, FALSE), covNamesNumbers = c(TRUE, FALSE))
	
	### Betas
	#########
	samples <- posteriors$Beta
	
	# make one caterpillar plot per model term
	model_terms <- attr(terms(model$XFormula), 'term.labels')
	model_terms <- model_terms[model_terms != 'rhizobiome_or_bulk']
	for (model_term in model_terms) {

		say(model_term)

		pars <- paste0('B[', model_term, ', ', taxa$taxon, ']')
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

		ggsave(caterpillars, file = paste0(out_dir, '/coefficient_estimates_', file_name, '.png'), width = 10, height = 8)

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
	
	taxas <- taxa$taxon
	
	omegas <- expand.grid(taxon1 = taxas, taxon2 = taxas)
	omegas$omega <- NA_real_
	n=0
	for (i in 1:nrow(omegas)) {
		
		taxon1 <- grepl(taxa1, pattern = omegas$taxon1[i])
		taxon2 <- grepl(taxa2, pattern = omegas$taxon2[i])
		this_col <- taxon1 & taxon2
		
		# NB taking the mean of the value obviates cases when 2 values are drawn... but is this correct? why >1 value in *some* cases?
		# omegas$omega[i] <- means[this_col]
		omegas$omega[i] <- mean(means[this_col])
		if (length(means[this_col]) > 1) n <- n + 1
	
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
	
	taxas <- taxa$taxon
	
	omegas <- expand.grid(taxon1 = taxas, taxon2 = taxas)
	omegas$omega <- NA_real_
	
	for (i in 1:nrow(omegas)) {
		
		taxon1 <- grepl(taxa1, pattern = omegas$taxon1[i])
		taxon2 <- grepl(taxa2, pattern = omegas$taxon2[i])
		this_col <- taxon1 & taxon2

		# NB taking the mean of the value obviates cases when 2 values are drawn... but is this correct? why >1 value in *some* cases?
		# omegas$omega[i] <- means[this_col]
		omegas$omega[i] <- mean(means[this_col])
	
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

	ggsave(omegas_plot, file = paste0(out_dir, '/omegas.png'), width = 18, height = 9, dpi = 600)

##############################
### predictions to present ###
##############################

	say('predictions to present')

	### save data frames of predictions for the mean predicted abundance, SD, and (optionally, 10th and 90th quantiles)

	### obtain rasters of climate and soil predictors
	
	# environmental rasters
	rasts_sq <- rast('./outputs_sonny/climatena_1961_2020_soilgrids_aggregated_8x.tif')

	names(rasts_sq)[names(rasts_sq) == 'ph_soilgrids'] <- 'ph_field'
	names(rasts_sq)[names(rasts_sq) == 'sand_soilgrids'] <- 'sand_field'
	names(rasts_sq)[names(rasts_sq) == 'silt_soilgrids'] <- 'silt_field'
	names(rasts_sq)[names(rasts_sq) == 'clay_soilgrids'] <- 'clay_field'
	names(rasts_sq)[names(rasts_sq) == 'soc_soilgrids_perc'] <- 'soc_field_perc'
	names(rasts_sq)[names(rasts_sq) == 'nitrogen_soilgrids_perc'] <- 'nitrogen_field'

	rast_predictors <- raw_predictors[raw_predictors %in% names(rasts_sq)]
	rasts_sq <- rasts_sq[[rast_predictors]]

	### sampling predictors
	sampling_ppt_mm <- median(env$sampling_ppt_mm)
	sampling_tmean_c <- median(env$sampling_tmean_c)

	sampling_ppt_mm_rast <- rasts_sq[[1]]
	sampling_ppt_mm_rast[] <- sampling_ppt_mm

	sampling_tmean_c_rast <- rasts_sq[[1]]
	sampling_tmean_c_rast[] <- sampling_tmean_c

	names(sampling_ppt_mm_rast) <- 'sampling_ppt_mm'
	names(sampling_tmean_c_rast) <- 'sampling_tmean_c'

	rasts_sq <- c(rasts_sq, sampling_ppt_mm_rast, sampling_tmean_c_rast)

	rast_env <- as.data.table(rasts_sq, cell = TRUE)
	rast_env <- rast_env[complete.cases(rast_env)]

	if (use_pc_axes_as_predictors) {

		cell <- rast_env[ , 'cell']
		rast_env <- predict(pca, rast_env)
		colnames(rast_env) <- tolower(colnames(rast_env))
		rast_env <- as.data.table(rast_env)
		rast_env[ , cell := cell]
		
	}

	rast_env[ , rhizobiome_or_bulk := 1]
	rast_env[ , rhizobiome_or_bulk := factor(rhizobiome_or_bulk, levels = 0:1)]

	# create data tables for predictions	
	predictions_mean <- predictions_sd <- rast_env[ , 'cell', drop = FALSE]

	### create one column per taxon in each data table for predictions
	for (i in 1:nrow(taxa)) {

		predictions_mean[ , DUMMY := NA_real_]
		predictions_sd[ , DUMMY := NA_real_]

		taxon <- as.character(taxa$taxon[i])
		colnames(predictions_mean)[ncol(predictions_mean)] <- taxon
		colnames(predictions_sd)[ncol(predictions_sd)] <- taxon

	}

	predictors <- if (use_pc_axes_as_predictors) {
		pc_predictors
	} else {
		raw_predictors
	}

	### make predictions in chunks to save on memory (HMSC crashes otherwise)
	taxas <- as.character(taxa$taxon)
	size <- round(nrow(rast_env) / 20)
	sets <- ceiling(nrow(rast_env) / size)
	for (set in seq_len(sets)) {

		# say(set, ' of ', sets, ' ', date())

		indices <- (1 + size * (set - 1)):min(nrow(rast_env), size * set)

		# select environmental data
		x_data <- rast_env[indices, ..predictors]
		x_data <- as.data.frame(x_data)
		n_cells <- nrow(x_data)

		# random effects
		site_code <- study_design$location
		site_code <- as.integer(as.factor(site_code))
		site_code <- unique(site_code[study_design$location == 'KS1'])

		pred_site_effect <- data.frame(site = rep(site_code, n_cells)) # assuming all places are Kansas
		pred_id_effect <- data.frame(id = rep(1, n_cells)) # assuming all samples are #40
		s_new <- list(site_effect = pred_site_effect, id_effect = pred_id_effect)

		# predict
		# creates a list with one element per posterior iteration, rows are cells of raster and columns are taxa
		grad <- prepareGradient(fit, XDataNew = x_data, sDataNew = s_new) # prepareGradient() is for new data
		preds <- predict(fit, Gradient = grad, expected = TRUE)

		# calculate mean, SD, quantiles across all MCMC samples
		preds <- abind(preds, along = 3)
		preds_mean <- apply(preds, c(1, 2), mean)
		preds_sd <- apply(preds, c(1, 2), sd)

		preds_mean <- as.data.table(preds_mean)
		preds_sd <- as.data.table(preds_sd)

		# remember predictions
		predictions_mean[indices, (taxas) := preds_mean]
		predictions_sd[indices, (taxas) := preds_sd]

	}

	saveRDS(predictions_mean, paste0(out_dir, '/predictions_to_raster_cells_1961_2020_mean.rds'))
	saveRDS(predictions_sd, paste0(out_dir, '/predictions_to_raster_cells_1961_2020_sd.rds'))

#############################
### predictions to future ###
#############################

	say('predictions to future')

	### save data frames of predictions for the mean predicted abundance, SD, and (optionally, 10th and 90th quantiles)

	futs <- c(
		'ssp245_2041_2070',
		'ssp245_2071_2100',
		'ssp370_2041_2070',
		'ssp370_2071_2100'
	)

	for (fut in futs) {

		### obtain rasters of climate and soil predictors
		
		# environmental rasters
		rasts <- rast(paste0('./outputs_sonny/climatena_ensemble_8GCMs_', fut, '_soilgrids_aggregated_8x.tif'))

		names(rasts)[names(rasts) == 'ph_soilgrids'] <- 'ph_field'
		names(rasts)[names(rasts) == 'sand_soilgrids'] <- 'sand_field'
		names(rasts)[names(rasts) == 'silt_soilgrids'] <- 'silt_field'
		names(rasts)[names(rasts) == 'clay_soilgrids'] <- 'clay_field'
		names(rasts)[names(rasts) == 'soc_soilgrids_perc'] <- 'soc_field_perc'
		names(rasts)[names(rasts) == 'nitrogen_soilgrids_perc'] <- 'nitrogen_field'

		rast_predictors <- raw_predictors[raw_predictors %in% names(rasts)]
		rasts <- rasts[[rast_predictors]]

		### sampling predictors
		sampling_ppt_mm <- median(env$sampling_ppt_mm)
		sampling_tmean_c <- median(env$sampling_tmean_c)

		sampling_ppt_mm_rast <- rasts[[1]]
		sampling_ppt_mm_rast[] <- sampling_ppt_mm

		sampling_tmean_c_rast <- rasts[[1]]
		sampling_tmean_c_rast[] <- sampling_tmean_c

		names(sampling_ppt_mm_rast) <- 'sampling_ppt_mm'
		names(sampling_tmean_c_rast) <- 'sampling_tmean_c'

		rasts <- c(rasts, sampling_ppt_mm_rast, sampling_tmean_c_rast)

		rast_env <- as.data.table(rasts, cell = TRUE)
		rast_env <- rast_env[complete.cases(rast_env)]

		if (use_pc_axes_as_predictors) {

			cell <- rast_env[ , 'cell']
			rast_env <- predict(pca, rast_env)
			colnames(rast_env) <- tolower(colnames(rast_env))
			rast_env <- as.data.table(rast_env)
			rast_env[ , cell := cell]
			
		}

		rast_env[ , rhizobiome_or_bulk := 1]
		rast_env[ , rhizobiome_or_bulk := factor(rhizobiome_or_bulk, levels = 0:1)]

		# create data tables for predictions	
		predictions_mean <- predictions_sd <- rast_env[ , 'cell', drop = FALSE]

		### create one column per taxon in each data table for predictions
		for (i in 1:nrow(taxa)) {

			predictions_mean[ , DUMMY := NA_real_]
			predictions_sd[ , DUMMY := NA_real_]

			taxon <- as.character(taxa$taxon[i])
			colnames(predictions_mean)[ncol(predictions_mean)] <- taxon
			colnames(predictions_sd)[ncol(predictions_sd)] <- taxon

		}

		predictors <- if (use_pc_axes_as_predictors) {
			pc_predictors
		} else {
			raw_predictors
		}

		### make predictions in chunks to save on memory (HMSC crashes otherwise)
		taxas <- as.character(taxa$taxon)
		size <- round(nrow(rast_env) / 20)
		sets <- ceiling(nrow(rast_env) / size)
		for (set in seq_len(sets)) {

			# say(set, ' of ', sets, ' ', date())

			indices <- (1 + size * (set - 1)):min(nrow(rast_env), size * set)

			# select environmental data
			x_data <- rast_env[indices, ..predictors]
			x_data <- as.data.frame(x_data)
			n_cells <- nrow(x_data)

			# random effects
			site_code <- study_design$location
			site_code <- as.integer(as.factor(site_code))
			site_code <- unique(site_code[study_design$location == 'KS1'])

			pred_site_effect <- data.frame(site = rep(site_code, n_cells)) # assuming all places are Kansas
			pred_id_effect <- data.frame(id = rep(1, n_cells)) # assuming all samples are #40
			s_new <- list(site_effect = pred_site_effect, id_effect = pred_id_effect)

			# predict
			# creates a list with one element per posterior iteration, rows are cells of raster and columns are taxa
			grad <- prepareGradient(fit, XDataNew = x_data, sDataNew = s_new) # prepareGradient() is for new data
			preds <- predict(fit, Gradient = grad, expected = TRUE)

			# calculate mean, SD, quantiles across all MCMC samples
			preds <- abind(preds, along = 3)
			preds_mean <- apply(preds, c(1, 2), mean)
			preds_sd <- apply(preds, c(1, 2), sd)

			preds_mean <- as.data.table(preds_mean)
			preds_sd <- as.data.table(preds_sd)

			# remember predictions
			predictions_mean[indices, (taxas) := preds_mean]
			predictions_sd[indices, (taxas) := preds_sd]

		}

		saveRDS(predictions_mean, paste0(out_dir, '/predictions_to_raster_cells_', fut, '_mean.rds'))
		saveRDS(predictions_sd, paste0(out_dir, '/predictions_to_raster_cells_', fut, '_sd.rds'))

	} # next future

#################################################
### create rasters of predictions to present  ###
#################################################

	say('create rasters of predictions to present')

	predictions_mean <- readRDS(paste0(out_dir, '/predictions_to_raster_cells_1961_2020.rds'))
	predictions_sd <-  readRDS(paste0(out_dir, '/predictions_to_raster_cells_1961_2020_sd.rds'))

	template <- rasts_sq[[1]]
	template[] <- NA_real_

	for (i in seq_along(taxa$taxon)) {

		taxon <- taxa$taxon[i]
		# say(taxon)

		### make maps of MEAN prediction
		# for infinite predictions to the maximum non-infinite value
		predictions_filtered <- predictions_mean[[taxon]]
		if (any(is.infinite(predictions_filtered))) {

			non_infinite_max <- max(predictions_filtered[!is.infinite(predictions_filtered)])
			predictions_filtered[is.infinite(predictions_filtered)] <- non_infinite_max

		}

		predictions_filtered <- log10(predictions_filtered)
		map_this_taxon_mean <- setValueByCell(template, val = predictions_filtered, cell = predictions_mean[['cell']])
		names(map_this_taxon_mean) <- taxon

		# ### make maps of SD of prediction
		# for infinite predictions to the maximum non-infinite value
		predictions_filtered <- predictions_sd[[taxon]]
		map_this_taxon_sd <- setValueByCell(template, val = predictions_filtered, cell = predictions_sd[['cell']])
		names(map_this_taxon_sd) <- taxon

		if (i == 1) {
			maps_mean <- map_this_taxon_mean
			maps_sd <- map_this_taxon_sd
		} else {
			maps_mean <- c(maps_mean, map_this_taxon_mean)
			maps_sd <- c(maps_sd, map_this_taxon_sd)
		}


	}

	maps_mean <- trim(maps_mean)
	maps_sd <- trim(maps_sd)

	writeRaster(maps_mean, paste0(out_dir, '/prediction_rasters_1961_2020_mean.tif'), overwrite = TRUE)
	writeRaster(maps_sd, paste0(out_dir, '/prediction_rasters_1961_2020_sd.tif'), overwrite = TRUE)

###############################################
### create rasters of predictions to future ###
###############################################

	say('create rasters of predictions to future')

	for (fut in futs) {

		template <- rasts_sq[[1]]
		template[] <- NA_real_

		predictions_mean <- readRDS(paste0(out_dir, '/predictions_to_raster_cells_', fut, '_mean.rds'))
		predictions_sd <- readRDS(paste0(out_dir, '/predictions_to_raster_cells_', fut, '_sd.rds'))

		for (i in seq_along(taxa$taxon)) {

			taxon <- taxa$taxon[i]
			# say(taxon)

			### make maps of MEAN prediction
			# for infinite predictions to the maximum non-infinite value
			predictions_filtered <- predictions_mean[[taxon]]
			if (any(is.infinite(predictions_filtered))) {

				non_infinite_max <- max(predictions_filtered[!is.infinite(predictions_filtered)])
				predictions_filtered[is.infinite(predictions_filtered)] <- non_infinite_max

			}

			predictions_filtered <- log10(predictions_filtered)
			map_this_taxon_mean <- setValueByCell(template, val = predictions_filtered, cell = predictions_mean[['cell']])
			names(map_this_taxon_mean) <- taxon

			# ### make maps of SD of prediction
			# for infinite predictions to the maximum non-infinite value
			predictions_filtered <- predictions_sd[[taxon]]
			map_this_taxon_sd <- setValueByCell(template, val = predictions_filtered, cell = predictions_sd[['cell']])
			names(map_this_taxon_sd) <- taxon

			if (i == 1) {
				maps_mean <- map_this_taxon_mean
				maps_sd <- map_this_taxon_sd
			} else {
				maps_mean <- c(maps_mean, map_this_taxon_mean)
				maps_sd <- c(maps_sd, map_this_taxon_sd)
			}


		}

		maps_mean <- trim(maps_mean)
		maps_sd <- trim(maps_sd)

		writeRaster(maps_mean, paste0(out_dir, '/prediction_rasters_', fut, '_mean.tif'), overwrite = TRUE)
		writeRaster(maps_sd, paste0(out_dir, '/prediction_rasters_', fut, '_sd.tif'), overwrite = TRUE)

	} # next future

##################################################################
### make maps of predicted abundance for each taxon to present ###
##################################################################

	say('make maps of predicted abundance for each taxon to present')

	# North American countries
	nam <- gadm(c('CAN', 'USA', 'MEX'), level = 1, path = 'C:/!scratch', resolution = 2)
	nam <- project(nam, rasts_sq)

	# add abundances
	sites_abundances <- cbind(site_xy, abund)
	sites_abundances <- vect(sites_abundances, geom = c('longitude', 'latitude'), crs = getCRS('WGS84'))
	sites_abundances <- project(sites_abundances, rasts_sq)

	# plotting extent
	maps_mean <- rast(paste0(out_dir, '/prediction_rasters_1961_2020_mean.tif'))

	extent <- c(-754404, 2156962, -2232910, 836303)
	extent <- ext(extent)
	extent <- as.polygons(extent, crs = crs(maps_mean))
	extent <- buffer(extent, 100000)
	maps_mean <- crop(maps_mean, extent)
	
	# ...but actually plot a smaller area (looks nicer)
	extent <- buffer(extent, -100000)
	extent <- ext(extent)
	extent_for_plot <- as.vector(extent)

	pdf(paste0(out_dir, '/prediction_maps_1961_2020_mean.pdf'), width = 8.5, height = 8)
	for (i in seq_along(taxa$taxon)) {
	
		taxon <- as.character(taxa$taxon[i])
		# say(taxon)

		# get map and project		
		this_map <- maps_mean[[taxon]]

		# to define colors, get min/max abundance across sites and the map
		this_abundance <- as.data.frame(sites_abundances[ , taxon])[ , 1, drop = TRUE]
		this_abundance <- log10(this_abundance + 1)
		map_min <- minmax(this_map)[1, ]
		map_max <- minmax(this_map)[2, ]
		min_abund <- min(map_min, this_abundance)
		
		# remove extreme abundances
		max_abund <- globalx(this_map, quantile, probs = 0.99, na.rm = TRUE)
		max_abund <- max(max_abund, this_abundance)
		this_map <- app(this_map, fun = function(x, max_abund) ifelse(x > max_abund, 0.999 * max_abund, x), max_abund = max_abund)

		if (is.infinite(min_abund)) {
		
			neg_inf_abund_masked <- app(this_map, fun = function(x) ifelse(is.infinite(x), max_abund, x))
			min_abund <- globalx(neg_inf_abund_masked, min)
			this_map <- app(this_map, fun = function(x, min_abund) ifelse(x < min_abund, min_abund, x), min_abund = min_abund)
			
		}

		abund_seq <- seq(min_abund, max_abund, length.out = 101)
		this_abundance_scaled <- (this_abundance - min_abund) / (max_abund - min_abund)
		this_abundance_scaled <- round(100 * this_abundance_scaled) + 1

		pallette <- viridis(n = 101)
		site_colors <- pallette[this_abundance_scaled]

		taxon_nice <- gsub(taxon, pattern = '_', replacement = ': ')
		title <- paste(taxon_nice, '(1961-2020)')

		par(mgp = c(3, 1, 0), mar = c(3, 4, 6, 4) + 0.1)
		plot(nam, col = 'gray30', ext = extent_for_plot, lwd = 0.1, main = title)
		# title(main = taxon_nice, sub = '1961-2020', xpd = NA)
		plot(this_map, range = c(min_abund, max_abund), add = TRUE)
		plot(nam, border = 'gray20', lwd = 0.1, add = TRUE)

		plot(sites_abundances, pch = 21, bg = site_colors, cex = 1.6, add = TRUE)
	
	}
	dev.off()

#######################################################################
### make delta maps of change in abundance for each taxon to future ###
#######################################################################

	say('make maps of change in predicted abundance for each taxon to future')

	# North American countries
	nam <- gadm(c('CAN', 'USA', 'MEX'), level = 1, path = 'C:/!scratch', resolution = 2)
	nam <- project(nam, rasts_sq)

	# add abundances
	sites_abundances <- cbind(site_xy, abund)
	sites_abundances <- vect(sites_abundances, geom = c('longitude', 'latitude'), crs = getCRS('WGS84'))
	sites_abundances <- project(sites_abundances, rasts_sq)

	extent <- c(-754404, 2156962, -2232910, 836303)
	extent <- ext(extent)
	extent <- as.polygons(extent, crs = crs(maps_mean))
	extent <- buffer(extent, 100000)
	
	# ...but actually plot a smaller area (looks nicer)
	extent <- buffer(extent, -100000)
	extent <- ext(extent)
	extent_for_plot <- as.vector(extent)

	for (fut in futs) {

		say(fut)

		# plotting extent
		maps_mean <- rast(paste0(out_dir, '/prediction_rasters_', fut, '_mean.tif'))
		# maps_mean <- crop(maps_mean, extent)

		pdf(paste0(out_dir, '/prediction_maps_', fut, '_mean.pdf'), width = 8.5, height = 8)
		for (i in seq_along(taxa$taxon)) {
		
			taxon <- as.character(taxa$taxon[i])
			# say(taxon)

			# get map
			this_map <- maps_mean[[taxon]]

			# to define colors, get min/max abundance across sites and the map
			this_abundance <- as.data.frame(sites_abundances[ , taxon])[ , 1, drop = TRUE]
			this_abundance <- log10(this_abundance + 1)
			map_min <- minmax(this_map)[1, ]
			map_max <- minmax(this_map)[2, ]
			min_abund <- min(map_min, this_abundance)
			
			# remove extreme abundances
			max_abund <- globalx(this_map, quantile, probs = 0.99, na.rm = TRUE)
			max_abund <- max(max_abund, this_abundance)
			this_map <- app(this_map, fun = function(x, max_abund) ifelse(x > max_abund, 0.999 * max_abund, x), max_abund = max_abund)

			if (is.infinite(min_abund)) {
			
				neg_inf_abund_masked <- app(this_map, fun = function(x) ifelse(is.infinite(x), max_abund, x))
				min_abund <- globalx(neg_inf_abund_masked, min)
				this_map <- app(this_map, fun = function(x, min_abund) ifelse(x < min_abund, min_abund, x), min_abund = min_abund)
				
			}

			abund_seq <- seq(min_abund, max_abund, length.out = 101)
			this_abundance_scaled <- (this_abundance - min_abund) / (max_abund - min_abund)
			this_abundance_scaled <- round(100 * this_abundance_scaled) + 1

			pallette <- viridis(n = 101)
			site_colors <- pallette[this_abundance_scaled]

			# title
			taxon_nice <- gsub(taxon, pattern = '_', replacement = ': ')
			fut_nice <- paste0('(SSP ', substr(fut, 4, 6), ': ', substr(fut, 8, 11), '-', substr(fut, 13, 16), ')')
			title <- paste(taxon_nice, fut_nice)

			par(mgp = c(3, 1, 0), mar = c(3, 4, 6, 4) + 0.1)
			plot(nam, col = 'gray30', ext = extent_for_plot, lwd = 0.1, main = title)
			plot(this_map, range = c(min_abund, max_abund), add = TRUE)
			plot(nam, border = 'gray20', lwd = 0.1, add = TRUE)

			plot(sites_abundances, pch = 21, bg = site_colors, cex = 1.6, add = TRUE)
		
		}
		dev.off()

	} # next future

say('DONE', deco = '!', level = 1)

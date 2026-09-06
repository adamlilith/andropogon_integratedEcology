### MODELING ANDROPOGON GERARDI DISTRIBUTION, PHENOTYPE, PHYSIOLOGY, GENOTYPE, and ASSOCIATED MICROBIAL COMMUNITIES
### Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2023-12
###
### This script constructs response plots for the integrated models of Andropogon gerardi occurrence, biomass, and non-biomass traits. It also constructs caterpillar plots comparing coefficients from univariate versus integrated models.
###
### source('C:/Kaji/R/andropogon_integratedEcology/sdm_pdm/sdm_pdm_99_postmodeling_anaysis.r')
###
### CONTENTS ###
### setup ###
### caterpillar plots of model coefficients comparing univariate versus integrated models ###
### pregenerate response curves for all facets vs all covariates ###
### response curves for supplement ###
### maps of psi and response curves plots for main text ###
### calculate correlations between facets within latent multivariate normal component ###
### 
#############
### setup ###
#############

	rm(list = ls())

	drive <- 'C:/Kaji/'

	setwd(paste0(drive, '/Research/Andropogon/Andropogon'))
	source(paste0(drive, '/R/andropogon_integratedEcology/sdm_pdm/sdm_pdm_00_shared_functions_and_variables.r'))

	library(viridis)

# say('#############################################################################################')
# say('### caterpillar plots of model coefficients comparing univariate versus integrated models ###')
# say('#############################################################################################')

# 		### tables of all models (except occurrence)
# 		top_models <- read_xlsx('./outputs_loretta/integrated_sdm_pdm/summary_of_ALL_top_models_manual_assessment.xlsx', sheet = 'summary_of_ALL_top_models')
# 		top_models <- as.data.table(top_models)

# 		### pre-compiled data
# 		nonbiomass_facet_names_short <- 'bla_can_hei'
# 		data_occs_counties_morph <- readRDS(paste0('./outputs_loretta/integrated_sdm_pdm/models_integrated/data_occs_counties_', nonbiomass_facet_names_short, '.rds'))
# 		data_nonbiomass_counties_morph <- readRDS(paste0('./outputs_loretta/integrated_sdm_pdm/models_integrated/data_nonbiomass_counties_', nonbiomass_facet_names_short, '.rds'))

# 	facets <- names(nonbiomass_facets)

# 	# form			formula for the facet, shared across models
# 	# filename_stub	text to put into filename relevant to the facet
# 	# chains_uni	MCMC chains for the univariate model (not subsetted or stacked)
# 	# chains_multi_morph MCMC chains for the integrated model with morphological traits (not subsetted or stacked) Leave as NULL if none.
# 	# chains_multi_phys	MCMC chains for the integrated model with physiological traits (not subsetted or stacked). Leave as NULL if none.
# 	# facet			label for the facet (R-friendly): 'occs', 'biomass', 'height', etc.
# 	# param_uni		basename in univariate model of the coefficient to plot
# 	# param_multi	basename in multivariate model of the coefficient to plot
# 	# title			title for entire plot
# 	plot_caterpillars <- function(facet, filename_stub, form, chains_uni, chains_multi_morph = NULL, chains_multi_phys = NULL, param_uni, param_multi, title) {

# 		terms <- terms(form)
# 		terms <- attr(terms, 'term.labels')
# 		terms <- toupper(terms)
# 		terms <- 'Intercept' %>% c(terms)
# 		terms <- gsub(terms, pattern = 'PH', replacement = 'pH')
# 		terms_order <- seq_along(terms)
# 		terms[terms == toupper('BIO1')] <- 'BIO01'
# 		terms[terms == toupper('bio12_log10p1')] <- expression(log[10](BIO12 + 1))
# 		terms[terms == toupper('I(BIO1^2)')] <- expression(BIO01^2)
# 		terms[terms == toupper('I(BIO12^2)')] <- expression(BIO12^2)
# 		terms[terms == toupper('I(BIO12_LOG10P1^2)')] <- expression(log[10](BIO12 + 1)^2)
# 		terms[terms == toupper('I(BIO15^2)')] <- expression(BIO15^2)
# 		betas <- paste0('beta_', facet, '[', 1:length(terms), ']')

# 		caters <- list() # stores all plots
# 		for (i in seq_along(terms)) {

# 			uni <- mc_subset(chains_uni, param_uni, j = i)
# 			if (!is.null(chains_multi_morph)) multi_morph <- mc_subset(chains_multi_morph, param_multi, j = i)
# 			if (!is.null(chains_multi_phys)) multi_phys <- mc_subset(chains_multi_phys, param_multi, j = i)

# 			uni <- mc_stack(uni)
# 			if (!is.null(chains_multi_morph)) multi_morph <- mc_stack(multi_morph)
# 			if (!is.null(chains_multi_phys)) multi_phys <- mc_stack(multi_phys)

# 			uni <- data.table(
# 				term = terms[i],
# 				model = 'Simple',
# 				estimate = c(uni)
# 			)

# 			if (!is.null(chains_multi_morph)) {
			
# 				multi_morph <- data.table(
# 					term = terms[i],
# 					model = 'Integrated Morph.',
# 					estimate = c(multi_morph)
# 				)
			
# 			} else { multi_morph <- data.table() }
			
# 			if (!is.null(chains_multi_phys)) {
			
# 				multi_phys <- data.table(
# 					term = terms[i],
# 					model = 'Integrated Phys.',
# 					estimate = c(multi_phys)
# 				)
			
# 			} else { multi_phys <- data.table() }
			
# 			estimates <- rbind(uni, multi_morph, multi_phys)

# 			caters[[length(caters) + 1]] <- ggplot(estimates, aes(x = model, y = estimate, color = model, fill = model)) +
# 				geom_violin(trim = TRUE, quantile.linetype = 'solid', quantiles = c(0.1, 0.5, 0.9)) +
# 				# scale_fill_manual(values = c('Simple' = 'khaki3', 'Integrated Morph.' = 'darkslategray3', 'Integrated Phys.' = 'indianred1')) +
# 				scale_color_manual(values = c('Simple' = 'khaki4', 'Integrated Morph.' = 'darkgreen', 'Integrated Phys.' = 'indianred4')) +
# 				scale_fill_manual(values = c('Simple' = 'khaki3', 'Integrated Morph.' = 'mediumseagreen', 'Integrated Phys.' = 'indianred1')) +
# 				ylab('Estimate') +
# 				labs(fill = 'Model') +
# 				ggtitle(terms[i]) +
# 				theme(
# 					axis.title.x = element_blank(),
# 					axis.text.x = element_text(angle = 45, hjust = 1),
# 					legend.position = 'none'
# 				)

# 		}

# 		if (length(betas) <= 4) {
# 			nrow <- 1
# 			height <- 4
# 		} else {
# 			nrow <- 2
# 			height <- 8
# 		}

# 		caters <- plot_grid(plotlist = caters, nrow = nrow, align = 'hv')
		
# 		plot_title <- ggdraw() +
# 			draw_label(title, x = 0,	hjust = 0) +
# 			theme(
# 				plot.margin = margin(0, 0, 0, 7)
# 			)
# 		caters <- plot_grid(plot_title, caters, ncol = 1, rel_heights = c(0.06, 1))

# 		ggsave(caters, filename = paste0('./outputs_loretta/integrated_sdm_pdm/coefficients_caterpillar_plots_', filename_stub, '_univar_vs_integrated.png'), width = 10, height = height, units = 'in', dpi = 600, bg = 'white')

# 		invisible(caters)

# 	} # EOF

# 	### occurrence coefficients
# 	###########################

# 		uni_dir <- './outputs_loretta/integrated_sdm_pdm/models_occurrence/[occs~hurdlepoisson(bio1^2_log(bio12)^2_bio15^2)]_[bias~logit(1)]'
# 		multi_morph_dir <- './outputs_loretta/integrated_sdm_pdm/models_integrated/non_centered_MVN([occs~hurdlepoisson_offset]_[biomass~hurdleln]_[bla_can_hei])_eta=2'
# 		multi_phys_dir <- './outputs_loretta/integrated_sdm_pdm/models_integrated/non_centered_MVN([occs~hurdlepoisson_offset]_[biomass~hurdleln]_[cnr_pho_spa_sto_tra])_eta=2'

# 		chains_uni <- readRDS(paste0(uni_dir, '/chains.rds'))
# 		chains_multi_morph <- readRDS(paste0(multi_morph_dir, '/chains.rds'))
# 		chains_multi_phys <- readRDS(paste0(multi_phys_dir, '/chains.rds'))

# 		niter <- nrow(chains_multi_morph$samples[[1]])
# 		keeps <- round(seq((niter / 2) + 1, niter, length.out = 1000))
# 		for (i in 1:mc_n_chains(chains_multi_morph)) {
# 			chains_multi_morph$samples[[i]] <- chains_multi_morph$samples[[i]][keeps, ]
# 		}

# 		niter <- nrow(chains_multi_phys$samples[[1]])
# 		keeps <- round(seq((niter / 2) + 1, niter, length.out = 1000))
# 		for (i in 1:mc_n_chains(chains_multi_phys)) {
# 			chains_multi_phys$samples[[i]] <- chains_multi_phys$samples[[i]][keeps, ]
# 		}

# 		uni_formulae <- readRDS(paste0(uni_dir, '/formulae.rds'))
# 		multi_morph_formulae <- readRDS(paste0(multi_morph_dir, '/formulae.rds'))
# 		multi_phys_formulae <- readRDS(paste0(multi_phys_dir, '/formulae.rds'))

# 		stopifnot(uni_formulae$formula_occs == multi_morph_formulae$formula_occs & uni_formulae$formula_occs == multi_phys_formulae$formula_occs)

# 		form <- uni_formulae$formula_occs

# 		caters_occs <- plot_caterpillars(facet = 'occs', filename_stub = 'occs', form = form, chains_uni = chains_uni, chains_multi_morph = chains_multi_morph, chains_multi_phys = chains_multi_phys, param_uni = 'beta_occs', param_multi = 'beta_occs', title = 'Abundance')

# 		caters_occs_psi <- plot_caterpillars(facet = 'psi', filename_stub = 'psi_occs', form = form, chains_uni = chains_uni, chains_multi_morph = chains_multi_morph, chains_multi_phys = chains_multi_phys, param_uni = 'beta_psi', param_multi = 'beta_psi', title = 'Abundance psi')

# 	### biomass coefficients
# 	########################

# 		uni_dir <- './outputs_loretta/integrated_sdm_pdm/models_biomass/[biomass~hurdleln_log(bio12)]'
# 		multi_morph_dir <- './outputs_loretta/integrated_sdm_pdm/models_integrated/non_centered_MVN([occs~hurdlepoisson_offset]_[biomass~hurdleln]_[bla_can_hei])_eta=2'
# 		multi_phys_dir <- './outputs_loretta/integrated_sdm_pdm/models_integrated/non_centered_MVN([occs~hurdlepoisson_offset]_[biomass~hurdleln]_[cnr_pho_spa_sto_tra])_eta=2'

# 		chains_uni <- readRDS(paste0(uni_dir, '/chains.rds'))
# 		chains_multi_morph <- readRDS(paste0(multi_morph_dir, '/chains.rds'))
# 		chains_multi_phys <- readRDS(paste0(multi_phys_dir, '/chains.rds'))

# 		niter <- nrow(chains_multi_morph$samples[[1]])
# 		keeps <- round(seq((niter / 2) + 1, niter, length.out = 1000))
# 		for (i in 1:mc_n_chains(chains_multi_morph)) {
# 			chains_multi_morph$samples[[i]] <- chains_multi_morph$samples[[i]][keeps, ]
# 		}

# 		niter <- nrow(chains_multi_phys$samples[[1]])
# 		keeps <- round(seq((niter / 2) + 1, niter, length.out = 1000))
# 		for (i in 1:mc_n_chains(chains_multi_phys)) {
# 			chains_multi_phys$samples[[i]] <- chains_multi_phys$samples[[i]][keeps, ]
# 		}

# 		uni_meta <- readRDS(paste0(uni_dir, '/!meta_biomass.rds'))
# 		uni_formulae <- uni_meta$formulae$formula_biomass
# 		multi_morph_formulae <- readRDS(paste0(multi_morph_dir, '/formulae.rds'))
# 		multi_phys_formulae <- readRDS(paste0(multi_phys_dir, '/formulae.rds'))

# 		stopifnot(uni_formulae == multi_morph_formulae$formula_biomass & uni_formulae == multi_phys_formulae$formula_biomass)

# 		form <- uni_formulae

# 		caters_biomass <- plot_caterpillars(facet = 'biomass', filename_stub = 'biomass', form = form, chains_uni = chains_uni, chains_multi_morph = chains_multi_morph, chains_multi_phys = chains_multi_phys, param_uni = 'beta_biomass', param_multi = 'beta_biomass', title = 'Biomass')

# 		caters_psi <- plot_caterpillars(facet = 'psi', filename_stub = 'psi_biomass', form = form, chains_uni = chains_uni, chains_multi_morph = chains_multi_morph, chains_multi_phys = chains_multi_phys, param_uni = 'beta_psi', param_multi = 'beta_psi', title = 'Biomass psi')

# 	### each non-biomass facet
# 	##########################

# 		for (i in seq_along(facets)) {

# 			# get selected model
# 		  	facet <- facets[i]
# 			say(facet)
# 			index <- which(top_models$facet == facet & top_models$selected)
# 			best <- top_models[index, ]

# 			# meta-data on this model
# 		  	resp_distrib <- best$resp_distrib
# 			filename_form <- nonbiomass_facets[[facet]]$filename
		  	
# 			# chains for each model
# 		  	univar_dir <- paste0('./outputs_loretta/integrated_sdm_pdm/models_', facet, '/[', facet, '~', tolower(resp_distrib), '(', filename_form, ')]')
		
# 			# predictions from multivariate model
# 			if (facet %in% c('blade_width', 'canopy_diameter', 'height')) {

# 				multi_dir <- paste0('./outputs_loretta/integrated_sdm_pdm/models_integrated/non_centered_MVN([occs~hurdlepoisson_offset]_[biomass~hurdleln]_[bla_can_hei])_eta=2')
				
# 			} else if (facet %in% c('cn_ratio', 'photosynthetic_rate', 'spad', 'stomatal_conductance', 'transpiration_rate')) {
				
# 				multi_dir <- paste0('./outputs_loretta/integrated_sdm_pdm/models_integrated/non_centered_MVN([occs~hurdlepoisson_offset]_[biomass~hurdleln]_[cnr_pho_spa_sto_tra])')

# 			}

# 			chains_uni <- readRDS(paste0(univar_dir, '/chains.rds'))
# 			chains_multi <- readRDS(paste0(multi_dir, '/chains.rds'))

# 			niter <- nrow(chains_multi$samples[[1]])
# 			keeps <- round(seq((niter / 2) + 1, niter, length.out = 1000))
# 			for (i in 1:mc_n_chains(chains_multi)) {
# 				chains_multi$samples[[i]] <- chains_multi$samples[[i]][keeps, ]
# 			}

# 			uni_formulae <- readRDS(paste0(univar_dir, '/formulae.rds'))
# 			uni_formulae <- uni_formulae$formula_facet

# 			# check that formulae match between models
# 			multi_formulae <- readRDS(paste0(multi_dir, '/formulae.rds'))
# 			stopifnot(uni_formulae == multi_formulae$nonbiomass_facets[[facet]]$formula)
# 			form <- uni_formulae

# 			# plotting data on this facet
# 			title <- get_nice_trait(facet)$short
# 			facet_codes <- fread(paste0(multi_dir, '/facet_codes.csv'))
# 			facet_num <- facet_codes$index[facet_codes$facet == facet] - 2
# 			param_multi <- paste0('beta_facet_', facet_num)

# 			# predictions from multivariate model
# 			if (facet %in% c('blade_width', 'canopy_diameter', 'height')) {

# 				caters_facet <- plot_caterpillars(facet = facet, filename_stub = facet, form = form, chains_uni = chains_uni, chains_multi_morph = chains_multi, param_uni = 'beta_facet', param_multi = param_multi, title = title)
				
# 				caters_facet_psi <- plot_caterpillars(facet = 'psi', filename_stub = paste0('psi_', facet), form = form, chains_uni = chains_uni, chains_multi_morph = chains_multi, param_uni = 'beta_psi', param_multi = 'beta_psi', title = paste0(get_nice_trait(facet)$short, ' psi'))

# 			} else if (facet %in% c('cn_ratio', 'photosynthetic_rate', 'spad', 'stomatal_conductance', 'transpiration_rate')) {
				
# 				caters_facet <- plot_caterpillars(facet = facet, filename_stub = facet, form = form, chains_uni = chains_uni, chains_multi_phys = chains_multi, param_uni = 'beta_facet', param_multi = param_multi, title = title)

# 				caters_facet_psi <- plot_caterpillars(facet = 'psi', filename_stub = paste0('psi_', facet), form = form, chains_uni = chains_uni, chains_multi_phys = chains_multi, param_uni = 'beta_psi', param_multi = param_multi, title = paste0(get_nice_trait(facet)$short, ' psi'))

# 			}
			
# 		} # next non-biomass facet

say('####################################################################')
say('### pregenerate response curves for all facets vs all covariates ###')
say('####################################################################')

	### user-defined
	################

	set.seed(1)

	n_sites <- 160 # number of sites along which to predict
	say('Using ', n_sites, ' sites for gradients!!!', level = 2, deco = '!')

	# use these quantiles to define range of BG sites and of occurrences
	bg_probs <- c(0.005, 0.995)
	occs_probs <- c(0.01, 0.99)

	user_input <- tolower(trimws(readline('Use abbreviated chains? (y/n): ')))
	if (user_input == 'y') {
		abbreviate <- TRUE # use just 1 chain for speed
		if (abbreviate) say('Using abbreviated chains!!!', level = 2, deco = '!')
	} else {
		abbreviate <- FALSE
	}

	# folders with morphological and physiological integrated models
	dir_multi_morph <- './outputs_loretta/integrated_sdm_pdm/models_integrated/non_centered_MVN([occs~hurdlepoisson_offset]_[biomass~hurdleln]_[bla_can_hei])_eta=2'

	dir_multi_phys <- './outputs_loretta/integrated_sdm_pdm/models_integrated/non_centered_MVN([occs~hurdlepoisson_offset]_[biomass~hurdleln]_[cnr_pho_spa_sto_tra])_eta=2'

	dir_occs_uni <- './outputs_loretta/integrated_sdm_pdm/models_occurrence/[occs~overdispersed_hurdlepoisson(bio1^2_log(bio12)^2_bio15^2)]_[bias~poaceae_offset_plus_exp_neg_1_WITH_zero_correction_covariate]'

	dir_biomass <- './outputs_loretta/integrated_sdm_pdm/models_biomass/[biomass~hurdleln_log(bio12)]'

	### centers and scales
	######################

	centers_scales <- fread('./outputs_loretta/integrated_sdm_pdm/centers_and_scales_for_covariates.csv')

	centers_counties <- centers_scales$center_occs
	scales_counties <- centers_scales$scale_occs

	centers_sites <- centers_scales$center_sites
	scales_sites <- centers_scales$scale_sites

	names(centers_counties) <- names(scales_counties) <- names(centers_sites) <- names(scales_sites) <- centers_scales$variable

	### tables of all models (except occurrence)
	############################################

	top_models <- read_xlsx('./outputs_loretta/integrated_sdm_pdm/summary_of_ALL_top_models_manual_assessment.xlsx', sheet = 'summary_of_ALL_top_models')
	top_models <- as.data.table(top_models)

	### pre-compiled data
	#####################

	say('pre-compiled data')

	nonbiomass_facet_names_short_morph <- 'bla_can_hei'
	nonbiomass_facet_names_short_phys <- 'cnr_pho_spa_sto_tra'
	
	data_occs_counties_morph <- readRDS(paste0('./outputs_loretta/integrated_sdm_pdm/models_integrated/data_occs_counties_', nonbiomass_facet_names_short_morph, '.rds'))

	data_biomass_counties_morph <- readRDS(paste0('./outputs_loretta/integrated_sdm_pdm/models_integrated/data_biomass_counties_', nonbiomass_facet_names_short_morph, '.rds'))
	
	data_biomass_sites_morph <- readRDS(paste0('./outputs_loretta/integrated_sdm_pdm/models_integrated/data_biomass_sites_', nonbiomass_facet_names_short_morph, '.rds'))
	
	data_nonbiomass_counties_morph <- readRDS(paste0('./outputs_loretta/integrated_sdm_pdm/models_integrated/data_nonbiomass_counties_', nonbiomass_facet_names_short_morph, '.rds'))
	
	data_nonbiomass_counties_phys <- readRDS(paste0('./outputs_loretta/integrated_sdm_pdm/models_integrated/data_nonbiomass_counties_', nonbiomass_facet_names_short_phys, '.rds'))

	data_nonbiomass_sites_morph <- readRDS(paste0('./outputs_loretta/integrated_sdm_pdm/models_integrated/data_nonbiomass_sites_', nonbiomass_facet_names_short_morph, '.rds'))

	data_nonbiomass_sites_phys <- readRDS(paste0('./outputs_loretta/integrated_sdm_pdm/models_integrated/data_nonbiomass_sites_', nonbiomass_facet_names_short_phys, '.rds'))

	### model metadata
	##################

	# formulae
	formulae_morph <- readRDS(paste0(dir_multi_morph, '/formulae.rds'))
	formulae_phys <- readRDS(paste0(dir_multi_phys, '/formulae.rds'))

	# facets for each integrated model
	facets_multi_morph <- fread(paste0(dir_multi_morph, '/facet_codes.csv'))
	facets_multi_phys <- fread(paste0(dir_multi_phys, '/facet_codes.csv'))

	facets_multi_morph <- facets_multi_morph$facet
	facets_multi_phys <- facets_multi_phys$facet

	facets_multi_morph <- facets_multi_morph[facets_multi_morph != 'occurrence' & facets_multi_morph != 'biomass']
	facets_multi_phys <- facets_multi_phys[facets_multi_phys != 'occurrence' & facets_multi_phys != 'biomass']

	### chains
	##########

	say('chains')

	### univariate occurrences chains
	chains_occs_uni <- readRDS(paste0(dir_occs_uni, '/chains.rds'))
	if (abbreviate) chains_occs_uni$samples <- chains_occs_uni$samples[[1]]

	chains_biomass_uni <- readRDS(file.path(dir_biomass, '/chains.rds'))
	if (abbreviate) chains_biomass_uni$samples <- chains_biomass_uni$samples[[1]]

	### integrated model chains

	# load chains, subset to fewest number of iterations done, then thin to 1000 each
	chains_multi_morph <- readRDS(paste0(dir_multi_morph, '/chains.rds'))
	chains_multi_phys <- readRDS(paste0(dir_multi_phys, '/chains.rds'))

	# make version of chains where correlation between facets is forced to 0
	chains_multi_morph_no_corr <- chains_multi_morph
	rows_cols <- expand.grid(row = 1:5, col = 1:5)
	nondiag <- rows_cols[[1]] != rows_cols[[2]]
	rows_cols <- rows_cols[nondiag, ]
	for (j in 1:nrow(rows_cols))  {
		u_star <- paste0('U_star[', rows_cols$row[j], ', ', rows_cols$col[j], ']')
		for (i in 1:4) {
			chains_multi_morph_no_corr$samples[[i]][ , u_star] <- 0
		}
	}

	chains_multi_phys_no_corr <- chains_multi_phys
	rows_cols <- expand.grid(row = 1:7, col = 1:7)
	nondiag <- rows_cols[[1]] != rows_cols[[2]]
	rows_cols <- rows_cols[nondiag, ]
	for (j in 1:nrow(rows_cols))  {
		u_star <- paste0('U_star[', rows_cols$row[j], ', ', rows_cols$col[j], ']')
		for (i in 1:4) {
			chains_multi_phys_no_corr$samples[[i]][ , u_star] <- 0
		}
	}

	if (abbreviate) chains_multi_morph$samples <- chains_multi_morph$samples[[1]]
	if (abbreviate) chains_multi_phys$samples <- chains_multi_phys$samples[[1]]

	### collate list of terms across facets
	form <- formulae_morph$formula_occs
	terms <- terms(form)
	terms <- attr(terms, 'term.labels')

	form <- formulae_morph$formula_biomass
	this_terms <- terms(form)
	this_terms <- attr(this_terms, 'term.labels')
	terms <- c(terms, this_terms)

	for (f in seq_along(nonbiomass_facets)) {

		form <- nonbiomass_facets[[f]]$formula
		this_terms <- terms(form)
		this_terms <- attr(this_terms, 'term.labels')
		terms <- c(terms, this_terms)

	}

	terms <- unique(terms)
	terms <- terms[terms != '1']
	terms <- terms[!grepl(terms, pattern = '\\^2')]
	terms <- terms[!grepl(terms, pattern = '\\:')]

	### create predictor matrices spanning the range of values of each variable across the species' known occupied range

	## get range of each variable
	nam_untrunc <- vect('./outputs_loretta/integrated_sdm_pdm/andropogon_gerardi_occurrences_with_environment_1961_2020_for_integration.gpkg')
	nam <- nam_untrunc[nam_untrunc$n_andropogon_gerardi > 0]

	nam_ssp370_2071_2100 <- vect('./data_from_adam_and_loretta/andropogon_gerardi_occurrences_with_environment_ensemble_8GCMs_ssp370_2071_2100_NEW.gpkg')

	nam$bio12_log10p1 <- log10(nam$bio12 + 1)

	nam_env <- as.data.table(nam)[ , ..terms]
	mins <- apply(nam_env, 2, min, na.rm = TRUE)
	maxs <- apply(nam_env, 2, max, na.rm = TRUE)

	ranges <- maxs - mins
	maxs <- maxs + 0.05 * ranges
	mins <- mins - 0.05 * ranges
	mins[names(mins) %in% paste0('bio', c(2, 3, 4, 12:19))] <- 0
	if (mins['insolation_2000_growing_season_kWh_per_m2'] < 0) mins['insolation_2000_growing_season_kWh_per_m2'] <- 0
	mins[names(mins) %in% 'bio12_log10p1'] <- 1E-6
	mins[names(mins) %in% c('sand', 'silt', 'clay')] <- 0
	mins[names(mins) %in% c('ph')] <- 3 # manual lower limit bc too low to show response
	maxs[names(maxs) %in% c('sand', 'silt', 'clay')] <- 1

	## for non-varying variables, use means across **sites**
	site_env <- data_biomass_sites_morph$site_data_raw[ , ..terms]
	means <- colMeans(site_env)

	## model matrix for occurrences and psi
	#######################################

		# make a list of model matrices for occurrences and psi
		# each MM allows one variable to increment across a gradient, while all others are held constant
		# also create a "fixed" version with just means (no variation across variables)

		say('model matrix for occurrences and psi')

		covars <- data_occs_counties_morph$covariates
		form <- formulae_morph$formula_occs

		fixed <- t(as.matrix(means))
		fixed <- fixed[ , covars, drop = FALSE]
		fixed <- fixed[rep(1, n_sites), , drop = FALSE]

		this_mins <- mins[covars]
		this_maxs <- maxs[covars]

		this_mins <- scale(matrix(this_mins, nrow = 1), center = centers_counties[covars], scale = scales_counties[covars])
		this_maxs <- scale(matrix(this_maxs, nrow = 1), center = centers_counties[covars], scale = scales_counties[covars])

		this_mins <- c(this_mins)
		this_maxs <- c(this_maxs)
		names(this_mins) <- covars
		names(this_maxs) <- covars

		fixed_scaled <- scale(fixed, center = centers_counties[covars], scale = scales_counties[covars])
		x_occs_fixed <- model.matrix(form, data = as.data.frame(fixed_scaled))
		x_occs <- list()
		for (i in seq_along(covars)) {

			covar <- covars[i]
			this_fixed_scaled <- fixed_scaled

			this_fixed_scaled[ , covar] <- seq(this_mins[covar], this_maxs[covar], length.out = n_sites)
			this_fixed_scaled <- as.data.frame(this_fixed_scaled)
			this_fixed_scaled <- model.matrix(form, data = this_fixed_scaled)

			x_occs[[covar]] <- this_fixed_scaled

		}
		names(x_occs) <- covars

	## model matrix for biomass
	###########################

		# make a list of model matrices for biomass
		# each MM allows one variable to increment across a gradient, while all others are held constant
		# also create a "fixed" version with just means (no variation across variables)

		say('model matrix for biomass')

		covars <- data_biomass_sites_morph$covariates
		form <- formulae_morph$formula_biomass

		this_mins <- mins[covars]
		this_maxs <- maxs[covars]

		this_mins <- scale(matrix(this_mins, nrow = 1), center = centers_sites[covars], scale = scales_sites[covars])
		this_maxs <- scale(matrix(this_maxs, nrow = 1), center = centers_sites[covars], scale = scales_sites[covars])

		this_mins <- c(this_mins)
		this_maxs <- c(this_maxs)
		names(this_mins) <- covars
		names(this_maxs) <- covars

		fixed <- t(as.matrix(means))
		fixed <- fixed[ , covars, drop = FALSE]
		fixed <- fixed[rep(1, n_sites), , drop = FALSE]

		fixed_scaled <- scale(fixed, center = centers_sites[covars], scale = scales_sites[covars])
		x_biomass_fixed <- model.matrix(form, data = as.data.frame(fixed_scaled))
		x_biomass <- list()
		for (i in seq_along(covars)) {

			covar <- covars[i]
			this_fixed_scaled <- fixed_scaled

			this_fixed_scaled[ , covar] <- seq(this_mins[covar], this_maxs[covar], length.out = n_sites)
			this_fixed_scaled <- as.data.frame(this_fixed_scaled)
			this_fixed_scaled <- model.matrix(form, data = this_fixed_scaled)

			x_biomass[[covar]] <- this_fixed_scaled

		}
		names(x_biomass) <- covars

	## model matrices for non-biomass facets
	########################################

		# make a list of model matrices for each non-biomass facet
		# the list has one element per facet, and each of these elements is itself list and contains one MM per variable with one variable incrementing across a gradient, while all others are held constant
		# also create list of a "fixed" MMs with just means (no variation across variables)

		say('model matrices for non-biomass facets')

		facets <- names(nonbiomass_facets)

		x_nonbiomass_fixed <- list()
		x_nonbiomass <- list()
		for (f in seq_along(facets)) {

			facet <- facets[f]

			if (facet %in% c('blade_width', 'canopy_diameter', 'height')) {

				covars <- data_nonbiomass_sites_morph[[facet]]$covariates
				form <- nonbiomass_facets[[facet]]$formula

			} else if (facet %in% c('cn_ratio', 'photosynthetic_rate', 'spad', 'stomatal_conductance', 'transpiration_rate')) {

				covars <- data_nonbiomass_sites_phys[[facet]]$covariates
				form <- nonbiomass_facets[[facet]]$formula

			}

			fixed <- t(as.matrix(means))
			fixed <- fixed[ , covars, drop = FALSE]
			fixed <- fixed[rep(1, n_sites), , drop = FALSE]

			fixed_scaled <- scale(fixed, center = centers_sites[covars], scale = scales_sites[covars])
			x_nonbiomass_fixed[[f]] <- model.matrix(form, data = as.data.frame(fixed_scaled))
			x_nonbiomass[[f]] <- list()
			for (i in seq_along(covars)) {

				covar <- covars[i]
				this_fixed_scaled <- fixed_scaled

				this_mins <- mins[covars]
				this_maxs <- maxs[covars]

				this_mins <- scale(matrix(this_mins, nrow = 1), center = centers_sites[covars], scale = scales_sites[covars])
				this_maxs <- scale(matrix(this_maxs, nrow = 1), center = centers_sites[covars], scale = scales_sites[covars])

				this_mins <- c(this_mins)
				this_maxs <- c(this_maxs)
				names(this_mins) <- covars
				names(this_maxs) <- covars

				this_fixed_scaled[ , covar] <- seq(this_mins[covar], this_maxs[covar], length.out = n_sites)
				this_fixed_scaled <- as.data.frame(this_fixed_scaled)
				this_fixed_scaled <- model.matrix(form, data = this_fixed_scaled)

				x_nonbiomass[[f]][[covar]] <- this_fixed_scaled

			}
			names(x_nonbiomass[[f]]) <- covars

		} # next non-biomass facet
		names(x_nonbiomass) <- facets
		names(x_nonbiomass_fixed) <- facets

	### for creating plots for main text
	resp_curves <- list() # list of lists, one sub-list per facet, each sublist has one panel per covariate

	### responses of ABUNDANCE
	###########################

	# biomass metadata
	meta <- readRDS(paste0(dir_biomass, '/!meta_biomass.rds'))
	resp_distrib_biomass <- meta$resp_distrib
	transform_biomass <- meta$transform

	covars <- data_occs_counties_morph$covariates
	responses <- ranges <- list() # to save pre-calculated response curves and ranges of BG, species, and sites
	for (covar in covars) {

		say('response curves for abundance vs ', covar)

		# values of x covariate (unscaled)
		x <- x_occs[[covar]][ , covar]
		x <- x * scales_counties[covar] + centers_counties[covar]
		if (covar == 'bio12_log10p1') x <- 10^x - 1

		# range of variable across baxckground
		this_covar <- if (covar == 'bio12_log10p1') { 'bio12' } else { covar }
		x_bg_empirical <- nam_untrunc[[this_covar]]
		x_bg_range_sq <- quantile(x_bg_empirical, probs = bg_probs, na.rm = TRUE)
		x_bg_empirical <- nam_ssp370_2071_2100[[this_covar]]
		x_bg_range_ssp370_2071_2100 <- quantile(x_bg_empirical, probs = bg_probs, na.rm = TRUE)

		# range of variable across known occurrences
		x_occs_empirical <- nam[nam$n_andropogon_gerardi > 0]
		this_covar <- if (covar == 'bio12_log10p1') { 'bio12' } else { covar }
		x_occs_empirical <- x_occs_empirical[[this_covar]]
		x_occs_range <- quantile(x_occs_empirical, probs = occs_probs, na.rm = TRUE)

		# range of variable across sites
		this_covar <- if (covar == 'bio12_log10p1') { 'bio12' } else { covar }
		x_sites_empirical <- unique(data_biomass_sites_morph$raw_data_biomass[[this_covar]])
		x_sites_range <- range(x_sites_empirical)

		# remember ranges of this covariate for later plotting
		ranges[[covar]] <- list(
			x_bg_range_sq = x_bg_range_sq,
			x_bg_range_ssp370_2071_2100 = x_bg_range_ssp370_2071_2100,
			x_occs_range = x_occs_range,
			x_sites_range = x_sites_range
		)

		# if biomass uses this covariate, use its incrementing model matrix... if not, used the fixed matrix
		x_biomass_varying <- x_biomass_fixed
		for (i in seq_along(x_biomass)) {

			if (any(names(x_biomass[[i]]) %in% covar)) {
				x_biomass_varying[[i]] <- x_biomass[[i]][[covar]]
			}

		}

		# if any non-biomass facets use this covariate, use its incrementing model matrix... if not, used the fixed matrix
		x_nonbiomass_varying <- x_nonbiomass_fixed
		for (i in seq_along(x_nonbiomass)) {

			if (any(names(x_nonbiomass[[i]]) %in% covar)) {
				x_nonbiomass_varying[[i]] <- x_nonbiomass[[i]][[covar]]
			}

		}

		# predict integrated model: for all facets that share this covariate, vary it
		# use integrated morphological model for these predictions
		say('   pred_multi_all_varying')
		pred_multi_all_varying <- predict_fully_integrated(
			chains = chains_multi_morph,
			nonbiomass_facets = nonbiomass_facets[c('blade_width', 'canopy_diameter', 'height')],
			resp_distrib_biomass = resp_distrib_biomass,
			transform_biomass = transform_biomass,
			x_occs = x_occs[[covar]],
			x_psi = x_occs[[covar]],
			x_biomass = x_biomass_varying,
			x_nonbiomass = x_nonbiomass_varying,
			w_occs = NULL,
			sampled = TRUE,
			type = 'samples',
			force_presence = FALSE
		)

		preds <- pred_multi_all_varying$preds_occs
		center_multi_all_varying <- apply(preds, 2, mean)
		upper_multi_all_varying <- apply(preds, 2, quantile, 0.9)
		lower_multi_all_varying <- apply(preds, 2, quantile, 0.1)

		# predict integrated model: for all facets that share this covariate, vary it
		# use integrated morphological model for these predictions
		say('   pred_multi_all_varying_force_presence')
		pred_multi_all_varying_force_presence <- predict_fully_integrated(
			chains = chains_multi_morph,
			nonbiomass_facets = nonbiomass_facets[c('blade_width', 'canopy_diameter', 'height')],
			resp_distrib_biomass = resp_distrib_biomass,
			transform_biomass = transform_biomass,
			x_occs = x_occs[[covar]],
			x_psi = x_occs[[covar]],
			x_biomass = x_biomass_varying,
			x_nonbiomass = x_nonbiomass_varying,
			w_occs = NULL,
			sampled = TRUE,
			type = 'samples',
			force_presence = TRUE
		)

		preds <- pred_multi_all_varying_force_presence$preds_occs
		center_multi_all_varying_force_presence <- apply(preds, 2, mean)
		upper_multi_all_varying_force_presence <- apply(preds, 2, quantile, 0.9)
		lower_multi_all_varying_force_presence <- apply(preds, 2, quantile, 0.1)

		# predict integrated model: for all facets that share this covariate, vary it but force facet correlatons to 0
		# use integrated morphological model for these predictions
		say('   pred_multi_all_varying_no_corr')
		pred_multi_all_varying_no_corr <- predict_fully_integrated(
			chains = chains_multi_morph_no_corr,
			nonbiomass_facets = nonbiomass_facets[c('blade_width', 'canopy_diameter', 'height')],
			resp_distrib_biomass = resp_distrib_biomass,
			transform_biomass = transform_biomass,
			x_occs = x_occs[[covar]],
			x_psi = x_occs[[covar]],
			x_biomass = x_biomass_varying,
			x_nonbiomass = x_nonbiomass_varying,
			w_occs = NULL,
			sampled = TRUE,
			type = 'samples',
			force_presence = FALSE
		)

		preds <- pred_multi_all_varying_no_corr$preds_occs
		center_multi_all_varying_no_corr <- apply(preds, 2, mean)
		upper_multi_all_varying_no_corr <- apply(preds, 2, quantile, 0.9)
		lower_multi_all_varying_no_corr <- apply(preds, 2, quantile, 0.1)

		# predict integrated model: only allow abundance to vary but force presence
		# use integrated morphological model for these predictions
		say('   pred_multi_abundance_varying')
		pred_multi_abundance_varying <- predict_fully_integrated(
			chains = chains_multi_morph,
			nonbiomass_facets = nonbiomass_facets[c('blade_width', 'canopy_diameter', 'height')],
			resp_distrib_biomass = resp_distrib_biomass,
			transform_biomass = transform_biomass,
			x_occs = x_occs[[covar]],
			x_psi = x_occs_fixed,
			x_biomass = x_biomass_fixed,
			x_nonbiomass = x_nonbiomass_fixed,
			w_occs = NULL,
			sampled = TRUE,
			type = 'samples',
			force_presence = FALSE
		)

		preds <- pred_multi_abundance_varying$preds_occs
		center_multi_abundance_varying <- apply(preds, 2, mean)
		upper_multi_abundance_varying <- apply(preds, 2, quantile, 0.9)
		lower_multi_abundance_varying <- apply(preds, 2, quantile, 0.1)

		# predict integrated model: allow covariate to vary only for occurrence and biomass
		# use integrated morphological model for these predictions
		say('   pred_multi_biomass_varying')
		pred_multi_biomass_varying <- predict_fully_integrated(
			chains = chains_multi_morph,
			nonbiomass_facets = nonbiomass_facets[c('blade_width', 'canopy_diameter', 'height')],
			resp_distrib_biomass = resp_distrib_biomass,
			transform_biomass =  transform_biomass,
			x_occs = x_occs[[covar]],
			x_psi = x_occs_fixed,
			x_biomass = x_biomass_varying,
			x_nonbiomass = x_nonbiomass_fixed,
			w_occs = NULL,
			sampled = TRUE,
			type = 'samples',
			force_presence = FALSE
		)

		preds <- pred_multi_biomass_varying$preds_occs
		center_multi_biomass_varying <- apply(preds, 2, mean)
		upper_multi_biomass_varying <- apply(preds, 2, quantile, 0.9)
		lower_multi_biomass_varying <- apply(preds, 2, quantile, 0.1)

		# predict integrated model: allow covariate to vary only for occurrence and psi
		# use integrated morphological model for these predictions
		say('   pred_multi_psi_varying')
		pred_multi_psi_varying <- predict_fully_integrated(
			chains = chains_multi_morph,
			nonbiomass_facets = nonbiomass_facets[c('blade_width', 'canopy_diameter', 'height')],
			resp_distrib_biomass = resp_distrib_biomass,
			transform_biomass = transform_biomass,
			x_occs = x_occs[[covar]],
			x_psi = x_occs[[covar]],
			x_biomass = x_biomass_fixed,
			x_nonbiomass = x_nonbiomass_fixed,
			w_occs = NULL,
			sampled = TRUE,
			type = 'samples',
			force_presence = FALSE
		)

		preds <- pred_multi_psi_varying$preds_occs
		center_multi_psi_varying <- apply(preds, 2, mean)
		upper_multi_psi_varying <- apply(preds, 2, quantile, 0.9)
		lower_multi_psi_varying <- apply(preds, 2, quantile, 0.1)

		# predict integrated model: allow covariate to vary only for occurrence and non-biomass facets
		# use integrated morphological model for these predictions
		say('   pred_multi_nonbiomass_facets_varying')
		pred_multi_nonbiomass_facets_varying <- predict_fully_integrated(
			chains = chains_multi_morph,
			nonbiomass_facets = nonbiomass_facets[c('blade_width', 'canopy_diameter', 'height')],
			resp_distrib_biomass = resp_distrib_biomass,
			transform_biomass = transform_biomass,
			x_occs = x_occs[[covar]],
			x_psi = x_occs_fixed,
			x_biomass = x_biomass_fixed,
			x_nonbiomass = x_nonbiomass_varying,
			w_occs = NULL,
			sampled = TRUE,
			type = 'samples',
			force_presence = FALSE
		)

		preds <- pred_multi_nonbiomass_facets_varying$preds_occs
		center_multi_nonbiomass_facets_varying <- apply(preds, 2, mean)
		upper_multi_nonbiomass_facets_varying <- apply(preds, 2, quantile, 0.9)
		lower_multi_nonbiomass_facets_varying <- apply(preds, 2, quantile, 0.1)

		# predict univariate model: allow covariate to vary for occurrences and psi
		say('   pred_uni_all_varying')
		pred_uni_all_varying <- predict_occs(
			chains = chains_occs_uni,
			x = x_occs[[covar]],
			x_psi = x_occs[[covar]],
			overdispersed = TRUE
		)

		center_uni_all_varying <- apply(pred_uni_all_varying, 2, mean)
		upper_uni_all_varying <- apply(pred_uni_all_varying, 2, quantile, 0.9)
		lower_uni_all_varying <- apply(pred_uni_all_varying, 2, quantile, 0.1)
	
		# predict univariate model: allow covariate to vary for occurrences and psi but force specis to be present
		say('   pred_uni_all_varying_force_presence')
		pred_uni_all_varying_force_presence <- predict_occs(
			chains = chains_occs_uni,
			x = x_occs[[covar]],
			x_psi = x_occs[[covar]],
			overdispersed = TRUE,
			force_presence = TRUE
		)

		center_uni_all_varying_force_presence <- apply(pred_uni_all_varying_force_presence, 2, mean)
		upper_uni_all_varying_force_presence <- apply(pred_uni_all_varying_force_presence, 2, quantile, 0.9)
		lower_uni_all_varying_force_presence <- apply(pred_uni_all_varying_force_presence, 2, quantile, 0.1)
	
		# predict univariate model: allow covariate to vary only for occurences, not psi
		say('   pred_uni_psi_constant')
		pred_uni_psi_constant <- predict_occs(
			chains = chains_occs_uni,
			x = x_occs[[covar]],
			x_psi = x_occs_fixed,
			overdispersed = TRUE
		)

		center_uni_psi_constant <- apply(pred_uni_psi_constant, 2, mean)
		upper_uni_psi_constant <- apply(pred_uni_psi_constant, 2, quantile, 0.9)
		lower_uni_psi_constant <- apply(pred_uni_psi_constant, 2, quantile, 0.1)

		# remember
		responses[[covar]] <- data.table(
			x = x,
			
			center_multi_all_varying = center_multi_all_varying,
			upper_multi_all_varying = upper_multi_all_varying,
			lower_multi_all_varying = lower_multi_all_varying,

			center_multi_all_varying_force_presence = center_multi_all_varying_force_presence,
			upper_multi_all_varying_force_presence = upper_multi_all_varying_force_presence,
			lower_multi_all_varying_force_presence = lower_multi_all_varying_force_presence,

			center_multi_all_varying_no_corr = center_multi_all_varying_no_corr,
			upper_multi_all_varying_no_corr = upper_multi_all_varying_no_corr,
			lower_multi_all_varying_no_corr = lower_multi_all_varying_no_corr,

			center_multi_abundance_varying = center_multi_abundance_varying,
			upper_multi_abundance_varying = upper_multi_abundance_varying,
			lower_multi_abundance_varying = lower_multi_abundance_varying,

			center_multi_biomass_varying = center_multi_biomass_varying,
			upper_multi_biomass_varying = upper_multi_biomass_varying,
			lower_multi_biomass_varying = lower_multi_biomass_varying,

			center_multi_psi_varying = center_multi_psi_varying,
			upper_multi_psi_varying = upper_multi_psi_varying,
			lower_multi_psi_varying = lower_multi_psi_varying,

			center_multi_nonbiomass_facets_varying = center_multi_nonbiomass_facets_varying,
			upper_multi_nonbiomass_facets_varying = upper_multi_nonbiomass_facets_varying,
			lower_multi_nonbiomass_facets_varying = lower_multi_nonbiomass_facets_varying,

			center_uni_all_varying = center_uni_all_varying,
			upper_uni_all_varying = upper_uni_all_varying,
			lower_uni_all_varying = lower_uni_all_varying,

			center_uni_all_varying_force_presence = center_uni_all_varying_force_presence,
			upper_uni_all_varying_force_presence = upper_uni_all_varying_force_presence,
			lower_uni_all_varying_force_presence = lower_uni_all_varying_force_presence,

			center_uni_psi_constant = center_uni_psi_constant,
			upper_uni_psi_constant = upper_uni_psi_constant,
			lower_uni_psi_constant = lower_uni_psi_constant

		)

	} # next covariate

	collated <- list(
		facet = 'abundance',
		ranges = ranges,
		responses = responses
	)
	saveRDS(collated, file = './outputs_loretta/integrated_sdm_pdm/response_curves_precalculated_abundance_univar_vs_integrated.rds')

	### responses of BIOMASS
	########################

	say('biomass', level = 2)

	meta <- readRDS(paste0(dir_biomass, '/!meta_biomass.rds'))
	resp_distrib_biomass <- meta$resp_distrib
	transform_biomass <- meta$transform

	form <- meta$formula$formula_biomass
	covars <- terms(form)
	covars <- attr(covars, 'term.labels')

	responses <- ranges <- list() # to save pre-calculated response curves and ranges of BG, species, and sites
	for (covar in covars) {

		say('biomass vs ', covar)

		# values of x covariate (unscaled)
		x <- x_occs[[covar]][ , covar]
		x <- x <- x * scales_counties[covar] + centers_counties[covar]
		if (grepl(covar, pattern = '_log10p1')) x <- 10^x - 1

		# range of variable across background
		this_covar <- if (covar == 'bio12_log10p1') { 'bio12' } else { covar }
		x_bg_empirical <- nam[[this_covar]]
		x_bg_range_sq <- quantile(x_bg_empirical, probs = bg_probs, na.rm = TRUE)
		x_bg_empirical <- nam_ssp370_2071_2100[[this_covar]]
		x_bg_range_ssp370_2071_2100 <- quantile(x_bg_empirical, probs = bg_probs, na.rm = TRUE)

		# range of variable across known occurrences
		x_occs_empirical <- nam[nam$n_andropogon_gerardi > 0]
		x_occs_empirical <- x_occs_empirical[[covar]]
		x_occs_range <- quantile(x_occs_empirical, probs = occs_probs, na.rm = TRUE)
		if (covar == 'bio12_log10p1') x_occs_range <- 10^x_occs_range - 1

		# range of variable across sites
		this_covar <- if (covar == 'bio12_log10p1') { 'bio12' } else { covar }
		x_sites_empirical <- unique(data_biomass_sites_morph$raw_data_biomass[[this_covar]])
		x_sites_range <- range(x_sites_empirical)

		# remember ranges of this covariate for later plotting
		ranges[[covar]] <- list(
			x_bg_range_sq = x_bg_range_sq,
			x_bg_range_ssp370_2071_2100 = x_bg_range_ssp370_2071_2100,
			x_occs_range = x_occs_range,
			x_sites_range = x_sites_range
		)

		# if any non-biomass facets use this covariate, use it's incrementing model matrix... if not, used the fixed matrix
		x_nonbiomass_varying <- x_nonbiomass_fixed
		for (i in seq_along(x_nonbiomass)) {

			if (any(names(x_nonbiomass[[i]]) %in% covar)) {
				x_nonbiomass_varying[[i]] <- x_nonbiomass[[i]][[covar]]
			}

		}

		# predict integrated model: for all facets that share this covariate, vary it
		# use integrated morphological model for these predictions
		say('   pred_multi_all_varying')
		pred_multi_all_varying <- predict_fully_integrated(
			chains = chains_multi_morph,
			nonbiomass_facets = nonbiomass_facets[c('blade_width', 'canopy_diameter', 'height')],
			resp_distrib_biomass = resp_distrib_biomass,
			transform_biomass = transform_biomass,
			x_occs = x_occs[[covar]],
			x_psi = x_occs[[covar]],
			x_biomass = x_biomass[[covar]],
			x_nonbiomass = x_nonbiomass_varying,
			w_occs = NULL,
			sampled = TRUE,
			type = 'samples'
		)

		preds <- pred_multi_all_varying$preds_biomass
		center_multi_all_varying <- apply(preds, 2, median)
		upper_multi_all_varying <- apply(preds, 2, quantile, 0.9)
		lower_multi_all_varying <- apply(preds, 2, quantile, 0.1)

		# predict integrated model: for all facets that share this covariate, vary it but force species to be present
		# use integrated morphological model for these predictions
		say('   pred_multi_all_varying_force_presence')
		pred_multi_all_varying_force_presence <- predict_fully_integrated(
			chains = chains_multi_morph,
			nonbiomass_facets = nonbiomass_facets[c('blade_width', 'canopy_diameter', 'height')],
			resp_distrib_biomass = resp_distrib_biomass,
			transform_biomass = transform_biomass,
			x_occs = x_occs[[covar]],
			x_psi = x_occs[[covar]],
			x_biomass = x_biomass[[covar]],
			x_nonbiomass = x_nonbiomass_varying,
			w_occs = NULL,
			sampled = TRUE,
			type = 'samples',
			force_presence = TRUE
		)

		preds <- pred_multi_all_varying_force_presence$preds_biomass
		center_multi_all_varying_force_presence <- apply(preds, 2, median)
		upper_multi_all_varying_force_presence <- apply(preds, 2, quantile, 0.9)
		lower_multi_all_varying_force_presence <- apply(preds, 2, quantile, 0.1)

		# predict integrated model: for all facets that share this covariate, vary it but force facet correlatons to 0
		# use integrated morphological model for these predictions
		say('   pred_multi_all_varying_no_corr')
		pred_multi_all_varying_no_corr <- predict_fully_integrated(
			chains = chains_multi_morph_no_corr,
			nonbiomass_facets = nonbiomass_facets[c('blade_width', 'canopy_diameter', 'height')],
			resp_distrib_biomass = resp_distrib_biomass,
			transform_biomass = transform_biomass,
			x_occs = x_occs[[covar]],
			x_psi = x_occs[[covar]],
			x_biomass = x_biomass_varying,
			x_nonbiomass = x_nonbiomass_varying,
			w_occs = NULL,
			sampled = TRUE,
			type = 'samples',
			force_presence = FALSE
		)

		preds <- pred_multi_all_varying_no_corr$preds_biomass
		center_multi_all_varying_no_corr <- apply(preds, 2, mean)
		upper_multi_all_varying_no_corr <- apply(preds, 2, quantile, 0.9)
		lower_multi_all_varying_no_corr <- apply(preds, 2, quantile, 0.1)

		# predict integrated model: only allow abundance to vary
		# use integrated morphological model for these predictions
		say('   pred_multi_abundance_varying')
		pred_multi_abundance_varying <- predict_fully_integrated(
			chains = chains_multi_morph,
			nonbiomass_facets = nonbiomass_facets[c('blade_width', 'canopy_diameter', 'height')],
			resp_distrib_biomass = resp_distrib_biomass,
			transform_biomass = transform_biomass,
			x_occs = x_occs[[covar]],
			x_psi = x_occs_fixed,
			x_biomass = x_biomass[[covar]],
			x_nonbiomass = x_nonbiomass_fixed,
			w_occs = NULL,
			sampled = TRUE,
			type = 'samples',
			force_presence = FALSE
		)

		preds <- pred_multi_abundance_varying$preds_biomass
		center_multi_abundance_varying <- apply(preds, 2, median)
		upper_multi_abundance_varying <- apply(preds, 2, quantile, 0.9)
		lower_multi_abundance_varying <- apply(preds, 2, quantile, 0.1)

		# predict integrated model: only vary biomass
		# use integrated morphological model for these predictions
		say('   pred_multi_biomass_varying')
		pred_multi_biomass_varying <- predict_fully_integrated(
			chains = chains_multi_morph,
			nonbiomass_facets = nonbiomass_facets[c('blade_width', 'canopy_diameter', 'height')],
			resp_distrib_biomass = resp_distrib_biomass,
			transform_biomass = transform_biomass,
			x_occs = x_occs[[covar]],
			x_psi = x_occs[[covar]],
			x_biomass = x_biomass[[covar]],
			x_nonbiomass = x_nonbiomass_fixed,
			w_occs = NULL,
			sampled = TRUE,
			type = 'samples'
		)

		preds <- pred_multi_biomass_varying$preds_biomass
		center_multi_biomass_varying <- apply(preds, 2, median)
		upper_multi_biomass_varying <- apply(preds, 2, quantile, 0.9)
		lower_multi_biomass_varying <- apply(preds, 2, quantile, 0.1)

		# predict integrated model: allow covariate to vary only for biomass and psi
		# use integrated morphological model for these predictions
		say('   pred_multi_psi_varying')
		pred_multi_psi_varying <- predict_fully_integrated(
			chains = chains_multi_morph,
			nonbiomass_facets = nonbiomass_facets[c('blade_width', 'canopy_diameter', 'height')],
			resp_distrib_biomass = resp_distrib_biomass,
			transform_biomass = transform_biomass,
			x_occs = x_occs_fixed,
			x_psi = x_occs[[covar]],
			x_biomass = x_biomass[[covar]],
			x_nonbiomass = x_nonbiomass_fixed,
			w_occs = NULL,
			sampled = TRUE,
			type = 'samples'
		)

		preds <- pred_multi_psi_varying$preds_biomass
		center_multi_psi_varying <- apply(preds, 2, median)
		upper_multi_psi_varying <- apply(preds, 2, quantile, 0.9)
		lower_multi_psi_varying <- apply(preds, 2, quantile, 0.1)

		# predict integrated model: allow covariate to vary only for biomass and non-biomass facets
		# use integrated morphological model for these predictions
		say('   pred_multi_nonbiomass_facets_varying')
		pred_multi_nonbiomass_facets_varying <- predict_fully_integrated(
			chains = chains_multi_morph,
			nonbiomass_facets = nonbiomass_facets[c('blade_width', 'canopy_diameter', 'height')],
			resp_distrib_biomass = resp_distrib_biomass,
			transform_biomass = transform_biomass,
			x_occs = x_occs_fixed,
			x_psi = x_occs_fixed,
			x_biomass = x_biomass[[covar]],
			x_nonbiomass = x_nonbiomass_varying,
			w_occs = NULL,
			sampled = TRUE,
			type = 'samples'
		)

		preds <- pred_multi_nonbiomass_facets_varying$preds_biomass
		center_multi_nonbiomass_facets_varying <- apply(preds, 2, median)
		upper_multi_nonbiomass_facets_varying <- apply(preds, 2, quantile, 0.9)
		lower_multi_nonbiomass_facets_varying <- apply(preds, 2, quantile, 0.1)

		# predict univariate model: allow covariate to vary for biomass and psi
		say('   pred_uni_all_varying')
		pred_uni_all_varying <- predict_biomass(
			chains = chains_biomass_uni,
			x = x_biomass[[covar]],
			x_psi = x_biomass[[covar]],
			resp_distrib = resp_distrib_biomass,
			transform = transform_biomass
		)

		center_uni_all_varying <- apply(pred_uni_all_varying, 2, median)
		upper_uni_all_varying <- apply(pred_uni_all_varying, 2, quantile, 0.9)
		lower_uni_all_varying <- apply(pred_uni_all_varying, 2, quantile, 0.1)
	
		# predict univariate model: allow covariate to vary for biomass and psi but force presence
		say('   pred_uni_all_varying_force_presence')
		pred_uni_all_varying_force_presence <- predict_biomass(
			chains = chains_biomass_uni,
			x = x_biomass[[covar]],
			x_psi = x_biomass[[covar]],
			resp_distrib = resp_distrib_biomass,
			transform = transform_biomass,
			force_presence = TRUE
		)

		center_uni_all_varying_force_presence <- apply(pred_uni_all_varying_force_presence, 2, median)
		upper_uni_all_varying_force_presence <- apply(pred_uni_all_varying_force_presence, 2, quantile, 0.9)
		lower_uni_all_varying_force_presence <- apply(pred_uni_all_varying_force_presence, 2, quantile, 0.1)
	
		# predict univariate model: allow covariate to vary only for biomass
		say('   pred_uni_psi_constant')
		pred_uni_psi_constant <- predict_biomass(
			chains = chains_biomass_uni,
			x = x_biomass[[covar]],
			x_psi = x_biomass_fixed,
			resp_distrib = resp_distrib_biomass,
			transform = transform_biomass
		)

		center_uni_psi_constant <- apply(pred_uni_psi_constant, 2, median)
		upper_uni_psi_constant <- apply(pred_uni_psi_constant, 2, quantile, 0.9)
		lower_uni_psi_constant <- apply(pred_uni_psi_constant, 2, quantile, 0.1)
	
		# remember
		responses[[covar]] <- data.table(
			x = x,
			
			center_multi_all_varying = center_multi_all_varying,
			upper_multi_all_varying = upper_multi_all_varying,
			lower_multi_all_varying = lower_multi_all_varying,

			center_multi_all_varying_force_presence = center_multi_all_varying_force_presence,
			upper_multi_all_varying_force_presence = upper_multi_all_varying_force_presence,
			lower_multi_all_varying_force_presence = lower_multi_all_varying_force_presence,

			center_multi_all_varying_no_corr = center_multi_all_varying_no_corr,
			upper_multi_all_varying_no_corr = upper_multi_all_varying_no_corr,
			lower_multi_all_varying_no_corr = lower_multi_all_varying_no_corr,

			center_multi_abundance_varying = center_multi_abundance_varying,
			upper_multi_abundance_varying = upper_multi_abundance_varying,
			lower_multi_abundance_varying = lower_multi_abundance_varying,

			center_multi_biomass_varying = center_multi_biomass_varying,
			upper_multi_biomass_varying = upper_multi_biomass_varying,
			lower_multi_biomass_varying = lower_multi_biomass_varying,

			center_multi_psi_varying = center_multi_psi_varying,
			upper_multi_psi_varying = upper_multi_psi_varying,
			lower_multi_psi_varying = lower_multi_psi_varying,

			center_multi_nonbiomass_facets_varying = center_multi_nonbiomass_facets_varying,
			upper_multi_nonbiomass_facets_varying = upper_multi_nonbiomass_facets_varying,
			lower_multi_nonbiomass_facets_varying = lower_multi_nonbiomass_facets_varying,

			center_uni_all_varying = center_uni_all_varying,
			upper_uni_all_varying = upper_uni_all_varying,
			lower_uni_all_varying = lower_uni_all_varying,

			center_uni_all_varying_force_presence = center_uni_all_varying_force_presence,
			upper_uni_all_varying_force_presence = upper_uni_all_varying_force_presence,
			lower_uni_all_varying_force_presence = lower_uni_all_varying_force_presence,

			center_uni_psi_constant = center_uni_psi_constant,
			upper_uni_psi_constant = upper_uni_psi_constant,
			lower_uni_psi_constant = lower_uni_psi_constant

		)

	} # next covariate

	collated <- list(
		facet = 'biomass',
		ranges = ranges,
		responses = responses
	)
	saveRDS(collated, file = './outputs_loretta/integrated_sdm_pdm/response_curves_precalculated_biomass_univar_vs_integrated.rds')

	### non-biomass traits
	###################################
	say('non-biomass traits', level = 2)

	# biomass metadata
	meta <- readRDS(paste0(dir_biomass, '/!meta_biomass.rds'))
	resp_distrib_biomass <- meta$resp_distrib
	transform_biomass <- meta$transform

	for (f in seq_along(nonbiomass_facets)) {

		facet <- names(nonbiomass_facets)[f]

		# meta-data on this model
		filename_form <- nonbiomass_facets[[facet]]$filename
		resp_distrib_focal_facet <- nonbiomass_facets[[facet]]$resp_distrib
		transform_focal_facet <- nonbiomass_facets[[facet]]$transform

		# predictions from univariate model
		univar_dir <- paste0('./outputs_loretta/integrated_sdm_pdm/models_', facet, '/[', facet, '~', tolower(resp_distrib_focal_facet), '(', filename_form, ')]')

		chains_facet_uni <- readRDS(paste0(univar_dir, '/chains.rds'))
		if (abbreviate) chains_facet_uni$samples <- chains_facet_uni$samples[1]
		covars <- names(x_nonbiomass[[facet]])

		if (facet %in% c('blade_width', 'canopy_diameter', 'height')) {
		
			this_data_nonbiomass_sites <- data_nonbiomass_sites_morph
			this_data_nonbiomass_counties <- data_nonbiomass_counties_morph
			this_nonbiomass_facets <- nonbiomass_facets[c('blade_width', 'canopy_diameter', 'height')]
			formulae <- formulae_morph
			chains_multi <- chains_multi_morph
			chains_multi_no_corr <- chains_multi_morph_no_corr
		
		} else if (facet %in% c('cn_ratio', 'photosynthetic_rate', 'spad', 'stomatal_conductance', 'transpiration_rate')) {

			this_data_nonbiomass_sites <- data_nonbiomass_sites_phys
			this_data_nonbiomass_counties <- data_nonbiomass_counties_phys
			this_nonbiomass_facets <- nonbiomass_facets[c('cn_ratio', 'photosynthetic_rate', 'spad', 'stomatal_conductance', 'transpiration_rate')]
			formulae <- formulae_phys
			chains_multi <- chains_multi_phys
			chains_multi_no_corr <- chains_multi_phys_no_corr
		
		}

		responses <- ranges <- list()
		for (covar in covars) {

			say(facet, ' vs ', covar)

			# values of x covariate (unscaled)
			x <- x_nonbiomass[[facet]][[covar]][ , covar]
			x <- x * scales_sites[covar] + centers_sites[covar]
			if (covar == 'bio12_log10p1') x <- 10^x - 1

			# range of variable across background
			x_bg_empirical <- nam[[covar]]
			if (covar == 'bio12_log10p1') x_bg_empirical <- 10^x_bg_empirical - 1
			x_bg_range_sq <- quantile(x_bg_empirical, probs = bg_probs, na.rm = TRUE)
			x_bg_empirical <- nam_ssp370_2071_2100[[this_covar]]
			x_bg_range_ssp370_2071_2100 <- quantile(x_bg_empirical, probs = bg_probs, na.rm = TRUE)

			# range of variable across known occurrences
			x_occs_empirical <- nam[nam$n_andropogon_gerardi > 0]
			x_occs_empirical <- x_occs_empirical[[covar]]
			if (covar == 'bio12_log10p1') x_occs_empirical <- 10^x_occs_empirical - 1
			x_occs_range <- quantile(x_occs_empirical, probs = occs_probs, na.rm = TRUE)

			# range of variable across sites
			this_covar <- if (covar == 'bio12_log10p1') { 'bio12' } else { covar }
			x_sites_empirical <- unique(this_data_nonbiomass_sites[[facet]]$site_data_raw[[this_covar]])
			x_sites_range <- range(x_sites_empirical)

			# remember ranges of this covariate for later plotting
			ranges[[covar]] <- list(
				x_bg_range_sq = x_bg_range_sq,
				x_bg_range_ssp370_2071_2100 = x_bg_range_ssp370_2071_2100,
				x_occs_range = x_occs_range,
				x_sites_range = x_sites_range
			)

			# create list with model matrices where for the focal variable, the covariate of interest changes but the covariate is constant for other, non-focal variables
			x_nonbiomass_only_focal_varying <- x_nonbiomass_fixed
			x_nonbiomass_only_focal_varying[[facet]] <- x_nonbiomass[[facet]][[covar]]

			# if any non-facet facets use this covariate, use it's incrementing model matrix... if not, used the fixed matrix
			x_nonbiomass_varying <- x_nonbiomass_fixed
			for (i in seq_along(x_nonbiomass)) {

				if (any(names(x_nonbiomass[[i]]) %in% covar)) {
					x_nonbiomass_varying[[i]] <- x_nonbiomass[[i]][[covar]]
				}

			}

			# if occurrences/psi use this covariate, use it's incrementing model matrix... if not, used the fixed matrix
			x_occs_covar <- if (any(names(x_occs) %in% covar)) {
				x_occs[[covar]]
			} else {
				x_occs_fixed
			}

			x_biomass_covar <- if (any(names(x_biomass) %in% covar)) {
				x_biomass[[covar]]
			} else {
				x_biomass_fixed
			}

			# predict integrated model: for all facets that share this covariate, vary it
			say('   pred_multi_all_varying')
			pred_multi_all_varying <- predict_fully_integrated(
				chains = chains_multi,
				nonbiomass_facets = this_nonbiomass_facets,
				resp_distrib_biomass = resp_distrib_biomass,
				transform_biomass = transform_biomass,
				x_occs = x_occs_covar,
				x_psi = x_occs_covar,
				x_biomass = x_biomass_covar,
				x_nonbiomass = x_nonbiomass_varying,
				w_occs = NULL,
				sampled = TRUE,
				type = 'samples'
			)

			preds <- pred_multi_all_varying$preds_nonbiomass[[facet]]
			center_multi_all_varying <- apply(preds, 2, median)
			upper_multi_all_varying <- apply(preds, 2, quantile, 0.9)
			lower_multi_all_varying <- apply(preds, 2, quantile, 0.1)

			# predict integrated model: for all facets that share this covariate, vary it but force presence
			say('   pred_multi_all_varying_force_presence')
			pred_multi_all_varying_force_presence <- predict_fully_integrated(
				chains = chains_multi,
				nonbiomass_facets = this_nonbiomass_facets,
				resp_distrib_biomass = resp_distrib_biomass,
				transform_biomass = transform_biomass,
				x_occs = x_occs_covar,
				x_psi = x_occs_covar,
				x_biomass = x_biomass_covar,
				x_nonbiomass = x_nonbiomass_varying,
				w_occs = NULL,
				sampled = TRUE,
				type = 'samples',
				force_presence = TRUE
			)

			preds <- pred_multi_all_varying_force_presence$preds_nonbiomass[[facet]]
			center_multi_all_varying_force_presence <- apply(preds, 2, median)
			upper_multi_all_varying_force_presence <- apply(preds, 2, quantile, 0.9)
			lower_multi_all_varying_force_presence <- apply(preds, 2, quantile, 0.1)

			# predict integrated model: for all facets that share this covariate, vary it but force facet correlatons to 0
			# use integrated morphological model for these predictions
			say('   pred_multi_all_varying_no_corr')
			pred_multi_all_varying_no_corr <- predict_fully_integrated(
				chains = chains_multi_no_corr,
				nonbiomass_facets = this_nonbiomass_facets,
				resp_distrib_biomass = resp_distrib_biomass,
				transform_biomass = transform_biomass,
				x_occs = x_occs_covar,
				x_psi = x_occs_covar,
				x_biomass = x_biomass_covar,
				x_nonbiomass = x_nonbiomass_varying,
				w_occs = NULL,
				sampled = TRUE,
				type = 'samples',
				force_presence = FALSE
			)

			preds <- pred_multi_all_varying_no_corr$preds_nonbiomass[[facet]]
			center_multi_all_varying_no_corr <- apply(preds, 2, mean)
			upper_multi_all_varying_no_corr <- apply(preds, 2, quantile, 0.9)
			lower_multi_all_varying_no_corr <- apply(preds, 2, quantile, 0.1)

			# predict integrated model: vary focal facet's covariate and occurrence covariate, hold all others at mean values
			say('   pred_multi_occs_varying')
			pred_multi_occs_varying <- predict_fully_integrated(
				chains = chains_multi,
				nonbiomass_facets = this_nonbiomass_facets,
				resp_distrib_biomass = resp_distrib_biomass,
				transform_biomass = transform_biomass,
				x_occs = x_occs_covar,
				x_psi = x_occs_fixed,
				x_biomass = x_biomass_fixed,
				x_nonbiomass = x_nonbiomass_only_focal_varying,
				w_occs = NULL,
				sampled = TRUE,
				type = 'samples'
			)

			preds <- pred_multi_occs_varying$preds_nonbiomass[[facet]]
			center_multi_abundance_varying <- apply(preds, 2, median)
			upper_multi_abundance_varying <- apply(preds, 2, quantile, 0.9)
			lower_multi_abundance_varying <- apply(preds, 2, quantile, 0.1)

			# predict integrated model: vary biomass and focal nonbiomass facet matrix, hold all others constant
			say('   pred_multi_biomass_varying')
			pred_multi_biomass_varying <- predict_fully_integrated(
				chains = chains_multi,
				nonbiomass_facets = this_nonbiomass_facets,
				resp_distrib_biomass = resp_distrib_biomass,
				transform_biomass = transform_biomass,
				x_occs = x_occs_fixed,
				x_psi = x_occs_fixed,
				x_biomass = x_biomass_covar,
				x_nonbiomass = x_nonbiomass_only_focal_varying,
				w_occs = NULL,
				sampled = TRUE,
				type = 'samples'
			)

			preds <- pred_multi_biomass_varying$preds_nonbiomass[[facet]]
			center_multi_biomass_varying <- apply(preds, 2, median)
			upper_multi_biomass_varying <- apply(preds, 2, quantile, 0.9)
			lower_multi_biomass_varying <- apply(preds, 2, quantile, 0.1)

			# predict integrated model: vary covariate of focal facet and psi, hold all others at mean values
			say('   pred_multi_psi_varying')
			pred_multi_psi_varying <- predict_fully_integrated(
				chains = chains_multi,
				nonbiomass_facets = this_nonbiomass_facets,
				resp_distrib_biomass = resp_distrib_biomass,
				transform_biomass = transform_biomass,
				x_occs = x_occs_fixed,
				x_psi = x_occs_covar,
				x_biomass = x_biomass_fixed,
				x_nonbiomass = x_nonbiomass_only_focal_varying,
				w_occs = NULL,
				sampled = TRUE,
				type = 'samples'
			)

			preds <- pred_multi_psi_varying$preds_nonbiomass[[facet]]
			center_multi_psi_varying <- apply(preds, 2, median)
			upper_multi_psi_varying <- apply(preds, 2, quantile, 0.9)
			lower_multi_psi_varying <- apply(preds, 2, quantile, 0.1)

			# predict integrated model: for focal and non-focal non-biomass facets, allow covariate to vary
			say('   pred_multi_nonbiomass_facets_varying')
			pred_multi_nonbiomass_facets_varying <- predict_fully_integrated(
				chains = chains_multi,
				nonbiomass_facets = this_nonbiomass_facets,
				resp_distrib_biomass = resp_distrib_biomass,
				transform_biomass = transform_biomass,
				x_occs = x_occs_fixed,
				x_psi = x_occs_fixed,
				x_biomass = x_biomass_fixed,
				x_nonbiomass = x_nonbiomass_varying,
				w_occs = NULL,
				sampled = TRUE,
				type = 'samples'
			)

			preds <- pred_multi_nonbiomass_facets_varying$preds_nonbiomass[[facet]]
			center_multi_nonbiomass_facets_varying <- apply(preds, 2, median)
			upper_multi_nonbiomass_facets_varying <- apply(preds, 2, quantile, 0.9)
			lower_multi_nonbiomass_facets_varying <- apply(preds, 2, quantile, 0.1)

			# predict univariate model: allow covariate for facet and psi to vary
			say('   pred_uni_all_varying')
			pred_uni_all_varying <- predict_nonbiomass_single_trait(
				chains = chains_facet_uni,
				x = x_nonbiomass[[facet]][[covar]],
				x_psi = x_nonbiomass[[facet]][[covar]],
				resp_distrib = resp_distrib_focal_facet,
				transform = transform_focal_facet
			)

			center_uni_all_varying <- apply(pred_uni_all_varying, 2, median)
			upper_uni_all_varying <- apply(pred_uni_all_varying, 2, quantile, 0.9)
			lower_uni_all_varying <- apply(pred_uni_all_varying, 2, quantile, 0.1)
		
			# predict univariate model: allow covariate for facet and psi to vary but force species to be present
			say('   pred_uni_all_varying_force_presence')
			pred_uni_all_varying_force_presence <- predict_nonbiomass_single_trait(
				chains = chains_facet_uni,
				x = x_nonbiomass[[facet]][[covar]],
				x_psi = x_nonbiomass[[facet]][[covar]],
				resp_distrib = resp_distrib_focal_facet,
				transform = transform_focal_facet,
				force_presence = TRUE
			)

			center_uni_all_varying_force_presence <- apply(pred_uni_all_varying_force_presence, 2, median)
			upper_uni_all_varying_force_presence <- apply(pred_uni_all_varying_force_presence, 2, quantile, 0.9)
			lower_uni_all_varying_force_presence <- apply(pred_uni_all_varying_force_presence, 2, quantile, 0.1)
		
			# predict univariate model: allow covariate to vary for focal facet, but not for psi
			say('   pred_uni_psi_constant')
			pred_uni_psi_constant <- predict_nonbiomass_single_trait(
				chains = chains_facet_uni,
				x = x_nonbiomass_varying[[facet]],
				x_psi = x_nonbiomass_fixed[[facet]],
				resp_distrib = resp_distrib_focal_facet,
				transform = transform_focal_facet
			)

			center_uni_psi_constant <- apply(pred_uni_psi_constant, 2, median)
			upper_uni_psi_constant <- apply(pred_uni_psi_constant, 2, quantile, 0.9)
			lower_uni_psi_constant <- apply(pred_uni_psi_constant, 2, quantile, 0.1)
		
			# remember
			responses[[covar]] <- data.table(
				x = x,
				
				center_multi_all_varying = center_multi_all_varying,
				upper_multi_all_varying = upper_multi_all_varying,
				lower_multi_all_varying = lower_multi_all_varying,

				center_multi_all_varying_force_presence = center_multi_all_varying_force_presence,
				upper_multi_all_varying_force_presence = upper_multi_all_varying_force_presence,
				lower_multi_all_varying_force_presence = lower_multi_all_varying_force_presence,

				center_multi_all_varying_no_corr = center_multi_all_varying_no_corr,
				upper_multi_all_varying_no_corr = upper_multi_all_varying_no_corr,
				lower_multi_all_varying_no_corr = lower_multi_all_varying_no_corr,

				center_multi_abundance_varying = center_multi_abundance_varying,
				upper_multi_abundance_varying = upper_multi_abundance_varying,
				lower_multi_abundance_varying = lower_multi_abundance_varying,

				center_multi_biomass_varying = center_multi_biomass_varying,
				upper_multi_biomass_varying = upper_multi_biomass_varying,
				lower_multi_biomass_varying = lower_multi_biomass_varying,
				
				center_multi_psi_varying = center_multi_psi_varying,
				upper_multi_psi_varying = upper_multi_psi_varying,
				lower_multi_psi_varying = lower_multi_psi_varying,

				center_multi_nonbiomass_facets_varying = center_multi_nonbiomass_facets_varying,
				upper_multi_nonbiomass_facets_varying = upper_multi_nonbiomass_facets_varying,
				lower_multi_nonbiomass_facets_varying = lower_multi_nonbiomass_facets_varying,

				center_uni_all_varying = center_uni_all_varying,
				upper_uni_all_varying = upper_uni_all_varying,
				lower_uni_all_varying = lower_uni_all_varying,

				center_uni_all_varying_force_presence = center_uni_all_varying_force_presence,
				upper_uni_all_varying_force_presence = upper_uni_all_varying_force_presence,
				lower_uni_all_varying_force_presence = lower_uni_all_varying_force_presence,

				center_uni_psi_constant = center_uni_psi_constant,
				upper_uni_psi_constant = upper_uni_psi_constant,
				lower_uni_psi_constant = lower_uni_psi_constant

			)

		} # next covariate

		collated <- list(
			facet = facet,
			ranges = ranges,
			responses = responses
		)
		saveRDS(collated, file = paste0('./outputs_loretta/integrated_sdm_pdm/response_curves_precalculated_', facet, '_univar_vs_integrated.rds'))

	} # next non-facet facet

# say('######################################')
# say('### response curves for supplement ###')
# say('######################################')

# 	# Create response curves for each facet. Each facet has one panel per predictor. Each panel has a response curve for:
# 	# * Single-facet model with the focal covariate changing for the focal facet and psi
# 	# * Single-facet model with the focal covariate changing for the focal facet but forcing presence
# 	# * Multi-facet model with the focal covariate changing for ALL facets
# 	# * Multi-facet model with the focal covariate changing for ALL facets but forcing presence
# 	# * Multi-facet model with the focal covariate changing for ALL facets but setting all inter-facet correlations to 0

# 	panel_title_size <- 10 # size of panel plot titles
# 	center_line_size <- 1 # thickness of lines for central tendency of each facet
# 	ymax_mult <- 1.05 # set ymax in plots equal to highest central tendency prediction times this value
# 	env_range_label_size <- 0.7 # size of text labeling environmental ranges

# 	### responses of ABUNDANCE
# 	##########################

# 	say('abundance')

# 	covars <- data_occs_counties_morph$covariates
# 	resp_curves_occs <- list()

# 	for (covar in covars) {

# 		say('response curves for abundance vs ', covar)

# 		# collate predictions for plotting
# 		n <- length(x)

# 		labels <- c(
# 			paste0('Simple (abundance & Ψ varying)'),
# 			paste0('Simple (only abundance varying)'),
# 			paste0('Integrated (all facets varying)'),
# 			paste0('Integrated (only abundance varying)'),
# 			paste0('Integrated (abundance & biomass varying)'),
# 			paste0('Integrated (abundance & Ψ varying)'),
# 			paste0('Integrated (abundance & non-biomass traits varying)')
# 		)
# 		linetypes <- c('solid', 'dashed', 'solid', 'twodash', 'dotdash', 'dotted', 'longdash')#, 'longdash'
# 		names(linetypes) <- labels

# 		pred_centers <- data.table(
# 			model = rep(labels, each = n),
# 			x = rep(x, 7),
# 			y = c(
# 				center_uni_all_varying,
# 				center_uni_psi_constant,
# 				center_multi_all_varying,
# 				center_multi_abundance_varying,
# 				center_multi_biomass_varying,
# 				center_multi_psi_varying,
# 				center_multi_nonbiomass_facets_varying
# 			)
# 		)

# 		pred_ci <- data.table(
# 			model = rep(labels, each = n),
# 			x = rep(x, 7),
# 			ymin = c(
# 				lower_uni_all_varying,
# 				lower_uni_psi_constant,
# 				lower_multi_all_varying,
# 				lower_multi_abundance_varying,
# 				lower_multi_biomass_varying,
# 				lower_multi_psi_varying,
# 				lower_multi_nonbiomass_facets_varying
# 			),
# 			ymax = c(
# 				upper_uni_all_varying,
# 				upper_uni_psi_constant,
# 				upper_multi_all_varying,
# 				upper_multi_abundance_varying,
# 				upper_multi_biomass_varying,
# 				upper_multi_psi_varying,
# 				upper_multi_nonbiomass_facets_varying
# 			)
# 		)

# 		pred_centers$model <- factor(pred_centers$model, levels = labels)
# 		pred_ci$model <- factor(pred_ci$model, levels = labels)

# 		# plot data
# 		nice <- get_nice_predictor(covar)
# 		title <- nice$short
# 		xlab <- nice$long

# 		xlim <- range(x)
# 		ylim <- c(0, ymax_mult * max(pred_centers$y[pred_centers$model != 'Simple (only abundance varying)' & !is.infinite(pred_centers$y)]))

# 		log10p1_trans <- function() {
# 			trans_new(
# 				'log10p1',
# 				transform = function(x) log10(x + 1),
# 				inverse   = function(x) 10^x - 1,
# 				domain    = c(0, Inf)
# 			)
# 		}

# 		if (covar != 'bio12_log10p1') {		
# 			bg_range_y <- ylim[1] + 0.005 * diff(ylim)
# 			species_range_y <- ylim[1] + 0.013 * diff(ylim)
# 			sites_range_y <- ylim[1] + 0.023 * diff(ylim)
# 		} else {
# 			bg_range_y <- ylim[1] + 0.1 * diff(log10(ylim + 1))
# 			species_range_y <- ylim[1] + 0.28 * diff(log10(ylim + 1))
# 			sites_range_y <- ylim[1] + 0.6 * diff(log10(ylim + 1))
# 		}

# 		resp_curves_occs[[length(resp_curves_occs) + 1]] <- ggplot() +
# 			annotate(
# 				'segment',
# 				x = x_bg_range_sq[1], xend = x_bg_range_sq[2], y = bg_range_y, yend = bg_range_y,
# 				size = 4, color = 'gray50', alpha = 0.7
# 			) +
# 			annotate(
# 				'segment',
# 				x = x_occs_range[1], xend = x_occs_range[2], y = species_range_y, yend = species_range_y,
# 				size = 4, color = 'chartreuse3', alpha = 0.7
# 			) +
# 			annotate(
# 				'segment',
# 				x = x_sites_range[1], xend = x_sites_range[2], y = sites_range_y, yend = sites_range_y,
# 				size = 4, color = 'yellow', alpha = 0.7
# 			) +
# 			# geom_ribbon(data = pred_ci, aes(x = x, ymin = ymin, ymax = ymax, color = model), fill = NA, size = 0.4, alpha = 0.5, linetype = 'dotted') +
# 			geom_line(data = pred_centers, aes(x = x, y = y, color = model, linetype = model), size = center_line_size) +
# 			scale_linetype_manual(values = linetypes, name = 'Model') +
# 			coord_cartesian(ylim = c(ylim[1], ylim[2] * 1.15), xlim = c(xlim[1], xlim[2])) +
# 			scale_y_continuous(trans = log10p1_trans()) +
# 			labs(color = 'Model', x = xlab, y = 'Relative abundance') +
# 			theme(
# 				legend.key.width = grid::unit(2, 'cm'),
# 				plot.title = element_text(size = panel_title_size),
# 				plot.margin = unit(c(1.5, 0.5, 0.5, 0.5), 'lines')
# 			)

# 	} # next covariate for ABUNDANCE

# 	names(resp_curves_occs) <- covars
# 	resp_curves[[length(resp_curves) + 1]] <- resp_curves_occs
# 	names(resp_curves)[length(resp_curves)] <- 'occs'

# 	# use 3 plots + 1 legend panel
# 	plot_panels <- resp_curves_occs[seq_len(min(3, length(resp_curves_occs)))]

# 	for (i in seq_along(plot_panels)) {
# 		plot_panels[[i]] <- plot_panels[[i]] + theme(legend.position = 'none')
# 	}

# 	legend_panel <- get_legend(
# 		resp_curves_occs[[1]] + theme(legend.position = 'right')
# 	)

# 	range_legend_panel <- get_legend(
# 		ggplot(
# 			data.frame(
# 			label = factor(c('Sites', 'Species', 'Background'), levels = c('Sites', 'Species', 'Background')),
# 			x = 1,
# 			y = 1:3
# 			),
# 			aes(x = x, y = y, fill = label)
# 		) +
# 		geom_tile() +
# 		scale_fill_manual(
# 			values = c(
# 				Sites = 'yellow',
# 				Species = 'chartreuse3',
# 				Background = 'grey50'
# 			),
# 			name = 'Environmental range'
# 		) +
# 		guides(fill = guide_legend(override.aes = list(color = NA))) +
# 		theme_void() +
# 		theme(legend.position = "right")
# 	)

# 	legend_stack <- plot_grid(
# 		legend_panel,
# 		range_legend_panel,
# 		# axis = 'trbl',
# 		ncol = 1,
# 		rel_heights = c(1, 0.8)
# 	)
	
# 	resps <- plot_grid(
# 		plotlist = c(plot_panels, list(legend_stack)),
# 		nrow = 1,
# 		ncol = 4,
# 		rel_widths = c(1, 1, 1, 1.2),
# 		align = 'hv'
# 	)

# 	plot_title <- ggdraw() +
# 	draw_label('Abundance', x = 0, hjust = 0) +
# 	theme(
# 		plot.margin = margin(0, 0, 0, 7)
# 	)

# 	resps <- plot_grid(plot_title, resps, ncol = 1, rel_heights = c(0.06, 1))
# 	ggsave(resps, filename = './outputs_loretta/integrated_sdm_pdm/response_curves_abundance_univar_vs_integrated.png', width = 14.6, height = 4, units = 'in', dpi = 600, bg = 'white')

# 	### responses of BIOMASS
# 	########################

# 	say('biomass', level = 2)

# 	covars <- data_biomass_sites_morph$covariates
# 	resp_curves_biomass <- list()

# 	for (covar in covars) {

# 		say('biomass vs ', covar)

# 		# collate predictions for plotting
# 		n <- length(x)

# 		labels <- c(
# 			'Simple (biomass & Ψ varying)',
# 			'Simple (only biomass varying)',
# 			'Integrated (all facets varying)',
# 			'Integrated (only biomass varying)',
# 			'Integrated (biomass & abundance varying)',
# 			'Integrated (biomass & Ψ varying)',
# 			'Integrated (biomass & non-biomass traits varying)'
# 		)
# 		linetypes <- c('solid', 'dashed', 'solid', 'twodash', 'dotdash', 'dotted', 'longdash')
# 		names(linetypes) <- labels

# 		pred_centers <- data.table(
# 			model = rep(labels, each = n),
# 			x = rep(x, 7),
# 			y = c(
# 				center_uni_all_varying,
# 				center_uni_psi_constant,
# 				center_multi_all_varying,
# 				center_multi_biomass_varying,
# 				center_multi_abundance_varying,
# 				center_multi_psi_varying,
# 				center_multi_nonbiomass_facets_varying
# 			)
# 		)

# 		pred_ci <- data.table(
# 			model = rep(labels, each = n),
# 			x = rep(x, 7),
# 			ymin = c(
# 				lower_uni_all_varying,
# 				lower_uni_psi_constant,
# 				lower_multi_all_varying,
# 				lower_multi_biomass_varying,
# 				lower_multi_abundance_varying,
# 				lower_multi_psi_varying,
# 				lower_multi_nonbiomass_facets_varying
# 			),
# 			ymax = c(
# 				upper_uni_all_varying,
# 				upper_uni_psi_constant,
# 				upper_multi_all_varying,
# 				upper_multi_biomass_varying,
# 				upper_multi_abundance_varying,
# 				upper_multi_psi_varying,
# 				upper_multi_nonbiomass_facets_varying
# 			)
# 		)

# 		pred_centers$model <- factor(pred_centers$model, levels = labels)
# 		pred_ci$model <- factor(pred_ci$model, levels = labels)

# 		ylim <- c(0, ymax_mult * max(pred_centers$y))
# 		xlim <- range(x)

# 		# plot data
# 		nice <- get_nice_predictor(covar)
# 		title <- nice$short
# 		xlab <- nice$long

# 		bg_range_y <- ylim[1] + 0.035 * diff(ylim)
# 		species_range_y <- ylim[1] + 0.09 * diff(ylim)
# 		sites_range_y <- ylim[1] + 0.16 * diff(ylim)

# 		resp_curves_biomass[[length(resp_curves_biomass) + 1]] <- ggplot() +
# 			annotate(
# 				'segment',
# 				x = x_bg_range_sq[1], xend = x_bg_range_sq[2], y = bg_range_y, yend = bg_range_y,
# 				size = 4, color = 'gray50', alpha = 0.7
# 			) +
# 			annotate(
# 				'segment',
# 				x = x_occs_range[1], xend = x_occs_range[2], y = species_range_y, yend = species_range_y,
# 				size = 4, color = 'chartreuse3', alpha = 0.7
# 			) +
# 			annotate(
# 				'segment',
# 				x = x_sites_range[1], xend = x_sites_range[2], y = sites_range_y, yend = sites_range_y,
# 				size = 4, color = 'yellow', alpha = 0.7
# 			) +
# 			# geom_ribbon(data = pred_ci, aes(x = x, ymin = ymin, ymax = ymax, color = model), fill = NA, size = 0.4, alpha = 0.5, linetype = 'dotted') +
# 			geom_line(data = pred_centers, aes(x = x, y = y, color = model, linetype = model), size = center_line_size) +
# 			scale_linetype_manual(values = linetypes, name = 'Model') +
# 			xlab(xlab) +
# 			ylab(expression('Biomass (g)')) +
# 			labs(color = 'Model') +
# 			coord_cartesian(ylim = c(ylim[1], ylim[2] * 1.15), xlim = c(xlim[1], xlim[2])) +
# 			theme(
# 				plot.title = element_text(size = panel_title_size),
# 				legend.key.width = grid::unit(1.8, 'cm'),
# 			)

# 	} # next covaraite for BIOMASS

# 	names(resp_curves_biomass) <- covars
# 	resp_curves[[length(resp_curves) + 1]] <- resp_curves_biomass
# 	names(resp_curves)[length(resp_curves)] <- 'biomass'

# 	# use plots + 1 legend panel
# 	plot_panels <- resp_curves_biomass
# 	for (i in seq_along(plot_panels)) {
# 		plot_panels[[i]] <- plot_panels[[i]] + theme(legend.position = 'none')
# 	}

# 	legend_panel <- get_legend(
# 		resp_curves_biomass[[1]] + theme(legend.position = 'right')
# 	)

# 	range_legend_panel <- get_legend(
# 		ggplot(
# 			data.frame(
# 			label = factor(c('Sites', 'Species', 'Background'), levels = c('Sites', 'Species', 'Background')),
# 			x = 1,
# 			y = 1:3
# 			),
# 			aes(x = x, y = y, fill = label)
# 		) +
# 		geom_tile() +
# 		scale_fill_manual(
# 			values = c(
# 				Sites = 'yellow',
# 				Species = 'chartreuse3',
# 				Background = 'grey50'
# 			),
# 			name = 'Environmental range'
# 		) +
# 		guides(fill = guide_legend(override.aes = list(color = NA))) +
# 		theme_void() +
# 		theme(legend.position = "right")
# 	)

# 	legend_stack <- plot_grid(
# 		legend_panel,
# 		range_legend_panel,
# 		# axis = 'trbl',
# 		ncol = 1,
# 		rel_heights = c(1, 0.5)
# 	)
	
# 	resps <- plot_grid(
# 		plotlist = c(plot_panels, list(legend_stack)),
# 		nrow = 1,
# 		ncol = length(plot_panels) + 1,
# 		rel_widths = c(rep(1, length(plot_panels)), 1.2),
# 		align = 'hv'
# 	)

# 	plot_title <- ggdraw() +
# 		draw_label('Biomass', x = 0, hjust = 0) +
# 		theme(
# 			plot.margin = margin(0, 0, 0, 7)
# 		)
# 	resps <- plot_grid(plot_title, resps, ncol = 1, rel_heights = c(0.06, 1))

# 	width <- 4.6 * length(covars) + 3.2

# 	ggsave(resps, filename = './outputs_loretta/integrated_sdm_pdm/response_curves_biomass_univar_vs_integrated.png', width = width, height = 4, units = 'in', dpi = 600, bg = 'white')

# 	### non-biomass traits
# 	###################################
# 	say('non-biomass traits', level = 2)

# 	for (f in seq_along(nonbiomass_facets)) {

# 		facet <- names(nonbiomass_facets)[f]
# 		covars <- names(x_nonbiomass[[facet]])

# 		resp_curves_facet <- list()
# 		for (covar in covars) {

# 			say(facet, ' vs ', covar)

# 			# collate predictions for plotting
# 			n <- length(x)

# 			facet_short <- tolower(get_nice_trait(facet)$short)
# 			labels <- c(
# 				paste0('Simple (', facet_short, ' & Ψ varying)'),
# 				paste0('Simple (only ', facet_short, ' varying)'),
# 				paste0('Integrated (all facets varying)'),
# 				paste0('Integrated (', facet_short, ' & biomass varying)'),
# 				paste0('Integrated (', facet_short, ' & abundance varying)'),
# 				paste0('Integrated (', facet_short, ' & Ψ varying)'),
# 				paste0('Integrated (non-biomass traits varying)')
# 			)
# 			linetypes <- c('solid', 'dashed', 'solid', 'twodash', 'dotdash', 'dotted')#, 'longdash'
# 			names(linetypes) <- labels

# 			pred_centers <- data.table(
# 				model = rep(labels, each = n),
# 				x = rep(x, 7),
# 				y = c(center_uni_all_varying, center_uni_psi_constant, center_multi_all_varying, center_multi_biomass_varying, center_multi_abundance_varying, center_multi_psi_varying, center_multi_nonbiomass_facets_varying)
# 			)

# 			pred_ci <- data.table(
# 				model = rep(labels, each = n),
# 				x = rep(x, 7),
# 				ymin = c(lower_uni_all_varying, lower_uni_psi_constant, lower_multi_all_varying, lower_multi_biomass_varying, lower_multi_abundance_varying, lower_multi_psi_varying, lower_multi_nonbiomass_facets_varying),
# 				ymax = c(upper_uni_all_varying, upper_uni_psi_constant, upper_multi_all_varying, upper_multi_biomass_varying, upper_multi_abundance_varying, upper_multi_psi_varying, upper_multi_nonbiomass_facets_varying)
# 			)

# 			pred_centers$model <- factor(pred_centers$model, levels = labels)
# 			pred_ci$model <- factor(pred_ci$model, levels = labels)

# 			ylim <- c(0, ymax_mult * max(pred_centers$y))
# 			xlim <- range(x)

# 			# plot data
# 			nice <- get_nice_predictor(covar)
# 			title <- nice$short
# 			xlab <- nice$long

# 			nice <- get_nice_trait(facet)
# 			ylab <- nice$long

# 			bg_range_y <- ylim[1] + 0.035 * diff(ylim)
# 			species_range_y <- ylim[1] + 0.09 * diff(ylim)
# 			sites_range_y <- ylim[1] + 0.16 * diff(ylim)

# 			resp_curves_facet[[length(resp_curves_facet) + 1]] <- ggplot() +
# 				annotate(
# 					'segment',
# 					x = x_bg_range_sq[1], xend = x_bg_range_sq[2], y = bg_range_y, yend = bg_range_y,
# 					size = 4, color = 'gray50', alpha = 0.7
# 				) +
# 				annotate(
# 					'segment',
# 					x = x_occs_range[1], xend = x_occs_range[2], y = species_range_y, yend = species_range_y,
# 					size = 4, color = 'chartreuse3', alpha = 0.7
# 				) +
# 				annotate(
# 					'segment',
# 					x = x_sites_range[1], xend = x_sites_range[2], y = sites_range_y, yend = sites_range_y,
# 					size = 4, color = 'yellow', alpha = 0.7
# 				) +
# 				# geom_ribbon(data = pred_ci, aes(x = x, ymin = ymin, ymax = ymax, color = model), fill = NA, size = 0.4, alpha = 0.5, linetype = 'dotted') +
# 				geom_line(data = pred_centers, aes(x = x, y = y, color = model, linetype = model), size = center_line_size) +
# 				scale_linetype_manual(values = linetypes, name = 'Model') +
# 				xlab(xlab) +
# 				ylab(ylab) +
# 				labs(color = 'Model') +
# 				# ggtitle(title) +
# 				coord_cartesian(ylim = c(ylim[1], ylim[2] * 1.15), xlim = c(xlim[1], xlim[2])) +
# 				theme(
# 					plot.title = element_text(size = panel_title_size),
# 					legend.key.width = grid::unit(1.8, 'cm')
# 				)

# 		} # next covariate
# 		names(resp_curves_facet) <- covars
# 		resp_curves[[length(resp_curves) + 1]] <- resp_curves_facet
# 		names(resp_curves)[length(resp_curves)] <- facet

# 		# use plots + 1 legend panel
# 		plot_panels <- resp_curves_facet
# 		for (i in seq_along(plot_panels)) {
# 			plot_panels[[i]] <- plot_panels[[i]] + theme(legend.position = 'none')
# 		}

# 		legend_panel <- get_legend(
# 			resp_curves_facet[[1]] + theme(legend.position = 'right')
# 		)

# 		range_legend_panel <- get_legend(
# 			ggplot(
# 				data.frame(
# 				label = factor(c('Sites', 'Species', 'Background'), levels = c('Sites', 'Species', 'Background')),
# 				x = 1,
# 				y = 1:3
# 				),
# 				aes(x = x, y = y, fill = label)
# 			) +
# 			geom_tile() +
# 			scale_fill_manual(
# 				values = c(
# 					Sites = 'yellow',
# 					Species = 'chartreuse3',
# 					Background = 'grey50'
# 				),
# 				name = 'Environmental range'
# 			) +
# 			guides(fill = guide_legend(override.aes = list(color = NA))) +
# 			theme_void() +
# 			theme(legend.position = "right")
# 		)

# 		legend_stack <- plot_grid(
# 			legend_panel,
# 			range_legend_panel,
# 			# axis = 'trbl',
# 			ncol = 1,
# 			rel_heights = c(1, 0.5)
# 		)
		
# 		resps <- plot_grid(
# 			plotlist = c(plot_panels, list(legend_stack)),
# 			nrow = 1,
# 			ncol = length(plot_panels) + 1,
# 			rel_widths = c(rep(1, length(plot_panels)), 1.2),
# 			align = 'hv'
# 		)

# 		plot_title <- ggdraw() +
# 			draw_label(get_nice_trait(facet)$short, x = 0, hjust = 0) +
# 			theme(
# 				plot.margin = margin(0, 0, 0, 7)
# 			)
# 		resps <- plot_grid(plot_title, resps, ncol = 1, rel_heights = c(0.06, 1))

# 		width <- 4.6 * length(covars) + 3.2

# 		ggsave(resps, filename = paste0('./outputs_loretta/integrated_sdm_pdm/response_curves_', facet, '_univar_vs_integrated.png'), width = width, height = 4, units = 'in', dpi = 600, bg = 'white')

# 	} # next non-facet facet

# 	saveRDS(resp_curves, file = './outputs_loretta/integrated_sdm_pdm/response_curves_univar_vs_integrated_all_plot.rds')

# say('###########################################################')
# say('### maps of psi and response curves plots for main text ###')
# say('###########################################################')

# 	# panel plot:
# 	# row 1: psi from integrated model | difference between biomass with presencew forced and without | difference between stomatal conductance (?) with presence forced and without
# 	# response curve of psi
# 	# response curves of abundance
# 	# response curves of biomass
# 	# response curves of stomatal conductance

# 	# # Mexico/US/Canadian states/provinces
# 	# nam <- vect('./data_from_gadm/gadm_4pt1_level_1_north_america_sans_alaska_lambert.gpkg')
# 	# nam <- simplifyGeom(nam, tolerance = 1000)

# 	# # predictions from integrated morphology model
# 	# pred_vect_morph_pres_only <- vect('./outputs_loretta/integrated_sdm_pdm/models_integrated/non_centered_MVN([occs~hurdlepoisson_offset]_[biomass~hurdleln]_[bla_can_hei])_eta=2/prediction_vector_nam_force_presence.gpkg')
# 	# pred_vect_morph_pres_abs <- vect('./outputs_loretta/integrated_sdm_pdm/models_integrated/non_centered_MVN([occs~hurdlepoisson_offset]_[biomass~hurdleln]_[bla_can_hei])_eta=2/prediction_vector_nam.gpkg')

# 	# # predictions from integrated physiology model
# 	# pred_vect_phys_pres_only <- vect('./outputs_loretta/integrated_sdm_pdm/models_integrated/non_centered_MVN([occs~hurdlepoisson_offset]_[biomass~hurdleln]_[cnr_pho_spa_sto_tra])_eta=2/prediction_vector_nam_force_presence.gpkg')
# 	# pred_vect_phys_pres_abs <- vect('./outputs_loretta/integrated_sdm_pdm/models_integrated/non_centered_MVN([occs~hurdlepoisson_offset]_[biomass~hurdleln]_[cnr_pho_spa_sto_tra])_eta=2/prediction_vector_nam.gpkg')

# 	# # plot extent
# 	# sites <- readRDS('./data_from_loretta/sdm_pdm_00_merged_site_data_with_climate/sites.rds')
# 	# sites <- vect(sites, geom = c('LONGITUDE', 'LATITUDE'), crs = getCRS('WGS84'))
# 	# sites <- project(sites, pred_vect_morph_pres_abs)

# 	# extent <- ext(sites)
# 	# extent <- as.polygons(extent, crs = pred_vect_morph_pres_abs)
# 	# extent <- buffer(extent, width = 200 * 1000) # nominal plot extent
# 	# extent_display <- buffer(extent, width = 300 * 1000) # larger than plot extent
# 	# extent <- ext(extent)
# 	# extent <- as.vector(extent)

# 	# nam <- crop(nam, extent_display)
# 	# pred_vect_morph_pres_only <- crop(pred_vect_morph_pres_only, extent_display)
# 	# pred_vect_morph_pres_abs <- crop(pred_vect_morph_pres_abs, extent_display)
# 	# pred_vect_phys_pres_only <- crop(pred_vect_phys_pres_only, extent_display)
# 	# pred_vect_phys_pres_abs <- crop(pred_vect_phys_pres_abs, extent_display)

# 	# ### map of psi from integrated model
# 	# ####################################

# 	# 	palette <- colorRampPalette(c('#b2182b', '#ef8a62', '#fddbc7', '#f7f7f7', '#d1e5f0', '#67a9cf', '#2166ac'))
# 	# 	palette <- palette(100)
# 	# 	this_min <- min(pred_vect_morph_pres_abs$psi_sq)
# 	# 	this_max <- max(pred_vect_morph_pres_abs$psi_sq)
# 	# 	palette <- palette[round(100 * this_min):round(100 * this_max)]

# 	# 	legend_title <- expression(italic('ψ'[county]^(occs)))
# 	# 	map_psi <- ggplot() +
# 	# 		layer_spatial(pred_vect_morph_pres_abs, aes(fill = psi_sq), color = NA) +
# 	# 		scale_fill_gradientn(name = legend_title, colors = palette, limits = c(this_min, this_max)) +
# 	# 		layer_spatial(nam, color = 'gray40', fill = NA, linewidth = 0.3) +
# 	# 		layer_spatial(sites, pch = 4, size = 4, fill = 'white') +
# 	# 		xlim(extent[1], extent[2]) + ylim(extent[3], extent[4]) +
# 	# 		ggtitle('A) Probability of occurrence') +
# 	# 		theme(
# 	# 			plot.title = element_text(size = 22),
# 	# 			plot.subtitle = element_text(size = 14),
# 	# 			legend.title = element_text(size = 20),
# 	# 			legend.text = element_text(size = 20)
# 	# 		)

# 	# ### map of difference between biomass predicted by integrated model wiith presence/absence allowed and presence enforced
# 	# ########################################################################################################################

# 	# 	pred_vect_morph_pres_abs$delta_biomass <- pred_vect_morph_pres_abs$biomass_mean_sq - pred_vect_morph_pres_only$biomass_mean_sq
# 	# 	palette <- c('#FFFDFD', '#FEE0D2', '#FC9272', '#DE2D26', '#67000D')

# 	# 	legend_title <- 'Difference\nin biomass (g)'
# 	# 	map_biomass <- ggplot() +
# 	# 		layer_spatial(pred_vect_morph_pres_abs, aes(fill = delta_biomass), color = NA) +
# 	# 		scale_fill_gradientn(name = legend_title, colors = palette) +
# 	# 		layer_spatial(nam, color = 'gray40', fill = NA, linewidth = 0.3) +
# 	# 		layer_spatial(sites, pch = 4, size = 4, fill = 'white') +
# 	# 		xlim(extent[1], extent[2]) + ylim(extent[3], extent[4]) +
# 	# 		ggtitle('B) Difference in biomass') +
# 	# 		theme(
# 	# 			plot.title = element_text(size = 22),
# 	# 			plot.subtitle = element_text(size = 14),
# 	# 			legend.title = element_text(size = 20),
# 	# 			legend.text = element_text(size = 20)
# 	# 		)

# 	# ### map of difference between stomatal conductance predicted by integrated model wiith presence/absence allowed and presence enforced
# 	# ########################################################################################################################

# 	# 	pred_vect_phys_pres_abs$delta_sc <- pred_vect_phys_pres_abs$stomatal_conductance_mean_sq - pred_vect_phys_pres_only$stomatal_conductance_mean_sq
# 	# 	palette <- c('#FFFDFD', '#FEE0D2', '#FC9272', '#DE2D26', '#67000D')

# 	# 	legend_title <- 'Difference\nin stom. con.\n(mol m⁻² s⁻¹)'
# 	# 	map_sc <- ggplot() +
# 	# 		layer_spatial(pred_vect_phys_pres_abs, aes(fill = delta_sc), color = NA) +
# 	# 		scale_fill_gradientn(name = legend_title, colors = palette) +
# 	# 		layer_spatial(nam, color = 'gray40', fill = NA, linewidth = 0.3) +
# 	# 		layer_spatial(sites, pch = 4, size = 4, fill = 'white') +
# 	# 		xlim(extent[1], extent[2]) + ylim(extent[3], extent[4]) +
# 	# 		ggtitle('B) Difference in stomatal conductance') +
# 	# 		theme(
# 	# 			plot.title = element_text(size = 22),
# 	# 			plot.subtitle = element_text(size = 14),
# 	# 			legend.title = element_text(size = 20),
# 	# 			legend.text = element_text(size = 20)
# 	# 		)

# 	# ### combine maps
# 	# ################

# 	# 	maps <- plot_grid(map_psi, map_biomass, map_sc, nrow = 1, align = 'hv')


# 	### response curves GENERIC
# 	###########################

# 		panel_title_size <- 9 # size of panel plot titles
# 		center_line_size <- 0.7 # thickness of lines for central tendency of each facet
# 		ymax_mult <- 1.05 # set ymax in plots equal to highest central tendency prediction times this value
# 		env_range_label_size <- 2 # size of text labeling environmental ranges
# 		axis_title_size <- 8
# 		axis_text_size <- 7
# 		legend_title_size <- 8
# 		legend_text_size <- 7
# 		x_range_width <- 3 # width of line segements showog range of sites, AG, BG

# 	### response curves for ABUNDANCE
# 	#################################

# 		precalc_resp <- readRDS('./outputs_loretta/integrated_sdm_pdm/response_curves_precalculated_abundance_univar_vs_integrated.rds')

# 		resp_curves_abund <- list()
# 		y_max_across_covars <- -Inf
# 		for (i in seq_along(precalc_resp$responses)) {

# 			covar <- names(precalc_resp$responses)[i]
# 			nice_pred <- get_nice_predictor(covar)

# 			# x-axis values
# 			x <- precalc_resp$responses[[covar]]$x

# 			# responses to plot
# 			center_multi_all_varying <- precalc_resp$response[[covar]]$center_multi_all_varying
# 			center_multi_all_varying_no_corr <- precalc_resp$response[[covar]]$center_multi_all_varying_no_corr
# 			center_multi_all_varying_force_presence <- precalc_resp$response[[covar]]$center_multi_all_varying_force_presence
			
# 			center_uni_all_varying <- precalc_resp$response[[covar]]$center_uni_all_varying
# 			center_uni_all_varying_force_presence <- precalc_resp$response[[covar]]$center_uni_all_varying_force_presence

# 			n <- length(x)
# 			df <- data.table(
# 				x = rep(x, 5),
# 				y = c(
# 					center_multi_all_varying,
# 					center_multi_all_varying_no_corr,
# 					center_multi_all_varying_force_presence,
# 					center_uni_all_varying,
# 					center_uni_all_varying_force_presence
# 				),
# 				experiment = c(
# 					rep('Integrated', n),
# 					rep('Integrated, no facet-facet correlations', n),
# 					rep('Integrated, forcing presence', n),
# 					rep('Abundance-only', n),
# 					rep('Abundance-only, forcing presence', n)
# 				)
# 			)

# 			xlim <- range(x)

# 			# calculate y limits... do not include an treatments that had infinite-predicted abundance
# 			if (any(is.infinite(df$y))) {

# 				is_inf <- which(is.infinite(df$y))
# 				inf_experiments <- df$experiment[is_inf]
# 				inf_experiments <- unique(inf_experiments)

# 				df_trunc <- df[df$experiment %notin% inf_experiments]
# 				y_max_across_covars <- max(y_max_across_covars, df_trunc$y)

# 			} else {
# 				y_max_across_covars <- max(y_max_across_covars, df$y)
# 			}

# 			resp_curves_abund[[covar]] <- ggplot(df, aes(x = x, y = y, color = experiment, linetype = experiment)) +
# 				geom_line(linewidth = center_line_size) +
# 				scale_color_manual(
# 					values = c(
# 						'Integrated' = 'black',
# 						'Integrated, no facet-facet correlations' = '#4C8CBE',
# 						'Integrated, forcing presence' = '#2B5C8F',
# 						'Abundance-only' = 'darkorange3',
# 						'Abundance-only, forcing presence' = 'darkorange'
# 					)
# 				) +
# 				scale_linetype_manual(
# 					values = c(
# 						'Integrated' = 'solid',
# 						'Integrated, no facet-facet correlations' = 'solid',
# 						'Integrated, forcing presence' = 'solid',
# 						'Abundance-only' = 'solid',
# 						'Abundance-only, forcing presence' = 'solid'
# 					)
# 				) +
# 				xlim(xlim[1], xlim[2]) +
# 				labs(
# 					color = 'Model/treatment',
# 					linetype = 'Model/treatment',
# 					x = nice_pred$long,
# 					y = 'Relative abundance'
# 				) +
# 				theme(
# 					legend.title = element_text(size = legend_title_size),
# 					legend.text = element_text(size = legend_text_size),
# 					legend.margin = margin(t = 0, r = 0, b = 10, l = 0)
# 				)

# 		} # next abundance covariate

# 		ylim <- c(0, y_max_across_covars)

# 		# use plots + 1 legend panel
# 		plot_panels <- resp_curves_abund
# 		for (i in seq_along(plot_panels)) {

# 			covar <- names(plot_panels)[i]

# 			# range of sites, occurrences, and BG
# 			x_sites_range <- precalc_resp$ranges[[covar]]$x_sites_range
# 			x_occs_range <- precalc_resp$ranges[[covar]]$x_occs_range
# 			x_bg_range_sq <- precalc_resp$ranges[[covar]]$x_bg_range_sq
# 			x_bg_range_fut = precalc_resp$ranges[[covar]]$x_bg_range_ssp370_2071_2100

# 			bg_sq_range_y <- ylim[1] + 0.005 * diff(ylim)
# 			bg_fut_range_y <- ylim[1] + 0.06 * diff(ylim)
# 			species_range_y <- ylim[1] + 0.12 * diff(ylim)
# 			sites_range_y <- ylim[1] + 0.18 * diff(ylim)

# 			plot_panels[[i]] <- plot_panels[[i]] +
# 				annotate(
# 					'segment',
# 					x = x_bg_range_sq[1], xend = x_bg_range_sq[2], y = bg_sq_range_y, yend = bg_sq_range_y,
# 					linewidth = x_range_width, color = 'gray50', alpha = 0.7
# 				) +
# 				annotate(
# 					'text',
# 					label = 'background - present',
# 					x = mean(x_bg_range_sq), y = bg_sq_range_y,
# 					color = 'black', size = env_range_label_size
# 				) +
# 				annotate(
# 					'segment',
# 					x = x_bg_range_fut[1], xend = x_bg_range_fut[2], y = bg_fut_range_y, yend = bg_fut_range_y,
# 					linewidth = x_range_width, color = 'gray50', alpha = 0.7
# 				) +
# 				annotate(
# 					'text',
# 					label = 'background - SSP370, 2080s',
# 					x = mean(x_bg_range_fut), y = bg_fut_range_y,
# 					color = 'black', size = env_range_label_size
# 				) +
# 				annotate(
# 					'segment',
# 					x = x_occs_range[1], xend = x_occs_range[2], y = species_range_y, yend = species_range_y,
# 					linewidth = x_range_width, color = 'chartreuse3', alpha = 0.7
# 				) +
# 				annotate(
# 					'text',
# 					label = 'occurrences',
# 					x = mean(x_occs_range), y = species_range_y,
# 					color = 'black', size = env_range_label_size
# 				) +
# 				annotate(
# 					'segment',
# 					x = x_sites_range[1], xend = x_sites_range[2], y = sites_range_y, yend = sites_range_y,
# 					linewidth = x_range_width, color = 'yellow', alpha = 0.7
# 				) +
# 				annotate(
# 					'text',
# 					label = 'sites',
# 					x = mean(x_sites_range), y = sites_range_y,
# 					color = 'black', size = env_range_label_size
# 				) +
# 				coord_cartesian(ylim = c(ylim[1], ylim[2])) +
# 				theme(
# 					legend.position = 'none',
# 					axis.title = element_text(size = axis_title_size),
# 					axis.text = element_text(size = axis_text_size)
# 				)

# 		}

# 		legend_panel <- get_legend(
# 			resp_curves_abund[[1]] + theme(legend.position = 'right')
# 		)

# 		range_legend_panel <- get_legend(
# 			ggplot(
# 				data.frame(
# 					label = factor(c('Sites', 'Species', 'Background (present or future)'), levels = c('Sites', 'Species', 'Background (present or future)')),
# 					x = 1,
# 					y = 1:3
# 				),
# 				aes(x = x, y = y, fill = label)
# 			) +
# 			geom_tile() +
# 			scale_fill_manual(
# 				values = c(
# 					Sites = 'yellow',
# 					Species = 'chartreuse3',
# 					Background = 'grey50'
# 				),
# 				name = 'Environmental range'
# 			) +
# 			guides(fill = guide_legend(override.aes = list(color = NA))) +
# 			theme_void() +
# 			theme(
# 				legend.position = 'right',
# 				legend.title = element_text(size = legend_title_size),
# 				legend.text = element_text(size = legend_text_size),
# 				legend.margin = margin(t = 10, r = 0, b = 0, l = 0)
# 			)
# 		)

# 		legend_stack <- plot_grid(
# 			legend_panel,
# 			range_legend_panel,
# 			# axis = 'trbl',
# 			ncol = 1,
# 			align = 'v',
# 			rel_heights = c(1, 0.7)
# 		)
		
# 		resps <- plot_grid(
# 			plotlist = c(plot_panels, list(legend_stack)),
# 			nrow = 1,
# 			ncol = length(plot_panels) + 1,
# 			rel_widths = c(rep(1, length(plot_panels)), 1.2),
# 			align = 'hv'
# 		)

# 		plot_title <- ggdraw() +
# 			draw_label('a) Abundance', x = 0, hjust = 0, size = panel_title_size) +
# 			theme(
# 				plot.margin = margin(0, 0, 0, 7)
# 			)
# 		resps <- plot_grid(plot_title, resps, ncol = 1, rel_heights = c(0.06, 1))

# 	ggsave(resps, filename = 'C:/!scratch/_test.png', bg = 'white', dpi = 600, width = 8, height = 2.8)


# say('#########################################################################################')
# say('### calculate correlations between facets within latent multivariate normal component ###')
# say('#########################################################################################')

# 	# calculate the correlation matrix between axes of the latent multivariate normal component in the integrated model

# 	model_dir <- './outputs_loretta/integrated_sdm_pdm/models_integrated/non_centered_MVN([occs~hurdlepoisson_offset]_[biomass~hurdleln]_[bla_can_hei])_eta=2'
# 	n_facets <- 5

# 	facets <- c('abundance', 'biomass', 'blade width', 'canopy diameter', 'height')

# 	chains <- readRDS(paste0(model_dir, '/chains.rds'))

# 	# extract Cholesky decomposition
# 	U_star_mean <- mc_extract(chains, 'U_star', j = TRUE, k = TRUE)
# 	U_star_mean <- matrix(U_star_mean, nrow = n_facets, ncol = n_facets)

# 	corr_mean <- t(U_star_mean) %*% U_star_mean
# 	colnames(corr_mean) <- rownames(corr_mean) <- facets

# 	# calculate correlations across iterations
# 	indices <- which(grepl(colnames(chains$samples[[1]]), pattern = 'U_star'))
# 	corr <- list()
# 	for (j in 1:4) {
# 		for (i in 1:1000) {

# 			U_star <- chains$samples[[j]][i, indices]
# 			U_star <- matrix(U_star, nrow = n_facets, ncol = n_facets)
# 			this_corr <- t(U_star) %*% U_star

# 			if (any(this_corr > 1 + eps()) | any(this_corr < -1 - eps())) stop()
# 			corr[[length(corr) + 1]] <- this_corr

# 		}
# 	}

# 	# mean correlation
# 	corr_mean <- apply(corr_array, c(1, 2), mean, na.rm = TRUE)
# 	colnames(corr_mean) <- rownames(corr_mean) <- facets

# 	# element-wise 10th and 90th quantiles across all matrices in corr
# 	corr_quantiles <- array(
# 		NA_real_,
# 		dim = c(n_facets, n_facets, 2),
# 		dimnames = list(facets, facets, c('q10', 'q90'))
# 	)
# 	corr_array <- simplify2array(corr)
# 	for (j in seq_len(n_facets)) {
# 		for (k in seq_len(n_facets)) {
# 			corr_quantiles[j, k, ] <- quantile(
# 				corr_array[j, k, ],
# 				probs = c(0.1, 0.9),
# 				na.rm = TRUE,
# 				names = FALSE
# 			)
# 		}
# 	}

# 	# calculate correlations between observed site-level means of each non-abundance facet
# 	obs <- read_xlsx('./data_from_loretta/!plant_sitelevel_data_24OCT2024 - Google Sheets [aggregated by Erica] before adding metadata 2024-12-28.xlsx', sheet = 'plantmaster_bysite_19NOV2024')
# 	obs <- as.data.table(obs)

# 	cols <- c('mBiomass', 'mBladeWidth', 'mCanopyDiam', 'mHeight')
# 	corr_obs <- cor(obs[ , ..cols], method = 'pearson')
	


say(date())
say('FINIS!', deco = '+', level = 1)

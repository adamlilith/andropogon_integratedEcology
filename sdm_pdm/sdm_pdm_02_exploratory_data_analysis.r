### MODELING ANDROPOGON GERARDI DISTRIBUTION, PHENOTYPE, PHYSIOLOGY, GENOTYPE, and ASSOCIATED MICROBIAL COMMUNITIES
### Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2023-12
###
### This script constructs a non-integrated model for AG geographic distribution.
###
### source('C:/Kaji/R/andropogon_integratedEcology/sdm_pdm/sdm_pdm_02_exploratory_data_analysis.r')
###
### CONTENTS ###
### setup ###
### trait histograms across sites ###
### trait histograms within sites ###
### correlations between environmental variables at phenotypic sample sites and occurrence sites ###
### univariate plots of each phenotypic trait by covariate ###
### pre-screen simple models for each trait ###
### pre-screen simple models of relationship between environment and occurrence ###
### plots of occurrence and traits in environmental space ###

#############
### setup ###
#############

	rm(list = ls())

	drive <- 'C:/Kaji/'
	
	setwd(paste0(drive, '/Research/Andropogon/Andropogon'))
	source(paste0(drive, '/R/andropogon_integratedEcology/sdm_pdm/sdm_pdm_00_shared_functions_and_variables.r'))

	library(AICcmodavg) # AICc
	library(cluster) # clustering
	library(cowplot) # combining ggplots
	library(DHARMa) # residuals analysis
	library(data.table) # fast data tables
	library(enmSdmX) # GIS & SDMing
	library(ggdendro) # dendrogram plots
	library(ggplot2) # plots
	library(ggspatial) # spatial plots
	library(grid) # tables
	library(gridExtra) # tables
	library(lme4) # mixed linear models
	library(MuMIn) # multi-model selection
	library(omnibus) # utilities
	library(patchwork) # combining ggplots
	library(predicts) # GIS & SDMing
	library(readxl) # Excel
	library(terra) # spatial objects
	library(tidyr) # data wrangling

# say('#####################################')
# say('### trait histograms across sites ###')
# say('#####################################')

# 	dirCreate('./outputs_loretta/integrated_sdm_pdm')

# 	sites <- readRDS('./data_from_loretta/sdm_pdm_00_merged_site_data_with_climate/sites.rds')
# 	biomass <- readRDS('./data_from_loretta/sdm_pdm_00_merged_site_data_with_climate/biomass.rds')
# 	morpho_phys <- readRDS('./data_from_loretta/sdm_pdm_00_merged_site_data_with_climate/morpho_phys.rds')

# 	### biomass variables
# 	#####################

# 	x <- biomass
# 	x <- x[ , c('SITE', 'VegBiomass', 'ReproBiomass', 'Biomass', 'VegRepro_ratio')]
# 	x <- x[complete.cases(x)]
	
	# long <- pivot_longer(
	# 	x,
	# 	cols = all_of(traits),
	# 	names_to = 'trait',
	# 	values_to = 'value'
	# )
	# long <- as.data.table(long)


# 	# plot
# 	plots <- ggplot(long, aes(x = value)) +
# 		geom_histogram(bins = 30) +
# 		facet_wrap(~ trait, scales = 'free') +
# 		theme_minimal() +
# 		theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
# 		labs(title = 'Biomass trait histograms')

# 	ggsave(plots, filename = './outputs_loretta/integrated_sdm_pdm/biomass_trait_histograms.png', width = 8, height = 8, dpi = 300, bg = 'white')

# 	### morphology/physiology variables
# 	###################################

# 	x <- morpho_phys
# 	x <- x[ , c('SITE', 'Delta13C', 'N_conc', 'CN_ratio', 'Height', 'BladeWidth', 'LeafThick', 'SPAD', 'CanopyDiam', 'WatPot', 'PhotoRate', 'StomCond', 'IntCO2', 'TranspRate')]

	# long <- pivot_longer(
	# 	x,
	# 	cols = all_of(traits),
	# 	names_to = 'trait',
	# 	values_to = 'value'
	# )
	# long <- as.data.table(long)

# 	# plot
# 	plots <- ggplot(long, aes(x = value)) +
# 		geom_histogram(bins = 30) +
# 		facet_wrap(~ trait, scales = 'free') +
# 		theme_minimal() +
# 		theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
# 		labs(title = 'Morphology/physiology trait histograms')

# 	ggsave(plots, filename = './outputs_loretta/integrated_sdm_pdm/morphology_trait_histograms.png', width = 16, height = 12, dpi = 300, bg = 'white')

# say('#####################################')
# say('### trait histograms within sites ###')
# say('#####################################')

	# dirCreate('./outputs_loretta/integrated_sdm_pdm')

	# sites <- readRDS('./data_from_loretta/sdm_pdm_00_merged_site_data_with_climate/sites.rds')
	# biomass <- readRDS('./data_from_loretta/sdm_pdm_00_merged_site_data_with_climate/biomass.rds')
	# morpho_phys <- readRDS('./data_from_loretta/sdm_pdm_00_merged_site_data_with_climate/morpho_phys.rds')

	# ### biomass variables
	# #####################

	# traits <- c('VegBiomass', 'ReproBiomass', 'Biomass', 'VegRepro_ratio')

	# x <- biomass
	# x <- x[ , c('SITE', ..traits)]
	# x <- x[complete.cases(x)]
	
	# long <- pivot_longer(
		# x,
		# cols = all_of(traits),
		# names_to = 'trait',
		# values_to = 'value'
	# )
	# long <- as.data.table(long)

	# for (trait in unique(long$trait)) {

		# keeps <- long$trait == trait
		# this_long <- long[keeps]
		# site_means <- this_long[ , .(mean_value = mean(value, na.rm = TRUE)), by = SITE]
		# site_means <- site_means[order(mean_value)]
		# this_long$SITE <- factor(this_long$SITE, levels = site_means$SITE)
		# this_long <- this_long[order(this_long$SITE)]

		# # plot
		# plots <- ggplot(this_long, aes(x = value, y = 0)) +
			# geom_point() +
			# facet_wrap(~ SITE) +
			# xlab(trait) +
			# theme(
				# axis.text.x = element_text(angle = 45, hjust = 1),
				# axis.title.y = element_blank(),
				# axis.text.y = element_blank(),
				# axis.ticks.y = element_blank()
			# )
			
		# ggsave(plots, filename = paste0('./outputs_loretta/integrated_sdm_pdm/trait_histograms_by_site_', tolower(trait), '.png'), width = 8, height = 8, dpi = 300, bg = 'white')

	# }

	# ### morphology/physiology variables
	# ###################################

	# traits <- c('Delta13C', 'N_conc', 'CN_ratio', 'Height', 'BladeWidth', 'LeafThick', 'SPAD', 'CanopyDiam', 'WatPot', 'PhotoRate', 'StomCond', 'IntCO2', 'TranspRate')

	# x <- morpho_phys
	# x <- x[ , c('SITE', ..traits)]
	
	# long <- pivot_longer(
		# x,
		# cols = all_of(traits),
		# names_to = 'trait',
		# values_to = 'value'
	# )
	# long <- as.data.table(long)

	# for (trait in traits) {

		# keeps <- long$trait == trait
		# this_long <- long[keeps]
		# site_means <- this_long[ , .(mean_value = mean(value, na.rm = TRUE)), by = SITE]
		# site_means <- site_means[order(mean_value)]
		# this_long$SITE <- factor(this_long$SITE, levels = site_means$SITE)
		# this_long <- this_long[order(this_long$SITE)]

		# # plot
		# plots <- ggplot(this_long, aes(x = value, y = 0)) +
			# geom_point() +
			# facet_wrap(~ SITE) +
			# xlab(trait) +
			# theme(
				# axis.text.x = element_text(angle = 45, hjust = 1),
				# axis.title.y = element_blank(),
				# axis.text.y = element_blank(),
				# axis.ticks.y = element_blank()
			# )

		# ggsave(plots, filename = paste0('./outputs_loretta/integrated_sdm_pdm/trait_histograms_by_site_', tolower(trait), '.png'), width = 8, height = 8, dpi = 300, bg = 'white')

	# }

say('####################################################################################################')
say('### correlations between environmental variables at phenotypic sample sites and occurrence sites ###')
say('####################################################################################################')

	dirCreate('./outputs_loretta/integrated_sdm_pdm')

	### sites

		sites <- readRDS('./data_from_loretta/sdm_pdm_00_merged_site_data_with_climate/sites.rds')
		vars <- c(paste0('bio', c(1, 5, 6, 7, 10, 12, 15, 18)), 'aridity', 'field_pH', 'SAND', 'SILT', 'CLAY')

		cors <- cor(sites[ , ..vars], method = 'spearman')
		# rownames(cors) <- colnames(cors) <- c(paste0('BIO', c(1, 5, 6, 7, 10, 12, 15, 18)), 'Aridity', 'pH', 'Sand', 'Silt', 'Clay')
		rownames(cors) <- colnames(cors) <- c('Mean Temp. (BIO1)', 'Max. Temp (BIO5)', 'Min. Temp (BIO6)', 'Temp. Range (BIO7)', 'Summer Temp (BIO10)', 'Annual Precip. (BIO12)', 'Precip. Variability (BIO15)', 'Summer Precip. (BIO18)', 'Aridity', 'pH', 'Sand', 'Silt', 'Clay')
		dists <- 1 - abs(cors)
		dists <- as.dist(dists)

		sites_cluster <- hclust(dists)

		dendro <- as.dendrogram(sites_cluster)
		sites_dendro <- dendro_data(dendro)

	### occurrences

		occs <- vect('./outputs_loretta/integrated_sdm_pdm/andropogon_gerardi_occurrences_with_environment_1961_2020_for_integration.gpkg')
		occs <- occs[occs$n_poaceae > 0]
		# occs <- occs[occs$n_andropogon_gerardi > 0]
		occs <- as.data.table(occs)

		vars <- c(paste0('bio', c(1, 5, 6, 7, 10, 12, 15, 18)), 'aridity', 'ph', 'sand', 'silt', 'clay')

		cors <- cor(occs[ , ..vars], method = 'spearman')
		# rownames(cors) <- colnames(cors) <- c(paste0('BIO', c(1, 5, 6, 7, 10, 12, 15, 18)), 'Aridity', 'pH', 'Sand', 'Silt', 'Clay')
		rownames(cors) <- colnames(cors) <- c('Mean Temp. (BIO1)', 'Max. Temp (BIO5)', 'Min. Temp (BIO6)', 'Temp. Range (BIO7)', 'Summer Temp (BIO10)', 'Annual Precip. (BIO12)', 'Precip. Variability (BIO15)', 'Summer Precip. (BIO18)', 'Aridity', 'pH', 'Sand', 'Silt', 'Clay')
		dists <- 1 - abs(cors)
		dists <- as.dist(dists)

		occs_cluster <- hclust(dists)

		dendro <- as.dendrogram(occs_cluster)
		occs_dendro <- dendro_data(dendro)

	### plot

        # Define consistent colors for labels
        label_names <- c('Mean Temp. (BIO1)', 'Max. Temp (BIO5)', 'Min. Temp (BIO6)', 'Temp. Range (BIO7)', 'Summer Temp (BIO10)', 'Annual Precip. (BIO12)', 'Precip. Variability (BIO15)', 'Summer Precip. (BIO18)', 'Aridity', 'pH', 'Sand', 'Silt', 'Clay')
        label_colors <- setNames(scales::hue_pal()(length(label_names)), label_names)

        sites_labels <- label(sites_dendro)
        sites_labels$color <- label_colors[as.character(sites_labels$label)]

        occs_labels <- label(occs_dendro)
        occs_labels$color <- label_colors[as.character(occs_labels$label)]

        sites_plot <- ggplot(segment(sites_dendro)) +
            geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
			geom_text(
				data = sites_labels,
				aes(x = x, y = y, label = label, color = label),
				hjust = 1.1, vjust = 0.5, angle = 90
			) +
            scale_color_manual(values = label_colors, guide = "none") +
            geom_hline(yintercept = 0.3, color = 'red', linetype = 'dashed') +
            coord_cartesian(clip = 'off') +
            theme_minimal() +
            theme(
                panel.grid = element_blank(),
                axis.text.x = element_blank(),
                plot.margin = margin(t = 5, r = 5, b = 100, l = 5)
            ) +
            labs(
                title = 'a) Correlations at field sites',
                y = "1 - |ρ|",
                x = NULL
            )

        occs_plot <- ggplot(segment(occs_dendro)) +
            geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
            geom_text(
                data = occs_labels,
                aes(x = x, y = y, label = label, color = label),
                hjust = 1.1, vjust = 0.5, size = 3, angle = 90
            ) +
            scale_color_manual(values = label_colors, guide = "none") +
            geom_hline(yintercept = 0.3, color = 'red', linetype = 'dashed') +
            coord_cartesian(clip = 'off') +
            theme_minimal() +
            theme(
                panel.grid = element_blank(),
                axis.text.x = element_blank(),
                plot.margin = margin(t = 5, r = 5, b = 100, l = 5)
            ) +
            labs(
                title = 'b) Correlations across SDM calibration sites',
                y = "1 - |ρ|",
                x = NULL
            )

		combo <- sites_plot + occs_plot
		ggsave(combo, filename = './outputs_loretta/integrated_sdm_pdm/correlations_between_covariates_sites_and_counties_with_poaceae.svg', height = 5, width = 12, dpi = 600)

# say('##############################################################')
# say('### univariate plots of each phenotypic trait by covariate ###')
# say('##############################################################')

# 	out_dir <- './outputs_loretta/integrated_sdm_pdm/univariate_plots_of_traits'
# 	dirCreate(out_dir)

# 	sites <- readRDS('./data_from_loretta/sdm_pdm_00_merged_site_data_with_climate/sites.rds')
# 	vars <- c(paste0('bio', c(1, 6, 7, 12, 15, 18)), 'aridity', 'field_pH', 'SAND', 'SILT', 'CLAY')

# 	biomass <- readRDS('./data_from_loretta/sdm_pdm_00_merged_site_data_with_climate/biomass.rds')
# 	morpho_phys <- readRDS('./data_from_loretta/sdm_pdm_00_merged_site_data_with_climate/morpho_phys.rds')

# 	biomass <- merge(sites[ , c('site_id', 'LONGITUDE', 'LATITUDE', 'field_pH', 'SAND', 'SILT', 'CLAY')], biomass, by.x = 'site_id', by.y = 'SITE')
# 	morpho_phys <- merge(sites[ , c('site_id', 'LONGITUDE', 'LATITUDE', 'field_pH', 'SAND', 'SILT', 'CLAY')], morpho_phys, by.x = 'site_id', by.y = 'SITE')

# 	traits <- c('Biomass', 'Delta13C', 'N_conc', 'CN_ratio', 'Height', 'BladeWidth', 'LeafThick', 'SPAD', 'CanopyDiam', 'WatPot', 'PhotoRate', 'StomCond', 'IntCO2', 'TranspRate')

# 	for (trait in traits) {

# 		say(trait)

# 		if (trait == 'Biomass') {
# 			x <- biomass
# 		} else {
# 			x <- morpho_phys
# 		}
	
# 		### plot univariate relationships between VALUES of each trait and environmental variables	
# 		x_tall <- data.table()
# 		for (i in seq_along(vars)) {
		
# 			var <- vars[i]
# 			x_tall <- rbind(
# 				x_tall,
# 				data.table(
# 					site_id = x$site_id,
# 					var = var,
# 					x = x[[var]],
# 					y = x[[trait]]
# 				)
# 			)

# 		}

# 		# Calculate R2 and P values for each variable
# 		stats_linear <- x_tall[ , .(
# 			r2 = summary(lm(y ~ x))$r.squared,
# 			p_value = summary(lm(y ~ x))$coefficients[2, 4],
# 			aicc = AICc(lm(y ~ x))
# 		), by = var]

# 		stats_quad <- x_tall[ , .(
# 			r2 = summary(lm(y ~ x + I(x^2)))$r.squared,
# 			p_value = summary(lm(y ~ x + I(x^2)))$coefficients[2, 4],
# 			aicc = AICc(lm(y ~ x + I(x^2)))
# 		), by = var]

# 		min_aicc <- min(stats_linear$aicc, stats_quad$aicc)
# 		stats_linear[ , delta_aicc := aicc - min_aicc]
# 		stats_quad[ , delta_aicc := aicc - min_aicc]

# 		raws <- ggplot(x_tall, aes(x = x, y = y, color = site_id)) +
# 			geom_point(pch = 21, fill = 'black') +
# 			geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = 'blue') +
# 			geom_smooth(method = 'lm', formula = y ~ x + I(x^2), se = FALSE, color = 'red') +
# 			facet_wrap(~ var, scales = 'free') +
# 			labs(title = trait, x = var, y = trait) +
# 			theme(
# 				legend.position = 'none'
# 			) +
# 			geom_text(
# 				data = stats_linear,
# 				aes(
# 					x = -Inf, y = Inf,
# 					label = paste0('∆AICc = ', round(delta_aicc, 2), '\nR² = ', round(r2, 2), '\nP = ', signif(p_value, 3))
# 				),
# 				hjust = -0.1, vjust = 1.1,
# 				color = 'blue', size = 3,
# 				inherit.aes = FALSE
# 			) +
# 			geom_text(
# 				data = stats_quad,
# 				aes(
# 					x = Inf, y = Inf,
# 					label = paste0('∆AICc = ', round(delta_aicc, 2), '\nR² = ', round(r2, 2), '\nP = ', signif(p_value, 3))
# 				),
# 				hjust = 1.1, vjust = 1.1,
# 				color = 'red', size = 3,
# 				inherit.aes = FALSE
# 			)

# 		### plot univariate relationships between SD of each trait and environmental variables	
# 		x_tall_sd <- x_tall[ , .(sd_y = sd(y, na.rm = TRUE), x = mean(x, na.rm = TRUE)), by = .(var, site_id)]

# 		# Calculate R2 and P values for each variable
# 		stats_linear <- x_tall_sd[ , .(
# 			r2 = summary(lm(sd_y ~ x))$r.squared,
# 			p_value = summary(lm(sd_y ~ x))$coefficients[2, 4],
# 			aicc = AICc(lm(sd_y ~ x))
# 		), by = var]

# 		stats_quad <- x_tall_sd[ , .(
# 			r2 = summary(lm(sd_y ~ x + I(x^2)))$r.squared,
# 			p_value = summary(lm(sd_y ~ x + I(x^2)))$coefficients[2, 4],
# 			aicc = AICc(lm(sd_y ~ x + I(x^2)))
# 		), by = var]

# 		min_aicc <- min(stats_linear$aicc, stats_quad$aicc)
# 		stats_linear[ , delta_aicc := aicc - min_aicc]
# 		stats_quad[ , delta_aicc := aicc - min_aicc]

# 		# Add R2 and P values as text annotations
# 		sds <- ggplot(x_tall_sd, aes(x = x, y = sd_y)) +
# 			geom_point() +
# 			geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = 'blue') +
# 			geom_smooth(method = 'lm', formula = y ~ x + I(x^2), se = FALSE, color = 'red') +
# 			facet_wrap(~ var, scales = 'free') +
# 			labs(title = paste0('SD of ', trait), x = var, y = paste0('SD of ', trait)) +
# 			theme(
# 				legend.position = 'none'
# 			) +
# 			geom_text(
# 				data = stats_linear,
# 				aes(
# 					x = -Inf, y = Inf,
# 					label = paste0('∆AICc = ', round(delta_aicc, 2), '\nR² = ', round(r2, 2), '\nP = ', signif(p_value, 3))
# 				),
# 				hjust = -0.1, vjust = 1.1,
# 				color = 'blue', size = 3,
# 				inherit.aes = FALSE
# 			) +
# 			geom_text(
# 				data = stats_quad,
# 				aes(
# 					x = Inf, y = Inf,
# 					label = paste0('∆AICc = ', round(delta_aicc, 2), '\nR² = ', round(r2, 2), '\nP = ', signif(p_value, 3))
# 				),
# 				hjust = 1.1, vjust = 1.1,
# 				color = 'red', size = 3,
# 				inherit.aes = FALSE
# 			)

# 		plots <- plot_grid(raws, sds, ncol = 2)

# 		ggsave(plots, filename = paste0(out_dir, '/', tolower(trait), '.png'), width = 16, height = 8, dpi = 300, bg = 'white')

	
# 	}

# say('###############################################')
# say('### pre-screen simple models for each trait ###')
# say('###############################################')

# 	# Evaluate univariate (linear and linear-plus-quadratic) and bivariate (linear and linear-plus-quadratic) models for each trait. Model fit and residuals are evaluated, and models ranked by AICc.

# 	# variables to evaluate
# 	vars <- c('bio1', 'bio12', 'bio18', 'aridity', 'field_pH', 'SAND', 'SILT', 'CLAY')
# 	twos <- combn(vars, 2, simplify = FALSE)

# 	### all possible combinations of two-variable models
# 	# remove combinations with bio12 and 18 bc thematically redundant
# 	bads <- integer()
# 	for (i in seq_along(twos)) {
# 		this_vars <- twos[[i]]
# 		# bad <- if (all(c('bio1', 'bio12') %in% this_vars) | all(c('bio12', 'bio18') %in% this_vars)) {
# 		bad <- if (all(c('bio12', 'bio18') %in% this_vars)) {
# 			i
# 		} else {
# 			NULL
# 		}
# 		bads <- c(bads, bad)
# 	}

# 	twos <- twos[-bads]

# 	forms <- vars
# 	forms <- c(forms, paste0(vars, ' + I(', vars, '^2)'))
# 	for (i in seq_along(twos)) {

# 		forms <- c(forms, paste0(twos[[i]], collapse = ' + '))	
# 		forms <- c(forms, paste0(c(twos[[i]], paste0('I(', twos[[i]][1], '^2)')), collapse = ' + '))
# 		forms <- c(forms, paste0(c(twos[[i]], paste0('I(', twos[[i]][2], '^2)')), collapse = ' + '))	
# 		forms <- c(forms, paste0(c(twos[[i]], paste0('I(', twos[[i]][1], '^2)'), paste0('I(', twos[[i]][2], '^2)')), collapse = ' + '))	

# 	}

# 	forms <- c(forms, '1')

# 	### setup
# 	out_dir <- './outputs_loretta/integrated_sdm_pdm'
# 	dirCreate(out_dir)

# 	# data
# 	sites <- readRDS('./data_from_loretta/sdm_pdm_00_merged_site_data_with_climate/sites.rds')
# 	vars <- c(paste0('bio', c(1, 6, 7, 12, 15, 18)), 'aridity', 'field_pH', 'SAND', 'SILT', 'CLAY')

# 	biomass <- readRDS('./data_from_loretta/sdm_pdm_00_merged_site_data_with_climate/biomass.rds')
# 	morpho_phys <- readRDS('./data_from_loretta/sdm_pdm_00_merged_site_data_with_climate/morpho_phys.rds')

# 	biomass <- merge(sites[ , c('site_id', 'LONGITUDE', 'LATITUDE', 'field_pH', 'SAND', 'SILT', 'CLAY')], biomass, by.x = 'site_id', by.y = 'SITE')
# 	morpho_phys <- merge(sites[ , c('site_id', 'LONGITUDE', 'LATITUDE', 'field_pH', 'SAND', 'SILT', 'CLAY')], morpho_phys, by.x = 'site_id', by.y = 'SITE')

# 	# traits <- c('Biomass', 'Delta13C', 'N_conc', 'CN_ratio', 'Height', 'BladeWidth', 'LeafThick', 'SPAD', 'CanopyDiam', 'WatPot', 'PhotoRate', 'StomCond', 'IntCO2', 'TranspRate')
# 	traits <- c('Biomass', 'Height', 'BladeWidth', 'LeafThick', 'SPAD', 'CanopyDiam', 'PhotoRate', 'StomCond', 'IntCO2', 'TranspRate')

# 	# CYCLE: trait, formula, response distribution (lognormal or gamma)
# 	models <- data.table()
# 	for (trait in traits) {
	
# 		say(trait)

# 		if (trait == 'Biomass') {
# 			x <- biomass
# 		} else {
# 			x <- morpho_phys
# 		}
	
# 		for (form in forms) {
		
# 			this_form <- form
# 			this_form <- paste('y ~', this_form)
# 			this_form_random <- paste0(this_form, ' + (1 | site_id)')
# 			this_form_as_form <- as.formula(this_form)

# 			terms <- terms(this_form_as_form)
# 			terms <- attr(terms, 'term.labels')
# 			terms <- terms[!grepl(terms, pattern = '\\^2')]

# 			cols <- c('site_id', trait, terms)
# 			data <- x[ , ..cols]
# 			colnames(data)[colnames(data) == trait] <- 'y'

# 			# Scale the columns in data that are named in terms
# 			data[ , (terms) := lapply(.SD, scale), .SDcols = terms]

# 			# Define the columns dynamically
# 			numeric_vars <- setdiff(names(data), c('site_id', 'y'))
# 			data_sd <- data[ , c(list(y = sd(y)), lapply(.SD, mean)), by = site_id, .SDcols = numeric_vars]
	
# 			# try different distributions for response
# 			for (distrib in c('lognormal', 'gamma')) {

# 				# model for raw values
# 				if (distrib == 'lognormal') {
				
# 					this_model_raw <- tryCatch({
# 						glmer(this_form_random, data = data, family = gaussian(link = 'log'))
# 					}, error = function(e) {
# 						message("Error fitting model: ", e$message)
# 						NULL
# 					})

# 					this_model_sd <- tryCatch({
# 						glm(this_form_as_form, data = data_sd, family = gaussian(link = 'log'))
# 					}, error = function(e) {
# 						message("Error fitting model: ", e$message)
# 						NULL
# 					})

# 				} else if (distrib == 'gamma') {

# 					this_model_raw <- tryCatch({
# 						glmer(this_form_random, data = data, family = Gamma(link = 'log'))
# 					}, error = function(e) {
# 						message("Error fitting model: ", e$message)
# 						NULL
# 					})

# 					this_model_sd <- tryCatch({
# 						glm(this_form_as_form, data = data_sd, family = Gamma(link = 'log'))
# 					}, error = function(e) {
# 						message("Error fitting model: ", e$message)
# 						NULL
# 					})

# 				}

# 				# remember
# 				if (!is.null(this_model_raw) & !isSingular(this_model_raw)) {

# 					this_model <- this_model_raw

# 					dharma <- tryCatch(
# 						simulateResiduals(fittedModel = this_model, plot = FALSE),
# 						 error = function(e) NULL
# 					)

# 					if (!is.null(dharma)) {
# 						uniformity <- testUniformity(dharma)$p.value
# 						outliers <- testOutliers(dharma)$p.value
# 						dispersion <- testDispersion(dharma)$p.value
# 						quants <- testQuantiles(dharma)$p.value
# 						if (is.null(quants)) quants <- NA
# 					} else {
# 						uniformity <- NA
# 						outliers <- NA
# 						dispersion <- NA
# 						quants <- NA
# 					}

# 					if (distrib == 'lognormal') {
# 						r2_marginal <- r.squaredGLMM(this_model)[1]
# 						r2_conditional <- r.squaredGLMM(this_model)[2]
# 					} else if (distrib == 'gamma') {
# 						r2_marginal <- r.squaredGLMM(this_model)[1, 1]
# 						r2_conditional <- r.squaredGLMM(this_model)[1, 2]
# 					}

# 					models <- rbind(
# 						models,
# 						data.table(
# 							trait = trait,
# 							response = 'raw',
# 							distrib = distrib,
# 							formula = form,
# 							r2_marginal = r2_marginal,
# 							r2_conditional = r2_conditional,
# 							aicc = AICc(this_model),
# 							uniformity = uniformity,
# 							dispersion = dispersion,
# 							outliers = outliers,
# 							quants = quants
# 						)
# 					)

# 				}

# 				if (!is.null(this_model_sd)) {

# 					this_model <- this_model_sd
# 					dharma <- tryCatch(
# 						simulateResiduals(fittedModel = this_model, plot = FALSE),
# 						 error = function(e) NULL
# 					)

# 					if (!is.null(dharma)) {
# 						uniformity <- testUniformity(dharma)$p.value
# 						outliers <- testOutliers(dharma)$p.value
# 						dispersion <- testDispersion(dharma)$p.value
# 						quants <- testQuantiles(dharma)$p.value
# 						if (is.null(quants)) quants <- NA
# 					} else {
# 						uniformity <- NA
# 						outliers <- NA
# 						dispersion <- NA
# 						quants <- NA
# 					}

# 					if (distrib == 'lognormal') {
# 						r2_marginal <- r.squaredGLMM(this_model)[1]
# 						r2_conditional <- r.squaredGLMM(this_model)[2]
# 					} else if (distrib == 'gamma') {
# 						r2_marginal <- r.squaredGLMM(this_model)[1, 1]
# 						r2_conditional <- r.squaredGLMM(this_model)[1, 2]
# 					}

# 					models <- rbind(
# 						models,
# 						data.table(
# 							trait = trait,
# 							response = 'sd',
# 							distrib = distrib,
# 							formula = form,
# 							r2_marginal = r2_marginal,
# 							r2_conditional = r2_conditional,
# 							aicc = AICc(this_model),
# 							uniformity = uniformity,
# 							dispersion = dispersion,
# 							outliers = outliers,
# 							quants = quants
# 						)
# 					)
# 				}

# 			} # next response distribution (lognormal/gamma)

# 		} # next model formula
	
# 	} # next trait

# 	# Create PDF output for subsets of the data.table
# 	pdf(file = paste0(out_dir, '/simple_trait_model_summary_tables.pdf'), width = 15, height = 20)

# 	for (trait in traits) {
# 		for (response in c('raw', 'sd')) {

# 			these <- which(models$trait == trait & models$response == response)
# 			subset_table <- models[these]
# 			subset_table[ , delta_aicc := aicc - min(aicc)]
# 			subset_table <- subset_table[order(subset_table$delta_aicc), ]

# 			subset_table <- subset_table[delta_aicc <= 10]

# 			rll <- exp(-0.5 * subset_table$aicc)
# 			sum_rll <- sum(rll)
# 			subset_table[ , weight := rll / sum_rll]

# 			subset_table[ , r2_marginal := round(r2_marginal, 2)]
# 			subset_table[ , r2_conditional := round(r2_conditional, 2)]
# 			subset_table[ , delta_aicc := round(delta_aicc, 2)]
# 			subset_table[ , aicc := round(aicc, 2)]
# 			subset_table[ , weight := round(weight, 2)]
			
# 			if (nrow(subset_table) > 0) {

# 				n <- min(nrow(subset_table), 60)
# 				subset_table <- subset_table[1:n]

# 				grid::grid.newpage()
# 				gridExtra::grid.table(
# 					subset_table,
# 					rows = NULL,
# 					theme = gridExtra::ttheme_default(
# 						core = list(fg_params = list(cex = 0.8)),
# 						colhead = list(fg_params = list(cex = 0.9, fontface = 'bold'))
# 					)
# 				)
# 				grid::grid.text(
# 					paste('Trait:', trait, '| Response:', response),
# 					x = 0.5, y = 0.95, gp = grid::gpar(fontsize = 14, fontface = 'bold')
# 				)
# 			} else {
# 				grid::grid.newpage()
# 				grid::grid.text(
# 					paste('No data for Trait:', trait, '| Response:', response),
# 					x = 0.5, y = 0.5, gp = grid::gpar(fontsize = 14, fontface = 'bold')
# 				)
# 			}
		
# 		}
# 	}

# 	dev.off()

# say('###################################################################################')
# say('### pre-screen simple models of relationship between environment and occurrence ###')
# say('###################################################################################')

# 	# Evaluate univariate (linear and linear-plus-quadratic) and bivariate (linear and linear-plus-quadratic) models for occurrence. Model fit and residuals are evaluated, and models ranked by AICc.

# 	# variables to evaluate
# 	vars <- c('bio1', 'bio12', 'bio18', 'aridity', 'ph', 'sand', 'silt', 'clay')
# 	twos <- combn(vars, 2, simplify = FALSE)

# 	### all possible combinations of two-variable models
# 	# remove combinations with bio12 and 18 bc thematically redundant
# 	bads <- integer()
# 	for (i in seq_along(twos)) {
# 		this_vars <- twos[[i]]
# 		# bad <- if (all(c('bio1', 'bio12') %in% this_vars) | all(c('bio12', 'bio18') %in% this_vars)) {
# 		bad <- if (all(c('bio12', 'bio18') %in% this_vars)) {
# 			i
# 		} else {
# 			NULL
# 		}
# 		bads <- c(bads, bad)
# 	}

# 	twos <- twos[-bads]

# 	forms <- vars
# 	forms <- c(forms, paste0(vars, ' + I(', vars, '^2)'))
# 	for (i in seq_along(twos)) {

# 		forms <- c(forms, paste0(twos[[i]], collapse = ' + '))	
# 		forms <- c(forms, paste0(c(twos[[i]], paste0('I(', twos[[i]][1], '^2)')), collapse = ' + '))
# 		forms <- c(forms, paste0(c(twos[[i]], paste0('I(', twos[[i]][2], '^2)')), collapse = ' + '))	
# 		forms <- c(forms, paste0(c(twos[[i]], paste0('I(', twos[[i]][1], '^2)'), paste0('I(', twos[[i]][2], '^2)')), collapse = ' + '))	

# 	}

# 	for (i in seq_along(twos)) {
# 		forms <- c(forms, paste0(twos[[i]][1], ' + ', twos[[i]][2], ' + ', twos[[i]][1], ':', twos[[i]][2]))
# 	}

# 	forms <- c(forms, '1')

# 	### setup
# 	out_dir <- './outputs_loretta/integrated_sdm_pdm'
# 	dirCreate(out_dir)

# 	# data
# 	occs <- vect('./outputs_loretta/integrated_sdm_pdm/andropogon_gerardi_occurrences_with_environment_1961_2020_for_integration.gpkg')
# 	occs <- as.data.table(occs)
# 	occs[ , log10_area_km2 := log10(area_km2)]
# 	occs[ , log10_num_poaceae_records := log10(n_poaceae + 1)]
# 	occs$log10_num_poaceae_records[is.na(occs$log10_num_poaceae_records)]  <- quantile(occs$log10_num_poaceae_records, 0.99, na.rm = TRUE)
# 	occs$n_andropogon_gerardi[is.na(occs$n_andropogon_gerardi)]  <- 0

# 	models <- data.table()

# 	for (form in forms) {
	
# 		this_form <- form
# 		this_form <- paste('~', this_form)
# 		this_form_as_form <- as.formula(this_form)

# 		terms <- terms(this_form_as_form)
# 		terms <- attr(terms, 'term.labels')
# 		terms <- terms[!grepl(terms, pattern = '\\^2')]
# 		terms <- terms[!grepl(terms, pattern = '\\:')]

# 		cols <- c('n_andropogon_gerardi', 'log10_area_km2', 'log10_num_poaceae_records', terms)
# 		data <- occs[ , ..cols]

# 		# Scale the columns in data that are named in terms
# 		data[ , (terms) := lapply(.SD, scale), .SDcols = terms]

# 		form1 <- if (length(terms) > 0) {
# 			paste('n_andropogon_gerardi ~ 1 + log10_area_km2 + log10_num_poaceae_records + ', paste0(terms, collapse = ' + '))
# 		} else {
# 			'n_andropogon_gerardi ~ 1 + log10_area_km2 + log10_num_poaceae_records'
# 		}
# 		form2 <- if (length(terms) > 0) {
# 			paste('n_andropogon_gerardi ~ 1 + log10_area_km2 + log10_num_poaceae_records + I(log10_num_poaceae_records^2) + ', paste0(terms, collapse = ' + '))
# 		} else {
# 			'n_andropogon_gerardi ~ 1 + log10_area_km2 + log10_num_poaceae_records + I(log10_num_poaceae_records^2)'
# 		}
# 		form1 <- as.formula(form1)
# 		form2 <- as.formula(form2)
# 		model1 <- glm(form1, data = data, family = poisson(link = 'log'))
# 		model2 <- glm(form2, data = data, family = poisson(link = 'log'))
# 		aicc1 <- AICc(model1)
# 		aicc2 <- AICc(model2)

# 		if (length(terms) == 0) terms <- 1
# 		models <- rbind(
# 			models,
# 			data.table(
# 				terms = paste(terms, collapse = ' '),
# 				env_formula = this_form,
# 				poaceae = c('linear', 'quadratic'),
# 				aicc = c(aicc1, aicc2)
# 			)
# 		)

# 	} # next model formula

# 	models[ , delta_aicc := aicc - min(aicc)]
# 	models <- models[order(delta_aicc)]

# 	# Create PDF output for subsets of the data.table
# 	pdf(file = paste0(out_dir, '/simple_occurrence_model_summary_table.pdf'), width = 8, height = 12)
# 		subset_models <- models[models$delta_aicc <= 500]

# 		grid::grid.newpage()
# 		gridExtra::grid.table(
# 			subset_models,
# 			rows = NULL,
# 			theme = gridExtra::ttheme_default(
# 				core = list(fg_params = list(cex = 0.8)),
# 				colhead = list(fg_params = list(cex = 0.9, fontface = 'bold'))
# 			)
# 		)
# 		grid::grid.text(
# 			paste('Top 1- or 2-covariate occurrence models'),
# 			x = 0.5, y = 0.95, gp = grid::gpar(fontsize = 14, fontface = 'bold')
# 		)

# 	dev.off()
	

# say('#############################################################')
# say('### plots of occurrence and traits in environmental space ###')
# say('#############################################################')

# 	### setup
# 	out_dir <- './outputs_loretta/integrated_sdm_pdm'
# 	dirCreate(out_dir)

# 	### data
# 	# occurrence data
# 	occs <- vect('./outputs_loretta/integrated_sdm_pdm/andropogon_gerardi_occurrences_with_environment_1961_2020_for_integration.gpkg')

# 	### plot of occurrence in environmental space
# 	#############################################

# 		occs_df <- as.data.table(occs)

# 		occs_bio1_bio12 <- ggplot() +
# 			geom_point(data = occs_df[occs$focal_region], mapping = aes(x = bio1, y = bio12), color = alpha('black', 0.2)) +
# 			geom_point(data = occs_df[occs_df$n_andropogon_gerardi > 0 & occs$focal_region], mapping = aes(x = bio1, y = bio12), color = 'blue') +
# 			xlab('Mean annual temperature (°C)') +
# 			ylab('Total annual precipitation (mm)')

# 		occs_bio1_bio18 <- ggplot() +
# 			geom_point(data = occs_df[occs$focal_region], mapping = aes(x = bio1, y = bio18), color = alpha('black', 0.2)) +
# 			geom_point(data = occs_df[occs_df$n_andropogon_gerardi > 0 & occs$focal_region], mapping = aes(x = bio1, y = bio18), color = 'green') +
# 			xlab('Mean annual temperature (°C)') +
# 			ylab('Summer precipitation (mm)')

# 		occs_sand_bio12 <- ggplot() +
# 			geom_point(data = occs_df[occs$focal_region], mapping = aes(x = sand, y = bio12), color = alpha('black', 0.2)) +
# 			geom_point(data = occs_df[occs_df$n_andropogon_gerardi > 0 & occs$focal_region], mapping = aes(x = sand, y = bio12), color = 'red') +
# 			xlab('Sand (SoilGrids)') +
# 			ylab('Total annual precipitation (mm)')

# 		occs_silt_bio12 <- ggplot() +
# 			geom_point(data = occs_df[occs$focal_region], mapping = aes(x = silt, y = bio12), color = alpha('black', 0.2)) +
# 			geom_point(data = occs_df[occs_df$n_andropogon_gerardi > 0 & occs$focal_region], mapping = aes(x = sand, y = bio12), color = 'orange') +
# 			# scale_color_continuous(name = 'Relative\ndensity', trans = 'log10', low = "blue", high = "firebrick1") +
# 			xlab('Silt (SoilGrids)') +
# 			ylab('Total annual precipitation (mm)')

# 		env_plots <- (occs_bio1_bio12 + occs_bio1_bio18) / (occs_sand_bio12 + occs_silt_bio12)

# 		ggsave(env_plots, filename = paste0(out_dir, '/occurrence_in_environmental_space.png'), width = 14, height = 12, bg = 'white')

# 	### maps of environment in geographic space
# 	###########################################

# 		# site data
# 		sites <- readRDS('./data_from_loretta/sdm_pdm_00_merged_site_data_with_climate/sites.rds')
# 		sites <- vect(sites, geom = c('LONGITUDE', 'LATITUDE'), crs = getCRS('WGS84'))
# 		sites <- project(sites, occs)

# 		extent <- ext(occs[occs$focal_region])
# 		extent <- as.vector(extent)

# 		preds <- c('bio1', 'bio12', 'bio15', 'bio18', 'aridity', 'sand', 'silt', 'ph')
# 		maps <- list()
# 		for (pred in preds) {
		
# 			say(pred)

# 			maps[[length(maps) + 1]] <- ggplot() +
# 				layer_spatial(data = occs, fill = 'gray', color = 'gray30') +
# 				layer_spatial(data = occs[occs$focal_region], mapping = aes(fill = .data[[pred]]), color = 0) +
# 				# scale_fill_continuous(trans = 'log10') +
# 				layer_spatial(sites, pch = 1, size = 3, color = 'yellow') +
# 				xlim(extent[1], extent[2]) + ylim(extent[3], extent[4]) +
# 				ggtitle(pred)
		
# 		}

# 	maps <- plot_grid(plotlist = maps, ncol = 4)
# 	ggsave(maps, filename = paste0(out_dir, '/maps_of_environmental_variables_in_focal_region.png'), width = 19.2, height = 10.8, bg = 'white')

say('FINIS!', level = 1)

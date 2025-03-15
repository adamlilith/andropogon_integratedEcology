### MODELING ANDROPOGON GERARDI DISTRIBUTION, PHENOTYPE, PHYSIOLOGY, GENOTYPE, and ASSOCIATED MICROBIAL COMMUNITIES
### Erica Newman | Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2023-12
###
### This script constructs a non-integrated model for AG geographic distribution.
###
### source('C:/Ecology/R/andropogon_integratedEcology/sdm_pdm/sdm_pdm_02_exploratory_data_analysis.r')
### source('C:/Subarashi/R/andropogon_integratedEcology/sdm_pdm/sdm_pdm_02_exploratory_data_analysis.r')
###
### CONTENTS ###
### setup ###
### trait histograms ###
### correlations between environmental variables at phenotypic sample sites ###

#############
### setup ###
#############

	rm(list = ls())

	# drive <- 'C:/Ecology/'
	drive <- 'C:/Subarashi/'

	setwd(paste0(drive, '/Research/Andropogon/Andropogon'))

	library(cluster) # clustering
	library(data.table) # fast data tables
	library(enmSdmX) # GIS & SDMing
	library(ggplot2) # plots
	library(omnibus) # utilities
	library(predicts) # GIS & SDMing
	library(readxl) # Excel
	library(terra) # spatial objects

# say('########################')
# say('### trait histograms ###')
# say('########################')

# 	dirCreate('./outputs_loretta/integrated_sdm_pdm')

# 	sites <- readRDS('./data_from_loretta/sdm_pdm_00_merged_site_data_with_climate/sites.rds')
# 	biomass <- readRDS('./data_from_loretta/sdm_pdm_00_merged_site_data_with_climate/biomass.rds')
# 	morpho_phys <- readRDS('./data_from_loretta/sdm_pdm_00_merged_site_data_with_climate/morpho_phys.rds')

# 	### biomass variables
# 	#####################

# 	x <- biomass
# 	x <- x[ , c('SITE', 'VegBiomass', 'ReproBiomass', 'Biomass', 'VegRepro_ratio')]
# 	x <- x[complete.cases(x)]
	
# 	long <- melt(x, id.vars = "SITE", variable.name = "trait", value.name = "value")

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
	
# 	long <- melt(x, id.vars = "SITE", variable.name = "trait", value.name = "value")

# 	# plot
# 	plots <- ggplot(long, aes(x = value)) +
# 		geom_histogram(bins = 30) +
# 		facet_wrap(~ trait, scales = 'free') +
# 		theme_minimal() +
# 		theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
# 		labs(title = 'Morphology/physiology trait histograms')

# 	ggsave(plots, filename = './outputs_loretta/integrated_sdm_pdm/morphology_trait_histograms.png', width = 16, height = 12, dpi = 300, bg = 'white')

say('###############################################################################')
say('### correlations between environmental variables at phenotypic sample sites ###')
say('###############################################################################')

	dirCreate('./outputs_loretta/integrated_sdm_pdm')

	sites <- readRDS('./data_from_loretta/sdm_pdm_00_merged_site_data_with_climate/sites.rds')
	vars <- c(paste0('bio', c(1, 6, 7, 12, 15, 18)), 'aridity')

	cors <- cor(sites[ , ..vars], method = 'spearman')
	dists <- 1 - abs(cors)
	dists <- as.dist(dists)

	cluster <- hclust(dists)

	png('./outputs_loretta/integrated_sdm_pdm/predictor_correlations_at_phenotyped_sites.png', width = 1200, height = 800, res = 150)

		plot(cluster, xlab = NULL, ylab = '1 - abs(correlation)')
		abline(h = 0.3, col = 'red', lty = 'dashed')

	dev.off()




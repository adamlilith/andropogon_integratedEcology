### MODELING ANDROPOGON GERARDI DISTRIBUTION, PHENOTYPE, PHYSIOLOGY, GENOTYPE, and ASSOCIATED MICROBIAL COMMUNITIES
### Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2023-12
###
### This script constructs a non-integrated model for AG geographic distribution.
###
### source('C:/Ecology/R/andropogon_integratedEcology/sdm/sdm_03_nonintegrated_model_n_maps.r')
### source('C:/Subarashi/R/andropogon_integratedEcology/sdm/sdm_03_nonintegrated_model_n_maps.r')
###
### CONTENTS ###
### setup ###
### map for publication of AG present-day distribution with sampling sites ###

#############
### setup ###
#############

	rm(list = ls())

	# drive <- 'C:/Ecology/'
	drive <- 'C:/Subarashi/'

	setwd(paste0(drive, '/Research/Andropogon/Andropogon'))

	library(cowplot) # combining ggplots
	library(enmSdmX) # GIS/SDMing
	library(ggplot2) # plotting
	library(ggspatial) # plotting spatial things
	library(omnibus) # helper functions
	library(readxl) # read Excel documents
	library(scales) # for plotting transparency
	library(terra) # spatial objects

say('####################################################################################')
say('### map for publication of AG sampling density with phenological sampling sites ###')
say('####################################################################################')

	out_dir <- './outputs_loretta/'

	### SDM
	ag_vect_sq <- vect('./data_from_adam_and_loretta/andropogon_gerardi_occurrences_with_environment_1961_2020_climatena.gpkg')

	# ratio of AG to Poaceae
	ag_vect_sq$ag_vs_poaceae_ratio <- ag_vect_sq$n_andropogon_gerardi / ag_vect_sq$n_poaceae

	# plot extent
	ag_vect_sq_pres <- ag_vect_sq[ag_vect_sq$n_andropogon_gerardi > 0]
	extent <- ext(ag_vect_sq_pres)
	extent <- as.vector(extent)
	x_range <- (extent[2] - extent[1])
	y_range <- (extent[4] - extent[3])
	extent[1] <- extent[1] + 0.15 * x_range
	extent[3] <- extent[3] + 0.125 * y_range
	extent[4] <- extent[4] - 0.2 * y_range

	### North America
	nam <- vect(paste0(drive, '/Research Data/GADM/Version 4.1/High Res North America Level 1 sans Great Lakes SpatVector WGS84.gpkg'))
	nam <- project(nam, ag_vect_sq)

	### sample site data
	pheno_sites <- read_xlsx('./data_from_loretta/!plant_sitelevel_data_11NOV2024 - Google Sheets [aggregated by Erica].xlsx', sheet = 'plantmaster_bysite_19NOV2024')

	pheno_sites <- vect(pheno_sites, geom = c('LONGITUDE', 'LATITUDE'), crs = getCRS('WGS84'))
	pheno_sites <- project(pheno_sites, ag_vect_sq)

	map <- ggplot() +
		# layer_spatial(ag_vect_sq, aes(fill = n_andropogon_gerardi), color = NA) +
		# layer_spatial(ag_vect_sq, aes(fill = ag_density), color = NA) +
		layer_spatial(ag_vect_sq, aes(fill = ag_vs_poaceae_ratio), color = NA) +
		layer_spatial(nam, color = 'gray20', fill = NA, linewidth = 0.3) +
		scale_fill_gradient(
			name = 'Sampling\nratio',
			guide = guide_legend(reverse = TRUE),
			label = function(x) sprintf("%.3f", x),
			trans = 'log2',
			low = '#c7e9c0',
			high = '#006d2c',
			na.value = 'gray80'
		) +
		layer_spatial(pheno_sites, pch = 21, fill = 'yellow', size = 4) +
		xlim(extent[1], extent[2]) + ylim(extent[3], extent[4]) +
		# ggtitle(expression('Present-day distribution of ' * italic('Andropogon gerardi')), subtitle = '1961-2020 | N-mixture model') +
		theme(
			plot.title = element_text(size = 16),
			plot.subtitle = element_text(size = 14),
			legend.title = element_text(size = 16),
			legend.text = element_text(size = 14),
			axis.text = element_text(size = 12)
		)

	ggsave(plot = map, filename = paste0(out_dir, '/ag_sampling_frequency_with_pheno_sample_sites.png'), width = 12, height = 10, dpi = 600)

say('####################################################################################')
say('### map for publication of AG present-day non-integrated SDM with sampling sites ###')
say('####################################################################################')

	in_out_dir <- './outputs_loretta/sdm_[nmixture]_[pseudoabsences_0.99]_[bio1^2_bio12^2_bio15^2]_[priors_ddnorm]/'

	### SDM
	ag_vect_sq <- vect(paste0(in_out_dir, '/sdm_nmixture_1961_2020_climate_focus.gpkg'))

	# plot extent
	ag_vect_sq_pres <- ag_vect_sq[ag_vect_sq$n_andropogon_gerardi > 0]
	extent <- ext(ag_vect_sq_pres)
	extent <- as.vector(extent)
	x_range <- (extent[2] - extent[1])
	y_range <- (extent[4] - extent[3])
	extent[1] <- extent[1] + 0.15 * x_range
	extent[3] <- extent[3] + 0.125 * y_range
	extent[4] <- extent[4] - 0.2 * y_range

	# centroids of counties with AG records
	cents_with_ag <- ag_vect_sq[ag_vect_sq$n_andropogon_gerardi > 0]
	cents_with_ag <- centroids(cents_with_ag)

	# color scheme for occurrence prediction
	# lambda_sq_quants <- quantile(ag_vect_sq$lambda_mean, c(0.25, 0.5, 0.75, 0.90, 0.95))
	# quant_labels <- c('[0 - 0.25)', '[0.25 - 0.50)', '[0.50 - 0.75)', '[0.75 - 0.90)', '[0.90 - 0.95)', '[0.95 - 1]')
	
	lambda_sq_quants <- quantile(ag_vect_sq$lambda_mean, c(0.5, 0.75, 0.90, 0.95))
	quant_labels <- c('[0 - 0.50)', '[0.50 - 0.75)', '[0.75 - 0.90)', '[0.90 - 0.95)', '[0.95 - 1]')
	
	# ag_vect_sq$quant_col <- NA
	# ag_vect_sq$quant_col[ag_vect_sq$lambda_mean < lambda_sq_quants[1]] <- quant_labels[1]
	# ag_vect_sq$quant_col[ag_vect_sq$lambda_mean >= lambda_sq_quants[1] & ag_vect_sq$lambda_mean < lambda_sq_quants[2]] <- quant_labels[2]
	# ag_vect_sq$quant_col[ag_vect_sq$lambda_mean >= lambda_sq_quants[2] & ag_vect_sq$lambda_mean < lambda_sq_quants[3]] <- quant_labels[3]
	# ag_vect_sq$quant_col[ag_vect_sq$lambda_mean >= lambda_sq_quants[3] & ag_vect_sq$lambda_mean < lambda_sq_quants[4]] <- quant_labels[4]
	# ag_vect_sq$quant_col[ag_vect_sq$lambda_mean >= lambda_sq_quants[4] & ag_vect_sq$lambda_mean < lambda_sq_quants[5]] <- quant_labels[5]
	# ag_vect_sq$quant_col[ag_vect_sq$lambda_mean >= lambda_sq_quants[5]] <- quant_labels[6]

	ag_vect_sq$quant_col <- NA
	ag_vect_sq$quant_col[ag_vect_sq$lambda_mean < lambda_sq_quants[1]] <- quant_labels[1]
	ag_vect_sq$quant_col[ag_vect_sq$lambda_mean >= lambda_sq_quants[1] & ag_vect_sq$lambda_mean < lambda_sq_quants[2]] <- quant_labels[2]
	ag_vect_sq$quant_col[ag_vect_sq$lambda_mean >= lambda_sq_quants[2] & ag_vect_sq$lambda_mean < lambda_sq_quants[3]] <- quant_labels[3]
	ag_vect_sq$quant_col[ag_vect_sq$lambda_mean >= lambda_sq_quants[3] & ag_vect_sq$lambda_mean < lambda_sq_quants[4]] <- quant_labels[4]
	ag_vect_sq$quant_col[ag_vect_sq$lambda_mean >= lambda_sq_quants[4]] <- quant_labels[5]

	# fill_scale <- c(
	# 	'[0 - 0.25)' = 'gray85',
	# 	'[0.25 - 0.5)' = alpha('forestgreen', 0.20),
	# 	'[0.50 - 0.75)' = alpha('forestgreen', 0.40),
	# 	'[0.75 - 0.90)' = alpha('forestgreen', 0.65),
	# 	'[0.90 - 0.95)' = alpha('forestgreen', 0.75),
	# 	'[0.95 - 1]' = 'forestgreen'
	# )

	fill_scale <- c(
		'[0 - 0.50)' = 'gray85',
		'[0.50 - 0.75)' = alpha('forestgreen', 0.40),
		'[0.75 - 0.90)' = alpha('forestgreen', 0.65),
		'[0.90 - 0.95)' = alpha('forestgreen', 0.75),
		'[0.95 - 1]' = 'forestgreen'
	)

	### North America
	nam <- vect(paste0(drive, '/Research Data/GADM/Version 4.1/High Res North America Level 1 sans Great Lakes SpatVector WGS84.gpkg'))
	nam <- project(nam, ag_vect_sq)

	### sample site data
	pheno_sites <- read_xlsx('./data_from_loretta/!plant_sitelevel_data_11NOV2024 - Google Sheets [aggregated by Erica].xlsx', sheet = 'plantmaster_bysite_19NOV2024')

	pheno_sites <- vect(pheno_sites, geom = c('LONGITUDE', 'LATITUDE'), crs = getCRS('WGS84'))
	pheno_sites <- project(pheno_sites, ag_vect_sq)

	map <- ggplot() +
		layer_spatial(ag_vect_sq, aes(fill = quant_col), color = NA) +
		layer_spatial(nam, color = 'gray20', fill = NA, linewidth = 0.3) +
		scale_fill_manual(
			name = 'Quantile\n of λ',
			values = fill_scale,
			guide = guide_legend(reverse = TRUE)
		) +
		layer_spatial(cents_with_ag, pch = 16, color = alpha('gray20', 0.5), size = 1) +
		layer_spatial(pheno_sites, pch = 21, fill = 'yellow', size = 4) +
		xlim(extent[1], extent[2]) + ylim(extent[3], extent[4]) +
		# ggtitle(expression('Present-day distribution of ' * italic('Andropogon gerardi')), subtitle = '1961-2020 | N-mixture model') +
		theme(
			plot.title = element_text(size = 16),
			plot.subtitle = element_text(size = 14),
			legend.title = element_text(size = 16),
			legend.text = element_text(size = 14),
			axis.text = element_text(size = 12)
		)

	ggsave(plot = map, filename = paste0(in_out_dir, '/sdm_nmixture_lambda_status_quo_with_pheno_sample_sites.png'), width = 12, height = 10, dpi = 600)





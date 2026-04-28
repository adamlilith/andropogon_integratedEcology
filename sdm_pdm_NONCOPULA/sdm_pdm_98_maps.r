### MODELING ANDROPOGON GERARDI DISTRIBUTION, PHENOTYPE, PHYSIOLOGY, GENOTYPE, and ASSOCIATED MICROBIAL COMMUNITIES
### Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2023-12
###
### This script constructs select maps of model inputs and outputs.
###
### source('C:/Kaji/R/andropogon_integratedEcology/sdm_pdm/sdm_pdm_98_maps.r')
###
### CONTENTS ###
### setup ###
### map of occurrence model input and output plus sample sites ###
### 
#############
### setup ###
#############

	rm(list = ls())

	drive <- 'C:/Kaji/'

	setwd(paste0(drive, '/Research/Andropogon/Andropogon'))
	source(paste0(drive, '/R/andropogon_integratedEcology/sdm_pdm/sdm_pdm_00_shared_functions_and_variables.r'))

	library(viridis)

###########################
### user-defined values ###
###########################

# say('##################################################################')
# say('### map of occurrence model input and output plus sample sites ###')
# say('##################################################################')

# 	occs <- vect('./data_from_adam_and_loretta/andropogon_gerardi_occurrences_with_environment_1961_2020.gpkg')
# 	preds <- vect('./outputs_loretta/integrated_sdm_pdm/models_occurrence/[occs_poisson~normal_homoscedastic_bio1^2_bio12^2_bio15^2_[bias~1]]/prediction_vector_nam.gpkg')
# 	sites <- readRDS('./data_from_loretta/sdm_pdm_00_merged_site_data_with_climate/sites.rds')
# 	sites <- vect(sites, geom = c('LONGITUDE', 'LATITUDE'), crs = getCRS('WGS84'))
# 	sites <- project(sites, occs)

# 	nam <- vect('./data_from_gadm/gadm_4pt1_level_1_north_america_sans_alaska_lambert.gpkg')

# 	# extent
# 	extent <- ext(sites)
# 	extent <- as.polygons(extent, crs = occs)
# 	extent <- buffer(extent, width = 200 * 1000) # nominal plot extent
# 	extent_display <- buffer(extent, width = 300 * 1000) # larger than plot extent
# 	extent <- ext(extent)
# 	extent <- as.vector(extent)

# 	pred_vect_display <- crop(preds, extent_display)

# 	pred_vect_display$color <- NA
# 	quants <- quantile(pred_vect_display$N_ag_county_sq, c(0.25, 0.5, 0.75, 0.95))
# 	pred_vect_display$color[pred_vect_display$N_ag_county_sq < quants[1]] <- 0
# 	pred_vect_display$color[pred_vect_display$N_ag_county_sq >= quants[1] & pred_vect_display$N_ag_county_sq < quants[2]] <- 1
# 	pred_vect_display$color[pred_vect_display$N_ag_county_sq >= quants[2] & pred_vect_display$N_ag_county_sq < quants[3]] <- 2
# 	pred_vect_display$color[pred_vect_display$N_ag_county_sq >= quants[3] & pred_vect_display$N_ag_county_sq < quants[4]] <- 3
# 	pred_vect_display$color[pred_vect_display$N_ag_county_sq >= quants[4]] <- 4
# 	pred_vect_display$color <- factor(pred_vect_display$color)

# 	viridis_palette <- viridis_pal(option = 'viridis')(256)
# 	occ_map <- ggplot() +
# 		layer_spatial(pred_vect_display, aes(fill = n_andropogon_gerardi), color = NA) +
# 		scale_fill_gradient(trans = 'log10', na.value = 'gray80', name = 'No. records', low = viridis_palette[1], high = viridis_palette[256]) +
# 		layer_spatial(nam, color = 'gray40', fill = NA, linewidth = 0.3) +
# 		layer_spatial(sites, pch = 21, size = 4, fill = 'white') +
# 		xlim(extent[1], extent[2]) + ylim(extent[3], extent[4]) +
# 		ggtitle('A) Occurrences') +
# 		theme(
# 			plot.title = element_text(size = 22),
# 			plot.subtitle = element_text(size = 14),
# 			legend.title = element_text(size = 20),
# 			legend.text = element_text(size = 20)
# 		)

# 	pred_map <- ggplot() +
# 		layer_spatial(pred_vect_display, aes(fill = color), color = NA) +
# 		scale_fill_manual(
# 			values = c('0' = '#f7fcf5', '1' = '#c7e9c0', '2' = '#74c476', '3' = '#238b45', '4' = '#00441b'),
# 			labels = c('≥0 to <0.25', '≥0.25 to <0.50', '≥0.50 to <0.75', '≥0.75 to <0.95', '≥0.95 to 1'),
# 			name = 'Relative\nabundance\nquantile'
# 		) +
# 		guides(fill = guide_legend(reverse = TRUE)) +
# 		layer_spatial(nam, color = 'gray40', fill = NA, linewidth = 0.3) +
# 		layer_spatial(sites, pch = 21, size = 4, fill = 'white') +
# 		xlim(extent[1], extent[2]) + ylim(extent[3], extent[4]) +
# 		ggtitle('B) Distribution model predictions') +
# 		theme(
# 			plot.title = element_text(size = 22),
# 			plot.subtitle = element_text(size = 14),
# 			legend.title = element_text(size = 20),
# 			legend.text = element_text(size = 20)
# 		)

# 	maps <- occ_map + pred_map

# 	filename <- paste0('./outputs_loretta/sem/occs_plus_sdm_from_sdm_pdm_bio1^2_bio12^2_bio15^2_bias~1.png')
# 	ggsave(plot = maps, filename = filename, width = 22, height = 10, dpi = 600)

say('##################################################################')
say('### maps of best models for each biomass and non-biomass facet ###')
say('##################################################################')

		top_models <- fread('./outputs_loretta/integrated_sdm_pdm/summary_of_ALL_top_models.csv')

		for (i in 1:nrow(top_models)) {
		
		  	facet <- top_models$facet[i]
		  	resp_distrib <- top_models$resp_distrib[i]
		  	
		  	out_dir <- paste0('./outputs_loretta/integrated_sdm_pdm/models_', facet, '/[', facet, '_', tolower(resp_distrib), '~normal~', preds_filename, ']', ifelse(log_precip, '_log_precip', ''), '/')
		
		}


say(date())
say('FINIS!', deco = '+', level = 1)

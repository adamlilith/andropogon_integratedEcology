#' Map biomass model predictions
#' 
#' source('C:/Kaji/R/andropogon_integratedEcology/sdm_pdm/functions/function_map_biomass_nonbiomass.r')
#'
#' facet			Name of facet ('biomass', 'n_concentration', etc.)
#' out_dir 			Folder in which to save map. Leave as `NULL` to not save.
#' filename_append String to append to file name. Will have "_" prefixed to it.
#' pred_vect_nam 	SpatVector with predictions.
#' response_var 	Any of: 'mu_<facet>_county_sq', 'mu_<facet>_county_ssp245_2041_2070', 'mu_<facet>_county_ssp245_2071_2100', 'mus_<facet>_ssp370_2041_2070', 'mu_<facet>_county_ssp370_2071_2100'
#' response_var_type 'mean', 'median', or 'inner_quant'
#' data_occs 		From prepare_occurrence_data().
#' data_biomass_nonbiomass 	From prepare_biomass_data() or prepare_nonbiomass_data().
#' title 			Character for plot title.
#' subtitle 		Character for plot subtitle.
#' legend_title 	Title for legend.
map_biomass_nonbiomass <- function(
	facet,
	out_dir,
	filename_append,
	pred_vect_nam,
	response_var,
	response_var_type,
	data_occs,
	data_biomass_nonbiomass,
	title = 'Biomass',
	subtitle = NULL,
	legend_title = 'Biomass (g)'
) {

	nam <- vect('./data_from_gadm/gadm_4pt1_level_1_north_america_sans_alaska_lambert.gpkg')
	nam <- simplifyGeom(nam, tolerance = 1000)

	# extent
	site_vect <- data_biomass_nonbiomass$site_vect
	site_vect <- project(site_vect, pred_vect_nam)
	extent <- ext(site_vect)
	extent <- as.polygons(extent, crs = pred_vect_nam)
	extent <- buffer(extent, width = 320 * 1000) # nominal plot extent
	extent <- ext(extent)
	extent_coords <- as.vector(extent)

	pred_vect_display <- crop(pred_vect_nam, extent)

	# delineate current range "core"
	range_core <- delineate_range_core(pred_vect_display, column = response_var, core_quant = core_quant)

	# get range of values for plotting
	vars <- c(
		paste0(facet, '_', response_var_type, '_sq'),
		paste0(facet, '_', response_var_type, '_ssp245_2041_2070'),
		paste0(facet, '_', response_var_type, '_ssp245_2071_2100'),
		paste0(facet, '_', response_var_type, '_ssp370_2041_2070'),
		paste0(facet, '_', response_var_type, '_ssp370_2071_2100')
	)

	# Calculate mean and standard deviation of biomass by site
	raw_data <- if (facet == 'biomass') { data_biomass_nonbiomass$raw_data_biomass } else { data_biomass_nonbiomass$raw_data_facet }
	raw_name <- get_raw_trait_name_from_rfriendly(facet)
	response_stats <- raw_data[ , .(response_mean = mean(get(raw_name)), response_sd = sd(get(raw_name))), by = SITE]

	min_val <- min(response_stats$response_mean)
	max_val <- max(response_stats$response_mean)
	for (var in vars) {

		x <- pred_vect_display[[var]]
		x <- unlist(x)
		x <- x[pred_vect_display$focal_region]
		x <- x[!is.infinite(x)]
		min_val <- min(min_val, x)
		max_val <- max(max_val, x)

	}
	min_val <- 0.99 * min_val
	max_val <- 1.01 * max_val

	resp_limits <- c(min_val, max_val)

	counties_with_ag <- pred_vect_display[pred_vect_display$n_andropogon_gerardi > 0]
	counties_with_ag <- centroids(counties_with_ag)

	center_col <- paste0(facet, '_', response_var_type)

	map <- ggplot() +
		layer_spatial(pred_vect_display, aes(fill = .data[[response_var]]), color = NA) +
		scale_fill_gradientn(
			name = legend_title,
			colors = c('#edf8e9', '#74c476', '#005a32'),
			limits = resp_limits
		) +
		layer_spatial(nam, color = 'gray40', fill = NA, linewidth = 0.3) +
		layer_spatial(counties_with_ag, pch = 16, color = alpha('gray20', 0.5), size = 0.35) +
		layer_spatial(site_vect, aes(fill = .data[[center_col]]), pch = 21, size = 4) +
		layer_spatial(range_core, color = 'orange2', fill = NA, linewidth = 1) +
		coord_sf(xlim = c(extent_coords[1], extent_coords[2]), ylim = c(extent_coords[3], extent_coords[4]), expand = FALSE) +
		ggtitle(title, subtitle) +
		theme(
			plot.title = element_text(size = 16),
			plot.subtitle = element_text(size = 14)
		)

	if (!is.null(out_dir)) ggsave(plot = map, filename = paste0(out_dir, '/map_', facet, '_', response_var_type, '_', filename_append, '.png'), width = 12, height = 10, dpi = 200)
	invisible(map)

}

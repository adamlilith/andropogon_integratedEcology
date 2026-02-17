#' Map probability a value is 0
#'
#' @param out_dir Folder in which to save map.
#' @param filename_append String to append to file name. Will have "_" prefixed to it.
#' @param pred_vect_nam SpatVector with predictions.
#' @param response_var Any of: 'psi_county_sq', 'psi_county_ssp245_2041_2070', 'psi_county_ssp245_2071_2100', 'psi_county_ssp370_2041_2070', or 'psi_county_ssp370_2071_2100'
#' @param title Character for plot title.
#' @param subtitle Character for plot subtitle.
map_psi <- function(
	out_dir,
	filename_append,
	pred_vect_nam,
	response_var,
	title = '<TRAIT>',
	subtitle = NULL
) {

	nam <- vect('./data_from_gadm/gadm_4pt1_level_1_north_america_sans_alaska_lambert.gpkg')
	data_biomass <- prepare_biomass_data(formula_biomass = ~ 1, log_precip = TRUE, n_response_curve_values = n_response_curve_values, calib = FALSE)
	site_vect <- data_biomass$site_vect_biomass

	# extent
	site_vect <- project(site_vect, pred_vect_nam)
	extent <- ext(site_vect)
	extent <- as.polygons(extent, crs = pred_vect_nam)
	extent <- buffer(extent, width = 320 * 1000) # nominal plot extent
	extent <- ext(extent)
	extent_coords <- as.vector(extent)

	pred_vect_display <- crop(pred_vect_nam, extent)

	# delineate current range "core"
	range_core_absent <- delineate_range_core(pred_vect_display, column = response_var, ag_core_quant = 0.1, rule = '<=')

	resp_limits <- c(0, 1)

	counties_with_ag <- pred_vect_display[pred_vect_display$n_andropogon_gerardi > 0]
	counties_with_ag <- centroids(counties_with_ag)

	map <- ggplot() +
		layer_spatial(pred_vect_display, aes(fill = .data[[response_var]]), color = NA) +
		scale_fill_continuous(
			name = 'Probability',
			type = 'viridis',
			limits = resp_limits
		) +
		layer_spatial(nam, color = 'gray40', fill = NA, linewidth = 0.3) +
		layer_spatial(counties_with_ag, pch = 16, color = alpha('black', 0.8), size = 0.35) +
		layer_spatial(range_core_absent, color = '#43145B', fill = NA, linewidth = 0.8) +
		layer_spatial(site_vect, pch = 4, size = 4, color = 'white') +
		coord_sf(xlim = c(extent_coords[1], extent_coords[2]), ylim = c(extent_coords[3], extent_coords[4]), expand = FALSE) +
		ggtitle(title, subtitle) +
		theme(
			plot.title = element_text(size = 16),
			plot.subtitle = element_text(size = 14)
		)

	ggsave(plot = map, filename = paste0(out_dir, '/map_psi_', filename_append, '.png'), width = 12, height = 10, dpi = 200)
	invisible(map)

}


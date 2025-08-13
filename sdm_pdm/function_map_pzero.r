#' Map probability a value is 0
#'
#' @param out_dir Folder in which to save map.
#' @param filename_append String to append to file name. Will have "_" prefixed to it.
#' @param pred_vect_nam SpatVector with predictions.
#' @param facet "R-friendly" name of facet: 'occs', 'biomass', 'height', 'canopy_diameter', etc.
#' @param response_var Any of: '<facet>_pzero_county_sq', '<facet>_pzero_county_ssp245_2041_2070', '<facet>_pzero_county_ssp245_2071_2100', '<facet>_pzero_county_ssp370_2041_2070', or '<facet>_pzero_county_ssp370_2071_2100'
#' @param data_traits From prepare_nonbiomass_traits().
#' @param title Character for plot title.
#' @param subtitle Character for plot subtitle.
#' @param plot_range_core If `TRUE`, extract from `pred_vect_nam` the range core based on `N_ag_county_sq` and plot it.
#' @param ag_core_quant Quantile used to delineate core from non-core.
map_pzero <- function(
	out_dir,
	filename_append,
	pred_vect_nam,
	facet,
	response_var,
	data_traits,
	title = '<TRAIT>',
	subtitle = NULL,
	plot_range_core = TRUE,
	ag_core_quant = 0.95
) {

	nam <- vect('./data_from_gadm/gadm_4pt1_level_1_north_america_sans_alaska_lambert.gpkg')

	# extent
	site_vect <- data_traits$site_vect_traits
	site_vect <- project(site_vect, pred_vect_nam)
	extent <- ext(site_vect)
	extent <- as.polygons(extent, crs = pred_vect_nam)
	extent <- buffer(extent, width = 200 * 1000)
	extent_display <- buffer(extent, width = 300 * 1000)
	extent <- ext(extent)
	extent <- as.vector(extent)

	pred_vect_display <- crop(pred_vect_nam, extent_display)

	# delineate current range "core"
	if (plot_range_core) range_core <- delineate_range_core(pred_vect_display, column = 'N_ag_county_sq', ag_core_quant = ag_core_quant)

	# get range of values for plotting
	vars <- c('_pzero_county_sq', '_pzero_county_ssp245_2041_2070', '_pzero_county_ssp245_2071_2100', '_pzero_county_ssp370_2041_2070', '_pzero_county_ssp370_2071_2100')
	vars <- paste0(facet, vars)

	# Calculate mean and standard deviation of biomass by site
	min_val <- max_val <- 0

	for (var in vars) {

		x <- pred_vect_display[[var]]
		x <- unlist(x)
		x <- x[pred_vect_display$focal_region]
		
		max_val <- max(max_val, min(1, 1.05 * max(x)))

	}

	resp_limits <- c(min_val, max_val)

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
		layer_spatial(counties_with_ag, pch = 16, color = alpha('gray20', 0.1), size = 0.35) +
		layer_spatial(data_traits$site_vect_traits, pch = 3, size = 4, color = 'orange') +
		xlim(extent[1], extent[2]) + ylim(extent[3], extent[4]) +
		ggtitle(title, subtitle) +
		theme(
			plot.title = element_text(size = 16),
			plot.subtitle = element_text(size = 14)
		)


	if (plot_range_core) map <- map + layer_spatial(range_core, color = 'cyan', fill = NA, size = 2)

	ggsave(plot = map, filename = paste0(out_dir, '/map_', facet, '_pzero_', filename_append, '.png'), width = 12, height = 10, dpi = 600)
	invisible(map)

}


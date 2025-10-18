#' Map change in probability of a zero value
#'
#' @param out_dir Folder in which to save map.
#' @param filename_append String to append to file name. Will have "_" prefixed to it.
#' @param fut String indicating future scenario, for example, 'ssp245_2041_2070'.
#' @param pred_vect_nam SpatVector with predictions.
#' @param facet R-friendly name of facet, like "occs", "biomass", "height", "canopy_diameter", etc.
#' @param response_var Column in `pred_vect_nam` representing non-present-day predictions.
#' @param response_var_sq Column in `pred_vect_nam` representing present-day predictions.
#' @param data_traits From function_prepare_nonbiomass_traits().
#' @param title Character for plot title.
#' @param subtitle Character for plot subtitle.
#' @param plot_range_core If `TRUE`, extract from `pred_vect_nam` the range core based on `N_ag_county_sq` and plot it.
#' @param ag_core_quant Quantile used to delineate core from non-core.
map_pzero_change <- function(
	out_dir,
	filename_append,
	fut,
	pred_vect_nam,
	facet,
	response_var,
	response_var_sq,
	data_traits,
	title = 'Change',
	subtitle = 'Change',
	plot_range_core = TRUE,
	ag_core_quant = 0.95
) {

	nam <- vect('./data_from_gadm/gadm_4pt1_level_1_north_america_sans_alaska_lambert.gpkg')

	# extent
	site_vect <- data_traits[['site_vect_traits']]
	site_vect <- project(site_vect, pred_vect_nam)
	extent <- ext(site_vect)
	extent <- as.polygons(extent, crs = pred_vect_nam)
	extent <- buffer(extent, width = 200 * 1000)
	extent_display <- buffer(extent, width = 300 * 1000) # extent for getting range of plotted values
	extent <- ext(extent)
	extent <- as.vector(extent)

	pred_vect_display <- crop(pred_vect_nam, extent_display)

	# delineate current and future range "core"
	if (plot_range_core) {
		
		column <- paste0('pzero_', facet, '_county_', fut)
		range_core_fut <- delineate_range_core(pred_vect_display, column = column, ag_core_quant)

		column <- paste0('pzero_', facet, '_county_sq')
		range_core_sq <- delineate_range_core(pred_vect_display, column = column, ag_core_quant)

	}

	# get range of values for plotting
	# get range of values for plotting
	vars <- c(
		paste0('pzero_', facet, '_county_sq'),
		paste0('pzero_', facet, '_county_ssp245_2041_2070'),
		paste0('pzero_', facet, '_county_ssp245_2071_2100'),
		paste0('pzero_', facet, '_county_ssp370_2041_2070'),
		paste0('pzero_', facet, '_county_ssp370_2071_2100')
	)

	min_val <- Inf
	max_val <- -Inf
	x_sq <- unlist(pred_vect_display[[response_var_sq]])
	for (var in vars) {

		x_fut <- pred_vect_display[[var]]
		x_fut <- unlist(x_fut)
		delta <- x_fut - x_sq
		min_val <- min(min_val, delta[!is.infinite(delta)])
		max_val <- max(max_val, delta[!is.infinite(delta)])

	}

	resp_limits <- c(min_val, max_val)

	counties_with_ag <- pred_vect_display[pred_vect_display$n_andropogon_gerardi > 0]
	counties_with_ag <- centroids(counties_with_ag)

	trait_nice <- tolower(get_nice_trait(facet)$short)

	sq_y <- unlist(pred_vect_display[[response_var_sq]])
	fut_y <- unlist(pred_vect_display[[response_var]])

	pred_vect_display$delta <- fut_y - sq_y

	map <- ggplot() +
		layer_spatial(pred_vect_display, aes(fill = delta), color = NA) +
		layer_spatial(nam, color = 'gray40', fill = NA, linewidth = 0.3) +
		layer_spatial(counties_with_ag, pch = 16, color = alpha('gray20', 0.5), size = 0.35) +
		layer_spatial(site_vect, pch = 3, size = 4) +
		scale_fill_gradient2(
			name = 'Change in\nprobability',
			low = '#2166ac',
			mid = 'beige',
			high = '#b2182b',
			midpoint = 0,
			limits = resp_limits
		) +
		xlim(extent[1], extent[2]) + ylim(extent[3], extent[4]) +
		ggtitle(title, subtitle) +
		theme(
			plot.title = element_text(size = 16),
			plot.subtitle = element_text(size = 14)
		)

	if (plot_range_core) {
		map <- map +
			layer_spatial(range_core_sq, color = 'gray40', fill = NA, linewidth = 1) +
			layer_spatial(range_core_fut, color = 'orange', fill = NA, linewidth = 1)
	}

	ggsave(plot = map, filename = paste0(out_dir, '/map_change_', facet, '_pzero_', filename_append, '.png'), width = 12, height = 10, dpi = 600)
	invisible(map)

}

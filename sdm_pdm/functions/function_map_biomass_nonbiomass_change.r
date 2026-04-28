#' Map biomass model change predictions
#'
#' facet			Name of facet
#' out_dir			Folder in which to save map. Leave as `NULL` to not save.
#' filename_append 	String to append to file name. Will have "_" prefixed to it.
#' fut 				String indicating future scenario, for example, 'ssp245_2041_2070'.
#' pred_vect_nam 	SpatVector with predictions.
#' response_var 	Any of: 'mu_<facet>_county_sq', 'mu_<facet>_county_ssp245_2041_2070', 'mu_<facet>_county_ssp245_2071_2100', 'mus_<facet>_ssp370_2041_2070', 'mu_<facet>_county_ssp370_2071_2100'
#' response_var_type 'mean', 'median', or 'inner_quant'
#' data_occs 		From prepare_occurrence_data().
#' data_biomass_nonbiomass 	From prepare_biomass_data().
#' title 			Character for plot title.
#' subtitle 		Character for plot subtitle.
map_biomass_nonbiomass_change <- function(
	facet,
	out_dir,
	filename_append,
	fut,
	pred_vect_nam,
	response_var,
	response_var_type,
	data_occs,
	data_biomass_nonbiomass,
	title = 'Title',
	subtitle = 'Subtitle',
	legend_title = 'Percent\nchange'
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

	# delineate current and future range "core"
	column <- paste0(facet, '_', response_var_type, '_', fut)
	range_core_fut <- delineate_range_core(pred_vect_display, column = column, core_quant)

	column <- paste0(facet, '_', response_var_type, '_sq')
	range_core_sq <- delineate_range_core(pred_vect_display, column = column, core_quant)

	# get range of values for plotting
	vars <- c(
		paste0(facet, '_', response_var_type, '_sq'),
		paste0(facet, '_', response_var_type, '_ssp245_2041_2070'),
		paste0(facet, '_', response_var_type, '_ssp245_2071_2100'),
		paste0(facet, '_', response_var_type, '_ssp370_2041_2070'),
		paste0(facet, '_', response_var_type, '_ssp370_2071_2100')
	)

	min_val <- Inf
	max_val <- -Inf
	x_sq <- pred_vect_display[[paste0(facet, '_', response_var_type, '_sq')]]
	x_sq <- unlist(x_sq)
	x_sq <- x_sq[pred_vect_display$focal_region]
	for (var in vars) {

		x_fut <- pred_vect_display[[var]]
		x_fut <- unlist(x_fut)
		x_fut <- x_fut[pred_vect_display$focal_region]
		delta <- (x_fut - x_sq) / x_sq
		min_val <- min(min_val, delta[!is.infinite(delta)])
		max_val <- max(max_val, delta[!is.infinite(delta)])

	}

	resp_limits <- c(min_val, max_val)

	counties_with_ag <- pred_vect_display[pred_vect_display$n_andropogon_gerardi > 0]
	counties_with_ag <- centroids(counties_with_ag)

	pred_vect_display$delta <- (pred_vect_display[[paste0(facet, '_', response_var_type, '_', fut)]] - pred_vect_display[[paste0(facet, '_', response_var_type, '_sq')]]) / pred_vect_display[[paste0(facet, '_', response_var_type, '_sq')]]

	map <- ggplot() +
		layer_spatial(pred_vect_display, aes(fill = delta), color = NA) +
		layer_spatial(nam, color = 'gray40', fill = NA, linewidth = 0.3) +
		layer_spatial(counties_with_ag, pch = 16, color = alpha('gray20', 0.5), size = 0.35) +
		layer_spatial(range_core_sq, color = 'gray40', fill = NA, linewidth = 1) +
		layer_spatial(range_core_fut, color = 'orange', fill = NA, linewidth = 1) +
		layer_spatial(site_vect, pch = 3, size = 4) +
		scale_fill_gradient2(
			name = legend_title,
			low = '#7b3294',
			mid = 'beige',
			high = '#008837',
			midpoint = 0,
			limits = resp_limits,
			labels = scales::label_percent(digits = 2)
		) +
		coord_sf(xlim = c(extent_coords[1], extent_coords[2]), ylim = c(extent_coords[3], extent_coords[4]), expand = FALSE) +
		ggtitle(title, subtitle) +
		theme(
			plot.title = element_text(size = 16),
			plot.subtitle = element_text(size = 14)
		)

	if (!is.null(out_dir)) ggsave(plot = map, filename = paste0(out_dir, '/map_change_', facet, '_', response_var_type, '_', filename_append, '.png'), width = 12, height = 10, dpi = 200)
	invisible(map)

}

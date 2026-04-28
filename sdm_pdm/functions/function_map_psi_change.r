#' Map change in probability of a zero value
#'
#' out_dir			Folder in which to save map. Leave as `NULL` to not save.
#' filename_append	String to append to file name. Will have "_" prefixed to it.
#' fut				String indicating future scenario, for example, 'ssp245_2041_2070'.
#' pred_vect_nam	SpatVector with predictions.
#' response_var		Column in `pred_vect_nam` representing non-present-day predictions.
#' response_var_sq	Column in `pred_vect_nam` representing present-day predictions.
#' title			Character for plot title.
#' subtitle			Character for plot subtitle.
map_psi_change <- function(
	out_dir,
	filename_append,
	fut,
	pred_vect_nam,
	response_var,
	response_var_sq,
	title = 'Change',
	subtitle = 'Change'
) {

	nam <- vect('./data_from_gadm/gadm_4pt1_level_1_north_america_sans_alaska_lambert.gpkg')
	nam <- simplifyGeom(nam, tolerance = 1000)

	data_biomass <- prepare_biomass_data(formula_biomass = ~ 1, n_response_curve_values = n_response_curve_values, calib = FALSE)
	site_vect <- data_biomass$site_vect

	# extent
	site_vect <- project(site_vect, pred_vect_nam)
	extent <- ext(site_vect)
	extent <- as.polygons(extent, crs = pred_vect_nam)
	extent <- buffer(extent, width = 320 * 1000) # nominal plot extent
	extent <- ext(extent)
	extent_coords <- as.vector(extent)

	pred_vect_display <- crop(pred_vect_nam, extent)

	# get range of values for plotting
	vars <- c(
		paste0('psi_sq'),
		paste0('psi_ssp245_2041_2070'),
		paste0('psi_ssp245_2071_2100'),
		paste0('psi_ssp370_2041_2070'),
		paste0('psi_ssp370_2071_2100')
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

	sq_y <- unlist(pred_vect_display[[response_var_sq]])
	fut_y <- unlist(pred_vect_display[[response_var]])

	pred_vect_display$delta <- fut_y - sq_y

	map <- ggplot() +
		layer_spatial(pred_vect_display, aes(fill = delta), color = NA) +
		layer_spatial(nam, color = 'gray40', fill = NA, linewidth = 0.3) +
		layer_spatial(counties_with_ag, pch = 16, color = alpha('gray20', 0.8), size = 0.35) +
		layer_spatial(site_vect, pch = 3, size = 4, color = 'black') +
		scale_fill_gradient2(
			name = 'Change in\nprobability\n(future -\npresent)',
			low = '#b2182b',
			mid = 'beige',
			high = '#2166ac',
			midpoint = 0,
			limits = resp_limits
		) +
		coord_sf(xlim = c(extent_coords[1], extent_coords[2]), ylim = c(extent_coords[3], extent_coords[4]), expand = FALSE) +
		ggtitle(title, subtitle) +
		theme(
			plot.title = element_text(size = 16),
			plot.subtitle = element_text(size = 14)
		)

	if (!is.null(out_dir)) ggsave(plot = map, filename = paste0(out_dir, '/map_psi_change_', filename_append, '.png'), width = 12, height = 10, dpi = 200)
	invisible(map)

}

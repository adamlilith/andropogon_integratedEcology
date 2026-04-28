#' Maps of change in abundance and probability of absence in 1930s
#'
#' source('C:/Kaji/R/andropogon_integratedEcology/sdm_pdm/functions/function_map_occurrence_change_1930s.r') 
#'
#' zero_inflated	`TRUE` or `FALSE`-- if `TRUE`, map for `psi` will be created, too.
#' pred_vect_nam	SpatVector with predictions for all of North America
#' pred_vect_1930s	SpatVector with predictions for CONUS
map_occurrence_change_1930s <- function(zero_inflated, pred_vect_nam, pred_vect_1930s) {

	data_traits <- prepare_nonbiomass_data(facet = 'height', formula = ~ 1, log_precip = FALSE, n_response_curve_values = n_response_curve_values, calib = calib)

	dust_bowl <- vect('./data_from_others/dust_bowl_counties_with_most_severe_wind_erosion.gpkg')
	dust_bowl <- aggregate(dust_bowl)

	# prepare spatial vector for mapping
	delta_vect <- pred_vect_nam[pred_vect_nam$state_province %in% c('Colorado', 'Nebraska', 'Kansas', 'Oklahoma', 'Texas', 'New Mexico')]

	# plot extent from Dust Bowl counties	
	extent <- ext(dust_bowl)
	extent <- as.polygons(extent, crs = delta_vect)
	extent <- buffer(extent, width = 20 * 1000) # nominal plot extent
	extent_display <- buffer(extent, width = 30 * 1000) # larger than plot extent
	extent <- ext(extent)
	extent <- as.vector(extent)

	nam <- vect('./data_from_gadm/gadm_4pt1_level_1_north_america_sans_alaska_lambert.gpkg')
	nam <- crop(nam, ext(extent_display))

	# identify 1930s range core from wider distribution than we will show in the map
	site_vect <- data_traits$site_vect_facet
	site_vect <- project(site_vect, pred_vect_1930s)
	extent_range_core <- ext(site_vect)
	extent_range_core <- as.polygons(extent_range_core, crs = pred_vect_1930s)
	extent_range_core <- buffer(extent_range_core, width = 500 * 1000) # nominal plot extent

	pred_vect_conus_restricted <- crop(pred_vect_1930s, extent_range_core)
	range_core <- delineate_range_core(pred_vect_conus_restricted, column = 'N_ag_county_mean_1930s', ag_core_quant = ag_core_quant)

	# change
	focal_n <- unlist(pred_vect_1930s[['N_ag_county_mean_1930s']])
	# delta_vect$delta_n <- focal_n / delta_vect$N_ag_county_mean_sq
	delta_vect$delta_n <- (focal_n - delta_vect$N_ag_county_mean_sq) / delta_vect$N_ag_county_mean_sq
	
	focal_psi <- unlist(pred_vect_1930s[['psi_county_1930s']])
	delta_vect$delta_psi <- focal_psi - delta_vect$psi_county_sq
	
	delta_vect_display <- crop(delta_vect, ext(extent_display))

	title <- paste0('Abundance')
	if (zero_inflated) title <- paste0('A) ', title)

	map_change <- ggplot() +
		layer_spatial(nam, fill = 'gainsboro', color = NA) +
		layer_spatial(
			delta_vect_display,
			aes(fill = delta_n),
			color = NA
		) +
		scale_fill_gradient2(
			name = 'Change (%)',
			low = '#c51b7d',
			mid = '#f7f7f7',
			high = '#4d9221',
			midpoint = 0,
			labels = scales::percent_format(accuracy = 1)
		) +
		layer_spatial(nam, color = 'gray40', fill = NA) +
		layer_spatial(dust_bowl, color = 'gray20', fill = NA, linewidth = 1) +
		layer_spatial(range_core, color = 'orange2', fill = NA, linewidth = 1) +
		layer_spatial(site_vect, pch = 3, size = 4) +
		annotation_scale() +
		coord_sf(xlim = c(extent[1], extent[2]), ylim = c(extent[3], extent[4]), expand = FALSE) +
		ggtitle(title) +
		theme(
			plot.title = element_text(size = 16),
			plot.subtitle = element_text(size = 14)
		)

	if (!zero_inflated) {
		maps <- map_change
	} else {

		title <- paste0('B) Change in probability of presence')

		map_change_psi <- ggplot() +
			layer_spatial(nam, fill = 'gainsboro', color = NA) +
			layer_spatial(
				delta_vect_display,
				aes(fill = delta_psi),
				color = NA
			) +
			scale_fill_gradient2(
				name = paste0('1930s\nminus\npresent'),
				low = '#b2182b',
				mid = 'beige',
				high = '#2166ac',
				midpoint = 0,
				# trans = 'log10',
				# limits = rep_limits_psi
			) +
			layer_spatial(nam, color = 'gray40', fill = NA) +
			layer_spatial(dust_bowl, color = 'gray20', fill = NA, linewidth = 1) +
			layer_spatial(range_core, color = 'orange2', fill = NA, linewidth = 1) +
			layer_spatial(site_vect, pch = 3, size = 4) +
			annotation_scale() +
			coord_sf(xlim = c(extent[1], extent[2]), ylim = c(extent[3], extent[4]), expand = FALSE) +
			ggtitle(title) +
			theme(
				plot.title = element_text(size = 16),
				plot.subtitle = element_text(size = 14)
			)

		maps <- plot_grid(map_change, map_change_psi, ncol = 2, align = 'hv')

	} # if zero-inflated

	filename <- paste0(out_dir, '/map_abundance_change_1930s.png')
	width <- if (zero_inflated) { 14 } else { 7 }
	ggsave(maps, filename = filename, width = width, height = 8, dpi = 200, bg = 'white')
	
	invisible(maps)

}

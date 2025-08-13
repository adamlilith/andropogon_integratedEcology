#' Maps of change in abundance and probability of absence in 1930s
#'
#' @param zero_inflated `TRUE` or `FALSE`-- if `TRUE`, map for `pzero` will be created, too.
#' @param pred_vect_nam SpatVector with predictions for all of North America
#' @param pred_vect_conus SpatVector with predictions for CONUS
#' pred_vect_nam
map_occurrence_change_1930s <- function(zero_inflated, pred_vect_nam, pred_vect_conus) {

	dust_bowl <- vect('./data_from_others/dust_bowl_counties_with_most_severe_wind_erosion.gpkg')
	dust_bowl <- aggregate(dust_bowl)

	# prepare spatial vector for mapping
	delta_vect <- pred_vect_nam[pred_vect_nam$country == 'United States']
	delta_vect <- delta_vect[delta_vect$state_province != 'Alaska']

	# plot extent from Dust Bowl counties
	
	extent <- ext(dust_bowl)
	extent <- as.polygons(extent, crs = delta_vect)
	extent <- buffer(extent, width = 40 * 1000) # nominal plot extent
	extent_display <- buffer(extent, width = 50 * 1000) # larger than plot extent
	extent <- ext(extent)
	extent <- as.vector(extent)

	nam <- vect('./data_from_gadm/gadm_4pt1_level_1_north_america_sans_alaska_lambert.gpkg')
	nam <- crop(nam, ext(extent_display))

	# identify 1930s range core from wider distribution than we will show in the map
	site_vect <- data_traits$site_vect_traits
	site_vect <- project(site_vect, pred_vect_conus)
	extent_range_core <- ext(site_vect)
	extent_range_core <- as.polygons(extent_range_core, crs = pred_vect_conus)
	extent_range_core <- buffer(extent_range_core, width = 500 * 1000) # nominal plot extent

	pred_vect_conus_restricted <- crop(pred_vect_conus, extent_range_core)
	range_core <- delineate_range_core(pred_vect_conus_restricted, column = 'N_ag_county_thirties', ag_core_quant = ag_core_quant)

	# change
	focal_n <- unlist(pred_vect_conus[['N_ag_county_thirties']])
	delta_vect$delta_n <- focal_n / delta_vect$N_ag_county_sq
	
	focal_pzero <- unlist(pred_vect_conus[['pzero_occs_county_thirties']])
	delta_vect$delta_pzero <- delta_vect$pzero_occs_county_sq - focal_pzero
	
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
			name = paste0('Ratio of\n1930s-to-present\nabundance'),
			low = '#c51b7d',
			mid = '#f7f7f7',
			high = '#4d9221',
			midpoint = 0,
			trans = 'log10',
			# limits = resp_limits_n
		) +
		layer_spatial(nam, color = 'gray40', fill = NA) +
		layer_spatial(dust_bowl, color = 'gray20', fill = NA, linewidth = 1) +
		layer_spatial(range_core, color = 'orange2', fill = NA, linewidth = 1) +
		layer_spatial(site_vect, pch = 3, size = 4) +
		annotation_scale() +
		xlim(extent[1], extent[2]) + ylim(extent[3], extent[4]) +
		ggtitle(title) +
		theme(
			plot.title = element_text(size = 16),
			plot.subtitle = element_text(size = 14)
		)

	if (!zero_inflated) {
		maps <- map_change
	} else {

		title <- paste0('B) Probability of absence')

		map_change_pzero <- ggplot() +
			layer_spatial(nam, fill = 'gainsboro', color = NA) +
			layer_spatial(
				delta_vect_display,
				aes(fill = delta_pzero),
				color = NA
			) +
			scale_fill_gradient2(
				name = paste0('Present\nminus\n', decade),
				low = '#b2182b',
				mid = 'beige',
				high = '#2166ac',
				midpoint = 0,
				# trans = 'log10',
				# limits = rep_limits_pzero
			) +
			layer_spatial(nam, color = 'gray40', fill = NA) +
			layer_spatial(dust_bowl, color = 'gray20', fill = NA, linewidth = 1) +
			layer_spatial(range_core, color = 'orange2', fill = NA, linewidth = 1) +
			layer_spatial(site_vect, pch = 3) +
			annotation_scale() +
			xlim(extent[1], extent[2]) + ylim(extent[3], extent[4]) +
			ggtitle(title) +
			theme(
				plot.title = element_text(size = 16),
				plot.subtitle = element_text(size = 14)
			)

		maps <- plot_grid(map_change, map_change_pzero, ncol = 2, align = 'h')

	} # if zero-inflated

	filename <- paste0(out_dir, '/map_abundance_change_1930s.png')
	width <- if (zero_inflated) { 14 } else { 7 }
	ggsave(maps, filename = filename, width = width, height = 8, dpi = 600, bg = 'white')
	
	invisible(maps)

}

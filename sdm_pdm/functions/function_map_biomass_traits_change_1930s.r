#' Maps of change in biomass or non-biomass traits in 1930s
#'
#' @param facet 'biomass', 'height', etc.
#' @param data_biomass_traits Either `data_biomass_mu` or `data_traits`
#' @param pred_vect_nam SpatVector with predictions for all of North America
#' @param pred_vect_conus SpatVector with predictions for CONUS
#' @param pred_vect_1930s,pred_vect_1950s SpatVectors with predictions from 1930s and 1950s
#' @param formula_mu,formula_sigma,formula_pzero `NULL` or formula for each
map_biomass_traits_change_1930s <- function(facet, data_biomass_traits, pred_vect_nam, pred_vect_1930s, formula_mu, formula_sigma, formula_pzero) {

	zero_inflated <- !is.null(formula_pzero)

	if (facet == 'biomass') {
		site_vect <- data_biomass_traits$site_vect_biomass
	} else {
		site_vect <- data_biomass_traits$site_vect_traits
	}

	dust_bowl <- vect('./data_from_others/dust_bowl_counties_with_most_severe_wind_erosion.gpkg')
	dust_bowl <- aggregate(dust_bowl)

	# prepare spatial vector for mapping
	delta_vect <- pred_vect_nam[pred_vect_nam$country == 'United States']
	delta_vect <- delta_vect[delta_vect$state_province != 'Alaska']

	# plot extent from Dust Bowl counties
	extent <- ext(dust_bowl)
	extent <- as.polygons(extent, crs = delta_vect)
	extent <- buffer(extent, width = 35 * 1000) # nominal plot extent
	extent_display <- buffer(extent, width = 50 * 1000) # larger than plot extent
	extent <- ext(extent)
	extent <- as.vector(extent)

	nam <- vect('./data_from_gadm/gadm_4pt1_level_1_north_america_sans_alaska_lambert.gpkg')
	nam <- crop(nam, ext(extent_display))

	# change
	focal_vals <- unlist(pred_vect_1930s[[paste0('mu_', facet, '_county_mean_1930s')]])
	sq_vals <- unlist(delta_vect[[paste0('mu_', facet, '_county_mean_sq')]])
	delta_vect$delta <- (focal_vals - sq_vals) / delta_vect$mu_biomass_county_mean_sq
	
	say('(1930s - SQ) / SQ')

	if (zero_inflated) {

		focal_vals <- unlist(pred_vect_1930s[[paste0('pzero_', facet, '_county_1930s')]])
		sq_vals <- unlist(delta_vect[[paste0('pzero_', facet, '_county_sq')]])
		delta_vect$delta_pzero <- focal_vals - sq_vals
		
	}
		
	delta_vect_display <- crop(delta_vect, ext(extent_display))

	nice_facet <- get_nice_trait(facet)
	title <- paste0(nice_facet$short, ' Change: Present vs. 1930s')
	if (zero_inflated) title <- paste0('A) ', title)
	legend_title <- paste0('Change (%)')

	map_change <- ggplot() +
		layer_spatial(nam, fill = 'gainsboro', color = NA) +
		layer_spatial(
			delta_vect_display,
			aes(fill = delta),
			color = NA
		) +
		scale_fill_gradient2(
			name = legend_title,
			low = '#c51b7d',
			mid = '#f7f7f7',
			high = '#4d9221',
			midpoint = 0,
			labels = scales::percent_format(accuracy = 1)
		) +
		layer_spatial(nam, color = 'gray40', fill = NA) +
		layer_spatial(dust_bowl, color = 'gray20', fill = NA, linewidth = 1) +
		# layer_spatial(range_core, color = 'orange2', fill = NA, linewidth = 1) +
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

		title <- paste0('B) Probability of Zero Value: 1930s minus Present')

		map_change_pzero <- ggplot() +
			layer_spatial(nam, fill = 'gainsboro', color = NA) +
			layer_spatial(
				delta_vect_display,
				aes(fill = delta_pzero),
				color = NA
			) +
			scale_fill_gradient2(
				name = paste0('Difference'),
				low = '#b2182b',
				mid = 'beige',
				high = '#2166ac',
				midpoint = 0,
				# trans = 'log10',
				# limits = rep_limits_pzero
			) +
			layer_spatial(nam, color = 'gray40', fill = NA) +
			layer_spatial(dust_bowl, color = 'gray20', fill = NA, linewidth = 1) +
			# layer_spatial(range_core, color = 'orange2', fill = NA, linewidth = 1) +
			layer_spatial(site_vect, pch = 3) +
			annotation_scale() +
			xlim(extent[1], extent[2]) + ylim(extent[3], extent[4]) +
			ggtitle(title) +
			theme(
				plot.title = element_text(size = 16),
				plot.subtitle = element_text(size = 14)
			)

		maps <- plot_grid(map_change, map_change_pzero, ncol = 2, align = 'h')

	}
	
	width <- if (zero_inflated) { 14 } else { 7 }
	filename <- paste0(out_dir, '/map_', facet, '_change_1930s.png')
	ggsave(maps, filename = filename, width = width, height = 8, dpi = 600, bg = 'white')
	
	invisible(maps)

}

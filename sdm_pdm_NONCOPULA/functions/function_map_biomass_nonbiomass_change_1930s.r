#' Maps of change in biomass or non-biomass traits in 1930s
#'
#' facet 'biomass', 'height', etc.
#' data_biomass_nonbiomass Either `data_biomass` or `data_traits`
#' pred_vect_nam SpatVector with predictions for all of North America
#' pred_vect_conus SpatVector with predictions for CONUS
#' pred_vect_1930s,pred_vect_1950s SpatVectors with predictions from 1930s and 1950s
#' formula_mu,formula_psi `NULL` or formula for each
map_biomass_nonbiomass_change_1930s <- function(facet, data_biomass_nonbiomass, pred_vect_nam, pred_vect_1930s, formula_mu, formula_psi) {

	zero_inflated <- !is.null(formula_psi)

	if (facet == 'biomass') {
		site_vect <- data_biomass_nonbiomass$site_vect_biomass
	} else {
		site_vect <- data_biomass_nonbiomass$site_vect_facet
	}

	dust_bowl <- vect('./data_from_others/dust_bowl_counties_with_most_severe_wind_erosion.gpkg')
	dust_bowl <- aggregate(dust_bowl)

	# prepare spatial vector for mapping
	delta_vect <- pred_vect_nam[pred_vect_nam$country == 'United States']
	delta_vect <- delta_vect[delta_vect$state_province %in% c('Colorado', 'Nebraska', 'Kansas', 'Oklahoma', 'Texas', 'New Mexico')]

	# plot extent from Dust Bowl counties
	extent <- ext(dust_bowl)
	extent <- as.polygons(extent, crs = delta_vect)
	extent <- buffer(extent, width = 20 * 1000) # nominal plot extent
	extent_display <- buffer(extent, width = 30 * 1000) # larger than plot extent
	extent <- ext(extent)
	extent <- as.vector(extent)

	nam <- vect('./data_from_gadm/gadm_4pt1_level_1_north_america_sans_alaska_lambert.gpkg')
	nam <- crop(nam, ext(extent_display))

	# change
	focal_vals <- unlist(pred_vect_1930s[[paste0('mu_', facet, '_county_mean_1930s')]])
	sq_vals <- unlist(delta_vect[[paste0('mu_', facet, '_county_mean_sq')]])
	
	delta_vect$delta <- (focal_vals - sq_vals) / sq_vals
	
	say('(1930s - SQ) / SQ')

	if (zero_inflated) {

		focal_vals <- unlist(pred_vect_1930s[[paste0('psi_county_1930s')]])
		sq_vals <- unlist(delta_vect[[paste0('psi_county_sq')]])
		delta_vect$delta_psi <- focal_vals - sq_vals
		
	}
		
	delta_vect_display <- crop(delta_vect, ext(extent_display))

	nice_facet <- get_nice_trait(facet)
	title <- paste0(nice_facet$short, ' Change: 1930s Relative to Present')
	if (zero_inflated) title <- paste0('A) ', title)

	map_change <- ggplot() +
		layer_spatial(nam, fill = 'gainsboro', color = NA) +
		layer_spatial(
			delta_vect_display,
			aes(fill = delta),
			color = NA
		) +
		scale_fill_gradient2(
			name = 'Change',
			low = '#c51b7d',
			mid = '#f7f7f7',
			high = '#4d9221',
			midpoint = 0,
			labels = scales::percent_format(accuracy = 1)
		) +
		layer_spatial(nam, color = 'gray40', fill = NA) +
		layer_spatial(dust_bowl, color = 'gray20', fill = NA, linewidth = 1) +
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

		title <- paste0('B) Probability of Presence: 1930s Minus Present')

		map_change_psi <- ggplot() +
			layer_spatial(nam, fill = 'gainsboro', color = NA) +
			layer_spatial(
				delta_vect_display,
				aes(fill = delta_psi),
				color = NA
			) +
			scale_fill_gradient2(
				name = paste0('Difference'),
				low = '#b2182b',
				mid = 'beige',
				high = '#2166ac',
				midpoint = 0,
				# trans = 'log10',
				# limits = rep_limits_psi
			) +
			layer_spatial(nam, color = 'gray40', fill = NA) +
			layer_spatial(dust_bowl, color = 'gray20', fill = NA, linewidth = 1) +
			# layer_spatial(range_core, color = 'orange2', fill = NA, linewidth = 1) +
			layer_spatial(site_vect, pch = 3, size = 4) +
			annotation_scale() +
			coord_sf(xlim = c(extent[1], extent[2]), ylim = c(extent[3], extent[4]), expand = FALSE) +
			ggtitle(title) +
			theme(
				plot.title = element_text(size = 16),
				plot.subtitle = element_text(size = 14)
			)

		maps <- plot_grid(map_change, map_change_psi, ncol = 2, align = 'h')

	}
	
	width <- if (zero_inflated) { 14 } else { 7 }
	filename <- paste0(out_dir, '/map_', facet, '_change_1930s.png')
	ggsave(maps, filename = filename, width = width, height = 8, dpi = 200, bg = 'white')
	
	invisible(maps)

}

#' Map biomass model change predictions
#'
#' @param out_dir Folder in which to save map.
#' @param filename_append String to append to file name. Will have "_" prefixed to it.
#' @param fut String indicating future scenario, for example, 'ssp245_2041_2070'.
#' @param pred_vect_display SpatVector with predictions.
#' @param data_occs From prepare_occurrences().
#' @param data_biomass From prepare_biomass().
#' @param title Character for plot title.
#' @param subtitle Character for plot subtitle.
#' @param plot_range_core If `TRUE`, extract from `pred_vect_display` the range core based on `N_ag_county_sq` and plot it.
#' @param ag_core_quant Quantile used to delineate core from non-core.
map_biomass_change <- function(
	out_dir,
	filename_append,
	fut,
	pred_vect_nam,
	response_var,
	data_occs,
	data_biomass,
	title = 'Biomass Change',
	subtitle = NULL,
	plot_range_core = TRUE,
	ag_core_quant = 0.95
) {

	nam <- vect('./data_from_gadm/gadm_4pt1_level_1_north_america_sans_alaska_lambert.gpkg')

	# extent
	site_vect <- data_biomass$site_vect_biomass
	site_vect <- project(site_vect, pred_vect_nam)
	extent <- ext(site_vect)
	extent <- as.polygons(extent, crs = pred_vect_nam)
	extent <- buffer(extent, width = 200 * 1000) # nominal plot extent
	# extent_display <- buffer(extent, width = 300 * 1000) # larger than plot extent
	extent_display <- buffer(extent, width = 200 * 1000) # larger than plot extent
	extent <- ext(extent)
	extent <- as.vector(extent)

	pred_vect_display <- crop(pred_vect_nam, extent_display)

	# delineate current and future range "core"
	if (plot_range_core) {
		
		column <- paste0('mu_biomass_county_mean_', fut)
		range_core_fut <- delineate_range_core(pred_vect_display, column = column, ag_core_quant)

		column <- paste0('mu_biomass_county_mean_sq')
		range_core_sq <- delineate_range_core(pred_vect_display, column = column, ag_core_quant)

	}

	# get range of values for plotting
	vars <- c('mu_biomass_county_mean_sq', 'mu_biomass_county_mean_ssp245_2041_2070', 'mu_biomass_county_mean_ssp245_2071_2100', 'mu_biomass_county_mean_ssp370_2041_2070', 'mu_biomass_county_mean_ssp370_2071_2100')

	min_val <- Inf
	max_val <- -Inf
	x_sq <- pred_vect_display$mu_biomass_county_mean_sq[pred_vect_display$focal_region]
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

	legend_title <- 'Percent\nbiomass\nchange'

	pred_vect_display$delta <- (pred_vect_display[[paste0('mu_biomass_county_mean_', fut)]] - pred_vect_display$mu_biomass_county_mean_sq) / pred_vect_display$mu_biomass_county_mean_sq

	map <- ggplot() +
		layer_spatial(pred_vect_display, aes(fill = delta), color = NA) +
		layer_spatial(nam, color = 'gray40', fill = NA, linewidth = 0.3) +
		layer_spatial(counties_with_ag, pch = 16, color = alpha('gray20', 0.5), size = 0.35) +
		layer_spatial(data_biomass$site_vect_biomass, pch = 3, size = 4) +
		scale_fill_gradient2(
			name = legend_title,
			low = '#7b3294',
			mid = 'beige',
			high = '#008837',
			midpoint = 0,
			limits = resp_limits,
			labels = scales::label_percent(digits = 2)
		) +
		xlim(extent[1], extent[2]) + ylim(extent[3], extent[4]) +
		ggtitle(title, subtitle) +
		theme(
			plot.title = element_text(size = 16),
			plot.subtitle = element_text(size = 14)
		)

	if (plot_range_core) {
		map <- map + layer_spatial(range_core_sq, color = 'gray40', fill = NA, linewidth = 1)
		map <- map + layer_spatial(range_core_fut, color = 'orange', fill = NA, linewidth = 1)
	}

	ggsave(plot = map, filename = paste0(out_dir, '/map_change_biomass_mean_', filename_append, '.png'), width = 12, height = 10, dpi = 600)
	invisible(map)

}

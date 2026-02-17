#' Map of change in abundance
#'
#' @param out_dir Folder in which to save image.
#' @param filename_append String to append to file name. Will be concatenated as per `paste0(out_dir, '/map_abundance_', filename_append, '.png')``.
#' @param pred_vect_nam SpatVector with predictions.
#' @param response_var Name of response variable: 'N_ag_county_sq', 'N_ag_county_ssp245_2041_2070', 'N_ag_county_ssp245_2071_2100', 'N_ag_county_ssp370_2041_2070', 'N_ag_county_ssp370_2071_2100'.
#' @param data_occs From prepare_occurrence_data().
#' @param title Character for plot title.
#' @param title Character for plot subtitle.
map_occurrence_change <- function(
	out_dir,
	filename_append,
	pred_vect_nam,
	response_var,
	data_occs,
	title = 'Change in AG Abundance',
	subtitle = NULL
) {

	nam <- vect('./data_from_gadm/gadm_4pt1_level_1_north_america_sans_alaska_lambert.gpkg')
	data_traits <- prepare_nonbiomass_data(facet = 'height', formula_facet = ~ 1, log_precip = FALSE)
	site_vect <- data_traits$site_vect_facet
	site_vect <- project(site_vect, pred_vect_nam)

	# extent
	extent <- ext(site_vect)
	extent <- as.polygons(extent, crs = pred_vect_nam)
	extent <- buffer(extent, width = 320 * 1000) # nominal plot extent
	extent <- ext(extent)
	extent_coords <- as.vector(extent)

	pred_vect_display <- crop(pred_vect_nam, extent)

	# delineate current and future range "core"
	column <- response_var
	range_core_fut <- delineate_range_core(pred_vect_display, column = column, ag_core_quant)

	column <- 'N_ag_county_mean_sq'
	range_core_sq <- delineate_range_core(pred_vect_display, column = column, ag_core_quant)

	# get range of values for plotting
	vars <- c('N_ag_county_mean_ssp245_2041_2070', 'N_ag_county_mean_ssp245_2071_2100', 'N_ag_county_mean_ssp370_2041_2070', 'N_ag_county_mean_ssp370_2071_2100')

	min_val <- Inf
	max_val <- -Inf
	x_sq <- pred_vect_display$N_ag_county_mean_sq[pred_vect_display$focal_region]
	for (var in vars) {
		
		x_fut <- pred_vect_display[[var]]
		x_fut <- unlist(x_fut)
		x_fut <- x_fut[pred_vect_display$focal_region]
		delta <- x_fut / x_sq
		min_val <- c(min_val, delta[!is.infinite(delta)], na.rm = TRUE)
		max_val <- max(max_val, delta[!is.infinite(delta)], na.rm = TRUE)

	}

	min_val <- min(min_val[!is.infinite(min_val) & min_val > 0], na.rm = TRUE)
	resp_limits <- log10(c(min_val, max_val))

	counties_with_ag <- pred_vect_display[pred_vect_display$n_andropogon_gerardi > 0]
	counties_with_ag <- centroids(counties_with_ag)

	legend_title <- 'Abundance\nchange ratio\nlog(future /\n  present)'

	this_pred_vect <- pred_vect_display
	this_pred_vect$delta <- log10(this_pred_vect[[response_var]] / this_pred_vect$N_ag_county_mean_sq)

	map <- ggplot() +
		layer_spatial(this_pred_vect, aes(fill = delta), color = NA) +
		layer_spatial(nam, color = 'gray40', fill = NA, linewidth = 0.3) +
		layer_spatial(counties_with_ag, pch = 16, color = alpha('gray20', 0.5), size = 0.35) +
		layer_spatial(data_traits$site_vect_facet, pch = 3, size = 3) +
		layer_spatial(range_core_sq, color = 'gray', fill = NA, linewidth = 1) +
		layer_spatial(range_core_fut, color = 'orange', fill = NA, linewidth = 1) +
		scale_fill_gradient2(
			name = legend_title,
			low = '#7b3294',
			mid = 'beige',
			high = '#008837',
			midpoint = 0,
			limits = resp_limits
			# labels = scales::label_percent(digits = 2)
		) +
		coord_sf(xlim = c(extent_coords[1], extent_coords[2]), ylim = c(extent_coords[3], extent_coords[4]), expand = FALSE) +
		ggtitle(title, subtitle = subtitle) +
		theme(
			plot.title = element_text(size = 16),
			plot.subtitle = element_text(size = 14)
		)

	filename <- paste0(out_dir, '/map_abundance_change_', filename_append, '.png')
	ggsave(plot = map, filename = filename, width = 12, height = 10, dpi = 300)
	invisible(map)

}

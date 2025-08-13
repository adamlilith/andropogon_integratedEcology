#' Map of abundance model predictions
#'
#' @param out_dir Folder in which to save image.
#' @param filename_append String to append to file name. Will be concatenated as per `paste0(out_dir, '/map_abundance_', filename_append, '.png')``.
#' @param pred_vect_nam SpatVector with predictions.
#' @param response_var Name of response variable: 'N_ag_county_sq', 'N_ag_county_ssp245_2041_2070', 'N_ag_county_ssp245_2071_2100', 'N_ag_county_ssp370_2041_2070', 'N_ag_county_ssp370_2071_2100'.
#' @param data_occs From prepare_occurrences().
#' @param title Character for plot title.
#' @param subtitle Character for plot subtitle.
#' @param ag_core_quant Quantile used to delineate core from non-core.
map_occurrence <- function(
	out_dir,
	filename_append,
	pred_vect_nam,
	response_var,
	data_occs,
	title = 'AG abundance',
	subtitle = NULL,
	ag_core_quant = 0.95
) {

	nam <- vect('./data_from_gadm/gadm_4pt1_level_1_north_america_sans_alaska_lambert.gpkg')
	data_traits <- prepare_nonbiomass_traits(trait = 'height', formula = ~ 1)

	# extent
	site_vect <- data_traits$site_vect_traits
	site_vect <- project(site_vect, pred_vect_nam)
	extent <- ext(site_vect)
	extent <- as.polygons(extent, crs = pred_vect_nam)
	extent <- buffer(extent, width = 200 * 1000) # nominal plot extent
	extent_display <- buffer(extent, width = 300 * 1000) # larger than plot extent
	extent <- ext(extent)
	extent <- as.vector(extent)

	pred_vect_display <- crop(pred_vect_nam, extent_display)

	# delineate current range "core"
	range_core <- delineate_range_core(pred_vect_display, column = response_var, ag_core_quant = ag_core_quant)

	# get range of values for plotting
	vars <- c('N_ag_county_sq', 'N_ag_county_ssp245_2041_2070', 'N_ag_county_ssp245_2071_2100', 'N_ag_county_ssp370_2041_2070', 'N_ag_county_ssp370_2071_2100')

	max_val <- -Inf
	for (var in vars) {

		x <- pred_vect_display[[var]]
		x <- unlist(x)
		x <- x[pred_vect_display$focal_region]
		max_val <- max(max_val, x[!is.infinite(x)])

	}

	# resp_limits <- c(0.1, max_val)
	resp_limits <- c(0, max_val)

	counties_with_ag <- pred_vect_display[pred_vect_display$n_andropogon_gerardi > 0]
	counties_with_ag <- centroids(counties_with_ag)

	legend_title <- 'Abundance'

	map <- ggplot() +
		layer_spatial(pred_vect_display, aes(fill = .data[[response_var]]), color = NA) +
		layer_spatial(nam, color = 'gray40', fill = NA, linewidth = 0.3) +
		layer_spatial(counties_with_ag, pch = 16, color = alpha('gray20', 0.5), size = 0.35) +
		layer_spatial(range_core, color = 'orange2', fill = NA, linewidth = 1) +
		layer_spatial(data_traits$site_vect_traits, pch = 3, size = 4) +
		# scale_fill_continuous(
		# 	name = legend_title,
		# 	type = 'viridis',
		# 	limits = resp_limits,
		# 	trans = 'log10'
		# ) +
		scale_fill_gradientn(
			name = legend_title,
			colors = c('#edf8e9', '#74c476', '#005a32'),
			limits = resp_limits
			# trans = 'log10'
		) +
		# scale_fill_viridis_c(
		# 	name = legend_title,
		# 	# type = 'viridis',
		# 	option = 'inferno',
		# 	limits = resp_limits,
		# 	trans = 'log10'
		# ) +
		xlim(extent[1], extent[2]) + ylim(extent[3], extent[4]) +
		ggtitle(title, subtitle) +
		theme(
			plot.title = element_text(size = 16),
			plot.subtitle = element_text(size = 14)
		)

	filename <- paste0(out_dir, '/map_abundance_', filename_append, '.png')
	ggsave(plot = map, filename = filename, width = 12, height = 10, dpi = 600)
	invisible(map)

}

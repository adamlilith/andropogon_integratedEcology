#' Map biomass model predictions
#'
#' @param out_dir Folder in which to save map.
#' @param filename_append String to append to file name. Will have "_" prefixed to it.
#' @param pred_vect_nam SpatVector with predictions.
#' @param response_var Any of: 'mu_biomass_county_sq', 'mu_biomass_county_ssp245_2041_2070', 'mu_biomass_county_ssp245_2071_2100', 'mu_biomass_county_ssp370_2041_2070', 'mu_biomass_county_ssp370_2071_2100' OR 'biomass_sigma_county_sq'
#' @param response_var_type Either 'mu' or 'sigma', depending on `response_var`.
#' @param data_occs From prepare_occurrences().
#' @param data_biomass From prepare_biomass().
#' @param title Character for plot title.
#' @param subtitle Character for plot subtitle.
#' @param plot_range_core If `TRUE`, extract from `pred_vect_nam` the range core based on `N_ag_county_sq` and plot it.
#' @param ag_core_quant Quantile used to delineate core from non-core.
map_biomass <- function(
	out_dir,
	filename_append,
	pred_vect_nam,
	response_var,
	response_var_type,
	data_occs,
	data_biomass,
	title = 'Biomass',
	subtitle = NULL,
	legend_title = 'Biomass',
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
	extent_display <- buffer(extent, width = 300 * 1000) # larger than plot extent
	extent <- ext(extent)
	extent <- as.vector(extent)

	pred_vect_display <- crop(pred_vect_nam, extent_display)

	# delineate current range "core"
	if (plot_range_core) range_core <- delineate_range_core(pred_vect_display, column = response_var, ag_core_quant = ag_core_quant)

	# get range of values for plotting
	if (response_var_type == 'mu') {
		vars <- c('mu_biomass_county_mean_sq', 'mu_biomass_county_mean_ssp245_2041_2070', 'mu_biomass_county_mean_ssp245_2071_2100', 'mu_biomass_county_mean_ssp370_2041_2070', 'mu_biomass_county_mean_ssp370_2071_2100')
	} else if (response_var_type == 'sigma') {
		vars <- 'biomass_sigma_county_mean_sq'
	}

	# Calculate mean and standard deviation of biomass by site
	biomass_stats <- data_biomass$raw_data_biomass[ , .(biomass_mean = mean(Biomass), biomass_sd = sd(Biomass)), by = SITE]

	if (response_var_type == 'mu') {
		max_val <- 1.05 * max(biomass_stats$biomass_mean)
	} else {
		max_val <- 1.05 * max(biomass_stats$biomass_sd)
	}
	for (var in vars) {

		x <- pred_vect_display[[var]]
		x <- unlist(x)
		x <- x[pred_vect_display$focal_region]
		
		# max_val <- max(max_val, 1.05 * quantile(x, 0.99, na.rm = TRUE))
		max_val <- max(max_val, 1 * quantile(x, 0.999, na.rm = TRUE))

	}

	resp_limits <- c(0, max_val)

	counties_with_ag <- pred_vect_display[pred_vect_display$n_andropogon_gerardi > 0]
	counties_with_ag <- centroids(counties_with_ag)

	if (response_var_type == 'mu') {

		map <- ggplot() +
			layer_spatial(pred_vect_display, aes(fill = .data[[response_var]]), color = NA) +
			# scale_fill_gradientn(
			# 	name = legend_title,
			# 	colors = c('#fcfbfd', '#6a51a3', '#3f007d'),
			# 	limits = resp_limits
			# 	# trans = 'log10'
			# ) +
			scale_fill_gradientn(
				name = legend_title,
				colors = c('#edf8e9', '#74c476', '#005a32'),
				limits = resp_limits
				# trans = 'log10'
			) +
			layer_spatial(nam, color = 'gray40', fill = NA, linewidth = 0.3) +
			layer_spatial(counties_with_ag, pch = 16, color = alpha('gray20', 0.5), size = 0.35) +
			layer_spatial(data_biomass$site_vect_biomass, aes(fill = biomass_mean), pch = 21, size = 4) +
			xlim(extent[1], extent[2]) + ylim(extent[3], extent[4]) +
			ggtitle(title, subtitle) +
			theme(
				plot.title = element_text(size = 16),
				plot.subtitle = element_text(size = 14)
			)

	} else if (response_var_type == 'sigma') {

		map <- ggplot() +
			layer_spatial(pred_vect_display, aes(fill = .data[[response_var]]), color = NA) +
			scale_fill_gradientn(
				name = legend_title,
				colors = c('#ffffcc', '#fd8d3c', '#800026'),
				limits = resp_limits
				# trans = 'log10'
			) +
			# scale_fill_continuous(
			# 	name = legend_title,
			# 	type = 'viridis',
			# 	trans = 'log10',
			# 	limits = resp_limits,
			# 	labels = scales::label_log(digits = 1),
			# 	direction = -1
			# ) +
			layer_spatial(nam, color = 'gray40', fill = NA, linewidth = 0.3) +
			layer_spatial(counties_with_ag, pch = 16, color = alpha('gray20', 0.5), size = 0.35) +
			layer_spatial(data_biomass$site_vect_biomass, aes(fill = biomass_sd), pch = 21, size = 4) +
			xlim(extent[1], extent[2]) + ylim(extent[3], extent[4]) +
			ggtitle(title, subtitle) +
			theme(
				plot.title = element_text(size = 16),
				plot.subtitle = element_text(size = 14)
			)
	
	}

	if (plot_range_core) map <- map + layer_spatial(range_core, color = 'cyan', fill = NA, linewidth = 1)

	ggsave(plot = map, filename = paste0(out_dir, '/map_biomass_', response_var_type, '_mean_', filename_append, '.png'), width = 12, height = 10, dpi = 600)
	invisible(map)

}

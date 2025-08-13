#' Map non-biomass model predictions
#'
#' @param out_dir Folder in which to save map.
#' @param filename_append String to append to file name. Will have "_" prefixed to it.
#' @param pred_vect_nam SpatVector with predictions.
#' @param trait "R-friendly" name of trait
#' @param response_var Any of: '<trait>_mu_county_sq', '<trait>_mu_county_ssp245_2041_2070', '<trait>_mu_county_ssp245_2071_2100', '<trait>_mu_county_ssp370_2041_2070', or '<trait>_mu_county_ssp370_2071_2100'
#' @param data_traits From prepare_biomass().
#' @param title Character for plot title.
#' @param subtitle Character for plot subtitle.
#' @param plot_range_core If `TRUE`, extract from `pred_vect_nam` the range core based on `N_ag_county_sq` and plot it.
#' @param ag_core_quant Quantile used to delineate core from non-core.
map_nonbiomass_trait <- function(
	out_dir,
	filename_append,
	pred_vect_nam,
	trait,
	response_var,
	data_traits,
	title = '<TRAIT>',
	subtitle = NULL,
	plot_range_core = TRUE,
	ag_core_quant = 0.95
) {

	nam <- vect('./data_from_gadm/gadm_4pt1_level_1_north_america_sans_alaska_lambert.gpkg')

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

	# # delineate current range "core"
	if (plot_range_core) range_core <- delineate_range_core(pred_vect_display, column = 'N_ag_county_sq', ag_core_quant = ag_core_quant)

	# get range of values for plotting
	vars <- c('_mu_county_sq', '_mu_county_ssp245_2041_2070', '_mu_county_ssp245_2071_2100', '_mu_county_ssp370_2041_2070', '_mu_county_ssp370_2071_2100')
	vars <- paste0(trait, vars)

	legend_title <- get_nice_trait(trait)$legend_title

	# Calculate mean and standard deviation of trait by site
	observed <- rep(NA_real_, data_traits$n_pheno_sites)
	sites <- unique(data_traits$raw_data_traits$SITE)
	for (i in seq_along(sites)) {
	
		site <- sites[i]
		trait_column_name <- get_raw_trait_name(trait)
		observed[i] <- mean(data_traits$raw_data_traits[[trait_column_name]][data_traits$raw_data_traits$SITE == site])

	}
	min_val <- 0.95 * min(observed)
	max_val <- 1.05 * max(observed)

	for (var in vars) {

		x <- pred_vect_display[[var]]
		x <- unlist(x)
		x <- x[pred_vect_display$focal_region]
		
		min_val <- min(min_val, 0.95 * min(x, na.rm = TRUE))
		# max_val <- max(max_val, 1.05 * quantile(x, 0.99, na.rm = TRUE))
		max_val <- max(max_val, 1.05 * max(x, na.rm = TRUE))

	}

	resp_limits <- c(min_val, max_val)

	counties_with_ag <- pred_vect_display[pred_vect_display$n_andropogon_gerardi > 0]
	counties_with_ag <- centroids(counties_with_ag)

	trait_mean_column <- paste0(trait, '_mean')

	map <- ggplot() +
		layer_spatial(pred_vect_display, aes(fill = .data[[response_var]]), color = NA) +
		scale_fill_continuous(
			name = legend_title,
			type = 'viridis',
			limits = resp_limits,
			trans = 'log10',
			labels = scales::label_log(digits = 2)
		) +
		layer_spatial(nam, color = 'gray40', fill = NA, linewidth = 0.3) +
		layer_spatial(counties_with_ag, pch = 16, color = alpha('gray20', 0.5), size = 0.35) +
		layer_spatial(data_traits$site_vect_traits, aes(fill = .data[[trait_mean_column]]), pch = 21, size = 4) +
		xlim(extent[1], extent[2]) + ylim(extent[3], extent[4]) +
		ggtitle(title, subtitle) +
		theme(
			plot.title = element_text(size = 16),
			plot.subtitle = element_text(size = 14)
		)


	if (plot_range_core) map <- map + layer_spatial(range_core, color = 'cyan', fill = NA, size = 2)

	ggsave(plot = map, filename = paste0(out_dir, '/map_', trait, '_mu_', filename_append, '.png'), width = 12, height = 10, dpi = 600)
	invisible(map)

}


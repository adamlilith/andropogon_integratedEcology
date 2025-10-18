#' Graph response of occurrence (abundance) to biomass
#'
#' For graphing response of abundance against the estimated biomass.
#' @param out_dir Folder in which to save graphs.
#' @param chains MCMC chains
#' @param data_response_curves_occs_vs_biomass Response curve array or matrix with biomass values for occurrence-vs-biomass responses. Output from create_response_curve_array_occ_vs_biomass().
#' @param data_biomass From prepare_biomass().
#' @param quant_threshold Used to define upper limit of the y axis based on range of values of the upper CI limits.
graph_response_curves_occurrence_vs_biomass <- function(
	out_dir,
	chains,
	data_response_curves_occs_vs_biomass,
	data_biomass = NULL,
	quant_threshold = 0.5
) {

	y_lab <- bquote('Expected abundance')
	nice_title <- bquote('Abundance versus biomass')

	# unscaled predictor value... note, we need to do this the same
	x <- data_response_curves_occs_vs_biomass$unscaled_biomass

	param <-'response_curve_occs_vs_biomass'
	response_mean <- hammer_extract(chains, param = param, j = 1:n_response_curve_values, stat = 'mean')
	response_lower <- hammer_extract(chains, param = param, j = 1:n_response_curve_values, stat = 'lower')
	response_upper <- hammer_extract(chains, param = param, j = 1:n_response_curve_values, stat = 'upper')

	# data frames to hold predictions in long format
	df_mean <- data.frame(
		x = x,
		response = response_mean
	)

	df_lower <- data.frame(
		x = x,
		response = response_lower
	)

	df_upper <- data.frame(
		x = x,
		response = response_upper
	)

	this_df_ci <- data.frame(
		x = c(x, rev(x)),
		response = c(df_upper$response, rev(df_lower$response))
	)

	max_val <- max(df_mean$response, quantile(df_upper$response[!is.infinite(df_upper$response)], quant_threshold))

	response <- ggplot() +
		geom_polygon(
			data = this_df_ci,
			mapping = aes(x = x, y = response),
			color = NA,
			fill = alpha('blue', 0.1)
		) +
		geom_line(
			data = df_mean,
			mapping = aes(x = x, y = response)
		) +
		coord_cartesian(ylim = c(0, max_val)) +
		xlab('Mean (μ) site-level biomass (g)') +
		ylab(paste0('Expected abundance (λ)')) +
		ggtitle(nice_title) +
		theme(
			legend.position = 'none',
			plot.title = element_text(size = 10)
		)

	ggsave(plot = response, filename = paste0(out_dir, '/response_curves_abundance_vs_biomass.png'), width = 12, height = 8, dpi = 600)
	invisible(response)

}


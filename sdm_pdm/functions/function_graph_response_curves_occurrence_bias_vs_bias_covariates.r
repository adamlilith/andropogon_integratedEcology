#' Graph of response curves for occurrence vs environment
#' 
#' @param out_dir Folder in which to save graphs.
#' @param chains MCMC chains
#' @param data_occs From prepare_occurrence_data().
graph_response_curves_occurrence_bias_vs_bias_covariates <- function(out_dir, chains, data_occs) {

	n_covariates <- data_occs$n_covariates_occs_bias
	covariates <- data_occs$covariates_occs_bias
	resp_curves_x <- data_occs$resp_curves_w_occs
	resp_curves_unscaled <- data_occs$resp_curves_w_occs_bias_unscaled
	centers <- data_occs$w_centers_occs_bias
	scales <- data_occs$w_scales_occs_bias
	y_lab <- bquote('Sampling Rate')

	# make a plot for how biomass mu and sigma respond to each predictor
	responses <- list()
	for (i in 1:n_covariates) {

		pred <- covariates[i]

		if (pred == 'area_km2_log10') {
			nice_title <- 'County area'
			nice_axis <- 'log₁₀(county area, km²)'
		} else if (pred == 'n_poaceae_log10p1') {
			nice_title <- 'Number of Poaceae'
			nice_axis <- 'log10(number of Poaceae + 1)'
		}
		nice_title <- bquote('Thinning rate' * ' versus ' * .(nice_title))

		# unscaled predictor value
		x <- resp_curves_unscaled[ , pred]
		
		param <-'response_curves_occs_bias'
		if (n_covariates == 1) {
			
			response_mean <- mc_extract(chains, param = param, j = TRUE, stat = 'mean')
			response_lower <- mc_extract(chains, param = param, j = TRUE, stat = 'lower')
			response_upper <- mc_extract(chains, param = param, j = TRUE, stat = 'upper')

		} else if (n_covariates > 1) {
		
			response_mean <- mc_extract(chains, param = param, j = TRUE, k = i, stat = 'mean')
			response_lower <- mc_extract(chains, param = param, j = TRUE, k = i, stat = 'lower')
			response_upper <- mc_extract(chains, param = param, j = TRUE, k = i, stat = 'upper')
	
		}

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
			xlab(nice_axis) +
			ylab(y_lab) +
			ylim(0, 1) +
			ggtitle(nice_title) +
			theme(
				legend.position = 'none',
				plot.title = element_text(size = 9)
			)

		responses[[i]] <- response

	} # next predictor

	if (n_covariates == 1) {
		nrow <- 1
		width <- 10
		height <- 8
	} else if (n_covariates <= 2) {
		nrow <- 1
		width <- 8
		height <- 4
	} else if (n_covariates <= 3) {
		nrow <- 1
		width <- 14
		height <- 4
	} else {
		nrow <- 2
		width <- 16
		height <- 10	
	}

	responses <- plot_grid(plotlist = responses, nrow = nrow)
	ggsave(plot = responses, filename = paste0(out_dir, '/response_curves_abundance_bias.png'), width = width, height = height, dpi = 600)
	invisible(responses)

}

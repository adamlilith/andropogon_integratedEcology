#' Graph of response curves for occurrence vs environment
#' 
#' out_dir			Folder in which to save graphs.
#' chains			MCMC chains
#' log_precip		If `TRUE`, bios 12-14 and 16-19 were logged
#' zero_inflated	If `TRUE`, model is zero-inflated
graph_response_curves_occurrence_vs_environment <- function(out_dir, chains, log_precip, zero_inflated) {

	y_lab <- 'Predicted Relative Abundance'

	data_occs <- prepare_occurrence_data(formula_occs = formula_occs, formula_occs_bias = formula_occs_bias, log_precip = log_precip, n_response_curve_values = n_response_curve_values, psa_quant = psa_quant, calib = calib)

	covariates <- data_occs$covariates_occs
	n_covariates <- length(covariates)
	x_array <- data_occs$resp_curves_x_occs

	y_max <- -Inf
	responses <- list()
	for (i in seq_along(covariates)) {

		covariate <- covariates[i]
		nice <- get_nice_predictor(covariate)
		nice_title <- nice$short
		nice_axis <- nice$long
		nice_title <- bquote('Abundance' * ' versus ' * .(nice_title))

		if (n_covariates == 1) {
			x <- x_array
		} else {
			x <- x_array[ , , i]
		}

		# predict multiple times to smooth over RNG
		times <- if (trial) { 3 } else { 30 }
		response_mean <- response_lower <- response_upper <- matrix(NA_real_, nrow = times, ncol = nrow(x))
		for (time in 1:times) {

			preds <- predict_occs(chains = chains, x = x, zero_inflated = zero_inflated, x_psi = x, times = 1)

			response_mean[time, ] <- colMeans(preds)
			response_lower[time, ] <- apply(preds, 2, quantile, 0.05)
			response_upper[time, ] <- apply(preds, 2, quantile, 0.95)

		}

		response_mean <- colMeans(response_mean)
		response_lower <- colMeans(response_lower)
		response_upper <- colMeans(response_upper)

		# unscaled predictor value
		x_plot <- data_occs$resp_curves_x_occs_unscaled[ , covariate]
		if (log_precip & covariate %in% paste0('bio', c(12:14, 16:19))) x_plot <- 10^x_plot - 1

		# data frames to hold predictions in long format
		df_mean <- data.frame(
			x = x_plot,
			response = response_mean
		)

		df_ci <- data.frame(
			x = c(x_plot, rev(x_plot)),
			response = c(response_upper, rev(response_lower))
		)

		responses[[i]] <- ggplot() +
			geom_polygon(
				data = df_ci,
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
			ggtitle(nice_title) +
			theme(
				legend.position = 'none',
				plot.title = element_text(size = 9)
			)

		y_max <- max(y_max, df_ci$response)

	}

	for (i in 1:n_covariates) {
		responses[[i]] <- responses[[i]] + coord_cartesian(ylim = c(0, y_max)) 
	}

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
	ggsave(plot = responses, filename = paste0(out_dir, '/response_curves_abundance_mean_inner_quant.png'), width = width, height = height, dpi = 120)
	invisible(responses)

}

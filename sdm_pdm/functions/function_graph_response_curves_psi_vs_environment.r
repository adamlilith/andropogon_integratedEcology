#' Graph of response curves for psi vs environment
#' 
#' out_dir			Folder in which to save graphs.
#' chains			MCMC chains
graph_response_curves_psi_vs_environment <- function(out_dir, chains) {

	y_lab <- 'Predicted Probability of Occurrence'

	data_psi <- prepare_occurrence_data(formula_occs = formula_psi, formula_bias = ~ 1, n_response_curve_values = n_response_curve_values, psa_quant = psa_quant, calib = calib)

	covariates <- data_psi$covariates
	n_covariates <- length(covariates)
	x_array <- data_psi$resp_curves_x

	responses <- list()
	for (i in seq_along(covariates)) {

		covariate <- covariates[i]
		nice <- get_nice_predictor(covariate)
		nice_title <- nice$short
		nice_axis <- nice$long
		nice_title <- bquote('Probability of Occurrence' * ' versus ' * .(nice_title))

		if (n_covariates == 1) {
			x <- x_array
		} else {
			x <- x_array[ , , i]
		}

		preds <- predict_psi(chains = chains, x = x)

		means <- colMeans(preds)
		hdi <- apply(preds, 2, function(x) hdi(x, credMass = 0.90))
		lowers <- hdi[1, ]
		uppers <- hdi[2, ]

		# unscaled predictor value
		x_plot <- data_psi$resp_curves_x_unscaled[ , covariate]
		if (covariate %in% paste0('bio', c(12:14, 16:19), '_log10p1')) x_plot <- 10^x_plot + 1

		# data frames to hold predictions in long format
		df_mean <- data.frame(
			x = x_plot,
			response = means
		)

		df_ci <- data.frame(
			x = c(x_plot, rev(x_plot)),
			response = c(uppers, rev(lowers))
		)

		responses[[i]] <- ggplot() +
			geom_polygon(
				data = df_ci,
				mapping = aes(x = x, y = response),
				color = NA,
				fill = 'lightblue',
				alpha = 0.5
			) +
			geom_line(
				data = df_mean,
				mapping = aes(x = x, y = response)
			) +
			coord_cartesian(ylim = c(0, 1)) +
			xlab(nice_axis) +
			ylab(y_lab) +
			ggtitle(nice_title) +
			theme(
				legend.position = 'none',
				plot.title = element_text(size = 9)
			)

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
	ggsave(plot = responses, filename = paste0(out_dir, '/response_curves_psi_mean_inner_quant.png'), width = width, height = height, dpi = 120)
	invisible(responses)

}

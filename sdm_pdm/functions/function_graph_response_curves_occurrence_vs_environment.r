#' Graph of response curves for occurrence vs environment
#' 
#' out_dir			Folder in which to save graphs. If NULL, do not save file.
#' chains			MCMC chains
#' zero_inflated	If `TRUE`, model is zero-inflated
#'
#' source('C:/Kaji/R/andropogon_integratedEcology/sdm_pdm/functions/function_graph_response_curves_occurrence_vs_environment.r')
graph_response_curves_occurrence_vs_environment <- function(out_dir, chains, zero_inflated, overdispersed = FALSE) {

	y_lab <- 'Relative Abundance'

	data_occs <- prepare_occurrence_data(formula_occs = formula_occs, formula_bias = ~1, n_response_curve_values = n_response_curve_values, psa_quant = psa_quant, calib = calib)

	covariates <- data_occs$covariates
	n_covariates <- length(covariates)
	x_occs_array <- data_occs$resp_curves_x

	y_max <- -Inf
	responses <- list()
	for (i in seq_along(covariates)) {

		covariate <- covariates[i]
		nice <- get_nice_predictor(covariate)
		nice_title <- nice$short
		nice_axis <- nice$long
		nice_title <- bquote('Abundance' * ' versus ' * .(nice_title))

		if (n_covariates == 1) {
			x <- x_occs_array
		} else {
			x <- x_occs_array[ , , i]
		}

		# predict multiple times to smooth over RNG
		preds <- predict_occs(chains = chains, x = x, x_psi = x, overdispersed = overdispersed)

		# centers <- apply(preds, 2, median)
		centers <- apply(preds, 2, mean)
		hdi <- apply(preds, 2, function(x) hdi(x, credMass = 0.90))
		lowers <- hdi[1, ]
		uppers <- hdi[2, ]

		# unscaled predictor value
		x_unscaled <- data_occs$resp_curves_x_unscaled[ , covariate]
		if (covariate %in% paste0('bio', c(12:14, 16:19), '_log10p1')) x_unscaled <- 10^x_unscaled + 1

		# data frames to hold predictions in long format
		df_center <- data.frame(
			x = x_unscaled,
			response = centers
		)

		df_ci <- data.frame(
			x = c(x_unscaled, rev(x_unscaled)),
			response = c(lowers, rev(uppers))
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
				data = df_center,
				mapping = aes(x = x, y = response)
			) +
			xlab(nice_axis) +
			ylab(y_lab) +
			ggtitle(nice_title) +
			theme(
				legend.position = 'none',
				plot.title = element_text(size = 10)
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
	if (!is.null(out_dir)) ggsave(plot = responses, filename = paste0(out_dir, '/response_curves_abundance_mean_hdpi.png'), width = width, height = height, dpi = 120)
	invisible(responses)

}

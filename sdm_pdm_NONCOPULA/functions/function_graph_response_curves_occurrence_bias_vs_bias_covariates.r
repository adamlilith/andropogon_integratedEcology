#' Graph of response curves for occurrence vs environment
#' 
#' @param out_dir Folder in which to save graphs.
#' @param chains MCMC chains
#' @param data_occs From prepare_occurrence_data().
graph_response_curves_occurrence_bias_vs_bias_covariates <- function(out_dir, chains) {

	data_occs <- prepare_occurrence_data(formula_occs = formula_occs, formula_occs_bias = formula_occs_bias, log_precip = log_precip, n_response_curve_values = n_response_curve_values, psa_quant = psa_quant, calib = calib)

	n_covariates <- data_occs$n_covariates_occs_bias
	covariates <- data_occs$covariates_occs_bias
	centers <- data_occs$w_centers_occs_bias
	scales <- data_occs$w_scales_occs_bias
	y_lab <- 'Sampling Rate'
	w_array <- data_occs$resp_curves_w_occs
	w_array_unscaled <- data_occs$resp_curves_w_occs_bias_unscaled

	# make a plot for how biomass mu and sigma respond to each predictor
	responses <- list()
	for (i in 1:n_covariates) {

		pred <- covariates[i]

		if (pred == 'area_km2_log10') {
			nice_title <- 'County Area'
			nice_axis <- 'log₁₀(County Area, km²)'
		} else if (pred == 'n_poaceae_log10p1') {
			nice_title <- 'Number of Poaceae'
			nice_axis <- 'log10(Number of Poaceae + 1)'
		}
		nice_title <- bquote('Thinning Rate' * ' versus ' * .(nice_title))

		# unscaled predictor value
		if (n_covariates == 1) {
			w <- w_array
		} else {
			w <- w_array[ , , i]
		}
		w_unscaled <- w_array_unscaled[ , pred]

		n_iters <- nrow(chains$samples[[1]])
		n_chains <- mc_n_chains(chains)
		preds <- matrix(NA_real_, nrow = n_iters * n_chains, ncol = nrow(w))
		alphas_occs <- mc_subset(chains, param = 'alpha_occs', j = TRUE)
		k <- 1
		for (chain in 1:n_chains) {
			for (n_iter in 1:n_iter) {

				this_alpha <- alphas_occs$samples[[chain]][n_iter, ]
				this_alpha <- cbind(this_alpha)

				preds_untrans <- w %*% this_alpha
				preds_untrans <- preds_untrans[ , 1]
				pred <- expit(preds_untrans)

				preds[k, ] <- pred

				k <- k + 1

			}
		}

		response_mean <- colMeans(preds)
		response_lower <- apply(preds, 2, quantile, 0.1)
		response_upper <- apply(preds, 2, quantile, 0.9)

		# data frames to hold predictions in long format
		df_mean <- data.frame(
			x = w_unscaled,
			response = response_mean
		)

		df_ci <- data.frame(
			x = c(w_unscaled, rev(w_unscaled)),
			response = c(response_upper, rev(response_lower))
		)

		response <- ggplot() +
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

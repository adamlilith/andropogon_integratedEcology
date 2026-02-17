#' Graph of response curves for occurrence vs environment
#' 
#' @param out_dir Folder in which to save graphs.
#' @param resp_type "mu" or "sigma" or "psi"
#' @param chains MCMC chains
#' @param data From prepare_occurrence_data().
#' @param log_precip If `TRUE`, bios 12-14 and 16-19 were logged
graph_response_curves_occurrence_vs_environment <- function(out_dir, resp_type, chains, data, log_precip) {

	n_covariates <- data$n_covariates_occs
	covariates <- data$covariates_occs
	resp_curves_x <- data$resp_curves_x_occs
	resp_curves_unscaled <- data$resp_curves_x_occs_unscaled
	centers <- data$x_centers_occs
	scales <- data$x_scales_occs
	
	if (resp_type == 'mu') {
		param <-'response_curves_occs_mu'
		y_lab <- bquote('Expected Abundance')
	} else if (resp_type == 'sigma') {
		param <-'response_curves_occs_sigma'
		y_lab <- bquote('Standard Deviation of Log Abundance')
	} else if (resp_type == 'psi') {
		param <-'response_curves_psi'
		y_lab <- bquote('Probability of Presence')
	}

	# make a plot for how biomass mu and sigma respond to each predictor
	responses <- list()
	# max_val <- if (resp_type == 'psi') { 1 } else { -Inf }
	max_val <- -Inf
	for (i in 1:n_covariates) {

		pred <- covariates[i]
		nice <- get_nice_predictor(pred)
		nice_title <- nice$short
		nice_axis <- nice$long

		if (resp_type == 'mu') {
			nice_title <- bquote('Abundance' * ' versus ' * .(nice_title))
		} else if (resp_type == 'sigma') {
			nice_title <- bquote('Standard Deviation in Abundance' * ' versus ' * .(nice_title))
		} else if (resp_type == 'psi') {
			nice_title <- bquote('Probability of Presence' * ' versus ' * .(nice_title))
		}

		# unscaled predictor value
		x <- resp_curves_unscaled[ , pred]

		if (log_precip & pred %in% paste0('bio', c(12:14, 16:19))) x <- 10^x - 1
		
		response_mean <- mc_extract(chains, param = param, j = 1:n_response_curve_values, k = i, stat = 'mean')
		response_lower <- mc_extract(chains, param = param, j = 1:n_response_curve_values, k = i, stat = 'lower')
		response_upper <- mc_extract(chains, param = param, j = 1:n_response_curve_values, k = i, stat = 'upper')

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
			ggtitle(nice_title) +
			theme(
				legend.position = 'none',
				plot.title = element_text(size = 9)
			)

			# if (resp_type != 'psi') {

				# mmax <- max(df_mean$response, na.rm = TRUE)
				# if (df_mean$response[1] == mmax | df_mean$response[nrow(df_mean)] %==na% mmax) {
					# mmax <- 1 * mmax
				# } else {
					mmax <- if (resp_type == 'mu') {
						quantile(df_upper$response, 0.9, na.rm = TRUE)
					} else if (resp_type == 'psi') {
						min(1, 1.05 * max(df_upper$response, na.rm = TRUE))
					} else if (resp_type == 'sigma') {
						1.05 * max(df_upper$response, na.rm = TRUE)
					}
				# }

				max_val <- max(max_val, mmax)

			# }

		responses[[i]] <- response

	} # next predictor

	for (i in 1:n_covariates) {
		responses[[i]] <- responses[[i]] + coord_cartesian(ylim = c(0, max_val)) 
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
	ggsave(plot = responses, filename = paste0(out_dir, '/response_curves_abundance_', resp_type, '.png'), width = width, height = height, dpi = 600)
	invisible(responses)

}

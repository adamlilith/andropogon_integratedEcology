#' Graph response of biomass (heteroscedastic gamma) to the environment
#'
#' For graphing response of biomass against the environment, assuming biomass ~ dgamma() where the shape and rate are functions of the environment.
#' @param out_dir Folder in which to save graphs.
#' @param chains MCMC chains
#' @param resp_type "mu" or "sigma" or "pzero"
#' @param data_biomass_mu From prepare_biomass().
#' @param quant_threshold Used to define upper limit of the y axis based on range of values of the upper CI limits.
graph_response_curves_biomass_vs_environment <- function(out_dir, chains, resp_type, data_biomass_mu, quant_threshold = 0.5) {

	n_covariates <- data_biomass_mu$n_covariates_biomass
	covariates <- data_biomass_mu$covariates_biomass
	resp_curves_x <- data_biomass_mu$resp_curves_x_biomass
	resp_curves_unscaled <- data_biomass_mu$resp_curve_x_biomass_unscaled
	centers <- data_biomass_mu$x_centers_biomass
	scales <- data_biomass_mu$x_scales_biomass

	if (resp_type == 'mu') {
		ylab <- bquote('Mean Biomass (g)')
		max_val <- max(data_biomass_mu$y_biomass)
		response_type_nice <- 'Mean (μ)'
		param <- paste0('response_curves_mu')
	} else if (resp_type == 'pzero') {
		ylab <- bquote('Probability of Zero (' * italic(p)[0] * ')')
		max_val <- -Inf
		response_type_nice <- bquote(italic(p)[0])
		param <- paste0('response_curves_pzero')
	}

	# make a plot for how biomass mu and sigma respond to each predictor
	responses <- list()
	for (i in 1:n_covariates) {

		pred <- covariates[i]

		nice_names <- get_nice_predictor(pred)
		nice_pred <- nice_names$short
		nice_axis <- nice_names$long

		nice_title <- bquote('Biomass ' * .(response_type_nice) * ' versus ' * .(nice_pred))

		# unscaled predictor value
		x <- resp_curves_unscaled[ , pred]

		# create data frame with SDM response
		if (n_covariates == 1) {

			response_mean <- hammer_extract(chains, param = param, j = 1:n_response_curve_values, stat = 'mean')
			response_lower <- hammer_extract(chains, param = param, j = 1:n_response_curve_values, stat = 'lower')
			response_upper <- hammer_extract(chains, param = param, j = 1:n_response_curve_values, stat = 'upper')

		} else {
		
			response_mean <- hammer_extract(chains, param = param, j = 1:n_response_curve_values, k = i, stat = 'mean')
			response_lower <- hammer_extract(chains, param = param, j = 1:n_response_curve_values, k = i, stat = 'lower')
			response_upper <- hammer_extract(chains, param = param, j = 1:n_response_curve_values, k = i, stat = 'upper')

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
			ylab(ylab) +
			ggtitle(nice_title) +
			theme(
				legend.position = 'none',
				plot.title = element_text(size = 10),
				axis.title = element_text(size = 8),
				axis.text = element_text(size = 8)
			)

		if (resp_type == 'mu') {

			data_TEMP_by_site <- data_biomass_mu$raw_data_biomass
			if (pred == 'ph') {
				names(data_TEMP_by_site)[names(data_TEMP_by_site) == 'site_soil_pH'] <- 'ph'
				x <- data_biomass_mu$raw_data_biomass[ , .(val = mean(.SD[['site_soil_pH']])), by = SITE][['val']]
			} else {
				x <- data_biomass_mu$raw_data_biomass[ , .(val = mean(.SD[[pred]])), by = SITE][['val']]				
			}
			observed <- data_biomass_mu$raw_data_biomass[ , .(val = mean(Biomass)), by = SITE][['val']]
			data_TEMP_by_site_collated <- data.table(
				x = x,
				y = observed,
				site = unique(data_TEMP_by_site$SITE)
			)

			response <- response + geom_point(
				data = data_TEMP_by_site_collated,
				mapping = aes(x = x, y = y, fill = site),
				pch = 21, size = 3
			)

			data_TEMP_by_plant <- data_biomass_mu$raw_data_biomass
			if (pred == 'ph') {
				data_TEMP_by_plant$x <- data_TEMP_by_plant[['site_soil_pH']]
			} else {
				data_TEMP_by_plant$x <- data_TEMP_by_plant[[pred]]
			}
			data_TEMP_by_plant$y <- data_TEMP_by_plant$Biomass
			data_TEMP_by_plant$site <- data_TEMP_by_plant$SITE

			response <- response + geom_point(
				data = data_TEMP_by_plant,
				mapping = aes(x = x, y = y, fill = site),
				pch = 21, size = 1
			)

		}

		max_val <- max(
			max_val,
			quantile(df_upper$response[!is.infinite(df_upper$response)], quant_threshold)
		)

		responses[[i]] <- response

	} # next predictor

	for (i in 1:n_covariates) {
		responses[[i]] <- responses[[i]] + coord_cartesian(ylim = c(0, max_val))
	}

	# ncol <- ceiling(sqrt(length(resp_type)))
	# width <- ncol * 6
	# height <- 5 * ceiling(length(resp_type) / ncol)

	if (data_biomass_mu$n_covariates == 2) {
		ncol <- 2
	} else {
		ncol <- ceiling(sqrt(length(data_biomass_mu$n_covariates)))
	}
	width <- 6 * ncol
	height <- 5 * ceiling(length(data_biomass_mu$n_covariates) / ncol)

	responses <- plot_grid(plotlist = responses, ncol = ncol)

	filename <- paste0(out_dir, '/response_curves_biomass_', resp_type, '.png')
	ggsave(plot = responses, filename = filename, width = width, height = height, dpi = 600)
	invisible(responses)
	
}


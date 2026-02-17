#' Graph response of non-biomass traits to the environment
#'
#' For graphing response of biomass against the environment, assuming biomass ~ dgamma() where the shape and rate are functions of the environment.
#' facet			Name of facet (R-friendly version)
#' out_dir 			Folder in which to save graphs.
#' chains 			MCMC chains
#' resp_type 		"mu" or "psi"
#' log_precip		If `TRUE`, log BIOs 12-14 and 16-29
#' data_facet 		From prepare_nonbiomass_data().
#' quant_threshold 	Used to define upper limit of the y axis based on range of values of the upper CI limits.
graph_response_curves_nonbiomass_vs_environment <- function(facet, out_dir, chains, resp_type, log_precip, data_facet, quant_threshold = 0.5) {

	facet_raw <- data_facet$facet_raw
	n_covariates <- data_facet$n_covariates_facet
	covariates <- data_facet$covariates_facet
	resp_curves_x <- data_facet$resp_curves_x_facet
	resp_curves_unscaled <- data_facet$resp_curve_x_facet_unscaled
	centers <- data_facet$x_centers_facet
	scales <- data_facet$x_scales_facet

	if (resp_type == 'mu') {
		facet_nice <- get_nice_trait(facet)
		title_facet <- facet_nice$short
		ylab <- facet_nice$long
		max_val <- max(data_facet$y_facet)
		response_type_nice <- 'Mean (μ)'
		param <- paste0('response_curves_facet_mu')
	} else if (resp_type == 'psi') {
		ylab <- bquote('Probability of presence (ψ)')
		title_facet <- 'Probability of Presence'
		max_val <- 1
		response_type_nice <- 'ψ'
		param <- paste0('response_curves_psi')
	}

	# make a plot for how biomass mu and sigma respond to each predictor
	responses <- list()
	for (i in 1:n_covariates) {

		pred <- covariates[i]

		nice_names <- get_nice_predictor(pred)
		nice_pred <- nice_names$short
		nice_axis <- nice_names$long

		nice_title <- bquote(.(title_facet) * ' ' * .(response_type_nice) * ' versus ' * .(nice_pred))

		# unscaled predictor value
		x <- resp_curves_unscaled[ , pred]
		if (pred %in% paste0('bio', c(12:14, 16:19)) & log_precip) x <- 10^x

		# create data frame with SDM response
		if (n_covariates == 1) {

			response_mean <- mc_extract(chains, param = param, j = 1:n_response_curve_values, stat = 'mean')
			response_lower <- mc_extract(chains, param = param, j = 1:n_response_curve_values, stat = 'lower')
			response_upper <- mc_extract(chains, param = param, j = 1:n_response_curve_values, stat = 'upper')

		} else {
		
			response_mean <- mc_extract(chains, param = param, j = 1:n_response_curve_values, k = i, stat = 'mean')
			response_lower <- mc_extract(chains, param = param, j = 1:n_response_curve_values, k = i, stat = 'lower')
			response_upper <- mc_extract(chains, param = param, j = 1:n_response_curve_values, k = i, stat = 'upper')

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

			data_TEMP_by_site <- data_facet$raw_data_facet
			if (pred == 'ph') {
				names(data_TEMP_by_site)[names(data_TEMP_by_site) == 'site_ph'] <- 'ph'
				x <- data_facet$raw_data_facet[ , .(val = mean(.SD[['site_ph']])), by = SITE][['val']]
			} else {
				x <- data_facet$raw_data_facet[ , .(val = mean(.SD[[pred]])), by = SITE][['val']]				
			}
			observed <- data_facet$raw_data_facet[ , .(val = mean(get(facet_raw))), by = SITE][['val']]
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

			data_TEMP_by_plant <- data_facet$raw_data_facet
			if (pred == 'ph') {
				data_TEMP_by_plant$x <- data_TEMP_by_plant[['site_ph']]
			} else {
				data_TEMP_by_plant$x <- data_TEMP_by_plant[[pred]]
			}
			data_TEMP_by_plant$y <- data_TEMP_by_plant[[facet_raw]]
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

	if (data_facet$n_covariates == 2) {
		ncol <- 2
	} else {
		ncol <- ceiling(sqrt(length(data_facet$n_covariates)))
	}
	width <- 6 * ncol
	height <- 5 * ceiling(length(data_facet$n_covariates) / ncol)

	responses <- plot_grid(plotlist = responses, ncol = ncol)

	filename <- paste0(out_dir, '/response_curves_', facet, '_', resp_type, '.png')
	ggsave(plot = responses, filename = filename, width = width, height = height, dpi = 600)
	invisible(responses)
	
}


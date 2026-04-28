#' Graph response of non-biomass facet to the environment
#'
#' source('C:/Kaji/R/andropogon_integratedEcology/sdm_pdm/functions/function_graph_response_curves_biomass_nonbiomass_vs_environment.r')
#' 
#' out_dir 			Folder in which to save graphs.
#' chains 			MCMC chains
#' resp_type 		"mu" or "psi"
#' data_biomass_nonbiomass 		From prepare_biomass_data() or prepare_nonbiomass_data().
#' quant_threshold 	Used to define upper limit of the y axis based on range of values of the upper CI limits.
graph_response_curves_biomass_nonbiomass_vs_environment <- function(
	facet,
	out_dir,
	chains,
	resp_type,
	data_biomass_nonbiomass,
	quant_threshold = 0.5
) {

	if (facet == 'biomass') {

		facet_raw <- get_raw_trait_name_from_rfriendly(facet)

		n_covariates <- data_biomass_nonbiomass$n_covariates
		covariates <- data_biomass_nonbiomass$covariates
		resp_curves_x <- data_biomass_nonbiomass$resp_curves_x
		resp_curves_unscaled <- data_biomass_nonbiomass$resp_curves_x_unscaled
		centers <- data_biomass_nonbiomass$x_centers
		scales <- data_biomass_nonbiomass$x_scales

	} else {
	
		facet_raw <- get_raw_trait_name_from_rfriendly(facet)

		n_covariates <- data_biomass_nonbiomass$n_covariates
		covariates <- data_biomass_nonbiomass$covariates
		resp_curves_x <- data_biomass_nonbiomass$resp_curves_x
		resp_curves_unscaled <- data_biomass_nonbiomass$resp_curves_x_unscaled
		centers <- data_biomass_nonbiomass$x_centers
		scales <- data_biomass_nonbiomass$x_scales

	}

	if (resp_type == 'mu' & facet == 'biomass') {
	
		facet_nice <- get_nice_trait(facet)$short
		
		ylab <- 'Mean Biomass (g)'
		max_val <- max(data_biomass_nonbiomass$y_biomass)
		param <- paste0('response_curves_biomass_mu')
	
	} else if (resp_type == 'mu') {
	
		facet_nice <- get_nice_trait(facet)$short

		ylab <- get_nice_trait(facet)$long
		max_val <- max(data_biomass_nonbiomass$y_facet)
		param <- 'response_curves_facet_mu'
	
	} else if (resp_type == 'psi') {

		facet_nice <- 'Probability of Presence'

		ylab <- bquote('Probability of Presence (ψ)')
		max_val <- 1
		param <- paste0('response_curves_psi')
	
	}

	# make a plot for how biomass mu and sigma respond to each predictor
	responses <- list()
	for (i in 1:n_covariates) {

		pred <- covariates[i]

		nice_names <- get_nice_predictor(pred)
		nice_pred <- nice_names$short
		nice_axis <- nice_names$long

		nice_title <- bquote(.(facet_nice) * ' ' * ' versus ' * .(nice_pred))

		# unscaled predictor value
		x <- resp_curves_unscaled[ , pred]

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

			data_TEMP_by_site <- if (facet == 'biomass') {
				data_biomass_nonbiomass$raw_data_biomass
			} else {
				data_biomass_nonbiomass$raw_data_facet
			}

			if (pred == 'ph') {
				names(data_TEMP_by_site)[names(data_TEMP_by_site) == 'site_ph'] <- 'ph'
				x <- if (facet == 'biomass') {
					data_biomass_nonbiomass$raw_data_biomass[ , .(val = mean(.SD[['site_ph']])), by = SITE][['val']]
				} else {
					data_biomass_nonbiomass$raw_data_facet[ , .(val = mean(.SD[['site_ph']])), by = SITE][['val']]
				}
			} else {
				if (facet == 'biomass') {
					x <- data_biomass_nonbiomass$raw_data_biomass[ , .(val = mean(.SD[[pred]])), by = SITE][['val']]				
				} else {
					x <- data_biomass_nonbiomass$raw_data_facet[ , .(val = mean(.SD[[pred]])), by = SITE][['val']]				
				}
			}
			
			if (facet == 'biomass') {
				observed <- data_biomass_nonbiomass$raw_data_biomass[ , .(val = mean(Biomass)), by = SITE][['val']]
			} else {
				observed <- data_biomass_nonbiomass$raw_data_facet[ , .(val = mean(get(facet_raw))), by = SITE][['val']]			
			}
			
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

			data_TEMP_by_plant <- if (facet == 'biomass') {
				data_biomass_nonbiomass$raw_data_biomass
			} else {
				data_biomass_nonbiomass$raw_data_facet
			}

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

	if (data_biomass_nonbiomass$n_covariates == 2) {
		ncol <- 2
	} else {
		ncol <- ceiling(sqrt(length(data_biomass_nonbiomass$n_covariates)))
	}
	width <- 6 * ncol
	height <- 5 * ceiling(length(data_biomass_nonbiomass$n_covariates) / ncol)

	responses <- plot_grid(plotlist = responses, ncol = ncol)

	filename <- paste0(out_dir, '/response_curves_', facet, '_', resp_type, '.png')
	ggsave(plot = responses, filename = filename, width = width, height = height, dpi = 600)
	invisible(responses)
	
}


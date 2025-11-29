#' Graph response of non-biomass trait to the environment
#'
#' For graphing response of a non-biomass trait against the environment.
#' @param trait "R-friendly" version of name of trait
#' @param quant_threshold Used to define upper limit of the y axis based on range of values of the upper CI limits.
graph_response_curves_nonbiomass_traits_vs_environment <- function(trait, quant_threshold = 0.5) {

	n_covariates <- data_traits$n_covariates_traits
	covariates <- data_traits$covariates_traits
	resp_curves_x <- data_traits$resp_curves_x_traits
	resp_curves_unscaled <- data_traits$resp_curve_x_traits_unscaled
	centers <- data_traits$x_centers_traits
	scales <- data_traits$x_scales_traits
	
	nice_trait_names <- get_nice_trait(trait)
	trait_raw <- get_raw_trait_name(trait)
	y_lab <- paste0('Site-level Mean ', nice_trait_names$long)

	# min/max values for axis limits
	y <- data$y_traits
	y_range <- nice_trait_names$range_fx(y)
	
	responses <- list()
	for (i in 1:n_covariates) {

		pred <- covariates[i]

		nice_pred_names <- get_nice_predictor(pred)
		x_lab <- nice_pred_names$long

		title <- paste0(nice_trait_names$short, ' vs ', nice_pred_names$short)

		# unscaled predictor value
		x <- resp_curves_unscaled[ , pred]
		xlim <- nice_pred_names$range_fx(x)

		# create data frame with SDM response
		param <- paste0('response_curves_traits_mu')
		if (n_covariates == 1) {

			response_mean <- hammer_extract(chains, param = param, j = TRUE, stat = 'mean')
			response_lower <- hammer_extract(chains, param = param, j = TRUE, stat = 'lower')
			response_upper <- hammer_extract(chains, param = param, j = TRUE, stat = 'upper')

		} else {
		
			response_mean <- hammer_extract(chains, param = param, j = TRUE, k = i, stat = 'mean')
			response_lower <- hammer_extract(chains, param = param, j = TRUE, k = i, stat = 'lower')
			response_upper <- hammer_extract(chains, param = param, j = TRUE, k = i, stat = 'upper')

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
			xlab(x_lab) +
			ylab(y_lab) +
			ggtitle(title) +
			coord_cartesian(xlim = c(xlim[1], xlim[2]), ylim = c(y_range[1], y_range[2])) +
			theme(
				legend.position = 'none',
				plot.title = element_text(size = 10),
				axis.title = element_text(size = 8),
				axis.text = element_text(size = 8)
			)

			### add observations
			data_TEMP_by_site <- data_traits$raw_data_traits
			if (pred == 'ph') {
				names(data_TEMP_by_site)[names(data_TEMP_by_site) == 'site_soil_pH'] <- 'ph'
				x <- data_traits$raw_data_traits[ , .(val = mean(.SD[['site_soil_pH']])), by = SITE][['val']]
			} else {
				x <- data_traits$raw_data_traits[ , .(val = mean(.SD[[pred]])), by = SITE][['val']]				
			}

			observed <- rep(NA_real_, data_traits$n_pheno_sites)
			sites <- unique(data_traits$raw_data_traits$SITE)
			for (count in seq_along(sites)) {
			
				site <- sites[count]
				trait_column_name <- get_raw_trait_name(trait)
				observed[count] <- mean(data_traits$raw_data_traits[[trait_column_name]][data_traits$raw_data_traits$SITE == site])

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

			data_TEMP_by_plant <- data_traits$raw_data_traits
			if (pred == 'ph') {
				data_TEMP_by_plant$x <- data_TEMP_by_plant[['site_soil_pH']]
			} else {
				data_TEMP_by_plant$x <- data_TEMP_by_plant[[pred]]
			}
			data_TEMP_by_plant$y <- data_TEMP_by_plant[[trait_raw]]
			data_TEMP_by_plant$site <- data_TEMP_by_plant$SITE

			response <- response + geom_point(
				data = data_TEMP_by_plant,
				mapping = aes(x = x, y = y, fill = site),
				pch = 21, size = 1
			)

		responses[[i]] <- response

	} # next predictor

	if (data_traits$n_covariates == 2) {
		ncol <- 2
	} else {
		ncol <- ceiling(sqrt(length(data_traits$n_covariates)))
	}
	width <- 6 * ncol
	height <- 5 * ceiling(length(data_traits$n_covariates) / ncol)

	responses <- plot_grid(plotlist = responses, ncol = ncol)

	ggsave(plot = responses, filename = paste0(out_dir, '/response_curves_', trait, '.png'), width = width, height = height, dpi = 600)
	invisible(responses)
	
}


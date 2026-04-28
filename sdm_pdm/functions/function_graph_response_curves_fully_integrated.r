#' Graph of response curves for occurrence vs environment
#' 
#' source('C:/Kaji/R/andropogon_integratedEcology/sdm_pdm/functions/function_graph_response_curves_fully_integrated.r')
#' 
#' out_dir			Folder in which to save graphs.
#' chains			MCMC chains
#' zero_inflated	If `TRUE`, model is zero-inflated
graph_response_curves_fully_integrated <- function(
	out_dir,
	chains,
	log_precip_occs,
	resp_distrib_biomass,
	log_precip_biomass,
	transform_biomass,
	nonbiomass_facets
) {

	### prepare metadata on each facet
	##################################

	# to plot each facet, we set the model matrices of the other facets at their centered (mean) values

	# occurrences
	covariates_occs <- data_occs_counties$covariates
	n_covariates_occs <- length(covariates_occs)
	n_terms_occs <- data_occs_counties$n_terms - 1 # do not count intercept
	x_occs_array <- data_occs_counties$resp_curves_x

	# occurrence model matrix
	if (n_terms_occs == 1) {
		x_occs_centered <- x_occs_array
	} else {
		x_occs_centered <- x_occs_array[ , , 1]
	}
	x_occs_centered[ , 2:ncol(x_occs_centered)] <- 0

	# biomass
	covariates_biomass <- data_biomass_counties$covariates
	n_covariates_biomass <- length(covariates_biomass)
	n_terms_biomass <- data_biomass_counties$n_terms - 1 # do not count intercept
	x_biomass_array <- data_biomass_counties$resp_curves_x

	# centered biomass model matrix
	if (n_terms_biomass == 1) {
		x_biomass_centered <- x_biomass_array
	} else {
		x_biomass_centered <- x_biomass_array[ , , 1]
	}
	x_biomass_centered[ , 2:ncol(x_biomass_centered)] <- 0

	# non-biomass facets
	covariates_nonbiomass <- x_nonbiomass_array <- list()
	n_terms_nonbiomass <- n_covariates_nonbiomass <- numeric()
	for (f in seq_along(nonbiomass_facets)) {

		facet <- names(nonbiomass_facets)[f]

		covariates_nonbiomass[[f]] <- attr(terms(nonbiomass_facets[[facet]]$formula), 'term.labels')
		n_covariates_nonbiomass[f] <- length(covariates_nonbiomass[[f]])
		n_terms_nonbiomass[f] <- data_nonbiomass_counties[[facet]]$n_terms - 1 # do not count intercept
		x_nonbiomass_array[[f]] <- data_nonbiomass_counties[[facet]]$resp_curves_x

	}
	names(covariates_nonbiomass) <- names(n_covariates_nonbiomass) <- names(n_terms_nonbiomass) <- names(x_nonbiomass_array) <- names(nonbiomass_facets)

	# centered non-biomass model matrices
	x_nonbiomass_centered <- list()
	for (f in seq_along((nonbiomass_facets))) {

		facet <- names(nonbiomass_facets)[f]
		x <- data_nonbiomass_counties[[facet]]$resp_curves_x
		if (n_terms_nonbiomass[facet] > 1) x <- x[ , , 1]
		x[ , 2:ncol(x)] <- 0
		x_nonbiomass_centered[[f]] <- x

	}
	names(x_nonbiomass_centered) <- names(nonbiomass_facets)

	# zero-inflation
	covariates_psi <- data_psi_counties$covariates
	n_covariates_psi <- length(covariates_psi)
	n_terms_psi <- data_psi_counties$n_terms - 1 # do not count intercept
	x_psi_array <- data_psi_counties$resp_curves_x

	# centered zero-inflation model matrix
	x_psi_centered <- data_psi_counties$resp_curves_x
	n_terms_psi <- data_psi_counties$n_terms
	if (n_terms_psi > 1) x_psi_centered <- x_psi_centered[ , , 1]
	x_psi_centered[ , 2:ncol(x_psi_centered)] <- 0

	### occurrence
	##############

		# fix all other variables at 0 (bc centered)

		say('   occurrence...')
		y_lab <- 'Relative Abundance'

		y_max <- -Inf
		responses <- list()
		
		for (i in seq_along(covariates_occs)) {

			covariate_occs <- covariates_occs[i]
			nice <- get_nice_predictor(covariate_occs)
			nice_title <- nice$short
			nice_axis <- nice$long
			nice_title <- bquote('Abundance' * ' versus ' * .(nice_title))

			# abundance model matrix
			if (n_covariates_occs == 1) {
				x_occs <- x_occs_array
			} else {
				x_occs <- x_occs_array[ , , i]
			}

			# predict
			preds <- predict_fully_integrated(
				chains = chains,
				nonbiomass_facets = nonbiomass_facets,
				resp_distrib_biomass = resp_distrib_biomass,
				transform_biomass = transform_biomass,
				x_occs = x_occs,
				x_psi = x_psi_centered,
				x_biomass = x_biomass_centered,
				x_nonbiomass = x_nonbiomass_centered
			)

			response_mean <- colMeans(preds$preds_occs)
			response_lower <- apply(preds$preds_occs, 2, quantile, 0.05)
			response_upper <- apply(preds$preds_occs, 2, quantile, 0.95)

			# unscaled predictor value
			x_plot <- data_occs_counties$resp_curves_x_unscaled[ , covariate_occs]
			if (log_precip_occs & covariate_occs %in% paste0('bio', c(12:14, 16:19))) x_plot <- 10^x_plot - 1

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

		for (i in 1:n_covariates_occs) {
			responses[[i]] <- responses[[i]] + coord_cartesian(ylim = c(0, y_max)) 
		}

		if (n_covariates_occs == 1) {
			nrow <- 1
			width <- 8
			height <- 6
		} else if (n_covariates_occs <= 2) {
			nrow <- 1
			width <- 10
			height <- 5
		} else if (n_covariates_occs <= 3) {
			nrow <- 1
			width <- 14
			height <- 4
		} else {
			nrow <- 2
			width <- 16
			height <- 10	
		}

		responses_occs <- plot_grid(plotlist = responses, nrow = nrow)
		ggsave(plot = responses_occs, filename = paste0(out_dir, '/response_curves_abundance_mean_inner_quant.png'), width = width, height = height, dpi = 120)

	### biomass
	###########

		# fix all other variables at 0 (bc centered)

		say('   biomass...')
		y_lab <- 'Biomass (g)'

		y_max <- -Inf
		responses <- list()
		
		for (i in seq_along(covariates_biomass)) {

			covariate_biomass <- covariates_biomass[i]
			nice <- get_nice_predictor(covariate_biomass)
			nice_title <- nice$short
			nice_axis <- nice$long
			nice_title <- bquote('Biomass' * ' versus ' * .(nice_title))

			# biomass model matrix
			if (n_covariates_biomass == 1) {
				x_biomass <- x_biomass_array
			} else {
				x_biomass <- x_biomass_array[ , , i]
			}

			# predict
			preds <- predict_fully_integrated(
				chains = chains,
				nonbiomass_facets = nonbiomass_facets,
				resp_distrib_biomass = resp_distrib_biomass,
				transform_biomass = transform_biomass,
				x_occs = x_occs_centered,
				x_psi = x_psi_centered,
				x_biomass = x_biomass,
				x_nonbiomass = x_nonbiomass_centered
			)

			response_mean <- colMeans(preds$preds_biomass)
			response_lower <- apply(preds$preds_biomass, 2, quantile, 0.05)
			response_upper <- apply(preds$preds_biomass, 2, quantile, 0.95)

			# unscaled predictor value
			x_plot <- data_biomass_counties$resp_curves_x_unscaled[ , covariate_biomass]
			if (log_precip_biomass & covariate_biomass %in% paste0('bio', c(12:14, 16:19))) x_plot <- 10^x_plot - 1

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

		for (i in 1:n_covariates_biomass) {
			responses[[i]] <- responses[[i]] + coord_cartesian(ylim = c(0, y_max)) 
		}

		if (n_covariates_biomass == 1) {
			nrow <- 1
			width <- 8
			height <- 6
		} else if (n_covariates_biomass <= 2) {
			nrow <- 1
			width <- 10
			height <- 5
		} else if (n_covariates_biomass <= 3) {
			nrow <- 1
			width <- 14
			height <- 4
		} else {
			nrow <- 2
			width <- 16
			height <- 10	
		}

		responses_biomass <- plot_grid(plotlist = responses, nrow = nrow)
		ggsave(plot = responses_biomass, filename = paste0(out_dir, '/response_curves_biomass_mean_inner_quant.png'), width = width, height = height, dpi = 120)

	### non-biomass
	###############

		for (f in seq_along(nonbiomass_facets)) {

			# fix all other variables at 0 (bc centered)
			facet <- names(nonbiomass_facets)[f]

			nice_facet <- get_nice_trait(facet)
			say('   ', facet, '...')

			y_max <- -Inf
			responses <- list()
			
			for (i in seq_along(covariates_nonbiomass[[facet]])) {

				covariate_nonbiomass <- covariates_nonbiomass[[facet]][i]
				nice <- get_nice_predictor(covariate_nonbiomass)
				nice_title <- capIt(nice$short)
				nice_axis <- nice$long
				nice_title <- bquote(.(nice_facet$short) * ' versus ' * .(nice_title))

				this_x_nonbiomass <- x_nonbiomass_centered
				if (n_covariates_nonbiomass[[facet]] == 1) {
					this_x_nonbiomass[[facet]] <- x_nonbiomass_array[[facet]]
				} else {
					this_x_nonbiomass[[facet]] <- x_nonbiomass_array[[facet]][ , , i]
				}

				# predict
				preds <- predict_fully_integrated(
					chains = chains,
					nonbiomass_facets = nonbiomass_facets,
					resp_distrib_biomass = resp_distrib_biomass,
					transform_biomass = transform_biomass,
					x_occs = x_occs_centered,
					x_psi = x_psi_centered,
					x_biomass = x_biomass_centered,
					x_nonbiomass = this_x_nonbiomass
				)

				response_mean <- colMeans(preds$preds_nonbiomass[[facet]])
				response_lower <- apply(preds$preds_nonbiomass[[facet]], 2, quantile, 0.05)
				response_upper <- apply(preds$preds_nonbiomass[[facet]], 2, quantile, 0.95)

				# unscaled predictor value
				x_plot <- data_nonbiomass_counties[[facet]]$resp_curves_x_unscaled[ , covariate_nonbiomass]

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
					ylab(nice_facet$long) +
					ggtitle(nice_title) +
					theme(
						legend.position = 'none',
						plot.title = element_text(size = 9)
					)

				y_max <- max(y_max, df_ci$response)

			}

			for (i in 1:n_covariates_nonbiomass[[facet]]) {
				responses[[i]] <- responses[[i]] + coord_cartesian(ylim = c(0, y_max)) 
			}

			if (n_covariates_nonbiomass[[facet]] == 1) {
				nrow <- 1
				width <- 8
				height <- 6
			} else if (n_covariates_nonbiomass[[facet]] <= 2) {
				nrow <- 1
				width <- 10
				height <- 5
			} else if (n_covariates_nonbiomass[[facet]] <= 3) {
				nrow <- 1
				width <- 14
				height <- 4
			} else {
				nrow <- 2
				width <- 16
				height <- 10	
			}

			responses_nonbiomass <- plot_grid(plotlist = responses, nrow = nrow)
			ggsave(plot = responses_nonbiomass, filename = paste0(out_dir, '/response_curves_', facet, '_mean_inner_quant.png'), width = width, height = height, dpi = 120)

		} # next non-biomass facet

	### zero-inflation
	##################

		# fix all other variables at 0 (bc centered)

		say('   zero-inflation...')
		y_lab <- 'Probability of Occurrence'

		responses <- list()
		for (i in seq_along(covariates_psi)) {

			covariate_psi <- covariates_psi[i]
			nice <- get_nice_predictor(covariate_psi)
			nice_title <- nice$short
			nice_axis <- nice$long
			nice_title <- bquote('Probability of Occurrence' * ' versus ' * .(nice_title))

			# abundance model matrix
			if (n_covariates_psi == 1) {
				x_psi <- x_psi_array
			} else {
				x_psi <- x_psi_array[ , , i]
			}

			# predict
			preds <- predict_fully_integrated(
				chains = chains,
				nonbiomass_facets = nonbiomass_facets,
				resp_distrib_biomass = resp_distrib_biomass,
				transform_biomass = transform_biomass,
				x_occs = x_occs_centered,
				x_psi = x_psi,
				x_biomass = x_biomass_centered,
				x_nonbiomass = x_nonbiomass_centered
			)

			response_mean <- colMeans(preds$preds_psi)
			response_lower <- apply(preds$preds_psi, 2, quantile, 0.05)
			response_upper <- apply(preds$preds_psi, 2, quantile, 0.95)

			# unscaled predictor value (use `log_precip_occs` for psi!)
			x_plot <- data_psi_counties$resp_curves_x_unscaled[ , covariate_psi]
			if (log_precip_occs & covariate_psi %in% paste0('bio', c(12:14, 16:19))) x_plot <- 10^x_plot - 1

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

		}

		for (i in 1:n_covariates_psi) {
			responses[[i]] <- responses[[i]] + coord_cartesian(ylim = c(0, 1)) 
		}

		if (n_covariates_psi == 1) {
			nrow <- 1
			width <- 8
			width <- 6
		} else if (n_covariates_psi <= 2) {
			nrow <- 1
			width <- 8
			height <- 4
		} else if (n_covariates_psi <= 3) {
			nrow <- 1
			width <- 14
			height <- 4
		} else {
			nrow <- 2
			width <- 16
			height <- 10	
		}

		responses_psi <- plot_grid(plotlist = responses, nrow = nrow)
		ggsave(plot = responses_psi, filename = paste0(out_dir, '/response_curves_psi_mean_inner_quant.png'), width = width, height = height, dpi = 120)

}

#' Post-modeling workflow for any model.
#'
#' @param facet Name of the facet being modeled (e.g., occurrence, biomass, plant height, etc.)
#' @param formulae `List` of model formula.
#' @param descrip Textual description of the model
#' @param homoscedastic `TRUE` if model is homoscedastic.
#' @param out_dir Folder in which to save results.
workflow_postmodeling_generic <- function(facet, formulae, descrip, homoscedastic, out_dir) {

	# ### model fit
	# #####################
	# say('model fit', level = 2)

	# 	### WAIC
	# 	########

	# 	sink(paste0(out_dir, '/waic.txt'), split = TRUE)
	# 	say('WAIC')
	# 	say(date(), post = 2)
	# 	print(chains$WAIC)
	# 	sink()

	### model convergence
	#####################
	say('model convergence', level = 2)

		### density/trace plots
		vars <- monitors_coeffs_not_indexed
		if (length(vars) > 0) {
			for (var in vars) {

				mcmc <- hammer_subset(chains, var)
				mcmc <- mcmc$samples
				mcmc <- ggs(mcmc)

				trace <- ggs_traceplot(mcmc, family = var)
				density <- ggs_density(mcmc, family = var, hpd = TRUE)
				combo <- trace + density
				filename <- paste0(out_dir, '/', var, '_density_trace.png')
				ggsave(combo, file = filename, width = 18, height = 12)
			
			}
		}

		vars <- monitors_coeffs_single_index
		if (length(vars) > 0) {
			for (var in vars) {

				mcmc <- hammer_subset(chains, var, j = TRUE)
				mcmc <- mcmc$samples
				mcmc <- ggs(mcmc)

				trace <- ggs_traceplot(mcmc, family = var)
				density <- ggs_density(mcmc, family = var, hpd = TRUE)
				combo <- trace + density
				filename <- paste0(out_dir, '/', var, '_density_trace.png')
				ggsave(combo, file = filename, width = 18, height = 12)
			
			}
		}

		vars <- monitors_coeffs_double_index
		if (length(vars) > 0) {
			for (var in vars) {

				mcmc <- hammer_subset(chains, var, j = TRUE, k = TRUE)
				mcmc <- mcmc$samples
				mcmc <- ggs(mcmc)

				trace <- ggs_traceplot(mcmc, family = var)
				density <- ggs_density(mcmc, family = var, hpd = TRUE)
				combo <- trace + density
				filename <- paste0(out_dir, '/', var, '_density_trace.png')
				ggsave(combo, file = filename, width = 18, height = 12)
			
			}
		}

		vars <- monitors_derived_not_indexed
		if (length(vars) > 0) {
			for (var in vars) {

				mcmc <- hammer_subset(chains, var)
				mcmc <- mcmc$samples
				mcmc <- ggs(mcmc)

				trace <- ggs_traceplot(mcmc, family = var)
				density <- ggs_density(mcmc, family = var, hpd = TRUE)
				combo <- trace + density
				filename <- paste0(out_dir, '/', var, '_density_trace.png')
				ggsave(combo, file = filename, width = 18, height = 12)
			
			}
		}

		### Gelman-Rubin statistic, effective sample size, and coefficients
		###################################################################

		if (exists('mcmc', inherits = FALSE)) rm(mcmc)
		if (exists('indices', inherits = FALSE)) rm(indices)
		params <- character()
		indices <- list()
		
		vars <- monitors_coeffs_not_indexed
		if (length(vars) > 0) {
			for (var in vars) {

				params <- c(params, var)
				indices[[length(indices) + 1]] <- list()

			}
		}

		vars <- monitors_coeffs_single_index
		if (length(vars) > 0) {
			for (var in vars) {

				params <- c(params, var)
				indices[[length(indices) + 1]] <- list(j = TRUE)

			}
		}

		vars <- monitors_coeffs_double_index
		if (length(vars) > 0) {
			for (var in vars) {

				params <- c(params, var)
				indices[[length(indices) + 1]] <- list(j = TRUE, k = TRUE)

			}
		}


		mcmc <- hammer_subset(chains, param = params, indices = indices)
		rhats <- tryCatch(
			gelman.diag(mcmc$samples, autoburnin = FALSE, multivariate = TRUE),
			error = function(cond) FALSE
		)

		ess <- effectiveSize(mcmc$samples)

		### coefficients
		coeffs_means <- hammer_extract(chains, param = params, indices = indices)
		coeffs_lowers <- hammer_extract(chains, param = params, indices = indices, stat = 'lower')
		coeffs_uppers <- hammer_extract(chains, param = params, indices = indices, stat = 'upper')

		### sampling autocorrelation
		############################
		say('sampling autocorrelation', level = 2)

		vars <- monitors_coeffs_not_indexed
		if (length(vars) > 0) {

			for (var in vars) {

				mcmc <- hammer_subset(chains, var)
				mcmc <- mcmc$samples
				ggs_mcmc <- ggs(mcmc)

				ac <- ggs_autocorrelation(ggs_mcmc, family = var)
				ggsave(ac, file = paste0(out_dir, '/', var, '_autocorrelation.png'), width = 19.2, height = 10.8, dpi = 300)

			}

		}

		vars <- monitors_coeffs_single_index
		if (length(vars) > 0) {

			for (var in vars) {

				mcmc <- hammer_subset(chains, var, j = TRUE)
				mcmc <- mcmc$samples
				ggs_mcmc <- ggs(mcmc)

				ac <- ggs_autocorrelation(ggs_mcmc, family = var)
				ggsave(ac, file = paste0(out_dir, '/', var, '_autocorrelation.png'), width = 19.2, height = 10.8, dpi = 300)

			}

		}

		vars <- monitors_coeffs_double_index
		if (length(vars) > 0) {

			for (var in vars) {

				mcmc <- hammer_subset(chains, var, j = TRUE, k = TRUE)
				mcmc <- mcmc$samples
				ggs_mcmc <- ggs(mcmc)

				ac <- ggs_autocorrelation(ggs_mcmc, family = var)
				ggsave(ac, file = paste0(out_dir, '/', var, '_autocorrelation.png'), width = 19.2, height = 10.8, dpi = 300)

			}

		}

	### parameter estimates
	#######################
	say('parameter estimates', level = 2)

		graphs <- list()
		
		vars <- monitors_coeffs_not_indexed
		if (length(vars) > 0) {
			for (var in vars) {

				mcmc <- hammer_subset(chains, var)
				mcmc <- ggs(mcmc$samples)

				graphs[[length(graphs) + 1]] <- 
					ggs_caterpillar(mcmc, family = var) +
					xlab('Estimated value') +
					ggtitle(var) +
					theme(
						plot.title = element_text(size = 12),
						axis.title.y = element_blank()
					)

			}
		}

		vars <- monitors_coeffs_single_index
		if (length(vars) > 0) {
			for (var in vars) {

				mcmc <- hammer_subset(chains, var, j = TRUE)
				mcmc <- ggs(mcmc$samples)

				graphs[[length(graphs) + 1]] <- 
					ggs_caterpillar(mcmc, family = var) +
					xlab('Estimated value') +
					ggtitle(var) +
					theme(
						plot.title = element_text(size = 12),
						axis.title.y = element_blank()
					)

			}
		}

		vars <- monitors_coeffs_double_index
		if (length(vars) > 0) {
			for (var in vars) {

				mcmc <- hammer_subset(chains, var, j = TRUE, k = TRUE)
				mcmc <- ggs(mcmc$samples)

				graphs[[length(graphs) + 1]] <- 
					ggs_caterpillar(mcmc, family = var) +
					xlab('Estimated value') +
					ggtitle(var) +
					theme(
						plot.title = element_text(size = 12),
						axis.title.y = element_blank()
					)

			}
		}

		vars <- monitors_derived_single_index
		if (length(vars) > 0) {
			for (var in vars) {

				mcmc <- hammer_subset(chains, var, j = TRUE)
				mcmc <- ggs(mcmc$samples)

				graphs[[length(graphs) + 1]] <- 
					ggs_caterpillar(mcmc, family = var) +
					xlab('Estimated value') +
					ggtitle(var) +
					theme(
						plot.title = element_text(size = 12),
						axis.title.y = element_blank()
					)

			}
		}

		vars <- monitors_derived_double_index
		if (length(vars) > 0) {
			for (var in vars) {

				mcmc <- hammer_subset(chains, var, j = TRUE, k = TRUE)
				mcmc <- ggs(mcmc$samples)

				graphs[[length(graphs) + 1]] <- 
					ggs_caterpillar(mcmc, family = var) +
					xlab('Estimated value') +
					ggtitle(var) +
					theme(
						plot.title = element_text(size = 12),
						axis.title.y = element_blank()
					)

			}
		}

		combo <- plot_grid(plotlist = graphs, ncol = 3, align = 'v', axis = 'l')
		ggsave(combo, file = paste0(out_dir, '/coefficients_parameters.png'), width = 19.2, height = 10.8, bg = 'white')

	### metadata
	############

		if (is.logical(rhats)) {
			rhat_means <- NA
			rhat_upper <- NA
		} else {
			rhat_means <- rhats$psrf[ , 1]
			rhat_upper <- rhats$psrf[ , 2]
		}

		coeffs <- data.table(
			coeff = names(coeffs_means),
			lower = coeffs_lowers,
			mean = coeffs_means,
			upper = coeffs_uppers,
			significant = NA,
			eff_sample_size = ess,
			rhat = rhat_means,
			rhat_upper = rhat_upper
			
		)

		insignificant <- coeffs_lowers <= 0 & coeffs_uppers >= 0
		insignificant <- ifelse(insignificant, 'ns', '*')

		coeffs$significant <- insignificant

		meta <- list(
			facet = facet,
			descrip = descrip,
			date = date(),
			formulae = formulae,
			waic = chains$WAIC,
			monitors = c(
				monitors_coeffs_not_indexed = monitors_coeffs_not_indexed,
				monitors_coeffs_single_index = monitors_coeffs_single_index,
				monitors_coeffs_double_index =monitors_coeffs_double_index,
				monitors_derived_not_indexed = monitors_derived_not_indexed,
				monitors_derived_single_index = monitors_derived_single_index,
				monitors_derived_double_index = monitors_derived_double_index,
				# monitors_geog_nam = monitors_geog_nam,
				# monitors_geog_conus = monitors_geog_conus,
				monitors_dharma = monitors_dharma,
				monitors_resp_curves = monitors_resp_curves
			),
			coeffs = coeffs
		)

		saveRDS(meta, paste0(out_dir, '/!meta_generic.rds'))
		sink(paste0(out_dir, '/!meta_generic.txt'), split = TRUE)
			print(meta)
		sink()

}

#' Post-modeling workflow for any model.
#' source('C:/Kaji/R/andropogon_integratedEcology/sdm_pdm/workflows/workflow_postmodeling_generic.r')
#'
#' facet		Name of the facet being modeled (e.g., occurrence, biomass, plant height, etc.)
#' formulae		`List` of model formula.
#' descrip		Textual description of the model
#' homoscedastic `TRUE` if model is homoscedastic.
#' out_dir		Folder in which to save results.
workflow_postmodeling_generic <- function(facet, formulae, descrip, homoscedastic, out_dir) {

	### model convergence
	#####################
	say('model convergence', level = 2)

		### density/trace plots
		vars <- monitors_coeffs_not_indexed
		if (length(vars) > 0) {
			for (var in vars) {

				mcmc <- mc_subset(chains, var)
				mcmc <- mcmc$samples
				mcmc <- ggs(mcmc)

				trace <- ggs_traceplot(mcmc, family = var)
				density <- ggs_density(mcmc, family = var, hpd = TRUE)
				combo <- trace + density
				filename <- paste0(out_dir, '/density_trace_', var, '.png')
				ggsave(combo, file = filename, width = 18, height = 9)
			
			}
		}

		vars <- monitors_coeffs_single_index
		if (length(vars) > 0) {
			for (var in vars) {

				mcmc <- mc_subset(chains, var, j = TRUE)
				mcmc <- mcmc$samples
				mcmc <- ggs(mcmc)

				trace <- ggs_traceplot(mcmc, family = var)
				density <- ggs_density(mcmc, family = var, hpd = TRUE)
				combo <- trace + density
				filename <- paste0(out_dir, '/density_trace_', var, '.png')
				ggsave(combo, file = filename, width = 18, height = 12)
			
			}
		}

		vars <- monitors_coeffs_double_index
		if (length(vars) > 0) {
			for (var in vars) {

				if (var == 'correlation') {

					rows <- matrix(1:n_facets, byrow = FALSE, nrow = n_facets, ncol = n_facets)
					cols <- matrix(1:n_facets, byrow = TRUE, nrow = n_facets, ncol = n_facets)

					upper_rows <- rows[upper.tri(rows, diag = FALSE)]
					upper_cols <- cols[upper.tri(cols, diag = FALSE)]

					params <- paste0('correlation[', upper_rows, ', ', upper_cols, ']')
					mcmc <- mc_subset(chains, params)

				} else {
					mcmc <- mc_subset(chains, var, j = TRUE, k = TRUE)
				}

				mcmc <- mcmc$samples
				mcmc <- ggs(mcmc)

				trace <- ggs_traceplot(mcmc, family = var)
				density <- ggs_density(mcmc, family = var, hpd = TRUE)
				combo <- trace + density
				filename <- paste0(out_dir, '/density_trace_', var, '.png')
				ggsave(combo, file = filename, width = 18, height = 12)
			
			}
		}

		vars <- monitors_derived_not_indexed
		if (length(vars) > 0) {
			for (var in vars) {

				mcmc <- mc_subset(chains, var)
				mcmc <- mcmc$samples
				mcmc <- ggs(mcmc)

				trace <- ggs_traceplot(mcmc, family = var)
				density <- ggs_density(mcmc, family = var, hpd = TRUE)
				combo <- trace + density
				filename <- paste0(out_dir, '/density_trace_', var, '.png')
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

				if (var == 'correlation') {

					extract <- mc_extract(chains, 'correlation', j = TRUE, k = TRUE)
					n <- round(sqrt(length(extract)))

					rows <- matrix(rep(1:n, n), n, n)
					cols <- matrix(rep(1:n, each = n), n, n)
					rows[lower.tri(rows, diag = TRUE)] <- NA
					cols[lower.tri(cols, diag = TRUE)] <- NA

					expand <- matrix(c(rows, cols), ncol = 2)
					expand <- expand[complete.cases(expand), , drop = FALSE]

					this <- character()
					for (i in 1:nrow(expand)) {
						this[i] <- paste0('correlation[', expand[i, 1], ', ', expand[i, 2], ']')
						indices[[length(indices) + 1]] <- list()
					}
					params <- c(params, this)

				} else {
					params <- c(params, var)
					indices[[length(indices) + 1]] <- list(j = TRUE, k = TRUE)
				}

			}
		}


		mcmc <- mc_subset(chains, param = params, indices = indices)
		rhats <- tryCatch(
			gelman.diag(mcmc$samples, autoburnin = FALSE, multivariate = TRUE),
			error = function(cond) FALSE
		)

		ess <- effectiveSize(mcmc$samples)

		### coefficients
		coeffs_means <- mc_extract(chains, param = params, indices = indices)
		coeffs_lowers <- mc_extract(chains, param = params, indices = indices, stat = 'lower')
		coeffs_uppers <- mc_extract(chains, param = params, indices = indices, stat = 'upper')

		### sampling autocorrelation
		############################
		say('sampling autocorrelation', level = 2)

		vars <- monitors_coeffs_not_indexed
		if (length(vars) > 0) {

			for (var in vars) {

				if (var == 'correlation[1, 2]') {
				
					var <- 'correlation'
					ggs_mcmc <- ggs(chains$samples)

				} else {

					mcmc <- mc_subset(chains, var)
					mcmc <- mcmc$samples
					mcmc <- ggs(mcmc)
					mcmc <- mc_subset(chains, var)
					mcmc <- mcmc$samples
					ggs_mcmc <- ggs(mcmc)

				}

				ac <- ggs_autocorrelation(ggs_mcmc, family = var)
				ggsave(ac, file = paste0(out_dir, '/autocorrelation_', var, '.png'), width = 19.2, height = 10.8, dpi = 300)

			}

		}

		vars <- monitors_coeffs_single_index
		if (length(vars) > 0) {

			for (var in vars) {

				mcmc <- mc_subset(chains, var, j = TRUE)
				mcmc <- mcmc$samples
				ggs_mcmc <- ggs(mcmc)

				ac <- ggs_autocorrelation(ggs_mcmc, family = var)
				ggsave(ac, file = paste0(out_dir, '/autocorrelation_', var, '.png'), width = 19.2, height = 10.8, dpi = 300)

			}

		}

		vars <- monitors_coeffs_double_index
		if (length(vars) > 0) {

			for (var in vars) {

				mcmc <- mc_subset(chains, var, j = TRUE, k = TRUE)
				mcmc <- mcmc$samples
				ggs_mcmc <- ggs(mcmc)

				ac <- ggs_autocorrelation(ggs_mcmc, family = var)
				ggsave(ac, file = paste0(out_dir, '/autocorrelation_', var, '.png'), width = 19.2, height = 10.8, dpi = 300)

			}

		}

		### correlation between parameters
		##################################
		say('correlation between parameters', level = 2)

			if (exists('param_stack', inherits = FALSE)) rm(param_stack)
			n_chains <- mc_n_chains(chains)
			vars <- monitors_coeffs_not_indexed
			if (length(vars) > 0) {
				for (var in vars) {

					this_param_stack <- mc_subset(chains, var)
					this_param_stack <- mc_stack(this_param_stack)
					if (exists('param_stack')) {
						param_stack <- cbind(param_stack, this_param_stack)
					} else {
						param_stack <- this_param_stack
					}

				}
			}

			vars <- monitors_coeffs_single_index
			if (length(vars) > 0) {
				for (var in vars) {

					this_param_stack <- mc_subset(chains, var, j = TRUE)
					this_param_stack <- mc_stack(this_param_stack)
					if (exists('param_stack')) {
						param_stack <- cbind(param_stack, this_param_stack)
					} else {
						param_stack <- this_param_stack
					}
				
				}
			}

			vars <- monitors_coeffs_double_index
			if (length(vars) > 0) {
				for (var in vars) {

					this_param_stack <- mc_subset(chains, var, j = TRUE, k = TRUE)
					this_param_stack <- mc_stack(this_param_stack)
					if (exists('param_stack')) {
						param_stack <- cbind(param_stack, this_param_stack)
					} else {
						param_stack <- this_param_stack
					}
				
				}
			}

			cors <- cor(param_stack, method = 'pearson', use = 'complete.obs')
			cors_df <- as.data.frame(as.table(cors))
			names(cors_df) <- c('Var1', 'Var2', 'Correlation')
			cors_df$high_cor <- abs(cors_df$Correlation) > 0.5 & cors_df$Var1 != cors_df$Var2

			cor_plot <- ggplot(cors_df, aes(x = Var1, y = Var2, fill = Correlation)) +
				geom_tile() +
				geom_point(data = subset(cors_df, high_cor), aes(x = Var1, y = Var2), 
					color = 'black', shape = 1, size = 3, stroke = 1.5) +
				scale_fill_gradient2(low = 'red', mid = 'white', high = 'blue', 
					midpoint = 0, limits = c(-1, 1)) +
				theme_minimal() +
				theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
				labs(x = '', y = '', title = 'Parameter Correlations')

			ggsave(cor_plot, file = paste0(out_dir, '/parameter_correlations.png'), width = 11, height = 10, dpi = 300, bg = 'white')


	### parameter estimates
	#######################
	say('parameter estimates', level = 2)

		graphs <- list()
		
		vars <- monitors_coeffs_not_indexed
		if (length(vars) > 0) {

			for (var in vars) {

				if (var == 'correlation[1, 2]') {
				
					var <- 'correlation'
					ggs_mcmc <- ggs(chains$samples)

				} else {

					mcmc <- mc_subset(chains, var)
					mcmc <- mcmc$samples
					ggs_mcmc <- ggs(mcmc)

				}

				graphs[[length(graphs) + 1]] <- 
					ggs_caterpillar(ggs_mcmc, family = var) +
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

				mcmc <- mc_subset(chains, var, j = TRUE)
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

				mcmc <- mc_subset(chains, var, j = TRUE, k = TRUE)
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

				mcmc <- mc_subset(chains, var, j = TRUE)
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

				mcmc <- mc_subset(chains, var, j = TRUE, k = TRUE)
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
		ggsave(combo, file = paste0(out_dir, '/coefficient_estimates.png'), width = 19.2, height = 10.8, bg = 'white')

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
			niter = niter,
			nburnin = nburnin,
			thin = thin,
			n_samples_per_chain = (niter - nburnin) / thin,
			nchains = nchains,
			date = date(),
			formulae = formulae,
			waic = chains$WAIC,
			monitors = c(
				monitors_coeffs_not_indexed = monitors_coeffs_not_indexed,
				monitors_coeffs_single_index = monitors_coeffs_single_index,
				monitors_coeffs_double_index =monitors_coeffs_double_index,
				monitors_derived_not_indexed = monitors_derived_not_indexed,
				monitors_derived_single_index = monitors_derived_single_index,
				monitors_derived_double_index = monitors_derived_double_index
			),
			coeffs = coeffs
		)

		if (exists('monitors_resp_curves', inherits = TRUE)) meta$monitors$monitors_resp_curves <- monitors_resp_curves
		if (any(monitors == 'log_lik')) {

			meta$log_lik = c(
				lower = unname(mc_extract(chains, 'log_lik', stat = 'lower')),
				ll_mean = unname(mc_extract(chains, 'log_lik')),
				upper = unname(mc_extract(chains, 'log_lik', stat = 'upper'))
			)

		}

		saveRDS(meta, paste0(out_dir, '/!meta_generic.rds'))
		sink(paste0(out_dir, '/!meta_generic.txt'), split = TRUE)
			print(meta)
		sink()

}

#' Post-modeling workflow for biomass-only models and models with a biomass component where it does not depend on other traits or abundance. Assumes that the number of predictor variables is one (but the model can have higher-order terms).
#'
#' source('C:/Kaji/R/andropogon_integratedEcology/sdm_pdm_for_sem/workflows/workflow_postmodeling_biomass.r')
#'
#' @param chains From NIMBLE.
#' @param descrip Model description.
#' @param formula_biomass_mu Biomass-environment formula for site-level mean.
#' @param formula_biomass_sigma Biomass-environment formula for site-level sd.
#' @param formula_pzero Biomass-environment probability of 0 biomass.
#' @param zero_inflated If `TRUE`, then model is zero-inflated.
#' @param crossvalidate If `TRUE`, do cross-validations.
#' @param out_dir Folder in which to save results.
workflow_postmodeling_biomass <- function(chains, descrip, formula_biomass_mu, formula_pzero, crossvalidate, out_dir) {

	zero_inflated <- !is.null(formula_pzero)

	### BIOMASS: plant-level residuals analysis
	###########################################
	say('BIOMASS: plant-level residuals analysis', level = 2)

	# NB this uses the mean predicted value of a site as a plant-level prediction
	sims_by_plant <- hammer_subset(chains, param = 'y_sim', j = TRUE)
	sims_by_plant <- hammer_rbind(sims_by_plant)
	sims_by_plant <- t(sims_by_plant)

	estimates <- hammer_subset(chains, param = 'mu', j = TRUE)
	estimates <- hammer_rbind(estimates)

	site_counts <- data_biomass_mu$raw_data_biomass[ , .N, by = SITE]
	fit_by_site <- apply(estimates, 2, median)
	fits <- numeric()
	for (i in seq_along(fit_by_site)) {
		fits <- c(fits, rep(fit_by_site[i], site_counts$N[i]))
	}

	dharma <- createDHARMa(simulatedResponse = sims_by_plant, observedResponse = data_biomass_mu$y_biomass, fittedPredictedResponse = fits, integerResponse = FALSE)

	dharma_quant_test <- testQuantiles(dharma, plot = FALSE)
	dharma_resid_test <- testResiduals(dharma, plot = FALSE) # uniformity, dispersion, outlier
	dev.off()

	file <- paste0(out_dir, '/dharma_by_plant_biomass.png')
	png(file, width = 1200, height = 800)
		plot(dharma)
	dev.off()
	
	file <- paste0(out_dir, '/dharma_by_plant_biomass_residuals.png')
	png(file, width = 1400, height = 1000)
		hist(dharma$scaledResiduals, main = 'DHARMa residuals for biomass by plant', xlab = 'Scaled residuals', breaks = 30)
	dev.off()

	### BIOMASS: site-level residuals analysis
	##########################################

		mean_est <- hammer_extract(chains, 'mu', j = TRUE)
		lower_est <- hammer_extract(chains, 'mu', j = TRUE, stat = 'lower')
		upper_est <- hammer_extract(chains, 'mu', j = TRUE, stat = 'upper')

		obs_est_biomass <- as.data.table(data_biomass_mu$site_vect_biomass)
		obs_est_biomass$mean_est <- mean_est
		obs_est_biomass$lower_est <- upper_est
		obs_est_biomass$upper_est <- lower_est

		vals <- c(obs_est_biomass$biomass_mean, obs_est_biomass$upper_est, obs_est_biomass$lower_est)
		lim <- c(0, 1.05 * max(vals))

		obs_vs_est <- ggplot() +
			geom_abline(intercept = 0, slope = 1) +
			geom_point(obs_est_biomass, mapping = aes(x = biomass_mean, y = mean_est), size = 3, pch = 1) +
			geom_errorbar(obs_est_biomass, mapping = aes(x = biomass_mean, ymin = lower_est, ymax = upper_est), width = 0) +
			coord_cartesian(xlim = lim, ylim = lim) +
			xlab('Observed mean site-level biomass (g)') +
			ylab('Estimated mean site-level biomass (g)') +
			ggtitle('Observed vs. predicted biomass') +
			theme(
				axis.title = element_text(size = 16),
				axis.text = element_text(size = 14)
			)

		ggsave(obs_vs_est, filename = paste0(out_dir, '/observed_vs_predicted_biomass_site_mu.png'), width = 9, height = 9, dpi = 600)

	### BIOMASS: map of residuals
	#############################
	say('BIOMASS: map of residuals', level = 2)

		# Using the mean DHARMa residual of a site for site-level residual

		residuals_by_plant <- dharma$scaledResiduals
		residuals_by_site <- rep(NA, constants$n_pheno_sites)
		site_ids <- unique(data_biomass_mu$raw_data_biomass$SITE)
		for (i in 1:constants$n_pheno_sites) residuals_by_site[i] <- mean(residuals_by_plant[data_biomass_mu$raw_data_biomass$SITE == site_ids[i]])

		site_vect_biomass_resid <- data_biomass_mu$site_vect_biomass
		site_vect_biomass_resid$residual <- residuals_by_site

		nam <- vect('./data_from_gadm/gadm_4pt1_level_1_north_america_sans_alaska_lambert.gpkg')

		# extent
		extent <- buffer(site_vect_biomass_resid, width = 200 * 1000)
		extent <- ext(extent)
		extent <- as.vector(extent)

		map <- ggplot() +
			layer_spatial(nam, color = 'gray30', fill = 'white', linewidth = 0.3) +
			layer_spatial(site_vect_biomass_resid, aes(fill = residual), pch = 21, size = 6) +
			scale_fill_gradient2(
				name = 'Mean of\nDHARMa\nresiduals',
				low = 'red',
				mid = 'beige',
				high = 'blue',
				midpoint = 0.5
			) +
			xlim(extent[1], extent[2]) + ylim(extent[3], extent[4]) +
			ggtitle(
				bquote('Residuals for ' * italic('Andropogon gerardi') * ' biomass '),
				subtitle = '1991-2020 | biomass-only model') +
			theme(
				plot.title = element_text(size = 16),
				plot.subtitle = element_text(size = 14)
			)

		ggsave(plot = map, filename = paste0(out_dir, '/map_residuals_biomass_site_mu.png'), width = 12, height = 9, dpi = 300)

	### BIOMASS: SAC in residuals
	#############################
	say('BIOMASS: SAC in residuals', level = 2)
	
		coords <- as.data.frame(crds(project(site_vect_biomass_resid, enmSdmX::getCRS('WGS84'))))

		# Compute Moran's I
		moran <- moran.test(residuals_by_site, nb2listw(knn2nb(knearneigh(coords, longlat = TRUE, k = 4))))

	### BIOMASS: plant-level residuals vs covariates
	################################################
	say('BIOMASS: plant-level residuals vs covariates', level = 2)

		resids <- dharma$scaledResiduals
		resids_vs_covariates <- list()
		resids_vs_covariates_results <- data.table()

		for (i in seq_len(data_biomass_mu$n_covariates)) {

			if (data_biomass_mu$covariates_biomass[i] == 'ph') {
				x <- data_biomass_mu$raw_data_biomass[['site_soil_pH']]
			} else {
				x <- data_biomass_mu$raw_data_biomass[[data_biomass_mu$covariates_biomass[i]]]
			}

			this_x <- data.frame(
				x = x,
				y = resids
			)

			resid_model <- mgcv::gam(logitAdj(y, epsilon = 0.0001) ~ s(x), data = this_x)

			F <- summary(resid_model)$s.table[1, 3]
			p <- summary(resid_model)$s.table[1, 4]
			sig <- ifelse(p < 0.05, '*', '-')
			
			resids_vs_covariates_results <- rbind(
				resids_vs_covariates_results,
				data.table(
					covariate = data_biomass_mu$covariates_biomass[i],
					F = F,
					p = p,
					significance = sig
				)
			)

			resids_vs_covariates[[i]] <- ggplot(this_x, aes(x = x, y = y)) +
				geom_point() +
				geom_smooth(method = loess, formula = y ~ x, se = FALSE) +
				xlab(data_biomass_mu$terms_biomass[i]) +
				ylab('Residual value') +
				ggtitle('DHARMa Residuals for Biomass Means')

		}

		ncol <- ceiling(sqrt(length(resids_vs_covariates)))
		width <- 6 * ncol
		height <- 5 * ceiling(length(resids_vs_covariates) / ncol)

		resids_vs_covariates <- plot_grid(plotlist = resids_vs_covariates, ncol = ncol)
		ggsave(plot = resids_vs_covariates, filename = paste0(out_dir, '/dharma_residuals_vs_covariates_biomass.png'), width = width, height = height, dpi = 600)

		### BIOMASS: pseudo-R2 (correlation between observed and estimate)
		##################################################################
		say('BIOMASS: model fit', level = 2)

		observed <- data_biomass_mu$raw_data_biomass[ , .(mean_biomass = mean(Biomass)), by = SITE][['mean_biomass']]
		estimated <- hammer_extract(chains, 'mu', j = TRUE)
		correl_mu <- cor(observed, estimated)

		rmse_mu <- unname(rmse_fx(observed, estimated))
		mean_error_mu <- mean((estimated - observed) / observed)
		mean_abs_error_mu <- mean(abs(estimated - observed) / observed)

	### BIOMASS: response curves
	############################
	say('BIOMASS: response curves', level = 2)
	
		responses <- graph_response_curves_biomass_vs_environment(
			out_dir = out_dir,
			chains = chains,
			resp_type = 'mu',
			data_biomass = data_biomass_mu,
			quant_threshold = 0.5
		)

		if (zero_inflated) {

			responses_pzero <- graph_response_curves_biomass_vs_environment(
				out_dir = out_dir,
				resp_type = 'pzero',
				chains = chains,
				data_biomass_mu = data_biomass_mu,
				quant_threshold = 1
			)

		}

	### BIOMASS: maps
	#################
	say('BIOMASS: burn prediction vectors', level = 2)

		pred_vect_nam <- burn_biomass_into_vector(chains = chains, formula_biomass_mu = formula_biomass_mu, formula_pzero = formula_pzero)

		writeVector(pred_vect_nam, paste0(out_dir, '/prediction_vector_nam.gpkg'), overwrite = TRUE)

	### BIOMASS: current map
	########################
	say('BIOMASS: current map', level = 2)

		mean_form <- paste(as.character(formula_biomass_mu), collapse = ' ')
		mean_form <- gsub(mean_form, pattern = 'I\\(', replacement = '')
		mean_form <- gsub(mean_form, pattern = '\\^2\\)', replacement = '²')
		mean_form <- gsub(mean_form, pattern = '*)', replacement = '×')

		subtitle <- paste0('Biomass μ ', mean_form, ' | 1991-2020')

		map <- map_biomass(
			out_dir = out_dir,
			filename_append = 'present',
			pred_vect_nam = pred_vect_nam,
			response_var = 'mu_biomass_county_mean_sq',
			response_var_type = 'mu',
			data_occs = data_occs,
			data_biomass = data_biomass_mu,
			title = bquote('Mean ' * italic('Andropogon gerardi') * ' biomass '),
			subtitle = subtitle,
			legend_title = 'Biomass (g)',
			plot_range_core = TRUE,
			ag_core_quant = 0.95
		)

		saveRDS(map, paste0(out_dir, '/map_biomass_mu_present.rds'))

		if (!is.null(formula_pzero)) {
			
			facet <- 'biomass'
			map <- map_pzero(
				out_dir = out_dir,
				filename_append = 'present',
				pred_vect_nam = pred_vect_nam,
				facet = facet,
				response_var = paste0('pzero_', facet, '_county_sq'),
				data_traits = data_biomass_mu,
				title = bquote('Probability of zero biomass '),
				subtitle = subtitle,
				plot_range_core = TRUE,
				ag_core_quant = 0.95
			)

		}

	# ### BIOMASS: future maps
	# ########################
	# say('BIOMASS: future maps', level = 2)

	# 	for (fut in futs) {

	# 		say(fut)

	# 		response_var <- paste0('mu_biomass_county_mean_', fut)

	# 		mean_form <- paste(as.character(formula_biomass_mu), collapse = ' ')
	# 		subtitle <- paste0('Biomass μ ', mean_form, '  | SSP ', substr(fut, 4, 6), ' ', substr(fut, 8, 11), '-', substr(fut, 13, 16))

	# 		map_biomass_fut <- map_biomass(
	# 			out_dir = out_dir,
	# 			filename_append = fut,
	# 			pred_vect_nam = pred_vect_nam,
	# 			response_var = response_var,
	# 			response_var_type = 'mu',
	# 			data_occs = data_occs,
	# 			data_biomass = data_biomass_mu,
	# 			title = bquote('Future distribution of mean ' * italic('Andropogon gerardi') * ' biomass '),
	# 			subtitle = subtitle,
	# 			legend_title = 'Biomass (g)',
	# 			plot_range_core = TRUE,
	# 			ag_core_quant = 0.95
	# 		)

	# 		if (!is.null(formula_pzero)) {

	# 			face <- 'biomass'				
	# 			map_fut_pzero <- map_pzero(
	# 				out_dir = out_dir,
	# 				filename_append = fut,
	# 				pred_vect_nam = pred_vect_nam,
	# 				facet = facet,
	# 				response_var = paste0('pzero_', facet, '_county_', fut),
	# 				data_traits = data_biomass_mu,
	# 				title = bquote('Future probability of zero biomass '),
	# 				subtitle = subtitle,
	# 				plot_range_core = TRUE,
	# 				ag_core_quant = 0.95
	# 			)

	# 		}

	# 	}

	# ### BIOMASS: future change maps
	# ###############################
	# say('BIOMASS: future change maps', level = 2)

	# 	maps_biomass_change <- list()
	# 	for (fut in futs) {

	# 		say(fut)

	# 		response_var <- paste0('mu_biomass_county_mean_', fut)
	# 		title <- bquote('Change in ' * italic('Andropogon gerardi') * ' biomass ')

	# 		mean_form <- paste(as.character(formula_biomass_mu), collapse = ' ')
	# 		subtitle <- paste0('Biomass μ ', mean_form, '  | SSP ', substr(fut, 4, 6), ' ', substr(fut, 8, 11), '-', substr(fut, 13, 16))

	# 		maps_biomass_change[[length(maps_biomass_change) + 1]] <- map_biomass_change(
	# 			out_dir = out_dir,
	# 			filename_append = fut,
	# 			fut = fut,
	# 			pred_vect_nam = pred_vect_nam,
	# 			response_var = response_var,
	# 			data_occs = data_occs,
	# 			data_biomass = data_biomass_mu,
	# 			title = title,
	# 			subtitle = subtitle,
	# 			plot_range_core = TRUE,
	# 			ag_core_quant = ag_core_quant
	# 		)

	# 		if (!is.null(formula_pzero)) {

	# 			title <- bquote('Change in probability of zero biomass ')

	# 			facet <- 'biomass'				
	# 			map_pzero_change <- map_pzero_change(
	# 				out_dir = out_dir,
	# 				filename_append = fut,
	# 				fut = fut,
	# 				pred_vect_nam = pred_vect_nam,
	# 				facet = facet,
	# 				response_var = paste0('pzero_', facet, '_county_', fut),
	# 				response_var_sq = paste0('pzero_', facet, '_county_sq'),
	# 				data_traits = data_biomass_mu,
	# 				title = title,
	# 				subtitle = subtitle,
	# 				plot_range_core = TRUE,
	# 				ag_core_quant = 0.95
	# 			)

	# 		}

	# 	}


	### occurrence: crossvalidation
	###############################

	if (crossvalidate) {
		
		say('BIOMASS: cross-validation', level = 2)

		max_k_folds <- if (trial) { 1 } else { k_folds }
		cv <- data.table()
		for (k in 1:max_k_folds) {

			say('fold: ', k, level = 3)
		
			fold_constants <- constants
			fold_data <- data
			fold_inits <- inits

			indices <- get(paste0('fold_', k, '_biomass'))

			# data
			fold_data$y_biomass <- fold_data$y_biomass[-indices]

			# constants: order matters in assigning these!
			test_sites <- sort(unique(fold_constants$site_index_biomass[indices]))
			train_sites <- sort(unique(fold_constants$site_index_biomass[-indices]))
			
			fold_constants$site_index_biomass <- renumSeq(fold_constants$site_index_biomass[fold_constants$site_index_biomass %in% train_sites])
			fold_constants$n_biomass <- fold_constants$n_biomass - length(indices)
			fold_constants$x_by_site_biomass <- fold_constants$x_by_site_biomass[train_sites, ]
			fold_constants$n_pheno_sites <- length(train_sites)
			fold_constants$n_counties <- 1 # faster!
			fold_constants$n_response_curve_values <- 1 # faster!
		
			# inits
			fold_inits$site_biomass_mu_log <- fold_inits$site_biomass_mu_log[train_sites]
			fold_inits$y_sim <- fold_inits$y_sim[-indices]

			fold_model <- nimbleModel(
				code = model_code,
				constants = fold_constants,
				data = fold_data,
				inits = fold_inits,
				check = FALSE,
				calculate = FALSE,
				buildDerivs = TRUE
			)

			fold_monitors <- c('beta_mu', 'mu', 'sigma')

			fold_conf <- configureMCMC(
				fold_model,
				monitors = fold_monitors,
				print = FALSE,
				enableWAIC = FALSE
			)

			# add no U-turn sampler (Hamiltonian Monte Carlo)
			# vars <- c('sigma', 'beta_mu', 'beta_pzero')
			vars <- c('sigma_log', monitors_coeffs_not_indexed, monitors_coeffs_single_index, monitors_coeffs_double_index)
			fold_conf$addSampler(target = vars, type = 'NUTS')
			say('NUTS sampler added to all continuous priors.')

			fold_build <- buildMCMC(fold_conf)
			fold_compiled <- compileNimble(fold_model, fold_build, showCompilerOutput = FALSE)

			fold_chains <- runMCMC(
				fold_compiled$fold_build,
				niter = 2 * niter,
				nburnin = 2 * nburnin,
				thin = 2 * thin,
				nchains = 1,
				inits = fold_inits,
				progressBar = TRUE,
				samplesAsCodaMCMC = TRUE,
				summary = TRUE,
				WAIC = FALSE,
				perChainWAIC = FALSE
			)

			test_x <- constants$x_by_site_biomass[test_sites, ]
			preds <- predict_biomass(chains = chains, x = test_x, zero_inflated = zero_inflated, type = 'mu')

			### evaluate predictions
			# observed site means
			site_means <- c()
			for (i in seq_along(test_sites)) {
				
				test_site <- test_sites[i]
				site_means <- c(site_means, mean(data$y[constants$site_index_biomass == test_site]))

			}

			# RMSE, mean abs prediction error, mean absolute error, correlation between observed and predicted
			rmse <- mape <- mae <- correl <- rep(NA_real_, nrow(preds))
			for (iter in 1:nrow(preds)) {
				correl[iter] <- cor(preds[iter, , drop = TRUE], site_means)
				mae[iter] <- mae_fx(preds[iter, , drop = TRUE], site_means)
				mape[iter] <- mean_abs_percent_error_fx(preds[iter, , drop = TRUE], site_means)
				rmse[iter] <- rmse_fx(preds[iter, , drop = TRUE], site_means)
			}

			# quantile of predictions in which observed values fall
			obs_quants <- rep(NA_real_, length(test_sites))
			n <- nrow(preds)
			for (i in seq_along(test_sites)) {
				obs_quants[i] <- sum(preds[ , i] < site_means[i]) / n
			}

			# proportion of observations within the inner 90th quantile of the distribution of predicted values
			obs_quants_prop_in_inner_90 <- sum(obs_quants >= 0.05 & obs_quants <= 0.95) / length(obs_quants)

			cv <- rbind(
				cv,
				data.table(
					k = k,

					n_train_sites = fold_constants$n_pheno_sites,
					n_test_sites = constants$n_pheno_sites - fold_constants$n_pheno_sites,

					n_train_plants = fold_constants$n_biomass,
					n_test_plants = constants$n_biomass - fold_constants$n_biomass,
					
					correl_lower = quantile(correl, 0.025),
					correl_mean = mean(correl),
					correl_upper = quantile(correl, 0.975),

					mae_lower = quantile(mae, 0.025),
					mae_mean = mean(mae),
					mae_upper = quantile(mae, 0.975),

					mape_lower = quantile(mape, 0.025),
					mape_mean = mean(mape),
					mape_upper = quantile(mape, 0.975),

					rmse_lower = quantile(rmse, 0.025),
					rmse_mean = mean(rmse),
					rmse_upper = quantile(rmse, 0.975),

					obs_quants_min = min(obs_quants),
					obs_quants_mean = mean(obs_quants),
					obs_quants_max = max(obs_quants),
					obs_quants_prop_in_inner_90 = obs_quants_prop_in_inner_90

				)
			)

		}

		cv <- rbind(
			cv,
			data.table(
					k = 'means',

					n_train_sites = mean(cv$n_train_sites),
					n_test_sites = mean(cv$n_test_sites),

					n_train_plants = mean(cv$n_train_plants),
					n_test_plants = mean(cv$n_test_plants),
					
					correl_lower = mean(cv$correl_lower),
					correl_mean = mean(cv$correl_mean),
					correl_upper = mean(cv$correl_upper),

					mae_lower = mean(cv$mae_lower),
					mae_mean = mean(cv$mae_mean),
					mae_upper = mean(cv$mae_upper),

					mape_lower = mean(cv$mape_lower),
					mape_mean = mean(cv$mape_mean),
					mape_upper = mean(cv$mape_upper),

					rmse_lower = mean(cv$rmse_lower),
					rmse_mean = mean(cv$rmse_mean),
					rmse_upper = mean(cv$rmse_upper),

					obs_quants_min = mean(cv$obs_quants_min),
					obs_quants_mean = mean(cv$obs_quants_mean),
					obs_quants_max = mean(cv$obs_quants_max),
					obs_quants_prop_in_inner_90 = mean(cv$obs_quants_prop_in_inner_90)

			)
		)

	} # if crossvalidating

	### remember
	############

		resid_p_values <- c(
			moran$estimate[1],
			dharma_resid_test$uniformity$p.value,
			dharma_resid_test$dispersion$p.value,
			dharma_resid_test$outliers$p.value,
			dharma_quant_test$p.value,
			summary(dharma_quant_test$qgamFits[[3]])$s.table[1, 4],
			summary(dharma_quant_test$qgamFits[[2]])$s.table[1, 4],
			summary(dharma_quant_test$qgamFits[[1]])$s.table[1, 4]
		)
		resid_p_values_sig <- ifelse(resid_p_values < 0.05, '*', '-')
		resid_test_statistic_values <- c(
			moran$statistic,
			dharma_resid_test$uniformity$statistic,
			dharma_resid_test$dispersion$statistic,
			dharma_resid_test$outliers$statistic,
			NA,
			summary(dharma_quant_test$qgamFits[[3]])$s.table[1, 3],
			summary(dharma_quant_test$qgamFits[[2]])$s.table[1, 3],
			summary(dharma_quant_test$qgamFits[[1]])$s.table[1, 3]
		)

		dharma_resids <- data.table(
			test = c('spatial autocorrelation', 'uniformity', 'dispersion', 'outliers', 'quantiles, overall', 'quantiles, upper', 'quantiles, middle', 'quantiles, lower'),
			p_value = resid_p_values,
			significant = resid_p_values_sig,
			test_statistic = c('Moran\'s I', names(dharma_resid_test$uniformity$statistic), names(dharma_resid_test$dispersion$statistic), 'exact binomial', NA, rep('chi squared', 3)),
			test_statistic_value = resid_test_statistic_values
		)

		meta_biomass <- list(
			facet = 'biomass',
			descrip = descrip,
			date = date(),
			formulae = formulae,
			dharma_resids = dharma_resids,
			dharma_resids_vs_covariates = resids_vs_covariates_results,
			fit_mu = c(correl_mu = correl_mu, rmse_mu = rmse_mu, mean_error_mu = mean_error_mu, mean_abs_error_mu = mean_abs_error_mu)
		)

		meta_biomass$crossvalidation <- if (crossvalidate) {
			cv
		} else {
			NA
		}

		saveRDS(meta_biomass, paste0(out_dir, '/!meta_biomass.rds'))
		sink(paste0(out_dir, '/!meta_biomass.txt'), split = TRUE)
			print(meta_biomass)
		sink()

}

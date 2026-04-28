#' Post-modeling workflow for biomass-only models and models with a biomass component where it does not depend on other traits or abundance.
#'
#' source('C:/Kaji/R/andropogon_integratedEcology/sdm_pdm/workflows/workflow_postmodeling_biomass.r')
#'
#' chains From NIMBLE.
#' descrip 					Model description.
#' formula_biomass 			Biomass-environment formula for site-level mean.
#' formula_biomass_sigma 	Biomass-environment formula for site-level sd.
#' formula_psi 				Biomass-environment probability of 0 biomass.
#' resp_distrib		 		Response distribution: 'gamma' or 'lognormal'
#' transform 				Either 'softplus' or 'exponential' or 'identity'.
#' crossvalidate 			If `TRUE`, do cross-validations.
#' out_dir 					Folder in which to save results.
workflow_postmodeling_biomass <- function(chains, descrip, formula_biomass, formula_biomass_sigma, formula_psi, resp_distrib, transform, crossvalidate, out_dir) {

	data_biomass <- prepare_biomass_data(formula_biomass = formula_biomass, n_response_curve_values = n_response_curve_values, calib = calib)

	### BIOMASS: plant-level residuals analysis
	###########################################
	say('BIOMASS: plant-level residuals analysis', level = 2)

		# NB this uses the mean predicted value of a site as a plant-level prediction

		# predicted values
		sims_by_plant <- mc_subset(chains, param = 'y_biomass_sim', j = TRUE)
		sims_by_plant <- mc_rbind(sims_by_plant)
		sims_by_plant <- t(sims_by_plant)

		preds <- predict_biomass(chains = chains, x = data_biomass$x_by_site, resp_distrib = resp_distrib, transform = transform)
		fit_by_site <- colMeans(preds)

		# estimates <- mc_subset(chains, param = 'biomass_site', j = TRUE)
		# estimates <- mc_rbind(estimates)

		site_counts <- data_biomass$raw_data_biomass[ , .N, by = SITE]
		# fit_by_site <- apply(estimates, 2, median)
		fits <- numeric()
		for (i in seq_along(fit_by_site)) {
			fits <- c(fits, rep(fit_by_site[i], site_counts$N[i]))
		}

		# dharma <- createDHARMa(simulatedResponse = sims_by_plant, observedResponse = data_biomass$y_biomass, fittedPredictedResponse = fits, integerResponse = FALSE)
		dharma <- createDHARMa(simulatedResponse = sims_by_plant, observedResponse = data_biomass$y_biomass)

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

	### BIOMASS: observed-versus-predicted
	######################################

		preds <- predict_biomass(chains = chains, x = data_biomass$x_by_site, resp_distrib = resp_distrib, transform = transform)
		mean_est <- colMeans(preds)
		lower_est <- apply(preds, 2, quantile, 0.025)
		upper_est <- apply(preds, 2, quantile, 0.975)

		# mean_est <- mc_extract(chains, 'biomass_site', j = TRUE)
		# lower_est <- mc_extract(chains, 'biomass_site', j = TRUE, stat = 'lower')
		# upper_est <- mc_extract(chains, 'biomass_site', j = TRUE, stat = 'upper')

		obs_est_biomass <- as.data.table(data_biomass$site_vect)
		obs_est_biomass$mean_est <- mean_est
		obs_est_biomass$lower_est <- lower_est
		obs_est_biomass$upper_est <- upper_est

		max_val <- max(obs_est_biomass$upper_est)
		lim <- c(0, 1.05 * max_val)

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

		ggsave(obs_vs_est, filename = paste0(out_dir, '/observed_vs_predicted_biomass_site_mu.png'), width = 9, height = 9, dpi = 200)

	### BIOMASS: map of residuals
	#############################
	say('BIOMASS: map of residuals', level = 2)

		# Using the mean DHARMa residual of a site for site-level residual

		residuals_by_plant <- dharma$scaledResiduals
		residuals_by_site <- rep(NA, constants$n_pheno_sites)
		site_ids <- unique(data_biomass$raw_data_biomass$SITE)
		for (i in 1:constants$n_pheno_sites) residuals_by_site[i] <- mean(residuals_by_plant[data_biomass$raw_data_biomass$SITE == site_ids[i]])

		site_vect_biomass_resid <- data_biomass$site_vect
		site_vect_biomass_resid$residual <- residuals_by_site

		coords <- as.data.frame(crds(project(site_vect_biomass_resid, enmSdmX::getCRS('WGS84'))))

		# Compute Moran's I
		moran <- moran.test(residuals_by_site, nb2listw(knn2nb(knearneigh(coords, longlat = TRUE, k = 4))))
		moran_p <- moran$p.value
		moran_p <- paste0('Moran P = ', sprintf('%.3f', round(moran_p, 3)))

		nam <- vect('./data_from_gadm/gadm_4pt1_level_1_north_america_sans_alaska_lambert.gpkg')
		nam <- simplifyGeom(nam, tolerance = 1000)

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
			annotate('text', x = extent[1], y = extent[3], label = moran_p, hjust = 0, vjust = 0, size = 6, color = 'red') +
			xlim(extent[1], extent[2]) + ylim(extent[3], extent[4]) +
			ggtitle(
				bquote('Residuals for ' * italic('Andropogon gerardi') * ' biomass '),
				subtitle = '1991-2020 | biomass-only model') +
			theme(
				plot.title = element_text(size = 16),
				plot.subtitle = element_text(size = 14)
			)

		ggsave(plot = map, filename = paste0(out_dir, '/map_residuals_biomass_site_mu.png'), width = 12, height = 9, dpi = 120)

	### BIOMASS: plant-level residuals vs covariates
	################################################
	say('BIOMASS: plant-level residuals vs covariates', level = 2)

		resids <- dharma$scaledResiduals
		resids_vs_covariates <- list()
		resids_vs_covariates_results <- data.table()

		for (i in seq_len(data_biomass$n_covariates)) {

			if (data_biomass$covariates[i] == 'ph') {
				x <- data_biomass$raw_data_biomass[['site_ph']]
			} else {
				x <- data_biomass$raw_data_biomass[[data_biomass$covariates[i]]]
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
					covariate = data_biomass$covariates[i],
					F = F,
					p = p,
					significance = sig
				)
			)

			resids_vs_covariates[[i]] <- ggplot(this_x, aes(x = x, y = y)) +
				geom_point() +
				geom_smooth(method = loess, formula = y ~ x, se = FALSE) +
				xlab(data_biomass$terms[i]) +
				ylab('Residual value') +
				ggtitle('DHARMa Residuals for Biomass Means')

		}

		ncol <- data_biomass$n_covariates
		width <- 6 * ncol
		height <- 5 * ceiling(length(resids_vs_covariates) / ncol)

		resids_vs_covariates <- plot_grid(plotlist = resids_vs_covariates, ncol = ncol)
		ggsave(plot = resids_vs_covariates, filename = paste0(out_dir, '/dharma_residuals_vs_covariates_biomass.png'), width = width, height = height, dpi = 200)

	### BIOMASS: pseudo-R2 (correlation between observed and estimate)
	##################################################################
	say('BIOMASS: model fit', level = 2)

		observed <- data_biomass$raw_data_biomass[ , .(mean_biomass = mean(Biomass)), by = SITE][['mean_biomass']]
		preds <- predict_biomass(chains = chains, x = data_biomass$x_by_site, resp_distrib = resp_distrib, transform = transform)
		estimated <- colMeans(preds)
		# estimated <- mc_extract(chains, 'biomass_site', j = TRUE)
		correl_mu <- cor(observed, estimated)

		rmse_mu <- unname(rmse_fx(observed, estimated))
		mean_error_mu <- mean((estimated - observed) / observed)
		mean_abs_error_mu <- mean(abs(estimated - observed) / observed)

	### BIOMASS: response curves
	############################
	say('BIOMASS: response curves', level = 2)
	
		n_covariates <- data_biomass$n_covariates

		responses <- list()
		for (i in seq_len(n_covariates)) {

			covariate <- data_biomass$covariates[i]

			x <- data_biomass$resp_curves_x
			if (n_covariates > 1) x <- x[ , , i]

			preds <- predict_biomass(chains = chains, x = x, resp_distrib = resp_distrib, transform = transform)
			
			hdi <- apply(preds, 2, function(x) hdi(x, credMass = 0.90))
			lowers <- hdi[1, ]
			uppers <- hdi[2, ]

			medians <- apply(preds, 2, median)

			unscaled <- data_biomass$resp_curves_x_unscaled
			x_unscaled <- unscaled[ , covariate, drop = TRUE]
			if (covariate %in% paste0('bio', c(12:14, 16:19), '_log10p1')) x_unscaled <- 10^x_unscaled - 1

			center <- data.table(
				x = x_unscaled,
				y = medians
			)

			cis <- data.table(
				x = c(x_unscaled, rev(x_unscaled)),
				ci = c(lowers, rev(uppers))
			)

			nice_covariate <- get_nice_predictor(covariate)
			title <- paste0('Biomass versus ', nice_covariate$short)
			xlab <- nice_covariate$long

			response <- ggplot() +
				geom_polygon(data = cis, mapping = aes(x = x, y = ci), fill = 'lightblue', alpha = 0.5) +
				geom_line(data = center, mapping = aes(x = x, y = y), color = 'black') +
				xlab(xlab) +
				ylab('Predicted Biomass (g)') +
				ggtitle(title) +
				theme(
					legend.position = 'none',
					plot.title = element_text(size = 14),
					axis.title = element_text(size = 12),
					axis.text = element_text(size = 12)
				)

			data_by_site <- data_biomass$raw_data_biomass

			### add measurments
			if (covariate == 'ph') {
				names(data_by_site)[names(data_by_site) == 'site_ph'] <- 'ph'
				x <- data_biomass$raw_data_biomass[ , .(val = mean(.SD[['site_ph']])), by = SITE][['val']]
			} else {
				x <- data_biomass$raw_data_biomass[ , .(val = mean(.SD[[covariate]])), by = SITE][['val']]
				if (covariate %in% paste0('bio', c(12:14, 16:19), '_log10p1')) x <- 10^x - 1
			}
			
			observed <- data_biomass$raw_data_biomass[ , .(val = median(Biomass)), by = SITE][['val']]
			
			data_by_site_collated <- data.table(
				x = x,
				y = observed,
				site = unique(data_by_site$SITE)
			)

			response <- response + geom_point(
				data = data_by_site_collated,
				mapping = aes(x = x, y = y, fill = site),
				pch = 21, size = 6
			)

			data_by_plant <- data_biomass$raw_data_biomass

			if (covariate == 'ph') {
				data_by_plant$x <- data_by_plant[['site_ph']]
			} else {
				data_by_plant$x <- data_by_plant[[covariate]]
			}
			if (covariate %in% paste0('bio', c(12:14, 16:19), '_log10p1')) data_by_plant$x <- 10^data_by_plant$x - 1
			facet_raw <- get_raw_trait_name_from_rfriendly('biomass')
			data_by_plant$y <- data_by_plant[[facet_raw]]
			data_by_plant$site <- data_by_plant$SITE

			response <- response + geom_point(
				data = data_by_plant,
				mapping = aes(x = x, y = y, fill = site),
				pch = 21, size = 2
			)

			responses[[i]] <- response

		}

		if (n_covariates == 2) {
			ncol <- 2
		} else {
			ncol <- ceiling(sqrt(n_covariates))
		}
		width <- 8 * ncol
		height <- 7 * ceiling(data_biomass$n_covariates / ncol)

		responses <- plot_grid(plotlist = responses, ncol = ncol)

		filename <- paste0(out_dir, '/response_curves_biomass_median_hdpi.png')
		ggsave(plot = responses, filename = filename, width = width, height = height, dpi = 120)

	### BIOMASS-PSI: response curves
	################################
	say('BIOMASS-PSI: response curves', level = 2)
	
		n_covariates <- data_biomass$n_covariates

		responses <- list()
		for (i in seq_len(n_covariates)) {

			covariate <- data_biomass$covariates[i]

			x <- data_biomass$resp_curves_x
			if (n_covariates > 1) x <- x[ , , i]

			preds <- predict_psi(chains = chains, x = x)
			
			hdi <- apply(preds, 2, function(x) hdi(x, credMass = 0.90))
			lowers <- hdi[1, ]
			uppers <- hdi[2, ]

			medians <- apply(preds, 2, median)

			unscaled <- data_biomass$resp_curves_x_unscaled
			x_unscaled <- unscaled[ , covariate, drop = TRUE]
			if (covariate %in% paste0('bio', c(12:14, 16:19), '_log10p1')) x_unscaled <- 10^x_unscaled - 1

			center <- data.table(
				x = x_unscaled,
				y = medians
			)

			cis <- data.table(
				x = c(x_unscaled, rev(x_unscaled)),
				ci = c(lowers, rev(uppers))
			)

			nice_covariate <- get_nice_predictor(covariate)
			title <- paste0('Probability of Occurrence versus ', nice_covariate$short)
			xlab <- nice_covariate$long
			ylab <- 'Probability of Occurrence'

			response <- ggplot() +
				geom_polygon(data = cis, mapping = aes(x = x, y = ci), fill = 'coral', alpha = 0.5) +
				geom_line(data = center, mapping = aes(x = x, y = y), color = 'darkred') +
				xlab(xlab) +
				ylab(ylab) +
				ylim(0, 1) +
				ggtitle(title, subtitle = facet_nice$short) +
				theme(
					legend.position = 'none',
					plot.title = element_text(size = 14),
					axis.title = element_text(size = 12),
					axis.text = element_text(size = 12)
				)

			responses[[i]] <- response

		}

		if (n_covariates == 2) {
			ncol <- 2
		} else {
			ncol <- ceiling(sqrt(n_covariates))
		}
		width <- 8 * ncol
		height <- 7 * ceiling(data_biomass$n_covariates / ncol)

		responses <- plot_grid(plotlist = responses, ncol = ncol)

		filename <- paste0(out_dir, '/response_curves_psi_biomass_median_hdpi.png')
		ggsave(plot = responses, filename = filename, width = width, height = height, dpi = 120)

	### BIOMASS: maps
	#################
	say('BIOMASS: burn prediction vectors', level = 2)

		pred_vect_nam <- burn_biomass_into_vector(demesne = 'nam', chains = chains, formula_biomass = formula_biomass, formula_psi = formula_psi, resp_distrib = resp_distrib, transform = transform)
		
		pred_vect_1930s <- burn_biomass_into_vector(demesne = '1930s', chains = chains, formula_biomass = formula_biomass, formula_psi = formula_psi, resp_distrib = resp_distrib, transform = transform)

		writeVector(pred_vect_nam, paste0(out_dir, '/prediction_vector_nam.gpkg'), overwrite = TRUE)
		writeVector(pred_vect_1930s, paste0(out_dir, '/prediction_vector_conus_1930s.gpkg'), overwrite = TRUE)

	### BIOMASS: current map
	########################
	say('BIOMASS: current map', level = 2)

		mean_form <- paste(as.character(formula_biomass), collapse = ' ')
		mean_form <- gsub(mean_form, pattern = 'I\\(', replacement = '')
		mean_form <- gsub(mean_form, pattern = '\\^2\\)', replacement = '²')
		mean_form <- gsub(mean_form, pattern = '*)', replacement = '×')
		mean_form <- gsub(mean_form, pattern = '\\:', replacement = '×')

		subtitle <- paste0('Biomass ', mean_form, ' | 1991-2020')

		map <- map_biomass_nonbiomass(
			facet = 'biomass',
			out_dir = out_dir,
			filename_append = 'present',
			pred_vect_nam = pred_vect_nam,
			response_var = 'biomass_mean_sq',
			response_var_type = 'mean',
			data_occs = data_occs,
			data_biomass_nonbiomass = data_biomass,
			title = bquote('Present-day Mean Biomass'),
			subtitle = subtitle,
			legend_title = 'Biomass\n(g)'
		)

		if (!is.null(formula_psi)) {
			
			map <- map_psi(
				out_dir = out_dir,
				filename_append = 'present',
				pred_vect_nam = pred_vect_nam,
				response_var = 'psi_sq',
				title = bquote('Present-day Probability of Occurrence'),
				subtitle = subtitle
			)

		}

	### BIOMASS: future maps
	########################
	say('BIOMASS: future maps', level = 2)

		for (fut in futs) {

			say('   ', fut)

			response_var <- paste0('biomass_mean_', fut)

			mean_form <- paste(as.character(formula_biomass), collapse = ' ')
			subtitle <- paste0('Biomass ', mean_form, '  | SSP ', substr(fut, 4, 6), ' ', substr(fut, 8, 11), '-', substr(fut, 13, 16))

			map_biomass_fut <- map_biomass_nonbiomass(
				facet = 'biomass',
				out_dir = out_dir,
				filename_append = fut,
				pred_vect_nam = pred_vect_nam,
				response_var = response_var,
				response_var_type = 'mean',
				data_occs = data_occs,
				data_biomass_nonbiomass = data_biomass,
				title = bquote('Future Mean Biomass'),
				subtitle = subtitle,
				legend_title = 'Biomass\n(g)'
			)

			if (!is.null(formula_psi)) {

				facet <- 'biomass'				
				map <- map_psi(
					out_dir = out_dir,
					filename_append = fut,
					pred_vect_nam = pred_vect_nam,
					response_var = paste0('psi_', fut),
					title = bquote('Future Probability of Occurrence'),
					subtitle = subtitle
				)

			}

		}

	### BIOMASS: future change maps
	###############################
	say('BIOMASS: future change maps', level = 2)

		for (fut in futs) {

			say('   ', fut)

			response_var <- paste0('biomass_mean_', fut)
			title <- bquote('Change in Mean Biomass ')

			mean_form <- paste(as.character(formula_biomass), collapse = ' ')
			subtitle <- paste0('Biomass ', mean_form, '  | SSP ', substr(fut, 4, 6), ' ', substr(fut, 8, 11), '-', substr(fut, 13, 16))

			map <- map_biomass_nonbiomass_change(
				facet = facet,
				out_dir = out_dir,
				filename_append = fut,
				fut = fut,
				pred_vect_nam = pred_vect_nam,
				response_var = response_var,
				response_var_type = 'mean',
				data_occs = data_occs,
				data_biomass_nonbiomass = data_biomass,
				title = title,
				subtitle = subtitle,
				legend_title = 'Percent\nchange'
			)

			if (!is.null(formula_psi)) {

				map <- map_psi_change(
					out_dir = out_dir,
					filename_append = fut,
					fut = fut,
					pred_vect_nam = pred_vect_nam,
					response_var = paste0('psi_', fut),
					response_var_sq = 'psi_sq',
					title = 'Change in Probability of Occurrence',
					subtitle = subtitle
				)

			}

		}

	### OCCURRENCE: 1930s change maps
	#################################
	say('BIOMASS: 1930s change maps', level = 2)

		map_thirties <- map_biomass_nonbiomass_change_1930s(
			facet = 'biomass',
			data_biomass = data_biomass,
			pred_vect_nam = pred_vect_nam,
			pred_vect_1930s = pred_vect_1930s,
			formula_psi = formula_psi
		)

	### occurrence: crossvalidation
	###############################

	if (crossvalidate) {
		
		say('BIOMASS: cross-validation', level = 2)

		max_k_folds <- if (trial) { 1 } else { k_folds }
		cv_df <- data.table()
		for (k in 1:max_k_folds) {

			say('fold: ', k, level = 3)
		
			fold_constants <- constants
			fold_data <- data
			fold_inits <- inits

			indices <- which(data_biomass$raw_data_biomass$geofold == k)

			# data
			fold_data$y_biomass <- fold_data$y_biomass[-indices]

			# constants: order matters in assigning these!
			test_sites <- sort(unique(fold_constants$site_index_biomass[indices]))
			train_sites <- sort(unique(fold_constants$site_index_biomass[-indices]))
			
			fold_constants$site_index_biomass <- renumSeq(fold_constants$site_index_biomass[fold_constants$site_index_biomass %in% train_sites])
			fold_constants$n_biomass <- fold_constants$n_biomass - length(indices)
			fold_constants$x_by_site_biomass <- fold_constants$x_by_site_biomass[train_sites, ]
			fold_constants$n_pheno_sites <- length(train_sites)
		
			# inits
			fold_inits$log_site_biomass_mu <- fold_inits$log_site_biomass_mu[train_sites]
			fold_inits$y_biomass_sim <- fold_inits$y_biomass_sim[-indices]
			fold_inits$beta_biomass <- mc_extract(chains, 'beta_biomass', j = TRUE)
			fold_inits$beta_psi <- mc_extract(chains, 'beta_psi', j = TRUE)
			fold_inits$log_sigma_biomass_within_sites <- log(mc_extract(chains, 'sigma_biomass_within_sites'))

			fold_model <- nimbleModel(
				code = model_code,
				constants = fold_constants,
				data = fold_data,
				inits = fold_inits,
				check = FALSE,
				calculate = FALSE,
				buildDerivs = TRUE
			)

			fold_monitors <- c(monitors_coeffs_not_indexed, monitors_coeffs_single_index, monitors_coeffs_double_index)

			fold_conf <- configureMCMC(
				fold_model,
				monitors = fold_monitors,
				print = FALSE,
				enableWAIC = FALSE
			)

			vars <- 'beta_biomass'
			fold_conf$removeSamplers(vars)
			fold_conf$addSampler(target = vars, type = 'AF_slice')
			say('AF slice sampler added to ', paste(vars, collapse = ' & '), '.')

			vars <- 'beta_psi'
			fold_conf$removeSamplers(vars)
			fold_conf$addSampler(target = vars, type = 'AF_slice')
			say('AF slice sampler added to ', paste(vars, collapse = ' & '), '.')

			fold_build <- buildMCMC(fold_conf)
			fold_compiled <- compileNimble(fold_model, fold_build, showCompilerOutput = FALSE)

			fold_chains <- runMCMC(
				fold_compiled$fold_build,
				niter = niter,
				nburnin = nburnin,
				thin = thin,
				nchains = 2,
				inits = fold_inits,
				progressBar = TRUE,
				samplesAsCodaMCMC = TRUE,
				summary = TRUE,
				WAIC = FALSE,
				perChainWAIC = FALSE
			)

			fold_x <- constants$x_by_site[test_sites, ]
			preds <- predict_biomass(chains = fold_chains, x = fold_x, resp_distrib = resp_distrib, transform = transform)

			### evaluate predictions
			# observed site means
			site_means <- c()
			for (i in seq_along(test_sites)) {
				
				test_site <- test_sites[i]
				site_means <- c(site_means, mean(data$y_biomass[constants$site_index_biomass == test_site]))

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

			cv_df <- rbind(
				cv_df,
				data.table(
					k = k,

					n_train_sites = fold_constants$n_pheno_sites,
					n_test_sites = constants$n_pheno_sites - fold_constants$n_pheno_sites,

					n_train_plants = fold_constants$n_biomass,
					n_test_plants = constants$n_biomass - fold_constants$n_biomass,
					
					correl_lower = quantile(correl, 0.025, na.rm = TRUE),
					correl_mean = mean(correl, na.rm = TRUE),
					correl_upper = quantile(correl, 0.975, na.rm = TRUE),

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

		} # next fold

		cv_df <- rbind(
			cv_df,
			data.table(
					k = 'means',

					n_train_sites = mean(cv_df$n_train_sites),
					n_test_sites = mean(cv_df$n_test_sites),

					n_train_plants = mean(cv_df$n_train_plants),
					n_test_plants = mean(cv_df$n_test_plants),
					
					correl_lower = mean(cv_df$correl_lower),
					correl_mean = mean(cv_df$correl_mean),
					correl_upper = mean(cv_df$correl_upper),

					mae_lower = mean(cv_df$mae_lower),
					mae_mean = mean(cv_df$mae_mean),
					mae_upper = mean(cv_df$mae_upper),

					mape_lower = mean(cv_df$mape_lower),
					mape_mean = mean(cv_df$mape_mean),
					mape_upper = mean(cv_df$mape_upper),

					rmse_lower = mean(cv_df$rmse_lower),
					rmse_mean = mean(cv_df$rmse_mean),
					rmse_upper = mean(cv_df$rmse_upper),

					obs_quants_min = mean(cv_df$obs_quants_min),
					obs_quants_mean = mean(cv_df$obs_quants_mean),
					obs_quants_max = mean(cv_df$obs_quants_max),
					obs_quants_prop_in_inner_90 = mean(cv_df$obs_quants_prop_in_inner_90)

			)
		)

	} # if crossvalidating

	### remember
	############

		resid_p_values <- c(
			moran$p.value,
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
			resp_distrib = resp_distrib,
			transform = transform,
			dharma_resids = dharma_resids,
			dharma_resids_vs_covariates = resids_vs_covariates_results,
			fit_mu = c(correl_mu = correl_mu, rmse_mu = rmse_mu, mean_error_mu = mean_error_mu, mean_abs_error_mu = mean_abs_error_mu)
		)

		meta_biomass$crossvalidation <- if (crossvalidate) {
			cv_df
		} else {
			NA
		}

		saveRDS(meta_biomass, paste0(out_dir, '/!meta_biomass.rds'))
		sink(paste0(out_dir, '/!meta_biomass.txt'), split = TRUE)
			print(meta_biomass)
		sink()

} # EOF

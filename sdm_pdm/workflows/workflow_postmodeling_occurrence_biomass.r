#' Post-modeling workflow for occurrence-only models and models with an occurrence component.
#'
#' source('C:/Kaji/R/andropogon_integratedEcology/sdm_pdm/workflows/workflow_postmodeling_occurrence_biomass.r')
#'
#' formula_... 		Formulae for occurrences/biomass mean, sigma, and probability of zero-inflation, and (for occurrences) bias.
#' resp_distrib 	Named vector of response distributions. For occurrence, this can be 'Poisson' or 'ZIP'. For biomass this can be 'gamma', 'hGamma' (zero-inflated gamma), 'lognormal', or 'hurdleLN' (zero-inflated lognormal)
#' transform		Named vector of transformations to translate MVN to mean occurrence intensity or biomass: 'identity', 'softplus' or 'exponential'.
#' out_dir 			Folder into which to save results.
workflow_postmodeling_occurrence_biomass <- function(
	formula_occs,
	formula_bias,
	formula_biomass,
	formula_psi,
	resp_distrib,
	transform,
	out_dir
) {

	zero_inflated <- !is.null(formula_psi)

	if (zero_inflated) data_traits <- prepare_nonbiomass_facets(trait = 'height', formula = ~ 1, n_response_curve_values = n_response_curve_values, calib = calib)

	### burn predictions into vector
	################################
	say('OCCURRENCE + BIOMASS: burn prediction vectors', level = 2)

		pred_vect_nam <- burn_occs_biomass_into_vector(
			demesne = 'nam',
			chains = chains,
			formula_occs = formula_occs,
			formula_biomass = formula_biomass,
			formula_psi = formula_psi,
			resp_distrib = resp_distrib,
			transform = transform
		)

		pred_vect_1930s <- burn_occs_biomass_into_vector(
			demesne = '1930s',
			chains = chains,
			formula_occs = formula_occs,
			formula_biomass = formula_biomass,
			formula_psi = formula_psi,
			resp_distrib = resp_distrib,
			transform = transform
		)

		writeVector(pred_vect_nam, paste0(out_dir, '/prediction_vector_nam.gpkg'), overwrite = TRUE)
		writeVector(pred_vect_1930s, paste0(out_dir, '/prediction_vector_conus_1930s.gpkg'), overwrite = TRUE)

	### response curves: OCCURRENCE vs environment
	##############################################
	say('OCCURRENCE + BIOMASS: response curves', level = 2)

		responses_mu <- graph_response_curves_occurrence_vs_environment(
			out_dir = out_dir,
			resp_type = 'mu',
			chains = chains,
			data = data_occs
		)

		if (zero_inflated) {

			responses_psi <- graph_response_curves_occurrence_vs_environment(
				out_dir = out_dir,
				resp_type = 'psi',
				chains = chains,
				data = data_county_psi
			)

		}

		if (data_occs$n_covariates_bias >= 1) {

			responses_bias <- graph_response_curves_occurrence_bias_vs_bias_covariates(
				out_dir = out_dir,
				chains = chains,
				data_occs = data_occs
			)

		}


	### OCCURRENCE & BIOMASS: facet vs facet responses
	##################################################
	say('OCCURRENCE & BIOMASS: facet vs facet responses', level = 2)

		x <- data_occs$ag_vect_sq$bio12[data_occs$ag_vect_sq$n_andropogon_gerardi > 0]
		x_min <- min(x)
		x_max <- max(x)

		for (fut in futs) {

			x <- data_occs[[paste0('counties_occs_', fut)]]
			x <- x$bio12[pred_vect_nam$n_andropogon_gerardi > 0]

			x_min <- min(x_min, x)
			x_max <- min(x_max, x)
		
		}

		mult <- 10
		covariates <- data_occs$covariates
		presences <- data_occs$ag_vect_sq$n_andropogon_gerardi > 0
		x <- data_occs$ag_vect_sq[presences, covariates]
		x <- as.data.table(x)
		x <- scale(x, center = data_occs$x_centers_occs, scale = data_occs$x_scales)
		means <- colMeans(x)
		x <- cbind(matrix(rep(means, mult * n_response_curve_values), byrow = TRUE, ncol = 3))
		colnames(x) <- covariates
		x_seq <- seq(x_min, x_max, length.out = mult * n_response_curve_values)
		x_seq <- (x_seq - data_occs$x_centers_occs['bio12']) / data_occs$x_scales['bio12']
		x[ , 'bio12'] <- x_seq
		x <- as.data.table(x)

		x_occs <- model.matrix(formula_occs, x)
		x_biomass <- model.matrix(formula_biomass, x)
		x_psi <- model.matrix(formula_psi, x)

		preds <- predict_occs_biomass(chains = chains, resp_distrib = resp_distrib, transform = transform, x_occs = x_occs, x_biomass = x_biomass, x_psi = x_psi)

		preds_occs <- colMeans(preds$preds_occs)
		preds_biomass <- colMeans(preds$preds_biomass)

		x_unscaled <- x_seq * data_occs$x_scales['bio12'] + data_occs$x_centers_occs['bio12']
		df_joint <- data.table(
			occurrence = preds_occs,
			biomass = preds_biomass,
			bio12 = x_unscaled
		)

		# hexagonal heatmap of predicted occurrence vs biomass,
		# filled by the mean bio12 within each hexagon
		facet_vs_facet <- ggplot(df_joint, aes(x = occurrence, y = biomass)) +
			stat_summary_hex(aes(z = bio12), fun = mean, bins = 30, color = NA) +
			scale_fill_gradient(name = 'BIO 12\n(mean)', low = 'yellow', high = 'blue') +
			xlab('Predicted abundance') +
			ylab('Predicted biomass (g)') +
			ggtitle('Occurrence vs biomass')

		ggsave(filename = paste0(out_dir, '/response_curves_occs_vs_biomass.png'), plot = facet_vs_facet, width = 8, height = 8, dpi = 600, bg = 'white')

	### OCCURRENCE: dharma residuals
	################################
	say('OCCURRENCE: dharma residuals', level = 2)

		sims <- mc_subset(chains, 'y_n_ag_sim', j = TRUE)
		sims <- mc_rbind(sims)
		sims <- sims[ , data_occs_counties$ag_vect_sq$focal_region]

		observed_y <- data_occs_counties$y_n_ag[data_occs_counties$ag_vect_sq$focal_region]

		nas_occs <- which(is.na(colSums(sims)))
		if (length(nas_occs) > 0) {
			sims <- sims[ , -nas_occs]
			observed_y <- observed_y[-nas_occs]
		}

		sims <- t(sims)
		dharma_occs <- createDHARMa(simulatedResponse = sims, observedResponse = observed_y, integerResponse = TRUE)

		dharma_quant_test_occs <- testQuantiles(dharma_occs, plot = FALSE)
		dharma_resid_test_occs <- testResiduals(dharma_occs, plot = FALSE) # uniformity, dispersion, outlier
		dev.off()

		file <- paste0(out_dir, '/dharma_n_ag.png')
		png(file, width = 1200, height = 800)
			plot(dharma_occs)
		dev.off()

		file <- paste0(out_dir, '/dharma_lambda_residuals_abundance_residuals.png')
		png(file, width = 1200, height = 800)
			hist(dharma_occs$scaledResiduals, main = 'DHARMa residuals for number of observed AG (y_n_ag)', xlab = 'Scaled residuals', breaks = 30)
		dev.off()

	### OCCURRENCE: spatial autocorrelation
	#######################################

		coords <- as.data.frame(crds(centroids(project(pred_vect_nam[pred_vect_nam$focal_region], enmSdmX::getCRS('WGS84')))))
		if (length(nas_occs) > 0) coords <- coords[-nas_occs, ]

		# Compute Moran's I
		moran_occs <- moran.test(dharma_occs$scaledResiduals, nb2listw(knn2nb(knearneigh(coords, longlat = TRUE, k = 4))))
		moran_p <- moran_occs$p.value
		moran_p <- paste0('Moran P = ', sprintf('%.3f', round(moran_p, 3)))

		residuals_vect <- pred_vect_nam[pred_vect_nam$focal_region]
		if (length(nas_occs) > 0) residuals_vect <- residuals_vect[-nas_occs]
		residuals_vect$residual <- dharma_occs$scaledResiduals

		extent <- ext(residuals_vect[residuals_vect$focal_region])
		extent <- as.vector(extent)

		map <- ggplot() +
			layer_spatial(pred_vect_nam, fill = 'gray50', color = 'gray40') +
			layer_spatial(residuals_vect, aes(fill = residual), color = NA) +
			scale_fill_distiller(palette = 'RdBu', limits = c(0, 1), na.value = 'grey80') +
			annotate('text', x = extent[1], y = extent[3], label = moran_p, hjust = 0, vjust = 0, size = 6, color = 'red') +
			xlim(extent[1], extent[2]) + ylim(extent[3], extent[4]) +
			ggtitle(bquote('Residuals for ' * italic('Andropogon gerardi') * ' occurrence')) +
			theme(
				plot.title = element_text(size = 16),
				plot.subtitle = element_text(size = 14)
			)

		ggsave(plot = map, filename = paste0(out_dir, '/map_residuals_occurrences.png'), width = 12, height = 9, dpi = 300)

	### OCCURRENCE: current map
	###########################
	say('OCCURRENCE: current map', level = 2)

		form <- paste(as.character(formula_occs), collapse = ' ')
		form <- gsub(form, pattern = 'I\\(', replacement = '')
		form <- gsub(form, pattern = '\\^2\\)', replacement = '²')
		form <- gsub(form, pattern = '*)', replacement = '×')

		form_bias <- paste(as.character(formula_bias), collapse = ' ')
		form_bias <- gsub(form_bias, pattern = 'I\\(', replacement = '')
		form_bias <- gsub(form_bias, pattern = '\\^2\\)', replacement = '²')
		form_bias <- gsub(form_bias, pattern = '*)', replacement = '×')

		if (!is.null(formula_psi)) {
			
			form_psi <- paste(as.character(formula_psi), collapse = ' ')
			form_psi <- gsub(form_psi, pattern = 'I\\(', replacement = '')
			form_psi <- gsub(form_psi, pattern = '\\^2\\)', replacement = '²')
			form_psi <- gsub(form_psi, pattern = '*)', replacement = '×')

			subtitle <- paste0('occ ', form, ' (bias ', form_bias, ')\nψ ', form_psi)

		} else {
		
			subtitle <- paste0('occ ', form, ' (bias ', form_bias, ')')

		}

		map <- map_occurrence(
			out_dir = out_dir,
			filename_append = 'present_day',
			pred_vect_nam = pred_vect_nam,
			response_var = 'N_ag_county_mean_sq',
			data_occs = data_occs,
			title = bquote('Present-day distribution of ' * italic('Andropogon gerardi') * ' abundance (1961-2020)'),
			subtitle = subtitle
		)

		if (zero_inflated) {

			form_psi <- paste(as.character(formula_psi), collapse = ' ')
			form_psi <- gsub(form_psi, pattern = 'I\\(', replacement = '')
			form_psi <- gsub(form_psi, pattern = '\\^2\\)', replacement = '²')
			form_psi <- gsub(form_psi, pattern = '*)', replacement = '×')

			psi_map <- map_psi(
				out_dir = out_dir,
				filename_append = 'present_day',
				pred_vect_nam = pred_vect_nam,
				response_var = 'psi_county_sq',
				title = 'Present-day probability of occurrence (1961-2020)',
				subtitle = subtitle
			)

		}

	### occurrence: future maps
	###########################
	say('OCCURRENCE: future maps', level = 2)

		maps_occs_fut <- list()
		for (fut in futs) {

			say(fut)

			form <- paste(as.character(formula_occs), collapse = ' ')
			form <- gsub(form, pattern = 'I\\(', replacement = '')
			form <- gsub(form, pattern = '\\^2\\)', replacement = '²')
			form <- gsub(form, pattern = '*)', replacement = '×')

			form_bias <- paste(as.character(formula_bias), collapse = ' ')
			form_bias <- gsub(form_bias, pattern = 'I\\(', replacement = '')
			form_bias <- gsub(form_bias, pattern = '\\^2\\)', replacement = '²')
			form_bias <- gsub(form_bias, pattern = '*)', replacement = '×')

			ssp <- paste0('SSP ', substr(fut, 4, 6), ' ', substr(fut, 8, 11), '-', substr(fut, 13, 16))
			title <- bquote('Future distribution of ' * italic('Andropogon gerardi') * ' abundance (' * ssp * ')')

			if (!is.null(formula_psi)) {
				
				form_psi <- paste(as.character(formula_psi), collapse = ' ')
				form_psi <- gsub(form_psi, pattern = 'I\\(', replacement = '')
				form_psi <- gsub(form_psi, pattern = '\\^2\\)', replacement = '²')
				form_psi <- gsub(form_psi, pattern = '*)', replacement = '×')

				subtitle <- paste0('occ ', form, ' (bias ', form_bias, ')\nψ ', form_psi)

			} else {
			
				subtitle <- paste0('occ ', form, ' (bias ', form_bias, ')')

			}

			response_var <- paste0('N_ag_county_mean_', fut)

			maps_occs_fut[[length(maps_occs_fut) + 1]] <- map_occurrence(
				out_dir = out_dir,
				filename_append = fut,
				pred_vect_nam = pred_vect_nam,
				response_var = response_var,
				data_occs = data_occs,
				title = title,
				subtitle = subtitle
			)

			if (zero_inflated) {

				title <- bquote('Future ' * italic('Andropogon gerardi') * ' probability of occurrence (' * ssp * ')')

				psi_map <- map_psi(
					out_dir = out_dir,
					filename_append = fut,
					pred_vect_nam = pred_vect_nam,
					response_var = paste0('psi_county_', fut),
					title = title,
					subtitle = subtitle
				)

			}

		}

	### occurrence: future change
	#############################
	say('OCCURRENCE: future change', level = 2)

		for (fut in futs) {

			say(fut)

			ssp <- paste0('SSP ', substr(fut, 4, 6), ' ', substr(fut, 8, 11), '-', substr(fut, 13, 16))
			title <- bquote('Change in distribution of ' * italic('Andropogon gerardi') * ' abundance (' * ssp * ')')

			response_var <- paste0('N_ag_county_mean_', fut)

			map <- map_occurrence_change(
				out_dir = out_dir,
				filename_append = fut,
				pred_vect_nam = pred_vect_nam,
				response_var = response_var,
				data_occs = data_occs,
				title = title,
				subtitle = subtitle
			)

			if (zero_inflated) {

				title <- bquote('Change in ' * italic('Andropogon gerardi') * ' probability of occurrence (' * ssp * ')')

				psi_map <- map_psi_change(
					out_dir = out_dir,
					filename_append = fut,
					fut = fut,
					pred_vect_nam = pred_vect_nam,
					response_var = paste0('psi_county_', fut),
					response_var_sq = 'psi_county_sq',
					title = title,
					subtitle = subtitle
				)

			}

		}

	### OCCURRENCE: 1930s
	#####################
	say('OCCURRENCE: 1930s change maps', level = 2)

		map_1930s <- map_occurrence_change_1930s(zero_inflated = zero_inflated, pred_vect_nam = pred_vect_nam, pred_vect_1930s = pred_vect_1930s)

	### summary
	###########
	
		if (length(dharma_quant_test_occs$qgamFits) > 0) {
			
			quant_test_upper <- summary(dharma_quant_test_occs$qgamFits[[3]])$s.table[1, 4]
			quant_test_middle <- summary(dharma_quant_test_occs$qgamFits[[2]])$s.table[1, 4]
			quant_test_lower <- summary(dharma_quant_test_occs$qgamFits[[1]])$s.table[1, 4]

		} else {
			quant_test_lower <- quant_test_middle <- quant_test_upper <- NA_real_
		}

		resid_p_values <- c(
			moran_occs$p.value,
			dharma_resid_test_occs$uniformity$p.value,
			dharma_resid_test_occs$dispersion$p.value,
			dharma_resid_test_occs$outliers$p.value,
			dharma_quant_test_occs$p.value,
			quant_test_upper,
			quant_test_middle,
			quant_test_lower
		)

		if (length(dharma_quant_test_occs$qgamFits) > 0) {
			
			quant_test_upper <- summary(dharma_quant_test_occs$qgamFits[[3]])$s.table[1, 4]
			quant_test_middle <- summary(dharma_quant_test_occs$qgamFits[[2]])$s.table[1, 3]
			quant_test_lower <- summary(dharma_quant_test_occs$qgamFits[[1]])$s.table[1, 3]

		} else {
			quant_test_lower <- quant_test_middle <- quant_test_upper <- NA_real_
		}

		resid_p_values_sig <- ifelse(resid_p_values < 0.05, '*', 'ns')
		resid_test_statistic_values <- c(
			moran_occs$statistic,
			dharma_resid_test_occs$uniformity$statistic,
			dharma_resid_test_occs$dispersion$statistic,
			dharma_resid_test_occs$outliers$statistic,
			NA,
			quant_test_upper,
			quant_test_middle,
			quant_test_lower
		)

		dharma_resids <- data.table(
			test = c('spatial autocorrelation', 'uniformity', 'dispersion', 'outliers', 'quantiles, overall', 'quantiles, upper', 'quantiles, middle', 'quantiles, lower'),
			p_value = resid_p_values,
			significant = resid_p_values_sig,
			test_statistic = c('Moran\'s I', names(dharma_resid_test_occs$uniformity$statistic), names(dharma_resid_test_occs$dispersion$statistic), 'exact binomial', NA, rep('chi squared', 3)),
			test_statistic_value = resid_test_statistic_values
		)

		# compile residuals analysis
		meta_occs <- list(
			facet = 'occurrence',
			date = date(),
			zero_inflated = zero_inflated,
			formulae = list(
				formula_occs = formula_occs,
				formula_bias = formula_bias,
				formula_psi = formula_psi
			),
			dharma_resids = dharma_resids
		)

		saveRDS(meta_occs, paste0(out_dir, '/!meta_occs.rds'))
		sink(paste0(out_dir, '/!meta_occs.txt'), split = TRUE)
			print(meta_occs)
		sink()

###############################################################
### BIOMASS BIOMASS BIOMASS BIOMASS BIOMASS BIOMASS BIOMASS ###
###############################################################

	### BIOMASS: plant-level residuals analysis
	###########################################
	say('BIOMASS: plant-level residuals analysis', level = 2)

		# NB this uses the mean predicted value of a site as a plant-level prediction
		sims_by_plant <- mc_subset(chains, param = 'y_biomass_sim', j = TRUE)
		sims_by_plant <- mc_rbind(sims_by_plant)
		sims_by_plant <- t(sims_by_plant)

		preds <- predict_biomass(chains = chains, x = data_biomass$x_by_site, resp_distrib = resp_distrib, transform = transform)
		estimates <- colMeans(preds)

		# estimates <- mc_subset(chains, param = 'biomass_site', j = TRUE)
		# estimates <- mc_rbind(estimates)

		site_counts <- data_biomass$raw_data_biomass[ , .N, by = SITE]
		fit_by_site <- apply(estimates, 2, median)
		fits <- numeric()
		for (i in seq_along(fit_by_site)) {
			fits <- c(fits, rep(fit_by_site[i], site_counts$N[i]))
		}

		dharma_biomass <- createDHARMa(simulatedResponse = sims_by_plant, observedResponse = data_biomass$y_biomass, fittedPredictedResponse = fits, integerResponse = FALSE)

		dharma_quant_test_biomass <- testQuantiles(dharma_biomass, plot = FALSE)
		dharma_resid_test_biomass <- testResiduals(dharma_biomass, plot = FALSE) # uniformity, dispersion, outlier
		dev.off()

		file <- paste0(out_dir, '/dharma_by_plant_biomass.png')
		png(file, width = 1200, height = 800)
			plot(dharma_biomass)
		dev.off()
		
		file <- paste0(out_dir, '/dharma_by_plant_biomass_residuals.png')
		png(file, width = 1400, height = 1000)
			hist(dharma_biomass$scaledResiduals, main = 'dharma_biomass residuals for biomass by plant', xlab = 'Scaled residuals', breaks = 30)
		dev.off()

	### BIOMASS: site-level residuals analysis
	##########################################

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

		# Using the mean dharma_biomass residual of a site for site-level residual
		residuals_by_plant <- dharma_biomass$scaledResiduals
		residuals_by_site <- rep(NA, constants$n_pheno_sites)
		site_ids <- unique(data_biomass$raw_data_biomass$SITE)
		for (i in 1:constants$n_pheno_sites) residuals_by_site[i] <- mean(residuals_by_plant[data_biomass$raw_data_biomass$SITE == site_ids[i]])

		site_vect_biomass_resid <- data_biomass$site_vect
		site_vect_biomass_resid$residual <- residuals_by_site

		nam <- vect('./data_from_gadm/gadm_4pt1_level_1_north_america_sans_alaska_lambert.gpkg')
		nam <- simplifyGeom(nam, tolerance = 1000)

		# extent
		extent <- buffer(site_vect_biomass_resid, width = 200 * 1000)
		extent <- ext(extent)
		extent <- as.vector(extent)

		# Compute Moran's I
		coords <- as.data.frame(crds(project(site_vect_biomass_resid, enmSdmX::getCRS('WGS84'))))
		moran_biomass <- moran.test(residuals_by_site, nb2listw(knn2nb(knearneigh(coords, longlat = TRUE, k = 4))))
		moran_p <- moran_biomass$p.value
		moran_p <- paste0('Moran P = ', sprintf('%.3f', round(moran_p, 3)))

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

		ggsave(plot = map, filename = paste0(out_dir, '/map_residuals_biomass_site_mu.png'), width = 12, height = 9, dpi = 300)

	### BIOMASS: plant-level residuals vs covariates
	################################################
	say('BIOMASS: plant-level residuals vs covariates', level = 2)

		resids <- dharma_biomass$scaledResiduals
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
				ggtitle('DHAMRa Residuals for Biomass Means')

		}

		if (data_biomass$n_covariates == 2) {
			ncol <- 2
		} else {
			ncol <- ceiling(sqrt(length(resids_vs_covariates)))
		}
		width <- 6 * ncol
		height <- 5 * ceiling(length(resids_vs_covariates) / ncol)

		resids_vs_covariates <- plot_grid(plotlist = resids_vs_covariates, ncol = ncol)
		ggsave(plot = resids_vs_covariates, filename = paste0(out_dir, '/dharma_residuals_vs_covariates_biomass.png'), width = width, height = height, dpi = 300)

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
	
		responses <- graph_response_curves_biomass_vs_environment(
			out_dir = out_dir,
			chains = chains,
			resp_type = 'mu',
			data_biomass = data_biomass,
			quant_threshold = 0.5
		)

	### BIOMASS: current map
	########################
	say('BIOMASS: current map', level = 2)

		mean_form <- paste(as.character(formula_biomass), collapse = ' ')
		mean_form <- gsub(mean_form, pattern = 'I\\(', replacement = '')
		mean_form <- gsub(mean_form, pattern = '\\^2\\)', replacement = '²')
		mean_form <- gsub(mean_form, pattern = '*)', replacement = '×')

		subtitle <- paste0('Biomass μ ', mean_form, ' | 1991-2020')

		map <- map_biomass_nonbiomass(
			out_dir = out_dir,
			filename_append = 'present',
			pred_vect_nam = pred_vect_nam,
			response_var = 'biomass_county_median_sq',
			response_var_type = 'mu',
			data_occs = data_occs,
			data_biomass_nonbiomass = data_biomass,
			title = bquote('Present-day distribution of mean ' * italic('Andropogon gerardi') * ' biomass '),
			subtitle = subtitle,
			legend_title = 'Biomass (g)'
		)

	### BIOMASS: future maps
	########################
	say('BIOMASS: future maps', level = 2)

		for (fut in futs) {

			say(fut)

			response_var <- paste0('biomass_county_median_', fut)

			mean_form <- paste(as.character(formula_biomass), collapse = ' ')
			subtitle <- paste0('Biomass μ ', mean_form, '  | SSP ', substr(fut, 4, 6), ' ', substr(fut, 8, 11), '-', substr(fut, 13, 16))
			
			map_biomass_fut <- map_biomass_nonbiomass(
				out_dir = out_dir,
				filename_append = fut,
				pred_vect_nam = pred_vect_nam,
				response_var = response_var,
				response_var_type = 'mu',
				data_occs = data_occs,
				data_biomass_nonbiomass = data_biomass,
				title = bquote('Future distribution of mean ' * italic('Andropogon gerardi') * ' biomass '),
				subtitle = subtitle,
				legend_title = 'Biomass (g)'
			)

		}

	### BIOMASS: future change maps
	###############################
	say('BIOMASS: future change maps', level = 2)

		site_vect <- data_biomass$site_vect

		maps_biomass_change <- list()
		for (fut in futs) {

			say(fut)

			response_var <- paste0('biomass_county_median_', fut)
			title <- bquote('Change in ' * italic('Andropogon gerardi') * ' biomass ')

			mean_form <- paste(as.character(formula_biomass), collapse = ' ')
			subtitle <- paste0('Biomass μ ', mean_form, '  | SSP ', substr(fut, 4, 6), ' ', substr(fut, 8, 11), '-', substr(fut, 13, 16))

			maps_biomass_change[[length(maps_biomass_change) + 1]] <- map_biomass_change(
				out_dir = out_dir,
				filename_append = fut,
				fut = fut,
				pred_vect_nam = pred_vect_nam,
				response_var = response_var,
				response_var_type = 'mean',
				data_occs = data_occs,
				data_biomass = data_biomass,
				title = title,
				subtitle = subtitle
			)

		}

	### BIOMASS: 1930s change maps
	##############################
	say('BIOMASS: 1930s change maps', level = 2)

		map_thirties <- map_biomass_nonbiomass_change_1930s(
			facet = 'biomass',
			data_biomass_traits = data_biomass,
			pred_vect_nam = pred_vect_nam,
			pred_vect_1930s = pred_vect_1930s,
			formula_psi = formula_psi
		)

	### summary
	###########
	
		if (length(dharma_quant_test_biomass$qgamFits) > 0) {
			
			quant_test_upper <- summary(dharma_quant_test_biomass$qgamFits[[3]])$s.table[1, 4]
			quant_test_middle <- summary(dharma_quant_test_biomass$qgamFits[[2]])$s.table[1, 4]
			quant_test_lower <- summary(dharma_quant_test_biomass$qgamFits[[1]])$s.table[1, 4]

		} else {
			quant_test_lower <- quant_test_middle <- quant_test_upper <- NA_real_
		}

		resid_p_values <- c(
			moran_occs$p.value,
			dharma_resid_test_biomass$uniformity$p.value,
			dharma_resid_test_biomass$dispersion$p.value,
			dharma_resid_test_biomass$outliers$p.value,
			dharma_quant_test_biomass$p.value,
			quant_test_upper,
			quant_test_middle,
			quant_test_lower
		)

		if (length(dharma_quant_test_biomass$qgamFits) > 0) {
			
			quant_test_upper <- summary(dharma_quant_test_biomass$qgamFits[[3]])$s.table[1, 4]
			quant_test_middle <- summary(dharma_quant_test_biomass$qgamFits[[2]])$s.table[1, 3]
			quant_test_lower <- summary(dharma_quant_test_biomass$qgamFits[[1]])$s.table[1, 3]

		} else {
			quant_test_lower <- quant_test_middle <- quant_test_upper <- NA_real_
		}

		resid_p_values_sig <- ifelse(resid_p_values < 0.05, '*', 'ns')
		resid_test_statistic_values <- c(
			moran_occs$statistic,
			dharma_resid_test_biomass$uniformity$statistic,
			dharma_resid_test_biomass$dispersion$statistic,
			dharma_resid_test_biomass$outliers$statistic,
			NA,
			quant_test_upper,
			quant_test_middle,
			quant_test_lower
		)

		dharma_resids <- data.table(
			test = c('spatial autocorrelation', 'uniformity', 'dispersion', 'outliers', 'quantiles, overall', 'quantiles, upper', 'quantiles, middle', 'quantiles, lower'),
			p_value = resid_p_values,
			significant = resid_p_values_sig,
			test_statistic = c('Moran\'s I', names(dharma_resid_test_biomass$uniformity$statistic), names(dharma_resid_test_biomass$dispersion$statistic), 'exact binomial', NA, rep('chi squared', 3)),
			test_statistic_value = resid_test_statistic_values
		)

		# compile residuals analysis
		meta_biomass <- list(
			facet = 'biomass',
			date = date(),
			zero_inflated = zero_inflated,
			formulae = list(
				formula_biomass = formula_biomass,
				formula_psi = formula_psi
			),
			dharma_resids = dharma_resids
		)

		saveRDS(meta_biomass, paste0(out_dir, '/!meta_biomass.rds'))
		sink(paste0(out_dir, '/!meta_biomass.txt'), split = TRUE)
			print(meta_biomass)
		sink()

}

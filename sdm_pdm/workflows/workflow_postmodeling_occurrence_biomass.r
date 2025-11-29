#' Post-modeling workflow for occurrence-only models and models with an occurrence component.
#'
#' source('C:/Kaji/R/andropogon_integratedEcology/sdm_pdm/workflows/workflow_postmodeling_occurrence_biomass.r')
#'
#' @param homoscedastic If `TRUE`, then do not analyze behavior of sigma
#' @param zero_inflated Logical.
#' @param formula_... Formulae for occurrences/biomass mean, sigma, and probability of zero-inflation, and (for occurrences) bias.
#' @param log_precip Logical. If `TRUE`, log bios 12-14 and 16-19
#' @param out_dir Folder into which to save results.
workflow_postmodeling_occurrence_biomass <- function(
	formula_occs,
	formula_occs_sigma,
	formula_occs_bias,
	formula_biomass_mu,
	formula_biomass_sigma,
	formula_pzero,
	log_precip,
	out_dir
) {

	homoscedastic_occs <- is.null(formula_occs_sigma)
	homoscedastic_biomass <- is.null(formula_biomass_sigma)
	zero_inflated <- !is.null(formula_pzero)

	if (zero_inflated) data_traits <- prepare_nonbiomass_traits(trait = 'height', formula = ~ 1, n_response_curve_values = n_response_curve_values, calib = calib)

	### burn predictions into vector
	################################
	say('OCCURRENCE + BIOMASS: burn prediction vectors', level = 2)

		pred_vect_nam <- burn_occs_biomass_into_vector(
			demesne = 'nam',
			chains = chains,
			formula_occs = formula_occs,
			formula_occs_sigma = formula_occs_sigma,
			formula_biomass_mu = formula_biomass_mu,
			formula_biomass_sigma = formula_biomass_sigma,
			formula_pzero = formula_pzero,
			log_precip = log_precip
		)

		pred_vect_1930s <- burn_occs_biomass_into_vector(
			demesne = '1930s',
			chains = chains,
			formula_occs = formula_occs,
			formula_occs_sigma = formula_occs_sigma,
			formula_biomass_mu = formula_biomass_mu,
			formula_biomass_sigma = formula_biomass_sigma,
			formula_pzero = formula_pzero,
			log_precip = log_precip
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
			data = data_occs,
			log_precip = log_precip
		)

		if (!homoscedastic_occs) {

			responses_sigma <- graph_response_curves_occurrence_vs_environment(
				out_dir = out_dir,
				resp_type = 'sigma',
				chains = chains,
				data = data_occs,
				log_precip = log_precip
			)

		}

		if (zero_inflated) {

			responses_pzero <- graph_response_curves_occurrence_vs_environment(
				out_dir = out_dir,
				resp_type = 'pzero',
				chains = chains,
				data = data_occs_pzero,
				log_precip = log_precip
			)

		}

		if (data_occs$n_covariates_occs_bias >= 1) {

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
		covariates <- data_occs$covariates_occs
		presences <- data_occs$ag_vect_sq$n_andropogon_gerardi > 0
		x <- data_occs$ag_vect_sq[presences, covariates]
		x <- as.data.table(x)
		x <- scale(x, center = data_occs$x_centers_occs, scale = data_occs$x_scales_occs)
		means <- colMeans(x)
		x <- cbind(matrix(rep(means, mult * n_response_curve_values), byrow = TRUE, ncol = 3))
		colnames(x) <- covariates
		x_seq <- seq(x_min, x_max, length.out = mult * n_response_curve_values)
		x_seq <- (x_seq - data_occs$x_centers_occs['bio12']) / data_occs$x_scales_occs['bio12']
		x[ , 'bio12'] <- x_seq
		x <- as.data.table(x)

		x_occs <- model.matrix(formula_occs, x)
		x_biomass <- model.matrix(formula_biomass_mu, x)

		preds <- predict_occs_biomass(chains, x_occs = x_occs, x_biomass = x_biomass, homoscedastic_occs = homoscedastic_occs, homoscedastic_biomass = homoscedastic_biomass, zero_inflated = zero_inflated, type = 'mu', x_pzero = NULL)

		preds_occs <- colMeans(preds$preds_occs)
		preds_biomass <- colMeans(preds$preds_biomass)

		# contour-style joint response of occurrence vs biomass
		x_unscaled <- x_seq * data_occs$x_scales_occs['bio12'] + data_occs$x_centers_occs['bio12']
		# x_unscaled$bio12 <- 10^x_unscaled + 1
		df_joint <- data.table(
			occurrence = preds_occs,
			biomass = preds_biomass,
			bio12 = x_unscaled
		)

		df_joint$bio12 <- 10^(df_joint$bio12) - 1

		# hexagonal heatmap of predicted occurrence vs biomass,
		# filled by the mean bio12 within each hexagon
		facet_vs_facet <- ggplot(df_joint, aes(x = occurrence, y = biomass)) +
			stat_summary_hex(aes(z = bio12), fun = mean, bins = 30, color = NA) +
			scale_fill_gradient(name = 'BIO 12\n(mean)', low = 'yellow', high = 'blue') +
			xlab('Predicted abundance') +
			ylab('Predicted biomass (g)') +
			ggtitle('Occurrence vs biomass')

		ggsave(filename = paste0(out_dir, '/facet_vs_facet_response_curve_occurrence_vs_biomass.png'), plot = facet_vs_facet, width = 8, height = 8, dpi = 600, bg = 'white')

	### OCCURRENCE: dharma residuals
	################################
	say('OCCURRENCE: dharma residuals', level = 2)

		sims <- hammer_subset(chains, param = 'y_n_ag_sim', j = TRUE, na.rm = TRUE)
		sims <- hammer_rbind(sims)
		sims <- sims[ , data_occs$ag_vect_sq$focal_region]
		sims <- t(sims)

		if (!anyNA(sims)) {

			fits <- hammer_extract(chains, param = 'lambda_mu_sq', j = TRUE, stat = 'mean')
			fits <- fits[data_occs$ag_vect_sq$focal_region]

			observed_y <- data_occs$y_n_ag[data_occs$ag_vect_sq$focal_region]

			dharma_occs <- createDHARMa(simulatedResponse = sims, observedResponse = observed_y, fittedPredictedResponse = fits, integerResponse = TRUE)

			dharma_quant_test <- testQuantiles(dharma_occs, plot = FALSE)
			dharma_resid_test <- testResiduals(dharma_occs, plot = FALSE) # uniformity, dispersion, outlier
			dev.off()

			file <- paste0(out_dir, '/dharma_n_ag.png')
			png(file, width = 1200, height = 800)
				plot(dharma_occs)
			dev.off()

			file <- paste0(out_dir, '/dharma_lambda_residuals_abundance_residuals.png')
			png(file, width = 1200, height = 800)
				hist(dharma_occs$scaledResiduals, main = 'DHARMa residuals for number of observed AG (y_n_ag)', xlab = 'Scaled residuals', breaks = 30)
			dev.off()

		}

	### OCCURRENCE: spatial autocorrelation
	#######################################

	if (!anyNA(sims)) {

		coords <- as.data.frame(crds(centroids(project(pred_vect_nam[pred_vect_nam$focal_region], enmSdmX::getCRS('WGS84')))))

		# Compute Moran's I
		moran_occs <- moran.test(dharma_occs$scaledResiduals, nb2listw(knn2nb(knearneigh(coords, longlat = TRUE, k = 4))))

		pred_vect_nam$residual <- NA_real_
		pred_vect_nam$residual[pred_vect_nam$focal_region] <- dharma_occs$scaledResiduals

		extent <- ext(pred_vect_nam[pred_vect_nam$focal_region])
		extent <- as.vector(extent)

		map <- ggplot() +
			layer_spatial(pred_vect_nam, fill = 'gray50', color = 'gray40') +
			layer_spatial(pred_vect_nam[pred_vect_nam$focal_region], aes(fill = residual), color = NA) +
			scale_fill_distiller(palette = 'RdBu', limits = c(0, 1), na.value = 'grey80') +
			xlim(extent[1], extent[2]) + ylim(extent[3], extent[4]) +
			ggtitle(bquote('Residuals for ' * italic('Andropogon gerardi') * ' occurrence')) +
			theme(
				plot.title = element_text(size = 16),
				plot.subtitle = element_text(size = 14)
			)

		ggsave(plot = map, filename = paste0(out_dir, '/map_residuals_occurrences.png'), width = 12, height = 9, dpi = 600)

	}

	### OCCURRENCE: current map
	###########################
	say('OCCURRENCE: current map', level = 2)

		form <- paste(as.character(formula_occs), collapse = ' ')
		form <- gsub(form, pattern = 'I\\(', replacement = '')
		form <- gsub(form, pattern = '\\^2\\)', replacement = '²')
		form <- gsub(form, pattern = '*)', replacement = '×')

		form_bias <- paste(as.character(formula_occs_bias), collapse = ' ')
		form_bias <- gsub(form_bias, pattern = 'I\\(', replacement = '')
		form_bias <- gsub(form_bias, pattern = '\\^2\\)', replacement = '²')
		form_bias <- gsub(form_bias, pattern = '*)', replacement = '×')

		map <- map_occurrence(
			out_dir = out_dir,
			filename_append = 'present_day',
			pred_vect_nam = pred_vect_nam,
			response_var = 'N_ag_county_mean_sq',
			data_occs = data_occs,
			title = bquote('Present-day distribution of ' * italic('Andropogon gerardi') * ' abundance'),
			subtitle = paste0('1961-2020 | occ ', form, ' (bias ', form_bias, ')'),
			ag_core_quant = ag_core_quant
		)

		if (zero_inflated) {

			pzero_map <- map_pzero(
				out_dir = out_dir,
				filename_append = 'present_day',
				pred_vect_nam = pred_vect_nam,
				facet = 'occs',
				response_var = 'pzero_occs_county_sq',
				data_traits = data_traits,
				title = 'Present-day distribution of probability of zero abundance',
				subtitle = paste0('1961-2020 | occ ', form, ' (bias ', form_bias, ')'),
				plot_range_core = TRUE,
				ag_core_quant = ag_core_quant
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

			form_bias <- paste(as.character(formula_occs_bias), collapse = ' ')
			form_bias <- gsub(form_bias, pattern = 'I\\(', replacement = '')
			form_bias <- gsub(form_bias, pattern = '\\^2\\)', replacement = '²')
			form_bias <- gsub(form_bias, pattern = '*)', replacement = '×')

			title <- bquote('Future distribution of ' * italic('Andropogon gerardi') * ' abundance')
			subtitle <- paste0('SSP ', substr(fut, 4, 6), ' ', substr(fut, 8, 11), '-', substr(fut, 13, 16), ' | occ ', form, ' (bias ', form_bias, ')')
			response_var <- paste0('N_ag_county_mean_', fut)

			maps_occs_fut[[length(maps_occs_fut) + 1]] <- map_occurrence(
				out_dir = out_dir,
				filename_append = fut,
				pred_vect_nam = pred_vect_nam,
				response_var = response_var,
				data_occs = data_occs,
				title = title,
				subtitle = subtitle,
				ag_core_quant = ag_core_quant
			)

			if (zero_inflated) {

				title <- bquote('Future distribution of ' * italic('Andropogon gerardi') * ' probability of zero abundance')

				pzero_map <- map_pzero(
					out_dir = out_dir,
					filename_append = fut,
					pred_vect_nam = pred_vect_nam,
					facet = 'occs',
					response_var = paste0('pzero_occs_county_', fut),
					data_traits = data_traits,
					title = title,
					subtitle = subtitle,
					plot_range_core = TRUE,
					ag_core_quant = ag_core_quant
				)

			}

		}

	### occurrence: future change
	#############################
	say('OCCURRENCE: future change', level = 2)

		for (fut in futs) {

			say(fut)

			response_var <- paste0('N_ag_county_mean_', fut)

			map <- map_occurrence_change(
				out_dir = out_dir,
				filename_append = fut,
				pred_vect_nam = pred_vect_nam,
				response_var = response_var,
				data_occs = data_occs,
				title = bquote('Change in ' * italic('Andropogon gerardi') * ' abundance'),
				subtitle = paste0('SSP ', substr(fut, 4, 6), ' ', substr(fut, 8, 11), '-', substr(fut, 13, 16)),
				ag_core_quant = 0.95
			)

			if (zero_inflated) {

				title <- bquote('Change in ' * italic('Andropogon gerardi') * ' probability of zero abundance')

				pzero_map <- map_pzero_change(
					out_dir = out_dir,
					filename_append = fut,
					fut = fut,
					pred_vect_nam = pred_vect_nam,
					facet = 'occs',
					response_var = paste0('pzero_occs_county_', fut),
					response_var_sq = 'pzero_occs_county_sq',
					data_traits = data_traits,
					title = title,
					subtitle = subtitle,
					plot_range_core = TRUE,
					ag_core_quant = 0.95
				)

			}

		}

	### OCCURRENCE: 1930s
	#####################
	say('OCCURRENCE: 1930s change maps', level = 2)

	map_1930s <- map_occurrence_change_1930s(zero_inflated = zero_inflated, pred_vect_nam = pred_vect_nam, pred_vect_1930s = pred_vect_1930s)

	### summary
	###########
	
	if (!anyNA(sims)) {
	
		resid_p_values <- c(
			moran_occs$estimate[1],
			dharma_resid_test$uniformity$p.value,
			dharma_resid_test$dispersion$p.value,
			dharma_resid_test$outliers$p.value,
			dharma_quant_test$p.value,
			summary(dharma_quant_test$qgamFits[[3]])$s.table[1, 4],
			summary(dharma_quant_test$qgamFits[[2]])$s.table[1, 4],
			summary(dharma_quant_test$qgamFits[[1]])$s.table[1, 4]
		)
		resid_p_values_sig <- ifelse(resid_p_values < 0.05, '*', 'ns')
		resid_test_statistic_values <- c(
			moran_occs$statistic,
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
	
	} else {
		dharma_resids <- NA
	}

	# compile residuals analysis
	meta_occs <- list(
		facet = 'occurrence',
		date = date(),
		homoscedastic_occs = homoscedastic_occs,
		homoscedastic_biomass = homoscedastic_biomass,
		zero_inflated = zero_inflated,
		formulae = list(
			formula_occs = formula_occs,
			formula_occs_bias = formula_occs_bias
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

	homoscedastic_biomass <- is.null(formula_biomass_sigma)

	### BIOMASS: plant-level residuals analysis
	###########################################
	say('BIOMASS: plant-level residuals analysis', level = 2)

	# NB this uses the mean predicted value of a site as a plant-level prediction
	sims_by_plant <- hammer_subset(chains, param = 'y_biomass_sim', j = TRUE)
	sims_by_plant <- hammer_rbind(sims_by_plant)
	sims_by_plant <- t(sims_by_plant)

	estimates <- hammer_subset(chains, param = 'mu_biomass_site', j = TRUE)
	estimates <- hammer_rbind(estimates)

	site_counts <- data_biomass_mu$raw_data_biomass[ , .N, by = SITE]
	fit_by_site <- apply(estimates, 2, median)
	fits <- numeric()
	for (i in seq_along(fit_by_site)) {
		fits <- c(fits, rep(fit_by_site[i], site_counts$N[i]))
	}

	dharma_biomass <- createDHARMa(simulatedResponse = sims_by_plant, observedResponse = data_biomass_mu$y_biomass, fittedPredictedResponse = fits, integerResponse = FALSE)

	dharma_quant_test <- testQuantiles(dharma_biomass, plot = FALSE)
	dharma_resid_test <- testResiduals(dharma_biomass, plot = FALSE) # uniformity, dispersion, outlier
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

		mean_est <- hammer_extract(chains, 'mu_biomass_site', j = TRUE)
		lower_est <- hammer_extract(chains, 'mu_biomass_site', j = TRUE, stat = 'lower')
		upper_est <- hammer_extract(chains, 'mu_biomass_site', j = TRUE, stat = 'upper')

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

		# Using the mean dharma_biomass residual of a site for site-level residual

		residuals_by_plant <- dharma_biomass$scaledResiduals
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

		ggsave(plot = map, filename = paste0(out_dir, '/map_residuals_biomass_site_mu.png'), width = 12, height = 9, dpi = 600)

	### BIOMASS: SAC in residuals
	#############################
	say('BIOMASS: SAC in residuals', level = 2)

		coords <- as.data.frame(crds(project(site_vect_biomass_resid, enmSdmX::getCRS('WGS84'))))

		# Compute Moran's I
		moran_biomass <- moran.test(residuals_by_site, nb2listw(knn2nb(knearneigh(coords, longlat = TRUE, k = 4))))

	### BIOMASS: plant-level residuals vs covariates
	################################################
	say('BIOMASS: plant-level residuals vs covariates', level = 2)

		resids <- dharma_biomass$scaledResiduals
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
				ggtitle('DHAMRa Residuals for Biomass Means')

		}

		# if (!homoscedastic_biomass) {
		
		# 	for (i in seq_len(data_biomass_mu$n_covariates)) {

		# 		if (data_biomass_mu$covariates_biomass[i] == 'ph') {
		# 			x <- data_biomass_mu$raw_data_biomass[['site_soil_pH']]
		# 		} else {
		# 			x <- data_biomass_mu$raw_data_biomass[[data_biomass_mu$covariates_biomass[i]]]
		# 		}

		# 		this_x <- data.frame(
		# 			x = x,
		# 			y = resids
		# 		)

		# 		resid_model <- mgcv::gam(logitAdj(y, epsilon = 0.0001) ~ s(x), data = this_x)

		# 		F <- summary(resid_model)$s.table[1, 3]
		# 		p <- summary(resid_model)$s.table[1, 4]
		# 		sig <- ifelse(p < 0.05, '*', '-')
				
		# 		resids_vs_covariates_results <- rbind(
		# 			resids_vs_covariates_results,
		# 			data.table(
		# 				covariate = data_biomass_mu$covariates_biomass[i],
		# 				F = F,
		# 				p = p,
		# 				significance = sig
		# 			)
		# 		)

		# 		resids_vs_covariates[[length(resids_vs_covariates) + 1]] <- ggplot(this_x, aes(x = x, y = y)) +
		# 			geom_point() +
		# 			geom_smooth(method = loess, formula = y ~ x, se = FALSE) +
		# 			xlab(data_biomass_mu$terms_biomass[i]) +
		# 			ylab('Residual value') +
		# 			ggtitle('dharma_biomass Residuals for Biomass Sigma')

		# 	}
		
		# }

		if (data_biomass_mu$n_covariates == 2 & homoscedastic_biomass) {
			ncol <- 2
		} else {
			ncol <- ceiling(sqrt(length(resids_vs_covariates)))
		}
		width <- 6 * ncol
		height <- 5 * ceiling(length(resids_vs_covariates) / ncol)

		resids_vs_covariates <- plot_grid(plotlist = resids_vs_covariates, ncol = ncol)
		ggsave(plot = resids_vs_covariates, filename = paste0(out_dir, '/dharma_residuals_vs_covariates_biomass.png'), width = width, height = height, dpi = 600)

	### BIOMASS: pseudo-R2 (correlation between observed and estimate)
	##################################################################
	say('BIOMASS: model fit', level = 2)

		observed <- data_biomass_mu$raw_data_biomass[ , .(mean_biomass = mean(Biomass)), by = SITE][['mean_biomass']]
		estimated <- hammer_extract(chains, 'mu_biomass_site', j = TRUE)
		correl_mu <- cor(observed, estimated)

		rmse_mu <- unname(rmse_fx(observed, estimated))
		mean_error_mu <- mean((estimated - observed) / observed)
		mean_abs_error_mu <- mean(abs(estimated - observed) / observed)

		if (!homoscedastic_biomass) {

			observed <- data_biomass_mu$raw_data_biomass[ , .(sd_biomass = sd(Biomass)), by = SITE][['sd_biomass']]
			estimated <- hammer_extract(chains, 'mu_biomass_site', j = TRUE)
			correl_sigma <- cor(observed, estimated)

			rmse_sigma <- unname(rmse_fx(observed, estimated))
			mean_error_sigma <- mean((estimated - observed) / observed)
			mean_abs_error_sigma <- mean(abs(estimated - observed) / observed)

		}
	
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

		if (!homoscedastic_biomass) {
			responses_sigma <- graph_response_curves_biomass_vs_environment(
				out_dir = out_dir,
				chains = chains,
				resp_type = 'sigma',
				data_biomass = data_biomass_mu,
				quant_threshold = 0.5
			)

		}

		if (zero_inflated) {

			responses_pzero <- graph_response_curves_biomass_vs_environment(
				out_dir = out_dir,
				resp_type = 'pzero',
				chains = chains,
				data_biomass_mu = data_biomass_mu
			)

		}


	### BIOMASS: current map
	########################
	say('BIOMASS: current map', level = 2)

		mean_form <- paste(as.character(formula_biomass_mu), collapse = ' ')
		mean_form <- gsub(mean_form, pattern = 'I\\(', replacement = '')
		mean_form <- gsub(mean_form, pattern = '\\^2\\)', replacement = '²')
		mean_form <- gsub(mean_form, pattern = '*)', replacement = '×')

		if (homoscedastic_biomass) {
			subtitle <- paste0('Biomass μ ', mean_form, ' | 1991-2020')
		} else {
			sd_form <- paste(as.character(formula_biomass_sigma), collapse = ' ')
			subtitle <- paste0('Biomass μ ', mean_form, ' & biomass σ ', sd_form, ' | 1991-2020')
		}

		map <- map_biomass(
			out_dir = out_dir,
			filename_append = 'present',
			pred_vect_nam = pred_vect_nam,
			response_var = 'mu_biomass_county_mean_sq',
			response_var_type = 'mu',
			data_occs = data_occs,
			data_biomass = data_biomass_mu,
			title = bquote('Present-day distribution of mean ' * italic('Andropogon gerardi') * ' biomass '),
			subtitle = subtitle,
			legend_title = 'Biomass (g)',
			plot_range_core = TRUE,
			ag_core_quant = 0.95
		)

		if (!is.null(formula_biomass_sigma)) {
			
			facet <- 'biomass'
			map_sigma <- map_biomass(
				out_dir = out_dir,
				filename_append = 'present',
				pred_vect_nam = pred_vect_nam,
				response_var = 'sigma_biomass_county_median_sq',
				response_var_type = 'sigma',
				data_occs = data_occs,
				data_biomass = data_biomass_mu,
				title = bquote('Present-day biomass among-site σ'),
				subtitle = subtitle,
				legend_title = 'Biomass σ (g)',
				plot_range_core = TRUE,
				ag_core_quant = 0.95
			)

		}

		if (!is.null(formula_pzero)) {
			
			facet <- 'biomass'
			map <- map_pzero(
				out_dir = out_dir,
				filename_append = 'present',
				pred_vect_nam = pred_vect_nam,
				facet = facet,
				response_var = paste0('pzero_', facet, '_county_sq'),
				data_traits = data_biomass_mu,
				title = bquote('Present-day probability of zero biomass '),
				subtitle = subtitle,
				plot_range_core = TRUE,
				ag_core_quant = 0.95
			)

		}

	### BIOMASS: future maps
	########################
	say('BIOMASS: future maps', level = 2)

		for (fut in futs) {

			say(fut)

			response_var <- paste0('mu_biomass_county_mean_', fut)

			mean_form <- paste(as.character(formula_biomass_mu), collapse = ' ')
			if (homoscedastic_biomass) {
				subtitle <- paste0('Biomass μ ', mean_form, '  | SSP ', substr(fut, 4, 6), ' ', substr(fut, 8, 11), '-', substr(fut, 13, 16))
			} else {

				sd_form <- paste(as.character(formula_biomass_sigma), collapse = ' ')
				subtitle <- paste0('Biomass μ ', mean_form, ' & biomass σ ', sd_form, ' | SSP ', substr(fut, 4, 6), ' ', substr(fut, 8, 11), '-', substr(fut, 13, 16))
			}

			map_biomass_fut <- map_biomass(
				out_dir = out_dir,
				filename_append = fut,
				pred_vect_nam = pred_vect_nam,
				response_var = response_var,
				response_var_type = 'mu',
				data_occs = data_occs,
				data_biomass = data_biomass_mu,
				title = bquote('Future distribution of mean ' * italic('Andropogon gerardi') * ' biomass '),
				subtitle = subtitle,
				legend_title = 'Biomass (g)',
				plot_range_core = TRUE,
				ag_core_quant = 0.95
			)

			if (!is.null(formula_pzero)) {

				face <- 'biomass'				
				map_fut_pzero <- map_pzero(
					out_dir = out_dir,
					filename_append = fut,
					pred_vect_nam = pred_vect_nam,
					facet = facet,
					response_var = paste0('pzero_', facet, '_county_', fut),
					data_traits = data_biomass_mu,
					title = bquote('Future probability of zero biomass '),
					subtitle = subtitle,
					plot_range_core = TRUE,
					ag_core_quant = 0.95
				)

			}

		}

	### BIOMASS: future change maps
	###############################
	say('BIOMASS: future change maps', level = 2)

		maps_biomass_change <- list()
		for (fut in futs) {

			say(fut)

			response_var <- paste0('mu_biomass_county_mean_', fut)
			title <- bquote('Change in ' * italic('Andropogon gerardi') * ' biomass ')

			mean_form <- paste(as.character(formula_biomass_mu), collapse = ' ')
			if (homoscedastic_biomass) {
				subtitle <- paste0('Biomass μ ', mean_form, '  | SSP ', substr(fut, 4, 6), ' ', substr(fut, 8, 11), '-', substr(fut, 13, 16))
			} else {

				sd_form <- paste(as.character(formula_biomass_sigma), collapse = ' ')
				subtitle <- paste0('Biomass μ ', mean_form, ' & biomass σ ', sd_form, ' | SSP ', substr(fut, 4, 6), ' ', substr(fut, 8, 11), '-', substr(fut, 13, 16))
			}

			maps_biomass_change[[length(maps_biomass_change) + 1]] <- map_biomass_change(
				out_dir = out_dir,
				filename_append = fut,
				fut = fut,
				pred_vect_nam = pred_vect_nam,
				response_var = response_var,
				data_occs = data_occs,
				data_biomass = data_biomass_mu,
				title = title,
				subtitle = subtitle,
				plot_range_core = TRUE,
				ag_core_quant = ag_core_quant
			)

			if (!is.null(formula_pzero)) {

				title <- bquote('Change in probability of zero biomass ')

				facet <- 'biomass'				
				map_pzero_change <- map_pzero_change(
					out_dir = out_dir,
					filename_append = fut,
					fut = fut,
					pred_vect_nam = pred_vect_nam,
					facet = facet,
					response_var = paste0('pzero_', facet, '_county_', fut),
					response_var_sq = paste0('pzero_', facet, '_county_sq'),
					data_traits = data_biomass_mu,
					title = title,
					subtitle = subtitle,
					plot_range_core = TRUE,
					ag_core_quant = 0.95
				)

			}

		}

	### BIOMASS: 1930s change maps
	##############################
	say('BIOMASS: 1930s change maps', level = 2)

	map_thirties <- map_biomass_traits_change_1930s(facet = 'biomass', data_biomass_traits = data_biomass_mu, pred_vect_nam = pred_vect_nam, pred_vect_1930s = pred_vect_1930s, formula_mu = formula_biomass_mu, formula_biomass_sigma = formula_biomass_sigma, formula_pzero = formula_pzero)

}

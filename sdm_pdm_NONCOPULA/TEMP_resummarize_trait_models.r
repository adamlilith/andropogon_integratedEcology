# source('C:/Kaji/R/andropogon_integratedEcology/sdm_pdm/!!!resummarize_trait_models.r')

	rm(list = ls())

	drive <- 'C:/Kaji/'

	setwd(paste0(drive, '/Research/Andropogon/Andropogon'))
	source(paste0(drive, '/R/andropogon_integratedEcology/sdm_pdm/sdm_pdm_00_shared_functions_and_variables.r'))

	traits <- c(
		'height',
		'blade_width',
		'leaf_thickness',
		'spad',
		'canopy_diameter',
		'photosynthetic_rate',
		'stomatal_conductance',
		'internal_co2',
		'transpiration_rate',
		'n_concentration',
		'cn_ratio'
	)

homoscedastic <- TRUE

for (trait in traits) {

	folders <- listFiles(paste0('C:/Kaji/Research/Andropogon/Andropogon/outputs_loretta/integrated_sdm_pdm/models_', trait), pattern = 'gamma\\~normal_homoscedastic')

	for (folder in folders) {
	
		say(folder, level = 1)
		chains <- readRDS(paste0(folder, '/chains.rds'))
		meta_old <- readRDS(paste0(folder, '/!meta_nonbiomass_trait.rds'))
		formulae <- readRDS(paste0(folder, '/formula.rds'))

		formula_traits <- formulae$formula_trait

		data_traits <- prepare_nonbiomass_facets(trait = trait, formula = formula_traits, n_response_curve_values = n_response_curve_values, calib = FALSE)


		### NON-BIOMASS TRAIT: plant-level residuals analysis
		#####################################################
		say('NON-BIOMASS TRAIT: plant-level residuals analysis', level = 2)

		# NB this uses the mean predicted value of a site as a plant-level prediction
		sims_by_plant <- mc_subset(chains, param = 'y_traits_sim', j = TRUE)
		sims_by_plant <- mc_rbind(sims_by_plant)
		sims_by_plant <- t(sims_by_plant)

		estimates <- mc_subset(chains, param = 'site_traits_mu', j = TRUE)
		estimates <- mc_rbind(estimates)

		site_counts <- data_traits$raw_data_traits[ , .N, by = SITE]
		fit_by_site <- apply(estimates, 2, median)
		fits <- numeric()
		for (i in seq_along(fit_by_site)) {
			fits <- c(fits, rep(fit_by_site[i], site_counts$N[i]))
		}

		dharma <- createDHARMa(simulatedResponse = sims_by_plant, observedResponse = data_traits$y_traits[ , 1], fittedPredictedResponse = fits, integerResponse = FALSE)

		dharma_quant_test <- testQuantiles(dharma, plot = FALSE)
		dharma_resid_test <- testResiduals(dharma, plot = FALSE) # uniformity, dispersion, outlier
		dev.off()

		# file <- paste0(out_dir, '/dharma_by_plant_biomass.png')
		# png(file, width = 1200, height = 800)
		# 	plot(dharma)
		# dev.off()
		
		# file <- paste0(out_dir, '/dharma_by_plant_traits_residuals.png')
		# png(file, width = 1400, height = 1000)
		# 	hist(dharma$scaledResiduals, main = 'DHARMa residuals for plant', xlab = 'Scaled residuals', breaks = 30)
		# dev.off()

		### NON-BIOMASS TRAIT: map of residuals
		#######################################
		say('NON-BIOMASS TRAIT: map of residuals', level = 2)

		# Using the mean DHARMa residual of a site for site-level residual

		residuals_by_plant <- dharma$scaledResiduals
		residuals_by_site <- rep(NA, data_traits$n_pheno_sites)
		site_ids <- unique(data_traits$raw_data_traits$SITE)
		for (i in 1:data_traits$n_pheno_sites) residuals_by_site[i] <- mean(residuals_by_plant[data_traits$raw_data_traits$SITE == site_ids[i]])

		site_vect_traits_resid <- data_traits$site_vect_traits
		site_vect_traits_resid$residual <- residuals_by_site

		# nam <- vect('./data_from_gadm/gadm_4pt1_level_1_north_america_sans_alaska_lambert.gpkg')
		# nam <- project(nam, site_vect_traits_resid)

		# # extent
		# extent <- buffer(site_vect_traits_resid, width = 200 * 1000)
		# extent <- ext(extent)
		# extent <- as.vector(extent)

		# map <- ggplot() +
		# 	layer_spatial(nam, color = 'gray30', fill = 'white', linewidth = 0.3) +
		# 	layer_spatial(site_vect_traits_resid, aes(fill = residual), pch = 21, size = 6) +
		# 	scale_fill_gradient2(
		# 		name = 'Mean of\nDHARMa\nresiduals',
		# 		low = 'red',
		# 		mid = 'beige',
		# 		high = 'blue',
		# 		midpoint = 0.5
		# 	) +
		# 	xlim(extent[1], extent[2]) + ylim(extent[3], extent[4]) +
		# 	ggtitle(
		# 		bquote('Residuals for ' * italic('Andropogon gerardi') * ' ' * trait),
		# 		subtitle = paste0('1991-2020 | ', get_nice_trait(trait)$long)
		# 	) +
		# 	theme(
		# 		plot.title = element_text(size = 16),
		# 		plot.subtitle = element_text(size = 14)
		# 	)

		# ggsave(plot = map, filename = paste0(out_dir, '/map_residuals_', trait, '_site_mu.png'), width = 12, height = 9, dpi = 600)

		### NON-BIOMASS TRAIT: SAC in residuals
		#######################################
		say('NON-BIOMASS TRAIT: SAC in residuals', level = 2)
		coords <- as.data.frame(crds(project(site_vect_traits_resid, enmSdmX::getCRS('WGS84'))))

		# Compute Moran's I
		moran <- moran.test(residuals_by_site, nb2listw(knn2nb(knearneigh(coords, longlat = TRUE, k = 4))))

		# # Save Moran's I results
		# file.remove(paste0(out_dir, '/dharma_residuals_spatial_autocorrelation_', trait, '.txt'))

		### NON-BIOMASS TRAIT: plant-level residuals vs covariates
		##########################################################
		say('NON-BIOMASS TRAIT: plant-level residuals vs covariates', level = 2)

		resids <- dharma$scaledResiduals
		resids_vs_covariates <- list()
		resids_vs_covariates_results <- data.table()
		for (i in seq_len(data_traits$n_covariates)) {

			if (data_traits$covariates_traits[i] == 'ph') {
				x <- data_traits$raw_data_traits[['site_ph']]
			} else {
				x <- data_traits$raw_data_traits[[data_traits$covariates_traits[i]]]
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
					covariate = data_traits$covariates_traits[i],
					F = F,
					p = p,
					significance = sig
				)
			)

		}

		# 	resids_vs_covariates[[i]] <- ggplot(this_x, aes(x = x, y = y)) +
		# 		geom_point() +
		# 		geom_smooth(method = loess, formula = y ~ x, se = FALSE) +
		# 		xlab(data_traits$terms_traits[i]) +
		# 		ylab('Residual value') +
		# 		ggtitle('DHARMa Residuals for ', trait)

		# }

		# if (data_traits$n_covariates == 2) {
		# 	ncol <- 2
		# } else {
		# 	ncol <- ceiling(sqrt(length(data_traits$n_covariates)))
		# }
		# width <- 6 * ncol
		# height <- 5 * ceiling(length(data_traits$n_covariates) / ncol)

		# resids_vs_covariates <- plot_grid(plotlist = resids_vs_covariates, ncol = ncol)
		# ggsave(plot = resids_vs_covariates, filename = paste0(out_dir, '/dharma_residuals_vs_covariates_biomass.png'), width = width, height = height, dpi = 600)

		### NON-BIOMASS TRAIT: pseudo-R2 (correlation between observed and estimate)
		############################################################################
		say('NON-BIOMASS TRAIT: pseudo-R2 (correlation between observed and estimate)', level = 2)

		# sink(paste0(out_dir, '/model_fit_', trait, '.txt'), split = TRUE)
		say('MODEL FIT', post = 2)

		observed <- rep(NA_real_, data_traits$n_pheno_sites)
		sites <- unique(data_traits$raw_data_traits$SITE)
		for (i in seq_along(sites)) {
		
			site <- sites[i]
			trait_column_name <- get_raw_trait_name_from_rfriendly(trait)
			observed[i] <- mean(data_traits$raw_data_traits[[trait_column_name]][data_traits$raw_data_traits$SITE == site])

		}
		estimated <- mc_extract(chains, 'site_traits_mu', j = TRUE)
		estimated <- unname(estimated)
		correl_mu <- cor(observed, estimated)
		say('Correlation between site-level mean observed and estimated mean trait value: ', correl_mu)

		rmse_mu <- rmse_fx(observed, estimated)
		say('RMSE between site-level mean observed and estimated mean trait value: ', rmse_mu)
		
		mean_error_mu <- mean((estimated - observed) / observed)
		mean_abs_error_mu <- mean(abs(estimated - observed) / observed)
		say('Mean error between estimated mean trait value and site-level mean observed trait value: ', mean_error_mu)
		say('Mean absolute error between estimated mean trait value and site-level mean observed trait value: ', mean_abs_error_mu)
		
		if (!homoscedastic) {
		
			observed <- data_traits$raw_data_biomass[ , .(sd_biomass = mean(Biomass)), by = SITE][['sd_biomass']]
			estimated <- mc_extract(chains, 'site_traits_sigma', j = TRUE)
			correl_sigma <- cor(observed, estimated)
			say('Correlation between site-level mean observed and estimated s.d. of taut values: ', correl)

			rmse_sigma <- rmse_fx(observed, estimated)
			say('RMSE between site-level mean observed and estimated mean trait value: ', rmse_sigma)
			
			mean_error_sigma <- mean((estimated - observed) / observed)
			mean_abs_error_sigma <- mean(abs(estimated - observed) / observed)
			say('Mean error between estimated trait value s.d. and site-level observed trait value s.d.: ', mean_error_sigma)
			say('Mean absolute error between estimated trait value s.d. and site-level observed trait value s.d.: ', mean_abs_error_sigma)

		
		}
	
		# sink()

	# ### NON-BIOMASS TRAIT: response curves
	# ######################################
	# say('NON-BIOMASS TRAIT: response curves', level = 2)

	# 	responses <- graph_response_curves_nonbiomass_traits_vs_environment(
	# 		trait = trait,
	# 		quant_threshold = 0.5
	# 	)


	# ### NON-BIOMASS TRAIT: cross validation
	# #######################################
	# say('NON-BIOMASS TRAIT: cross validation', level = 2)

	# 	if (!trial) {

	# 		loss_fx <- mean_abs_percent_error_fx

	# 		say('Cross-validation loss function:')
	# 		print(loss_fx)
	# 		say('')

	# 		folds_fx <- folds_for_single_nonbiomass_trait # change according to response data type we're using

	# 		cv <- runCrossValidate(
	# 			MCMCconfiguration = conf,
	# 			k = k_folds, # universal setting
	# 			foldFunction = folds_fx,
	# 			lossFunction = loss_fx,
	# 			MCMCcontrol = list(niter = niter, nburnin = nburnin),
	# 			returnSamples = FALSE,
	# 			nCores = 1,
	# 			nBootReps = 200,
	# 			silent = FALSE
	# 		)

	# 		mean_form <- paste(as.character(formula_traits), collapse = ' ')
	# 		mean_form <- gsub(mean_form, pattern = 'I\\(', replacement = '')
	# 		mean_form <- gsub(mean_form, pattern = '\\^2\\)', replacement = '²')
	# 		mean_form <- gsub(mean_form, pattern = '*)', replacement = '×')
	# 		form <- paste0(trait, ' ', mean_form)

	# 		cv_eval <- data.table(
	# 			model = form,
	# 			k = 'summary',
	# 			cv_value = cv$CVvalue,
	# 			cv_value_se = cv$CVstandardError
	# 		)

	# 		for (k in 1:k_folds) {
			
	# 			cv_eval <- rbind(
	# 				cv_eval,
	# 				data.table(
	# 					model = form,
	# 					k = k,
	# 					cv_value = cv$foldCVinfo[[k]][1],
	# 					cv_value_se = cv$foldCVinfo[[k]][2]
	# 				)
	# 			)

	# 		}

	# 		# # write.csv(cv_eval, paste0(out_dir, '/crossvalidation_', trait, '.csv'), row.names = FALSE)

	# 	}

	### summaries
	#############

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

		descrip <- gsub(meta_old$descrip, pattern = ' lognormal', replacement = '~dnorm')

		meta <- list(
			facet = trait,
			descrip = descrip,
			formulae = list(formula_traits = formula_traits),
			date = date(),
			dharma_resids = data.table(
				test = c('spatial autocorrelation', 'uniformity', 'dispersion', 'outliers', 'quantiles, overall', 'quantiles, upper', 'quantiles, middle', 'quantiles, lower'),
				p_value = resid_p_values,
				significant = resid_p_values_sig,
				test_statistic = c('Moran\'s I', names(dharma_resid_test$uniformity$statistic), names(dharma_resid_test$dispersion$statistic), 'exact binomial', NA, rep('chi squared', 3)),
				test_statistic_value = resid_test_statistic_values
			),
			dharma_resids_vs_covariates = resids_vs_covariates_results,
			fit_mu = c(correl_mu = correl_mu, rmse_mu = rmse_mu, mean_error_mu = mean_error_mu, mean_abs_error_mu = mean_abs_error_mu)
		)

		if (homoscedastic) {
			meta$fit_sd <- NA
		} else {
			meta$fit_sigma = c(correl_sigma = correl_sigma, rmse_sigma = rmse_sigma, mean_error_sigma = mean_error_sigma, mean_abs_error_sigma = mean_abs_error_sigma)
		}

		meta$crossvalidation <- meta_old$crossvalidation

		saveRDS(meta, paste0(folder, '/!meta_nonbiomass_trait.rds'))
		sink(paste0(folder, '/!meta_nonbiomass_trait.txt'), split = TRUE)
			print(meta)
		sink()

	} # next folder


}
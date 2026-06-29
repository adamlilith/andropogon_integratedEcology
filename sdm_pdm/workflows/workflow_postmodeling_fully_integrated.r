#' Post-modeling workflow for occurrence-only models and models with an occurrence component.
#'
#' source('C:/Kaji/R/andropogon_integratedEcology/sdm_pdm/workflows/workflow_postmodeling_occurrence_biomass.r')
#'
#' formula_occs			Formula for occurrences
#' formula_bias			Formula for bias in sampling occurrences
#' formula_psi			Formula for presence/absence
#' formula_biomass		Formula for biomass
#' nonbiomass_facets	Named list of non-biomass facets
#' 
#' out_dir 				Folder into which to save results.
workflow_postmodeling_fully_integrated <- function(
	formula_occs,
	formula_bias,
	formula_psi,
	formula_biomass,
	nonbiomass_facets,
	out_dir
) {
	
	### collate data
	################

		# data for OCCURRENCES at counties
		say('preparing data for occurrences')
		data_occs_counties <- prepare_occurrence_data(formula_occs = formula_occs, formula_bias = ~ 1, n_response_curve_values = n_response_curve_values, psa_quant = psa_quant, calib = calib)

		# data for OCCURRENCES using site-level environment
		data_occs_sites <- prepare_biomass_data(formula_biomass = formula_occs, n_response_curve_values = n_response_curve_values, calib = calib)

		say('preparing data for psi')

		# data for ZERO-INFLATION at counties
		data_psi_counties <- prepare_occurrence_data(formula_occs = formula_psi, formula_bias = ~ 1, n_response_curve_values = n_response_curve_values, psa_quant = psa_quant, calib = calib)

		# data for psi at sites
		data_psi_sites <- prepare_biomass_data(formula_biomass = formula_psi, n_response_curve_values = n_response_curve_values, calib = calib)

		say('preparing data for biomass')
		# data for BIOMASS at sites
		data_biomass_sites <- prepare_biomass_data(formula_biomass = formula_biomass, n_response_curve_values = n_response_curve_values, calib = calib)

		# data for BIOMASS using county-level environment
		data_biomass_counties <- prepare_occurrence_data(formula_occs = formula_biomass, formula_bias = ~ 1, n_response_curve_values = n_response_curve_values, psa_quant = psa_quant, calib = calib)

		# data for NON-BIOMASS FACETS at sites and counties
		data_nonbiomass_sites <- data_nonbiomass_counties <- list()
		for (f in seq_along(nonbiomass_facets)) {

			facet <- names(nonbiomass_facets)[f]
			say('preparing data for ', facet)
			
			formula_facet <- nonbiomass_facets[[facet]]$formula

			data_nonbiomass_sites[[f]] <- prepare_nonbiomass_data(facet = facet, formula_facet = formula_facet, n_response_curve_values = n_response_curve_values, calib = calib)

			data_nonbiomass_counties[[f]] <- prepare_occurrence_data(formula_occs = formula_facet, formula_bias = ~ 1, n_response_curve_values = n_response_curve_values, psa_quant = psa_quant, calib = calib)

		}

		names(data_nonbiomass_sites) <- names(data_nonbiomass_counties) <- names(nonbiomass_facets)

say('ASSUMING WE ARE USING CHAINS WITH DIFFERENT # OF ITERATIONS.\nREFORMATING INTO A SINGLE CHAIN.\nREMOVE THIS WHEN ALL ITERATIONS ARE DONE!!!', level = 1)
n_chains <- length(chains$samples)
n_iter_1 <- nrow(chains$samples[[1]])
diff_iter <- FALSE
if (n_chains > 1) {
	for (i in 2:n_chains) {
		if (n_iter_1 != nrow(chains$samples[[i]])) diff_iter <- TRUE
	}
}

if (diff_iter & n_chains > 1) {
	for (i in 2:n_chains) {
		chains$samples[[1]] <- rbind(chains$samples[[1]], chains$samples[[i]])
	}
	chains <- list(samples = as.mcmc(chains$samples[[1]]))
	chains$samples <- as.mcmc.list(chains$samples, start = 1, send = nrow(chains$samples[[1]]), thin = 1)
}


	resp_distrib_biomass <- formulae$meta_biomass$resp_distrib_biomass
	transform_biomass <- formulae$meta_biomass$transform_biomass
	log_precip_biomass <- formulae$meta_biomass$log_precip_biomass

	### burn predictions into vector
	################################
	say('OCCURRENCE + BIOMASS + NON-BIOMASS FACETS: burn prediction vectors', level = 2)

		pred_vect_nam <- burn_fully_integrated_into_vector(
			demesne = 'nam',
			chains = chains,
			formula_occs = formula_occs,
			formula_bias = formula_bias,
			formula_psi = formula_psi,
			formula_biomass = formula_biomass,
			resp_distrib_biomass = resp_distrib_biomass,
			transform_biomass = transform_biomass,
			log_precip_biomass = log_precip_biomass,
			nonbiomass_facets = nonbiomass_facets
		)

		pred_vect_1930s <- burn_fully_integrated_into_vector(
			demesne = '1930s',
			chains = chains,
			formula_occs = formula_occs,
			formula_bias = formula_bias,
			formula_psi = formula_psi,
			formula_biomass = formula_biomass,
			nonbiomass_facets = nonbiomass_facets
		)

		writeVector(pred_vect_nam, paste0(out_dir, '/prediction_vector_nam.gpkg'), overwrite = TRUE)
		writeVector(pred_vect_1930s, paste0(out_dir, '/prediction_vector_conus_1930s.gpkg'), overwrite = TRUE)

	### response curves
	###################
	say('RESPONSE CURVES', level = 2)

		graph_response_curves_fully_integrated(
			out_dir = out_dir,
			chains = chains,
			resp_distrib_biomass = resp_distrib_biomass,
			transform_biomass = transform_biomass,
			nonbiomass_facets = nonbiomass_facets
		)

		if (data_occs_counties$n_covariates_bias >= 1) {

			responses_bias <- graph_response_curves_occurrence_bias_vs_bias_covariates(
				out_dir = out_dir,
				chains = chains
			)

		}

	### DHARMA RESIDUALS: occurrences
	#################################
	say('DHARMA RESIDUALS: occurrences', level = 2)

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

		file <- paste0(out_dir, '/dharma_qq_quantiles_abundance.png')
		png(file, width = 1200, height = 800)
			plot(dharma_occs)
		dev.off()

		file <- paste0(out_dir, '/dharma_distribution_abundance_residuals.png')
		png(file, width = 1200, height = 800)
			hist(dharma_occs$scaledResiduals, main = 'DHARMa residuals for number of observed AG (y_n_ag)', xlab = 'Scaled residuals', breaks = 30)
		dev.off()

	### BIOMASS GENERAL DHARMa
	##########################

		say('   biomass general DHARMA...')

		sims <- mc_subset(chains, 'y_biomass_sim', j = TRUE)
		sims <- mc_rbind(sims)
		sims <- t(sims)

		observed_y <- data_biomass_sites$y_biomass

		dharma_biomass <- createDHARMa(simulatedResponse = sims, observedResponse = observed_y)

		dharma_quant_test_biomass <- testQuantiles(dharma_biomass, plot = FALSE)
		dharma_resid_test_biomass <- testResiduals(dharma_biomass, plot = FALSE) # uniformity, dispersion, outlier
		dev.off()

		file <- paste0(out_dir, '/dharma_qq_quantiles_biomass.png')
		png(file, width = 1200, height = 800)
			plot(dharma_biomass)
		dev.off()

		file <- paste0(out_dir, '/dharma_distribution_biomass_residuals.png')
		png(file, width = 1200, height = 800)
			hist(dharma_biomass$scaledResiduals, main = 'DHARMa residuals for site-level mean biomass', xlab = 'Scaled residuals', breaks = 30)
		dev.off()
	
	### BIOMASS DHAMRa vs COVARIATES
	################################

		# stores results for all "residuals vs covariates" tests for biomass and non-biomass facets
		resids_vs_covariates_table <- data.table()

		say('   biomass DHARMa vs covariates...')

		resids <- dharma_biomass$scaledResiduals
		resids_vs_covariates <- list()

		for (i in seq_len(data_biomass_sites$n_covariates)) {

			covariate <- data_biomass_sites$covariates[i]

			if (data_biomass_sites$covariates[i] == 'ph') {
				x <- data_biomass_sites$raw_data_biomass[['site_ph']]
			} else {
				x <- data_biomass_sites$raw_data_biomass[[covariate]]
			}

			this_x <- data.frame(
				x = x,
				y = resids
			)

			resid_model <- mgcv::gam(logitAdj(y, epsilon = 0.0001) ~ s(x), data = this_x)

			F <- summary(resid_model)$s.table[1, 3]
			p <- summary(resid_model)$s.table[1, 4]
			sig <- ifelse(p < 0.05, '*', '-')
			
			resids_vs_covariates_table <- rbind(
				resids_vs_covariates_table,
				data.table(
					facet = 'biomass',
					covariate = covariate,
					F = F,
					p = p,
					significance = sig
				)
			)

			resids_vs_covariates[[i]] <- ggplot(this_x, aes(x = x, y = y)) +
				geom_point() +
				geom_smooth(method = loess, formula = y ~ x, se = FALSE) +
				annotate('text', x = -Inf, y = -Inf, label = paste0('P = ', round(p, 3)), hjust = -0.1, vjust = -1.2, size = 4) +
				xlab(data_biomass_sites$terms[i]) +
				ylab('Residual value') +
				ggtitle('DHAMRa Residuals for Biomass Means')

		}

		if (data_biomass_sites$n_covariates == 2) {
			ncol <- 2
		} else {
			ncol <- ceiling(sqrt(length(resids_vs_covariates)))
		}
		width <- 6 * ncol
		height <- 5 * ceiling(length(resids_vs_covariates) / ncol)

		resids_vs_covariates <- plot_grid(plotlist = resids_vs_covariates, ncol = ncol)
		ggsave(plot = resids_vs_covariates, filename = paste0(out_dir, '/dharma_residuals_vs_covariates_biomass.png'), width = width, height = height, dpi = 120)

	### DHARMa RESIDUALS: non-biomass
	#################################

		dharma_nonbiomass <- list()
		for (f in seq_along(nonbiomass_facets)) {

			facet <- names(nonbiomass_facets)[f]
			
			### GENERAL DHARMa
			##################

			say('   ', facet, ' general DHAMRa...')

			sims <- mc_subset(chains, paste0('y_facet_sim_', f), j = TRUE)
			sims <- mc_rbind(sims)
			sims <- t(sims)

			observed_y <- data_nonbiomass_sites[[facet]]$y_facet

			dharma_nonbiomass[[f]] <- createDHARMa(simulatedResponse = sims, observedResponse = observed_y)

			dharma_quant_test_nonbiomass <- testQuantiles(dharma_nonbiomass[[f]], plot = FALSE)
			dharma_resid_test_nonbiomass <- testResiduals(dharma_nonbiomass[[f]], plot = FALSE) # uniformity, dispersion, outlier
			dev.off()

			file <- paste0(out_dir, '/dharma_qq_quantiles_', facet, '.png')
			png(file, width = 1200, height = 800)
				plot(dharma_nonbiomass[[f]])
			dev.off()

			file <- paste0(out_dir, '/dharma_distribution_', facet, '_residuals.png')
			png(file, width = 1200, height = 800)
				hist(dharma_nonbiomass[[f]]$scaledResiduals, main = paste('DHARMa residuals for site-level mean', facet), xlab = 'Scaled residuals', breaks = 30)
			dev.off()

			assign(paste0('dharma_quant_test_', facet), dharma_quant_test_nonbiomass)
			assign(paste0('dharma_resid_test_', facet), dharma_resid_test_nonbiomass)

		} # next facet
		names(dharma_nonbiomass) <- names(nonbiomass_facets)

	### NON-BIOMASS DHAMRa vs COVARIATES
	####################################

		for (f in seq_along(nonbiomass_facets)) {

			facet <- names(nonbiomass_facets)[f]
			say('   ', facet, ' DHARMa vs covariates...')

			resids <- dharma_nonbiomass[[facet]]
			resids <- resids$scaledResiduals
			resids_vs_covariates <- list()
		
			for (i in seq_len(data_nonbiomass_sites[[facet]]$n_covariates)) {

				covariate <- data_nonbiomass_sites[[facet]]$covariates[i]

				if (data_nonbiomass_sites[[facet]]$covariates[i] == 'ph') {
					x <- data_nonbiomass_sites[[facet]]$raw_data_facet[['site_ph']]
				} else {
					x <- data_nonbiomass_sites[[facet]]$raw_data_facet[[covariate]]
				}

				this_x <- data.frame(
					x = x,
					y = resids
				)

				resid_model <- mgcv::gam(logitAdj(y, epsilon = 0.0001) ~ s(x), data = this_x)

				F <- summary(resid_model)$s.table[1, 3]
				p <- summary(resid_model)$s.table[1, 4]
				sig <- ifelse(p < 0.05, '*', '-')
				
				resids_vs_covariates_table <- rbind(
					resids_vs_covariates_table,
					data.table(
						facet = facet,
						covariate = covariate,
						F = F,
						p = p,
						significance = sig
					)
				)

				resids_vs_covariates[[i]] <- ggplot(this_x, aes(x = x, y = y)) +
					geom_point() +
					geom_smooth(method = loess, formula = y ~ x, se = FALSE) +
					annotate('text', x = -Inf, y = -Inf, label = paste0('P = ', round(p, 3)), hjust = -0.1, vjust = -1.2, size = 4) +
					xlab(data_nonbiomass_sites[[facet]]$terms_nonbiomass[i]) +
					ylab('Residual value') +
					ggtitle('DHAMRa Residuals for Biomass Means')

			}

			if (data_nonbiomass_sites[[facet]]$n_covariates == 2) {
				ncol <- 2
			} else {
				ncol <- ceiling(sqrt(length(resids_vs_covariates)))
			}
			width <- 6 * ncol
			height <- 5 * ceiling(length(resids_vs_covariates) / ncol)

			resids_vs_covariates <- plot_grid(plotlist = resids_vs_covariates, ncol = ncol)
			ggsave(plot = resids_vs_covariates, filename = paste0(out_dir, '/dharma_residuals_vs_covariates_', facet, '.png'), width = width, height = height, dpi = 120)

		} # next non-biomass facet

	### OCCURRENCE: spatial autocorrelation
	#######################################
	say('OCCURRENCE: spatial autocorrelation', level = 2)

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

		ggsave(plot = map, filename = paste0(out_dir, '/map_residuals_occurrences.png'), width = 12, height = 9, dpi = 120)

	### BIOMASS: spatial autocorrelation
	####################################
	say('BIOMASS: spatial autocorrelation', level = 2)

		coords <- as.data.frame(crds(project(data_biomass_sites$site_vect, enmSdmX::getCRS('WGS84'))))

		# Compute Moran's I
		resids_raw <- dharma_biomass$scaledResiduals
		resids <- numeric()
		index <- data_biomass_sites$site_index_biomass
		sites <- unique(index)
		for (i in sites) resids <- c(resids, mean(resids_raw[index == i]))
		
		moran_biomass <- moran.test(resids, nb2listw(knn2nb(knearneigh(coords, longlat = TRUE, k = 4))))
		moran_p <- moran_biomass$p.value
		moran_p <- paste0('Moran P = ', sprintf('%.3f', round(moran_p, 3)))

		residuals_vect <- data_biomass_sites$site_vect
		residuals_vect$residual <- resids

		extent <- ext(residuals_vect)
		extent <- as.vector(extent)

		nam <- vect('./data_from_gadm/gadm_4pt1_level_1_north_america_sans_alaska_lambert.gpkg')
		nam <- simplifyGeom(nam, tolerance = 1000)
		map <- ggplot() +
			layer_spatial(nam, fill = 'gray95', color = 'gray40') +
			layer_spatial(residuals_vect, aes(fill = residual), pch = 21, size = 5) +
			scale_fill_distiller(palette = 'RdBu', limits = c(0, 1), na.value = 'grey80') +
			annotate('text', x = extent[1], y = extent[3], label = moran_p, hjust = 0, vjust = 0, size = 6, color = 'black') +
			xlim(extent[1], extent[2]) + ylim(extent[3], extent[4]) +
			ggtitle(bquote('Residuals for ' * italic('Andropogon gerardi') * ' biomass')) +
			theme(
				plot.title = element_text(size = 16),
				plot.subtitle = element_text(size = 14)
			)

		ggsave(plot = map, filename = paste0(out_dir, '/map_residuals_biomass.png'), width = 12, height = 9, dpi = 120)

	### NON-BIOMASS FACETS: spatial autocorrelation
	###############################################

		nam <- vect('./data_from_gadm/gadm_4pt1_level_1_north_america_sans_alaska_lambert.gpkg')
		nam <- simplifyGeom(nam, tolerance = 1000)
		for (f in seq_along(nonbiomass_facets)) {

			facet <- names(nonbiomass_facets)[f]

			say(toupper(facet), ': spatial autocorrelation', level = 2)

			coords <- as.data.frame(crds(project(data_nonbiomass_sites[[facet]]$site_vect, enmSdmX::getCRS('WGS84'))))

			# Compute Moran's I
			resids_raw <- dharma_nonbiomass[[facet]]$scaledResiduals
			resids <- numeric()
			index <- data_nonbiomass_sites[[facet]]$site_index_facet
			sites <- unique(index)
			for (i in sites) resids <- c(resids, mean(resids_raw[index == i]))

			moran_nonbiomass <- moran.test(resids, nb2listw(knn2nb(knearneigh(coords, longlat = TRUE, k = 4))))
			moran_p <- moran_nonbiomass$p.value
			moran_p <- paste0('Moran P = ', sprintf('%.3f', round(moran_p, 3)))

			residuals_vect <- data_nonbiomass_sites[[facet]]$site_vect
			residuals_vect$residual <- resids

			extent <- ext(residuals_vect)
			extent <- as.vector(extent)

			map <- ggplot() +
				layer_spatial(nam, fill = 'gray95', color = 'gray40') +
				layer_spatial(residuals_vect, aes(fill = residual), pch = 21, size = 5) +
				scale_fill_distiller(palette = 'RdBu', limits = c(0, 1), na.value = 'grey80') +
				annotate('text', x = extent[1], y = extent[3], label = moran_p, hjust = 0, vjust = 0, size = 6, color = 'black') +
				xlim(extent[1], extent[2]) + ylim(extent[3], extent[4]) +
				ggtitle(bquote('Residuals for ' * italic('Andropogon gerardi') * ' ' * .(facet))) +
				theme(
					plot.title = element_text(size = 16),
					plot.subtitle = element_text(size = 14)
				)

			ggsave(plot = map, filename = paste0(out_dir, '/map_residuals_', facet, '.png'), width = 12, height = 9, dpi = 120)

			assign(paste0('moran_', facet), moran_nonbiomass)

		}

	### CURRENT MAPS
	################
	say('CURRENT MAPS', level = 2)

		### occurrences
		say('     occurrence...')
		form <- paste(as.character(formula_occs), collapse = ' ')
		form <- gsub(form, pattern = 'I\\(', replacement = '')
		form <- gsub(form, pattern = '\\^2\\)', replacement = '²')
		form <- gsub(form, pattern = '*)', replacement = '×')

		form_bias <- paste(as.character(formula_bias), collapse = ' ')
		form_bias <- gsub(form_bias, pattern = 'I\\(', replacement = '')
		form_bias <- gsub(form_bias, pattern = '\\^2\\)', replacement = '²')
		form_bias <- gsub(form_bias, pattern = '*)', replacement = '×')

		form_psi <- paste(as.character(formula_psi), collapse = ' ')
		form_psi <- gsub(form_psi, pattern = 'I\\(', replacement = '')
		form_psi <- gsub(form_psi, pattern = '\\^2\\)', replacement = '²')
		form_psi <- gsub(form_psi, pattern = '*)', replacement = '×')

		subtitle <- paste0('occ ', form, ' (bias ', form_bias, ')\nψ ', form_psi)

		map <- map_occurrence(
			out_dir = out_dir,
			filename_append = 'present_day',
			pred_vect_nam = pred_vect_nam,
			response_var = 'N_ag_county_mean_sq',
			data_occs = data_occs,
			title = bquote('Present-day ' * italic('Andropogon gerardi') * ' relative abundance (1961-2020)'),
			subtitle = subtitle
		)

		### zero-inflation
		say('     zero-inflation...')

		form_psi <- paste(as.character(formula_psi), collapse = ' ')
		form_psi <- gsub(form_psi, pattern = 'I\\(', replacement = '')
		form_psi <- gsub(form_psi, pattern = '\\^2\\)', replacement = '²')
		form_psi <- gsub(form_psi, pattern = '*)', replacement = '×')

		psi_map <- map_psi(
			out_dir = out_dir,
			filename_append = 'present_day',
			pred_vect_nam = pred_vect_nam,
			response_var = 'psi_county_sq',
			title = 'Present-day Probability of Occurrence (1961-2020)',
			subtitle = subtitle
		)

		### biomass
		say('     biomass...')
		form <- paste(as.character(formula_biomass), collapse = ' ')
		form <- gsub(form, pattern = 'I\\(', replacement = '')
		form <- gsub(form, pattern = '\\^2\\)', replacement = '²')
		form <- gsub(form, pattern = '*)', replacement = '×')

		subtitle <- paste0('biomass ', form, '\nψ ', form_psi)

		map <- map_biomass_nonbiomass(
			facet = 'biomass',
			out_dir = out_dir,
			filename_append = 'present_day',
			pred_vect_nam = pred_vect_nam,
			response_var = 'biomass_county_mean_sq',
			response_var_type = 'mean',
			data_occs = data_occs_counties,
			data_biomass_nonbiomass = data_biomass_sites,
			title = 'Present-day Biomass (1991-2020)',
			subtitle = subtitle,
			legend_title = 'Biomass\n(g)'
		)

		### non-biomass
		for (f in seq_along(nonbiomass_facets)) {

			facet <- names(nonbiomass_facets)[f]
			say('     ', facet, '...')
			facet_nice <- get_nice_trait(facet)

			form <- nonbiomass_facets[[facet]]$formula
			form <- paste(as.character(form), collapse = ' ')
			form <- gsub(form, pattern = 'I\\(', replacement = '')
			form <- gsub(form, pattern = '\\^2\\)', replacement = '²')
			form <- gsub(form, pattern = '*)', replacement = '×')

			subtitle <- paste0(facet_nice$short, ' ', form, '\nψ ', form_psi)

			map <- map_biomass_nonbiomass(
				facet = facet,
				out_dir = out_dir,
				filename_append = 'present_day',
				pred_vect_nam = pred_vect_nam,
				response_var = paste0(facet, '_county_mean_sq'),
				response_var_type = 'mean',
				data_occs = data_occs_counties,
				data_biomass_nonbiomass = data_nonbiomass_sites[[facet]],
				title = paste0('Present-day ', capIt(facet_nice$short), ' (1991-2020)'),
				subtitle = subtitle,
				legend_title = facet_nice$legend_title
			)

		}

	### FUTURE MAPS
	###############
	say('FUTURE MAPS', level = 2)

		for (fut in futs) {

			say(fut)

			### occurrence
			say('    occurrence...')
			form <- paste(as.character(formula_occs), collapse = ' ')
			form <- gsub(form, pattern = 'I\\(', replacement = '')
			form <- gsub(form, pattern = '\\^2\\)', replacement = '²')
			form <- gsub(form, pattern = '*)', replacement = '×')

			form_bias <- paste(as.character(formula_bias), collapse = ' ')
			form_bias <- gsub(form_bias, pattern = 'I\\(', replacement = '')
			form_bias <- gsub(form_bias, pattern = '\\^2\\)', replacement = '²')
			form_bias <- gsub(form_bias, pattern = '*)', replacement = '×')

			ssp <- paste0('SSP ', substr(fut, 4, 6), ' ', substr(fut, 8, 11), '-', substr(fut, 13, 16))
			title <- bquote('Future ' * italic('Andropogon gerardi') * ' Relative Abundance (' * .(ssp) * ')')

			form_psi <- paste(as.character(formula_psi), collapse = ' ')
			form_psi <- gsub(form_psi, pattern = 'I\\(', replacement = '')
			form_psi <- gsub(form_psi, pattern = '\\^2\\)', replacement = '²')
			form_psi <- gsub(form_psi, pattern = '*)', replacement = '×')

			subtitle <- paste0('occ ', form, ' (bias ', form_bias, ')\nψ ', form_psi)
			response_var <- paste0('N_ag_county_mean_', fut)

			map_occs_fut <- map_occurrence(
				out_dir = out_dir,
				filename_append = fut,
				pred_vect_nam = pred_vect_nam,
				response_var = response_var,
				data_occs = data_occs,
				title = title,
				subtitle = subtitle
			)

			### zero-inflation
			say('    zero-inflation...')
			title <- bquote('Future ' * italic('Andropogon gerardi') * ' Probability of Occurrence (' * .(ssp) * ')')

			psi_map <- map_psi(
				out_dir = out_dir,
				filename_append = fut,
				pred_vect_nam = pred_vect_nam,
				response_var = paste0('psi_', fut),
				title = title,
				subtitle = subtitle
			)

			### biomass
			say('    biomass...')
			form <- paste(as.character(formula_biomass), collapse = ' ')
			form <- gsub(form, pattern = 'I\\(', replacement = '')
			form <- gsub(form, pattern = '\\^2\\)', replacement = '²')
			form <- gsub(form, pattern = '*)', replacement = '×')

			title <- bquote('Future ' * italic('Andropogon gerardi') * ' Biomass (' * .(ssp) * ')')
			subtitle <- paste0('biomass ', form, '\nψ ', form_psi)

			map <- map_biomass_nonbiomass(
				facet = 'biomass',
				out_dir = out_dir,
				filename_append = fut,
				pred_vect_nam = pred_vect_nam,
				response_var = paste0('biomass_county_mean_', fut),
				response_var_type = 'mean',
				data_occs = data_occs_counties,
				data_biomass_nonbiomass = data_biomass_sites,
				title = title,
				subtitle = subtitle,
				legend_title = 'Biomass\n(g)'
			)

			### non-biomass
			for (f in seq_along(nonbiomass_facets)) {

				facet <- names(nonbiomass_facets)[f]
				say('    ', facet, '...')
				facet_nice <- get_nice_trait(facet)

				form <- nonbiomass_facets[[facet]]$formula
				form <- paste(as.character(form), collapse = ' ')
				form <- gsub(form, pattern = 'I\\(', replacement = '')
				form <- gsub(form, pattern = '\\^2\\)', replacement = '²')
				form <- gsub(form, pattern = '*)', replacement = '×')

				title <- bquote('Future ' * italic('Andropogon gerardi') * ' ' * .(capIt(facet_nice$short)) * ' (' * .(ssp) * ')')
				subtitle <- paste0(facet_nice$short, ' ', form, '\nψ ', form_psi)

				map <- map_biomass_nonbiomass(
					facet = facet,
					out_dir = out_dir,
					filename_append = fut,
					pred_vect_nam = pred_vect_nam,
					response_var = paste0(facet, '_county_mean_', fut),
					response_var_type = 'mean',
					data_occs = data_occs_counties,
					data_biomass_nonbiomass = data_nonbiomass_sites[[facet]],
					title = title,
					subtitle = subtitle,
					legend_title = facet_nice$legend_title
				)

			}

		}

	### FUTURE CHANGE
	#################
	say('FUTURE CHANGE', level = 2)

		for (fut in futs) {

			say(fut)
			say('     occurrence...')

			ssp <- paste0('SSP ', substr(fut, 4, 6), ' ', substr(fut, 8, 11), '-', substr(fut, 13, 16))

			form <- paste(as.character(formula_occs), collapse = ' ')
			form <- gsub(form, pattern = 'I\\(', replacement = '')
			form <- gsub(form, pattern = '\\^2\\)', replacement = '²')
			form <- gsub(form, pattern = '*)', replacement = '×')

			form_bias <- paste(as.character(formula_bias), collapse = ' ')
			form_bias <- gsub(form_bias, pattern = 'I\\(', replacement = '')
			form_bias <- gsub(form_bias, pattern = '\\^2\\)', replacement = '²')
			form_bias <- gsub(form_bias, pattern = '*)', replacement = '×')

			ssp <- paste0('SSP ', substr(fut, 4, 6), ' ', substr(fut, 8, 11), '-', substr(fut, 13, 16))
			title <- bquote('Change in ' * italic('Andropogon gerardi') * ' Relative Abundance (' * .(ssp) * ')')

			form_psi <- paste(as.character(formula_psi), collapse = ' ')
			form_psi <- gsub(form_psi, pattern = 'I\\(', replacement = '')
			form_psi <- gsub(form_psi, pattern = '\\^2\\)', replacement = '²')
			form_psi <- gsub(form_psi, pattern = '*)', replacement = '×')

			subtitle <- paste0('occ ', form, ' (bias ', form_bias, ')\nψ ', form_psi)

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

			### zero-inflation
			say('     zero-inflation...')
			title <- bquote('Change in ' * italic('Andropogon gerardi') * ' Probability of Occurrence (' * .(ssp) * ')')

			psi_map <- map_psi_change(
				out_dir = out_dir,
				filename_append = fut,
				fut = fut,
				pred_vect_nam = pred_vect_nam,
				response_var = paste0('psi_', fut),
				response_var_sq = 'psi_county_sq',
				title = title,
				subtitle = subtitle
			)

			### biomass
			say('     biomass...')
			form <- paste(as.character(formula_biomass), collapse = ' ')
			form <- gsub(form, pattern = 'I\\(', replacement = '')
			form <- gsub(form, pattern = '\\^2\\)', replacement = '²')
			form <- gsub(form, pattern = '*)', replacement = '×')

			title <- bquote('Future ' * italic('Andropogon gerardi') * ' Biomass (' * .(ssp) * ')')
			subtitle <- paste0('biomass ', form, '\nψ ', form_psi)

			response_var <- paste0('biomass_county_mean_', fut)

			map <- map_biomass_nonbiomass_change(
				facet = 'biomass',
				out_dir = out_dir,
				filename_append = fut,
				fut = fut,
				pred_vect_nam = pred_vect_nam,
				response_var = response_var,
				response_var_type = 'mean',
				data_occs = data_occs_sites,
				data_biomass_nonbiomass = data_biomass_sites,
				title = title,
				subtitle = subtitle,
				legend_title = 'Percent\nchange'
			)

			### non-biomass
			for (f in seq_along(nonbiomass_facets)) {

				facet <- names(nonbiomass_facets)[f]
				say('     ', facet, '...')
				facet_nice <- get_nice_trait(facet)

				form <- nonbiomass_facets[[facet]]$formula
				form <- paste(as.character(form), collapse = ' ')
				form <- gsub(form, pattern = 'I\\(', replacement = '')
				form <- gsub(form, pattern = '\\^2\\)', replacement = '²')
				form <- gsub(form, pattern = '*)', replacement = '×')

				title <- bquote('Change in ' * italic('Andropogon gerardi') * ' ' * .(capIt(facet_nice$short)) * ' (' * .(ssp) * ')')
				subtitle <- paste0(facet_nice$short, ' ', form, '\nψ ', form_psi)

				response_var <- paste0(facet, '_county_mean_', fut)

				map <- map_biomass_nonbiomass_change(
					facet = facet,
					out_dir = out_dir,
					filename_append = fut,
					fut = fut,
					pred_vect_nam = pred_vect_nam,
					response_var = response_var,
					response_var_type = 'mean',
					data_occs = data_occs_sites,
					data_biomass_nonbiomass = data_nonbiomass_sites[[facet]],
					title = title,
					subtitle = subtitle,
					legend_title = 'Percent\nchange'
				)

			}

		}

	### 1930s CHANGE MAPS
	#####################
	say('1930S CHANGE MAPS', level = 2)

		say('     occuurrence...')
		map <- map_occurrence_change_1930s(
			zero_inflated = TRUE,
			pred_vect_nam = pred_vect_nam,
			pred_vect_1930s = pred_vect_1930s
		)

		say('     biomass...')
		map <- map_biomass_nonbiomass_change_1930s(
			facet = 'biomass',
			data_biomass_nonbiomass = data_biomass_sites,
			pred_vect_nam = pred_vect_nam,
			pred_vect_1930s = pred_vect_1930s,
			formula_psi = formula_psi
		)

		for (f in seq_along(nonbiomass_facets)) {

			facet <- names(nonbiomass_facets)[f]
			say('     ', facet, '...')

			map <- map_biomass_nonbiomass_change_1930s(
				facet = facet,
				data_biomass_nonbiomass = data_nonbiomass_sites[[facet]],
				pred_vect_nam = pred_vect_nam,
				pred_vect_1930s = pred_vect_1930s,
				formula_psi = formula_psi
			)

		}

	### SUMMARY
	###########
	say('SUMMARY', level = 2)
	
		### DHARMa residuals summary
		dharma_resids <- data.table()
		for (facet in c('occs', 'biomass', names(nonbiomass_facets))) {

			quant <- get(paste0('dharma_quant_test_', facet))
			resid <- get(paste0('dharma_resid_test_', facet))
			moran <- get(paste0('moran_', facet))

			if (length(quant$qgamFits) > 0) {
				
				quant_test_upper <- summary(quant$qgamFits[[3]])$s.table[1, 4]
				quant_test_middle <- summary(quant$qgamFits[[2]])$s.table[1, 4]
				quant_test_lower <- summary(quant$qgamFits[[1]])$s.table[1, 4]

			} else {
				quant_test_lower <- quant_test_middle <- quant_test_upper <- NA_real_
			}

			resid_p_values <- c(
				moran$p.value,
				resid$uniformity$p.value,
				resid$dispersion$p.value,
				resid$outliers$p.value,
				quant$p.value,
				quant_test_upper,
				quant_test_middle,
				quant_test_lower
			)

			if (length(quant$qgamFits) > 0) {
				
				quant_test_upper <- summary(quant$qgamFits[[3]])$s.table[1, 4]
				quant_test_middle <- summary(quant$qgamFits[[2]])$s.table[1, 3]
				quant_test_lower <- summary(quant$qgamFits[[1]])$s.table[1, 3]

			} else {
				quant_test_lower <- quant_test_middle <- quant_test_upper <- NA_real_
			}

			resid_p_values_sig <- ifelse(resid_p_values < 0.05, '*', 'ns')
			resid_test_statistic_values <- c(
				moran$statistic,
				resid$uniformity$statistic,
				resid$dispersion$statistic,
				resid$outliers$statistic,
				NA,
				quant_test_upper,
				quant_test_middle,
				quant_test_lower
			)

			this_dharma_resids <- data.table(
				facet = facet,
				test = c('spatial autocorrelation', 'uniformity', 'dispersion', 'outliers', 'quantiles, overall', 'quantiles, upper', 'quantiles, middle', 'quantiles, lower'),
				p_value = resid_p_values,
				significant = resid_p_values_sig,
				test_statistic = c('Moran\'s I', names(resid$uniformity$statistic), names(resid$dispersion$statistic), 'exact binomial', NA, rep('chi squared', 3)),
				test_statistic_value = resid_test_statistic_values
			)

			dharma_resids <- rbind(dharma_resids, this_dharma_resids)

		}

		# compile residuals analysis
		meta_fully_integrated <- list(
			facet = c('occurrence', 'biomass', names(nonbiomass_facets)),
			date = date(),
			zero_inflated = TRUE,
			formulae = formulae,
			nonbiomass_facets = nonbiomass_facets,
			resids_vs_covariates = resids_vs_covariates_table,
			dharma_resids = dharma_resids
		)

		saveRDS(meta_fully_integrated, paste0(out_dir, '/!meta_fully_integrated.rds'))
		sink(paste0(out_dir, '/!meta_fully_integrated.txt'), split = TRUE)
			print(meta_fully_integrated)
		sink()

}

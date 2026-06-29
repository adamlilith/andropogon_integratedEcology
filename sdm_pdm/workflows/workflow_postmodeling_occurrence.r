#' Post-modeling workflow for occurrence-only models and models with an occurrence component.
#'
#' source('C:/Kaji/R/andropogon_integratedEcology/sdm_pdm/workflows/workflow_postmodeling_occurrence.r')
#'
#' formula_occs,formula_psi,formula_bias Formulae for occurrences, sigma, probability of zero-inflation, and  occurrence bias
#' out_dir Folder into which to save results.
workflow_postmodeling_occurrence <- function(formula_occs, formula_psi, formula_bias, out_dir, overdispersed = FALSE) {

	zero_inflated <- !is.null(formula_psi)

	### burn predictions into vector
	################################
	say('OCCURRENCE: burn prediction vectors', level = 2)

		pred_vect_nam <- burn_occs_into_vector(demesne = 'nam', chains = chains, formula_occs = formula_occs, formula_psi = formula_psi, overdispersed = overdispersed)
		pred_vect_1930s <- burn_occs_into_vector(demesne = '1930s', chains = chains, formula_occs = formula_occs, formula_psi = formula_psi, overdispersed = overdispersed)

		writeVector(pred_vect_nam, paste0(out_dir, '/prediction_vector_nam.gpkg'), overwrite = TRUE)
		writeVector(pred_vect_1930s, paste0(out_dir, '/prediction_vector_conus_1930s.gpkg'), overwrite = TRUE)

	### response curves: OCCURRENCE vs environment
	##############################################
	say('OCCURRENCE: response curves', level = 2)

		responses <- graph_response_curves_occurrence_vs_environment(
			out_dir = out_dir,
			chains = chains,
			zero_inflated = zero_inflated,
			overdispersed = overdispersed
		)

		if (zero_inflated) {

			responses_psi <- graph_response_curves_psi_vs_environment(
				out_dir = out_dir,
				chains = chains
			)

		}

		if (data_occs$n_covariates_bias >= 1) {

			responses_bias <- graph_response_curves_bias_vs_bias_covariates(
				out_dir = out_dir,
				chains = chains
			)

		}

	### occurrence: DHARMa residuals
	################################
	say('OCCURRENCE: DHARMa residuals', level = 2)
		
		sims <- mc_subset(chains, 'y_n_ag_sim', j = TRUE)
		sims <- mc_rbind(sims)
		sims <- sims[ , data_occs$ag_vect_sq$focal_region]

		observed_y <- data_occs$y_n_ag[data_occs$ag_vect_sq$focal_region]

		nas_occs <- which(is.na(colSums(sims)))
		if (length(nas_occs) > 0) {
			sims <- sims[ , -nas_occs]
			observed_y <- observed_y[-nas_occs]
		}

		sims <- t(sims)
		dharma <- createDHARMa(simulatedResponse = sims, observedResponse = observed_y, integerResponse = TRUE)

		dharma_quant_test <- testQuantiles(dharma, plot = FALSE)
		dharma_resid_test <- testResiduals(dharma, plot = FALSE) # uniformity, dispersion, outlier
		dev.off()

		file <- paste0(out_dir, '/dharma_n_ag.png')
		png(file, width = 1200, height = 800)
			plot(dharma)
		dev.off()

		file <- paste0(out_dir, '/dharma_lambda_residuals_n_ag.png')
		png(file, width = 1200, height = 800)
			hist(dharma$scaledResiduals, main = 'DHARMa residuals for number of observed AG (y_n_ag)', xlab = 'Scaled residuals', breaks = 30)
		dev.off()

	### OCCURRENCE: spatial autocorrelation
	#######################################

		coords <- as.data.frame(crds(centroids(project(pred_vect_nam[pred_vect_nam$focal_region], enmSdmX::getCRS('WGS84')))))
		if (length(nas_occs) > 0) coords <- coords[-nas_occs, ]

		# Compute Moran's I
		moran <- moran.test(dharma$scaledResiduals, nb2listw(knn2nb(knearneigh(coords, longlat = TRUE, k = 4))))
		moran_p <- moran$p.value
		moran_p <- paste0('Moran P = ', sprintf('%.3f', round(moran_p, 3)))

		residuals_vect <- pred_vect_nam[pred_vect_nam$focal_region]
		if (length(nas_occs) > 0) residuals_vect <- residuals_vect[-nas_occs]
		residuals_vect$residual <- dharma$scaledResiduals

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

	### OCCURRENCE: current map
	###########################
	say('OCCURRENCE: current map', level = 2)

		form <- paste(as.character(formula_occs), collapse = ' ')
		form <- gsub(form, pattern = 'I\\(', replacement = '')
		form <- gsub(form, pattern = '\\^2\\)', replacement = '²')
		form <- gsub(form, pattern = '*)', replacement = '×')

		# form_bias <- paste(as.character(formula_bias), collapse = ' ')
		# form_bias <- gsub(form_bias, pattern = 'I\\(', replacement = '')
		# form_bias <- gsub(form_bias, pattern = '\\^2\\)', replacement = '²')
		# form_bias <- gsub(form_bias, pattern = '*)', replacement = '×')

		map <- map_occurrence(
			out_dir = out_dir,
			filename_append = 'present_day',
			pred_vect_nam = pred_vect_nam,
			response_var = 'N_ag_mean_sq',
			data_occs = data_occs,
			title = bquote('Present-day Relative Abundance'),
			subtitle = paste0('1961-2020 | occ ', form, ' (bias ~ offset)')
		)

		if (zero_inflated) {

			form_psi <- paste(as.character(formula_psi), collapse = ' ')
			form_psi <- gsub(form_psi, pattern = 'I\\(', replacement = '')
			form_psi <- gsub(form_psi, pattern = '\\^2\\)', replacement = '²')
			form_psi <- gsub(form_psi, pattern = '*)', replacement = '×')

			map <- map_psi(
				out_dir = out_dir,
				filename_append = 'present',
				pred_vect_nam = pred_vect_nam,
				response_var = 'psi_sq',
				title = 'Present-day Probability Occurrence',
				subtitle = paste0('1961-2020 | occ ', form, ' (ψ ', form_psi, ')')
			)

		}

	### occurrence: future maps
	###########################
	say('OCCURRENCE: future maps', level = 2)

		for (fut in futs) {

			say(fut)

			form <- paste(as.character(formula_occs), collapse = ' ')
			form <- gsub(form, pattern = 'I\\(', replacement = '')
			form <- gsub(form, pattern = '\\^2\\)', replacement = '²')
			form <- gsub(form, pattern = '*)', replacement = '×')

			title <- bquote('Future Relative Abundance')
			subtitle <- paste0('SSP ', substr(fut, 4, 6), ' ', substr(fut, 8, 11), '-', substr(fut, 13, 16))
			response_var <- paste0('N_ag_mean_', fut)

			map <- map_occurrence(
				out_dir = out_dir,
				filename_append = fut,
				pred_vect_nam = pred_vect_nam,
				response_var = response_var,
				data_occs = data_occs,
				title = title,
				subtitle = subtitle
			)

			if (zero_inflated) {

				title <- bquote('Future Probability of Occurrence')

				form_psi <- paste(as.character(formula_psi), collapse = ' ')
				form_psi <- gsub(form_psi, pattern = 'I\\(', replacement = '')
				form_psi <- gsub(form_psi, pattern = '\\^2\\)', replacement = '²')
				form_psi <- gsub(form_psi, pattern = '*)', replacement = '×')

				subtitle <- paste0('SSP ', substr(fut, 4, 6), ' ', substr(fut, 8, 11), '-', substr(fut, 13, 16))

				map <- map_psi(
					out_dir = out_dir,
					filename_append = fut,
					pred_vect_nam = pred_vect_nam,
					response_var = paste0('psi_', fut),
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

			response_var <- paste0('N_ag_mean_', fut)

			map <- map_occurrence_change(
				out_dir = out_dir,
				filename_append = fut,
				pred_vect_nam = pred_vect_nam,
				response_var = response_var,
				data_occs = data_occs,
				title = bquote('Change in Relative Abundance'),
				subtitle = paste0('SSP ', substr(fut, 4, 6), ' ', substr(fut, 8, 11), '-', substr(fut, 13, 16))
			)

			if (zero_inflated) {

				title <- bquote('Change in Probability of Occurrence')

				map <- map_psi_change(
					out_dir = out_dir,
					filename_append = fut,
					fut = fut,
					pred_vect_nam = pred_vect_nam,
					response_var = paste0('psi_', fut),
					response_var_sq = 'psi_sq',
					title = title,
					subtitle = paste0('SSP ', substr(fut, 4, 6), ' ', substr(fut, 8, 11), '-', substr(fut, 13, 16))
				)

			}

		}

	### OCCURRENCE: 1930s
	#####################
	say('OCCURRENCE: 1930s change maps', level = 2)

		map_1930s <- map_occurrence_change_1930s(
			zero_inflated = zero_inflated,
			pred_vect_nam = pred_vect_nam,
			pred_vect_1930s = pred_vect_1930s
		)

	### summary
	###########
	say('SUMMARY', level = 2)
	
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
	resid_p_values_sig <- ifelse(resid_p_values < 0.05, '*', 'ns')
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

	# compile residuals analysis
	meta_occs <- list(
		facet = 'occurrence',
		date = date(),
		zero_inflated = zero_inflated,
		formulae = list(
			formula_occs = formula_occs,
			formula_bias = NA,
			formula_psi = formula_psi
		),
		dharma_resids = dharma_resids
	)

	saveRDS(meta_occs, paste0(out_dir, '/!meta_occs.rds'))
	sink(paste0(out_dir, '/!meta_occs.txt'), split = TRUE)
		print(meta_occs)
	sink()

}

#' Post-modeling workflow for occurrence-only models and models with an occurrence component.
#'
#' @param laplace Output of `runLaplace()`
#' @param homoscedastic If `TRUE`, then do not analyze behavior of sigma
#' @param zero_inflated Logical.
#' @param formula_occs,formula_bias Formulae for occurrences and for occurrence bias
#' @param data_occs Output of `prepare_occurrence_data()`
#' @param out_dir Folder into which to save results.
workflow_postmodeling_occurrence_laplace <- function(laplace, homoscedastic, zero_inflated, formula_occs, formula_bias, data_occs, out_dir) {

	### response curves: occurrence vs environment
	##############################################

	say('OCCURRENCE: response curves: occurrence vs environment', level = 2)

	# linear terms of formula
	terms <- terms(formula_occs)
	terms <- attr(terms, 'term.labels')

	terms_linear <- terms[!grepl(terms, pattern = '\\^2\\)')]
	terms_linear <- terms_linear[!grepl(terms_linear, pattern = ':')]

	# mean values across counties with occurrences
	x <- data_occs$ag_vect_sq
	x <- as.data.table(x)
	x <- x[x$focal_region, ..terms_linear]
	mean_values <- sapply(x, mean, na.rm = TRUE)

	# min/max values across area of interest
	x <- data_occs$ag_vect_sq
	x <- as.data.table(x)
	keeps <- which(x$focal_region)
	x <- x[x$focal_region, ..terms_linear]

	min_values <- sapply(x, min, na.rm = TRUE)
	max_values <- sapply(x, max, na.rm = TRUE)

	for (fut in futs) {
	
		x <- data_occs[[paste0('counties_occs_', fut)]]
		x <- as.data.table(x)
		x <- x[keeps, ..terms_linear]
	
		min_values_fut <- sapply(x, min, na.rm = TRUE)
		max_values_fut <- sapply(x, max, na.rm = TRUE)
	
		min_values <- pmin(min_values, min_values_fut)
		max_values <- pmax(max_values, max_values_fut)

	}

	# construct array with all covariates held constant except for focal one
	# each "page" pertains to a different covariate's gradient
	n_rows <- 200
	n_terms_linear <- length(terms_linear)
	env_array <- array(NA, dim = c(n_rows, n_terms_linear, n_terms_linear))
	for (i in seq_along(terms_linear)) {
	
		term_linear <- terms_linear[i]
		env_array[ , , i] <- rep(mean_values, each = n_rows)
		env_array[ , i, i] <- seq(min_values[term_linear], max_values[term_linear], length.out = n_rows)
		colnames(env_array) <- terms_linear
	
	}

	# construct scaled matrix
	x_centers_occs <- data_occs$x_centers_occs
	x_scales_occs <- data_occs$x_scales
	
	x_centers_occs <- x_centers_occs[terms_linear]
	x_scales_occs <- x_scales_occs[terms_linear]

	scaled_array <- env_array
	for (i in seq_along(terms_linear)) {
		scaled_array[ , , i] <- scale(scaled_array[ , , i], center = x_centers_occs[terms_linear], scale = x_scales_occs[terms_linear])
	}

	# construct model matrices
	n_terms <- length(terms)
	mm_array <- array(NA, dim = c(n_rows, 1 + n_terms, n_terms_linear))
	for (i in seq_along(terms_linear)) {
		scaled_matrix <- as.data.frame(scaled_array[ , , i])
		mm_array[ , , i] <- model.matrix(formula_occs, scaled_matrix)
		colnames(mm_array) <- c('(Intercept)', terms)
	}

	# make predictions by simulating from posterior and plot
	n_sims <- 10000 # number of simulations

	coeffs <- laplace$summary$params
	coeffs <- coeffs[grepl(rownames(coeffs), pattern = 'beta_occs'), 'estimate', drop = FALSE]
	coeffs <- as.matrix(coeffs)

	if (homoscedastic) {

		log_lambda_sigma <- laplace$summary$params['log_lambda_sigma', 'estimate']
		log_lambda_sigma_sd <- laplace$summary$params['log_lambda_sigma', 'stdError']
		lambda_sigma <- rlnorm(1, log_lambda_sigma, log_lambda_sigma_sd)

	}

	responses <- list()
	max_val <- -Inf
	preds <- matrix(NA_real_, nrow = n_rows, ncol = n_sims)

	for (i in seq_along(terms_linear)) {
		
		# mean prediction
		this_pred_mean <- (mm_array[ , , i] %*% coeffs)[ , 1]

		if (homoscedastic) {
			
			for (j in 1:n_rows) {

				for (k in 1:n_sims) {

					lambda_sigma_sim <- rlnorm(1, log_lambda_sigma, log_lambda_sigma_sd)

					log_lambda_mu <- rnorm(1, this_pred_mean[j], sd = lambda_sigma_sim)
					lambda_mu <- exp(log_lambda_mu)
		
					preds[j, k] <- rpois(1, lambda_mu)
		
				}

			}

		} else {
			stop('DIDNT DO THIS PART YET')
		}
	
		pred_lower <- apply(preds, 1, quantile, 0.1, na.rm = TRUE)
		pred_upper <- apply(preds, 1, quantile, 0.9, na.rm = TRUE)
		pred_mean <- apply(preds, 1, mean, na.rm = TRUE)
		pred_median <- apply(preds, 1, median, na.rm = TRUE)

		# plot
		pred <- terms_linear[i]

		bounds <- data.frame(
			x = c(env_array[ , i, i], rev(env_array[ , i, i])),
			y = c(pred_lower, rev(pred_upper))
		)

		mean_median <- data.frame(
			x = rep(env_array[ , i, i], 2),
			y = c(pred_mean, pred_median),
			type = c(
				rep('mean', n_rows),
				rep('median', n_rows)
			)
		)

		# max_val <- max(max_val, quantile(bounds$y, 0.9))
		max_val <- max(max_val, bounds$y)

		pred_nice <- get_nice_predictor(pred)
		xlab <- pred_nice$long
		title <- paste0(LETTERS[i], ') ', pred_nice$short)
	
		# line type for each type
		line_types <- c('mean' = 'solid', 'median' = 'solid')
		line_widths <- c('mean' = 1.6, 'median' = 0.8)
		line_colors <- c('mean' = 'black', 'median' = 'gray50')
		mean_median$type <- factor(mean_median$type, levels = names(line_types))

		responses[[i]] <- ggplot() +
			geom_polygon(
				data = bounds,
				mapping = aes(x = x, y = y),
				color = NA,
				fill = alpha('blue', 0.1)
			) +
			geom_line(
				data = mean_median,
				mapping = aes(x = x, y = y, linetype = type, linewidth = type, color = type)
			) +
			scale_linetype_manual(name = NULL, values = line_types) +
			scale_linewidth_manual(name = NULL, values = line_widths) +
			scale_color_manual(name = NULL, values = line_colors) +
			xlab(xlab) + ylab('Relative abundance') +
			ggtitle(title) +
			# scale_y_log10() +
			theme(
				legend.position = 'inside',
				legend.position.inside = c(0.98, 0.98),
				legend.justification = c('right', 'top'),
				legend.background = element_blank()
			)

	} # next covariate

	for (i in seq_along(terms_linear)) {
		responses[[i]] <- responses[[i]] + coord_cartesian(ylim = c(0, max_val)) 
	}

	nrow <- floor(sqrt(n_terms_linear))
	ncol <- ceiling(n_terms_linear / nrow)
	responses <- plot_grid(plotlist = responses, nrow = nrow, ncol = ncol, align = 'h')

	width <- 6 * ncol
	height <- 6 * nrow

	ggsave(responses, filename = paste0(out_dir, '/response_curves_occs.png'), width = width, height = height, dpi = 400, bg = 'white')


	# ### OCCURRENCE: current map
	# ###########################

	# say('OCCURRENCE: current map', level = 2)
	
	# # linear terms of formula
	# terms <- terms(formula_occs)
	# terms <- attr(terms, 'term.labels')

	# terms_linear <- terms[!grepl(terms, pattern = '\\^2\\)')]
	# terms_linear <- terms_linear[!grepl(terms_linear, pattern = ':')]

	# # model matrix
	# x <- data_occs$ag_vect_sq
	# x <- as.data.table(x)
	# x <- x[ , ..terms_linear]

	# # scale
	# x_centers_occs <- data_occs$x_centers_occs
	# x_scales_occs <- data_occs$x_scales
	
	# x_centers_occs <- x_centers_occs[terms_linear]
	# x_scales_occs <- x_scales_occs[terms_linear]

	# x <- scale(x, center = x_centers_occs, scale = x_scales_occs)
	# x <- as.data.table(x)

	# mm <- model.matrix(formula_occs, x)

	# coeffs <- laplace$summary$params
	# coeffs <- coeffs[grepl(rownames(coeffs), pattern = 'beta_occs'), 'estimate', drop = FALSE]
	# coeffs <- as.matrix(coeffs)

	# pred <- mm %*% coeffs
	# pred <- pred[ , 1]
	# pred <- exp(pred)

	# x_vect <- data_occs$ag_vect_sq
	# x_vect$pred <- pred

	# x_vect <- x_vect[x_vect$focal_region]
	# pred_clip <- x_vect$pred

	# quant_class <- rep(NA_integer_, nrow(x_vect))
	# q1 <- quantile(pred_clip, 0.5)
	# q2 <- quantile(pred_clip, 0.75)
	# q3 <- quantile(pred_clip, 0.90)
	# q4 <- quantile(pred_clip, 0.95)
	# for (i in seq_along(pred_clip)) {
	# 	quant_class[i] <- if (pred_clip[i] <= q1) {
	# 		0
	# 	} else if (pred_clip[i] <= q2) {
	# 		1
	# 	} else if (pred_clip[i] <= q3) {
	# 		2
	# 	} else if (pred_clip[i] <= q4) {
	# 		3
	# 	} else {
	# 		4
	# 	}
	# }

	# ggplot() +
	# 	# layer_spatial(x_vect, aes(fill = pred), color = NA) #+
	# 	# layer_spatial(x_vect, aes(fill = q), color = NA) #+
	# 	layer_spatial(x_vect, aes(fill = quant_class), color = NA) +
	# 	# scale_fill_continuous(trans = 'log10') +
	# 	layer_spatial(nam, fill = NA, color = 'gray')









	# ### occurrence: DHARMa residuals
	# ################################
	# say('OCCURRENCE: DHARMa residuals', level = 2)

	# 	sims <- mc_subset(chains, param = 'y_n_ag_sim', j = TRUE)
	# 	sims <- mc_rbind(sims)
	# 	sims <- sims[ , data_occs$ag_vect_sq$focal_region]
	# 	sims <- t(sims)

	# 	fits <- mc_extract(chains, param = 'lambda_mu_sq', j = TRUE, stat = 'mean')
	# 	fits <- fits[data_occs$ag_vect_sq$focal_region]

	# 	observed_y <- data_occs$y_n_ag[data_occs$ag_vect_sq$focal_region]

	# 	dharma <- createDHARMa(simulatedResponse = sims, observedResponse = observed_y, fittedPredictedResponse = fits, integerResponse = TRUE)

	# 	dharma_quant_test <- testQuantiles(dharma, plot = FALSE)
	# 	dharma_resid_test <- testResiduals(dharma, plot = FALSE) # uniformity, dispersion, outlier
	# 	dev.off()

	# 	file <- paste0(out_dir, '/dharma_n_ag.png')
	# 	png(file, width = 1200, height = 800)
	# 		plot(dharma)
	# 	dev.off()

	# 	file <- paste0(out_dir, '/y_n_ag_dharma_lambda_residuals_n_ag.png')
	# 	png(file, width = 1200, height = 800)
	# 		hist(dharma$scaledResiduals, main = 'DHARMa residuals for number of observed AG (y_n_ag)', xlab = 'Scaled residuals', breaks = 30)
	# 	dev.off()

	# ### OCCURRENCE: spatial autocorrelation
	# #######################################

	# coords <- as.data.frame(crds(centroids(project(pred_vect_nam[pred_vect_nam$focal_region], enmSdmX::getCRS('WGS84')))))

	# # Compute Moran's I
	# moran <- moran.test(dharma$scaledResiduals, nb2listw(knn2nb(knearneigh(coords, longlat = TRUE, k = 4))))

	# ### OCCURRENCE: current map
	# ###########################
	# say('OCCURRENCE: current map', level = 2)

	# 	form <- paste(as.character(formula_occs), collapse = ' ')
	# 	form <- gsub(form, pattern = 'I\\(', replacement = '')
	# 	form <- gsub(form, pattern = '\\^2\\)', replacement = '²')
	# 	form <- gsub(form, pattern = '*)', replacement = '×')

	# 	form_bias <- paste(as.character(formula_bias), collapse = ' ')
	# 	form_bias <- gsub(form_bias, pattern = 'I\\(', replacement = '')
	# 	form_bias <- gsub(form_bias, pattern = '\\^2\\)', replacement = '²')
	# 	form_bias <- gsub(form_bias, pattern = '*)', replacement = '×')

	# 	map <- map_occurrence(
	# 		out_dir = out_dir,
	# 		filename_append = 'present_day',
	# 		pred_vect_nam = pred_vect_nam,
	# 		response_var = 'N_ag_county_sq',
	# 		data_occs = data_occs,
	# 		title = bquote('Present-day distribution of ' * italic('Andropogon gerardi') * ' abundance'),
	# 		subtitle = paste0('1961-2020 | occ ', form, '(bias ', form_bias, ')'),
	# 		ag_core_quant = ag_core_quant
	# 	)

	# 	if (zero_inflated) {

	# 		map_psi <- map_psi(
	# 			out_dir = out_dir,
	# 			filename_append = 'present_day',
	# 			pred_vect_nam = pred_vect_nam,
	# 			facet = 'occs',
	# 			response_var = 'psi_county_sq',
	# 			data_traits = data_traits,
	# 			title = 'Present-day distribution of probability of zero abundance abundance',
	# 			subtitle = paste0('1961-2020 | occ ', form, '(bias ', form_bias, ')'),
	# 			plot_range_core = TRUE,
	# 			ag_core_quant = ag_core_quant
	# 		)

	# 	}

	# ### occurrence: future maps
	# ###########################
	# say('OCCURRENCE: future maps', level = 2)

	# 	maps_occs_fut <- list()
	# 	for (fut in futs) {

	# 		say(fut)

	# 		form <- paste(as.character(formula_occs), collapse = ' ')
	# 		form <- gsub(form, pattern = 'I\\(', replacement = '')
	# 		form <- gsub(form, pattern = '\\^2\\)', replacement = '²')
	# 		form <- gsub(form, pattern = '*)', replacement = '×')

	# 		form_bias <- paste(as.character(formula_bias), collapse = ' ')
	# 		form_bias <- gsub(form_bias, pattern = 'I\\(', replacement = '')
	# 		form_bias <- gsub(form_bias, pattern = '\\^2\\)', replacement = '²')
	# 		form_bias <- gsub(form_bias, pattern = '*)', replacement = '×')

	# 		title <- bquote('Future distribution of ' * italic('Andropogon gerardi') * ' abundance')
	# 		subtitle <- paste0('SSP ', substr(fut, 4, 6), ' ', substr(fut, 8, 11), '-', substr(fut, 13, 16), ' | occ ', form, ' (bias ', form_bias, ')')
	# 		response_var <- paste0('N_ag_county_', fut)

	# 		maps_occs_fut[[length(maps_occs_fut) + 1]] <- map_occurrence(
	# 			out_dir = out_dir,
	# 			filename_append = fut,
	# 			pred_vect_nam = pred_vect_nam,
	# 			response_var = response_var,
	# 			data_occs = data_occs,
	# 			title = title,
	# 			subtitle = subtitle,
	# 			ag_core_quant = ag_core_quant
	# 		)

	# 	}

	# ### occurrence: future change
	# #############################
	# say('OCCURRENCE: future change', level = 2)

	# 	for (fut in futs) {

	# 		say(fut)

	# 		response_var <- paste0('N_ag_county_', fut)

	# 		map <- map_occurrence_change(
	# 			out_dir = out_dir,
	# 			filename_append = fut,
	# 			pred_vect_nam = pred_vect_nam,
	# 			response_var = response_var,
	# 			data_occs = data_occs,
	# 			title = bquote('Change in ' * italic('Andropogon gerardi') * ' abundance'),
	# 			subtitle = paste0('SSP ', substr(fut, 4, 6), ' ', substr(fut, 8, 11), '-', substr(fut, 13, 16)),
	# 			ag_core_quant = 0.95
	# 		)

	# 	}

	# ### OCCURRENCE: 1930s
	# #####################
	# say('OCCURRENCE: 1930s change maps', level = 2)

	# map_thirties <- map_occurrence_change_1930s(zero_inflated = zero_inflated, pred_vect_nam = pred_vect_nam, pred_vect_conus = pred_vect_conus)


	# # compile residuals analysis
	# resid_p_values <- c(
	# 	moran$estimate[1],
	# 	dharma_resid_test$uniformity$p.value,
	# 	dharma_resid_test$dispersion$p.value,
	# 	dharma_resid_test$outliers$p.value,
	# 	dharma_quant_test$p.value,
	# 	summary(dharma_quant_test$qgamFits[[3]])$s.table[1, 4],
	# 	summary(dharma_quant_test$qgamFits[[2]])$s.table[1, 4],
	# 	summary(dharma_quant_test$qgamFits[[1]])$s.table[1, 4]
	# )
	# resid_p_values_sig <- ifelse(resid_p_values < 0.05, '*', 'ns')
	# resid_test_statistic_values <- c(
	# 	moran$statistic,
	# 	dharma_resid_test$uniformity$statistic,
	# 	dharma_resid_test$dispersion$statistic,
	# 	dharma_resid_test$outliers$statistic,
	# 	NA,
	# 	summary(dharma_quant_test$qgamFits[[3]])$s.table[1, 3],
	# 	summary(dharma_quant_test$qgamFits[[2]])$s.table[1, 3],
	# 	summary(dharma_quant_test$qgamFits[[1]])$s.table[1, 3]
	# )

	# meta_occs <- list(
	# 	facet = 'occurrence',
	# 	date = date(),
	# 	homoscedastic = homoscedastic,
	# 	zero_inflated = zero_inflated,
	# 	formulae = list(
	# 		formula_occs = formula_occs,
	# 		formula_bias = formula_bias
	# 	),
	# 	dharma_resids = data.table(
	# 		test = c('spatial autocorrelation', 'uniformity', 'dispersion', 'outliers', 'quantiles, overall', 'quantiles, upper', 'quantiles, middle', 'quantiles, lower'),
	# 		p_value = resid_p_values,
	# 		significant = resid_p_values_sig,
	# 		test_statistic = c('Moran\'s I', names(dharma_resid_test$uniformity$statistic), names(dharma_resid_test$dispersion$statistic), 'exact binomial', NA, rep('chi squared', 3)),
	# 		test_statistic_value = resid_test_statistic_values
	# 	)
	# )

	# saveRDS(meta_occs, paste0(out_dir, '/!meta_occs.rds'))
	# sink(paste0(out_dir, '/!meta_occs.txt'), split = TRUE)
	# 	print(meta_occs)
	# sink()

}

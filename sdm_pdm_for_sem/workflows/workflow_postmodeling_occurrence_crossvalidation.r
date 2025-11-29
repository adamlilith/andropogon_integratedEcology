#' Post-modeling crossvalidation for occurrence-only models and models with an occurrence component
#'
#' @param homoscedastic If `TRUE`, then do not analyze behavior of sigma
#' @param zero_inflated Logical.
#' @param formula_occs,formula_occs_bias Formulae for occurrences and for occurrence bias
#' @param constants `constants` `list` from model preparation.
#' @param out_dir Folder into which to save results.
#'
workflow_postmodeling_occurrence_crossvalidation <- function(formula_occs, formula_occs_sigma,formula_pzero, formula_occs_bias, constants, out_dir = out_dir) {

	say('OCCURRENCE: cross-validation', level = 2)

	homoscedastic <- is.null(formula_occs_sigma)
	zero_inflated <- !is.null(formula_pzero)

	max_k_folds <- if (trial) { 1 } else { k_folds }
	crossvalidation <- data.table()
	for (k in 1:max_k_folds) {

		say('fold: ', k, level = 3)
	
		fold_constants <- constants
		fold_data <- data
		fold_inits <- inits

		indices <- get(paste0('fold_', k, '_occs'))

		# data
		fold_data$y_n_ag <- fold_data$y_n_ag[-indices]

		# constants
		n_counties_test <- length(indices)
		fold_constants$n_counties_occs_calib <- fold_constants$n_counties_occs_calib - n_counties_test
		fold_constants$counties_x_occs_calib_sq <- fold_constants$counties_x_occs_calib_sq[-indices, ]
		if (zero_inflated) fold_constants$counties_x_occs_pzero_calib_sq <- fold_constants$counties_x_occs_pzero_calib_sq[-indices, ]
		fold_constants$w_occs_bias <- fold_constants$w_occs_bias[-indices, ]
	
		test_x <- constants$counties_x_occs_calib_sq[indices, ]
		if (zero_inflated) test_x_pzero <- constants$counties_x_occs_pzero_calib_sq[indices, ]

		# inits
		fold_inits$y_n_ag_sim <- fold_inits$y_n_ag_sim[-indices]
		fold_inits$N <- fold_inits$N[-indices]

		fold_model <- nimbleModel(
			code = model_code,
			constants = fold_constants,
			data = fold_data,
			inits = fold_inits,
			check = FALSE,
			calculate = FALSE,
			buildDerivs = TRUE
		)

		fold_monitors <- if (homoscedastic) {
			c('beta_occs_mu', 'lambda_sigma')
		} else {

		}
		if (zero_inflated) fold_monitors <- c(fold_monitors, 'beta_pzero')

		fold_conf <- configureMCMC(
			fold_model,
			monitors = fold_monitors,
			print = FALSE,
			enableWAIC = FALSE
		)

		# vars <- c('alpha_occs', 'beta_occs_mu')
		# if (!homoscedastic) vars <- c(vars, 'beta_occs_sigma')
		# if (zero_inflated) vars <- c(vars, 'beta_pzero')
		# fold_conf$addSampler(target = vars, type = 'NUTS')

		# # AF_slice samplers for spearmanated parameters
		# vars <- c(monitors_coeffs_not_indexed, monitors_coeffs_single_index, monitors_coeffs_double_index)
		# for (var in vars) {
		# 	conf$removeSamplers(var)
		# }
		# conf$addSampler(target = vars, type = 'AF_slice')

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

		# evaluate predictions
		preds <- predict_occs(chains = fold_chains, x = test_x, homoscedastic = homoscedastic, zero_inflated = zero_inflated, type = 'mu', x_pzero = test_x_pzero)
		if (zero_inflated) preds_pzero <- predict_occs(chains = fold_chains, x = test_x_pzero, homoscedastic = homoscedastic, zero_inflated = zero_inflated, type = 'pzero', x_pzero = NULL)
		
		test_n <- data$y_n_ag[indices]
		test_binary <- as.numeric(test_n > 0)

		pres_indices <- which(test_binary == 1)
		bg_indices <- which(test_binary == 0)

		preds_pres <- preds[ , pres_indices]
		preds_bg <- preds[ , bg_indices]

		# calculate accuracy statistics
		cbi <- mae <- auc <- biserial <- spearman <- pzero_acc <- rep(NA_real_, nrow(preds))
		for (iter in 1:nrow(preds)) {

			this_preds <- preds[iter, , drop = TRUE]
			if (zero_inflated) this_preds_pzero <- preds_pzero[iter, , drop = TRUE]

			biserial[iter] <- cor(this_preds, test_binary)
			auc[iter] <- enmSdmX::evalAUC(preds_pres[iter, , drop = TRUE], preds_bg[iter, , drop = TRUE])
			cbi[iter] <- enmSdmX::evalContBoyce(this_preds, preds_bg[iter, , drop = TRUE])
			mae[iter] <- mae_fx(this_preds, test_n)
			spearman[iter] <- cor(this_preds, test_n, method = 'spearman')
			if (zero_inflated) pzero_acc[iter] <- cor(-1 * test_binary, this_preds_pzero)
		
		}

		crossvalidation <- rbind(
			crossvalidation,
			data.table(
				k = k,
				n_train = fold_constants$n_counties_occs_calib,
				n_test = length(indices),

				n_test_pres = ncol(preds_pres),
				n_test_bg = ncol(preds_bg),

				point_biserial_lower = quantile(biserial, 0.025, na.rm = TRUE),
				point_biserial_mean = mean(biserial, na.rm = TRUE),
				point_biserial_upper = quantile(biserial, 0.975, na.rm = TRUE),

				auc_lower = quantile(auc, 0.25, na.rm = TRUE),
				auc_mean = mean(auc, na.rm = TRUE),
				auc_upper = quantile(auc, 0.975, na.rm = TRUE),

				cbi_lower = quantile(cbi, 0.025, na.rm = TRUE),
				cbi_mean = mean(cbi, na.rm = TRUE),
				cbi_upper = quantile(cbi, 0.975, na.rm = TRUE),

				spearman_lower = quantile(spearman, 0.025, na.rm = TRUE),
				spearman_mean = mean(spearman, na.rm = TRUE),
				spearman_upper = quantile(spearman, 0.975, na.rm = TRUE),

				mae_lower = quantile(mae, 0.025, na.rm = TRUE),
				mae_mean = mean(mae, na.rm = TRUE),
				mae_upper = quantile(mae, 0.975, na.rm = TRUE),

				pzero_acc_lower = quantile(pzero_acc, 0.025, na.rm = TRUE),
				pzero_acc_mean = mean(pzero_acc, na.rm = TRUE),
				pzero_acc_upper = quantile(pzero_acc, 0.975, na.rm = TRUE)


			)
		)

	}

	crossvalidation <- rbind(
		crossvalidation,
		data.table(
				k = 'means',
				n_train = mean(crossvalidation$n_train),
				n_test = mean(crossvalidation$n_test),

				n_test_pres = mean(crossvalidation$n_test_pres),
				n_test_bg = mean(crossvalidation$n_test_bg),

				point_biserial_lower = mean(crossvalidation$point_biserial_lower),
				point_biserial_mean = mean(crossvalidation$point_biserial_mean),
				point_biserial_upper = mean(crossvalidation$point_biserial_upper),

				auc_lower = mean(crossvalidation$auc_lower),
				auc_mean = mean(crossvalidation$auc_mean),
				auc_upper = mean(crossvalidation$auc_upper),

				cbi_lower = mean(crossvalidation$cbi_lower),
				cbi_mean = mean(crossvalidation$cbi_mean),
				cbi_upper = mean(crossvalidation$cbi_upper),

				spearman_lower = mean(crossvalidation$spearman_lower),
				spearman_mean = mean(crossvalidation$spearman_mean),
				spearman_upper = mean(crossvalidation$spearman_upper),

				mae_lower = mean(crossvalidation$mae_lower),
				mae_mean = mean(crossvalidation$mae_mean),
				mae_upper = mean(crossvalidation$mae_upper),

				pzero_acc_lower = mean(crossvalidation$pzero_acc_lower),
				pzero_acc_mean = mean(crossvalidation$pzero_acc_mean),
				pzero_acc_upper = mean(crossvalidation$pzero_acc_upper)
		)
	)

	### nimble crossvalidation method
	# if (!trial) {

	# 	loss_fx <- mae_fx

	# 	say('Cross-validation loss function:')
	# 	print(loss_fx)
	# 	say('')

	# 	folds_fx <- folds_for_occurrences # change according to response data type we're using

	# 	cv <- runCrossValidate(
	# 		MCMCconfiguration = conf,
	# 		k = k_folds, # universal setting
	# 		foldFunction = folds_fx,
	# 		lossFunction = loss_fx,
	# 		MCMCcontrol = list(niter = niter, nburnin = nburnin),
	# 		returnSamples = FALSE,
	# 		nCores = 1,
	# 		nBootReps = 200,
	# 		silent = FALSE
	# 	)

	# 	mean_form <- paste(as.character(formula_occs), collapse = ' ')
	# 	mean_form <- gsub(mean_form, pattern = 'I\\(', replacement = '')
	# 	mean_form <- gsub(mean_form, pattern = '\\^2\\)', replacement = '²')
	# 	mean_form <- gsub(mean_form, pattern = '*)', replacement = '×')
	# 	form <- paste0(trait, ' ', mean_form)

	# 	crossvalidation <- data.table(
	# 		model = form,
	# 		k = 'summary',
	# 		cv_value = cv$CVvalue,
	# 		cv_value_se = cv$CVstandardError
	# 	)

	# 	for (k in 1:k_folds) {
		
	# 		crossvalidation <- rbind(
	# 			crossvalidation,
	# 			data.table(
	# 				model = form,
	# 				k = k,
	# 				cv_value = cv$foldCVinfo[[k]][1],
	# 				cv_value_se = cv$foldCVinfo[[k]][2]
	# 			)
	# 		)

	# 	}

	# 	# sink(paste0(out_dir, '/cross_validation_occurrence.txt'), split = TRUE)
	# 	# say('GEO-FOLD CROSS VALIDATION')
	# 	# say(date(), post = 2)
	# 	# say('Number of geo-folds: ', k_folds, post = 2)
	# 	# print(cv)
	# 	# sink()

	# }


	meta_crossvalidation <- list(
		facet = 'occurrence',
		date = date(),
		homoscedastic = homoscedastic,
		zero_inflated = zero_inflated,
		formulae = list(
			formula_occs = formula_occs,
			formula_occs_bias = formula_occs_bias
		),
		crossvalidation = crossvalidation
	)

	saveRDS(meta_crossvalidation, paste0(out_dir, '/!meta_crossvalidation.rds'))
	sink(paste0(out_dir, '/!meta_crossvalidation.txt'), split = TRUE)
		print(meta_crossvalidation)
	sink()

}

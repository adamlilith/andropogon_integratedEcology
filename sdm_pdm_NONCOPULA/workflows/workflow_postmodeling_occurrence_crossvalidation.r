#' Post-modeling crossvalidation for occurrence-only models and models with an occurrence component
#'
#' formula_occs, formula_psi, formula_occs_bias Formulae for occurrences, probability of presence, and for occurrence bias
#' zero_inflated Logical.
#' constants `constants` `list` from model preparation.
#' out_dir Folder into which to save results.
#'
workflow_postmodeling_occurrence_crossvalidation <- function(formula_occs, formula_psi = NULL, formula_occs_bias, constants, out_dir) {

	say('OCCURRENCE: cross-validation', level = 2)

	zero_inflated <- !is.null(formula_psi)

	ag_vect <- vect('./outputs_loretta/integrated_sdm_pdm/andropogon_gerardi_occurrences_with_environment_1961_2020_for_integration.gpkg')

	max_k_folds <- if (trial) { 1 } else { k_folds }
	crossvalidation <- data.table()
	for (k in 1:max_k_folds) {

		say('fold: ', k, level = 3)
	
		fold_constants <- constants
		fold_data <- data
		fold_inits <- inits

		indices <- which(data_occs$ag_vect$geofold == k)

		# data
		fold_data$y_n_ag <- fold_data$y_n_ag[-indices]

		# constants
		n_counties_test <- length(indices)
		fold_constants$n_counties_occs_calib <- fold_constants$n_counties_occs_calib - n_counties_test
		fold_constants$counties_x_occs_calib_sq <- fold_constants$counties_x_occs_calib_sq[-indices, ]
		if (zero_inflated) fold_constants$counties_x_occs_psi_calib_sq <- fold_constants$counties_x_occs_psi_calib_sq[-indices, ]
		fold_constants$w_occs_bias <- fold_constants$w_occs_bias[-indices, ]
	
		test_x <- constants$counties_x_occs_calib_sq[indices, ]
		if (zero_inflated) test_x_psi <- constants$counties_x_occs_psi_calib_sq[indices, ]

		# inits
		fold_inits$y_n_ag_sim <- fold_inits$y_n_ag_sim[-indices]
		fold_inits$N <- fold_inits$N[-indices]
		fold_inits$beta_occs <- mc_extract(chains, 'beta_occs', j = TRUE)
		if (zero_inflated) {
			fold_inits$beta_psi <- mc_extract(chains, 'beta_psi', j = TRUE)
			fold_inits$z_county <- fold_inits$z_county[-indices]
		}

		fold_model <- nimbleModel(
			code = model_code,
			constants = fold_constants,
			data = fold_data,
			inits = fold_inits,
			check = FALSE,
			calculate = FALSE,
			buildDerivs = TRUE
		)

		fold_monitors <- c('beta_occs', 'lambda_sigma')
		if (zero_inflated) fold_monitors <- c(fold_monitors, 'beta_psi')

		fold_conf <- configureMCMC(
			fold_model,
			monitors = fold_monitors,
			print = FALSE,
			enableWAIC = FALSE
		)

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
		preds <- predict_occs(chains = fold_chains, x = test_x, zero_inflated = zero_inflated, x_psi = test_x_psi)
		if (zero_inflated) preds_psi <- predict_psi(chains = fold_chains, x = test_x_psi)
		
		test_n <- data$y_n_ag[indices]
		test_binary <- as.numeric(test_n > 0)

		pres_indices <- which(test_binary == 1)
		bg_indices <- which(test_binary == 0)

		preds_pres <- preds[ , pres_indices]
		preds_bg <- preds[ , bg_indices]

		# calculate accuracy statistics
		cbi <- mae <- auc <- biserial <- spearman <- psi_acc <- rep(NA_real_, nrow(preds))
		for (iter in 1:nrow(preds)) {

			this_preds <- preds[iter, , drop = TRUE]
			if (zero_inflated) this_preds_psi <- preds_psi[iter, , drop = TRUE]

			biserial[iter] <- cor(this_preds, test_binary)
			auc[iter] <- enmSdmX::evalAUC(preds_pres[iter, , drop = TRUE], preds_bg[iter, , drop = TRUE])
			cbi[iter] <- enmSdmX::evalContBoyce(this_preds, preds_bg[iter, , drop = TRUE])
			mae[iter] <- mae_fx(this_preds, test_n)
			spearman[iter] <- cor(this_preds, test_n, method = 'spearman')
			if (zero_inflated) psi_acc[iter] <- cor(test_binary, this_preds_psi)
		
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

				psi_acc_lower = quantile(psi_acc, 0.025, na.rm = TRUE),
				psi_acc_mean = mean(psi_acc, na.rm = TRUE),
				psi_acc_upper = quantile(psi_acc, 0.975, na.rm = TRUE)


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

				psi_acc_lower = mean(crossvalidation$psi_acc_lower),
				psi_acc_mean = mean(crossvalidation$psi_acc_mean),
				psi_acc_upper = mean(crossvalidation$psi_acc_upper)
		)
	)

	# remember
	meta_crossvalidation <- list(
		facet = 'occurrence',
		date = date(),
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

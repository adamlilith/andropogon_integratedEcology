#' Post-modeling cv_occs for occurrence-only models and models with an occurrence component
#'
#' @param constants `constants` `list` from model preparation.
#' @param inits `inits` `list` from model preparation.
#' @param formula_occs,formula_occs_bias Formulae for occurrences and for occurrence bias
#' @param formula_biomass_mu,formula_biomass_sigma Biomass formulae.
#' @param formula_pzero Formula for probability of zero abundance and biomass.
#' @param out_dir Folder into which to save results.
#'
workflow_postmodeling_occurrence_biomass_crossvalidation <- function(chains, constants, inits, formula_occs, formula_occs_sigma, formula_occs_bias, formula_biomass_mu, formula_biomass_sigma, formula_pzero, out_dir) {

	say('OCCURRENCE + BIOMASS: cross-validation', level = 2)

	homoscedastic <- is.null(formula_occs_sigma)
	zero_inflated <- !is.null(formula_pzero)

	max_k_folds <- if (trial) { 1 } else { k_folds }
	cv_occs <- cv_biomass <- data.table()
	for (k in 1:max_k_folds) {

		say('fold: ', k, level = 3)
	
		fold_constants <- constants
		fold_data <- data
		fold_inits <- inits

		# county test/train sites
		test_counties <- get(paste0('fold_', k, '_occs'))

		# biomass test/train sites
		test_plants <- get(paste0('fold_', k, '_biomass'))
		test_sites <- sort(unique(fold_constants$site_index_biomass[test_plants]))
		train_sites <- sort(unique(fold_constants$site_index_biomass[-test_plants]))

		# data
		fold_data$y_n_ag <- fold_data$y_n_ag[-test_counties]
		fold_data$y_biomass <- fold_data$y_biomass[-test_plants]

		# constants: occurrences
		n_counties_test <- length(test_counties)
		fold_constants$n_counties_occs_calib <- fold_constants$n_counties_occs_calib - n_counties_test
		
		fold_constants$x_by_county_occs <- fold_constants$x_by_county_occs[-test_counties, ]
		fold_constants$x_by_county_biomass <- fold_constants$x_by_county_biomass[-test_counties, ]
		
		fold_constants$x_by_site_occs <- fold_constants$x_by_site_occs[train_sites, ]
		fold_constants$x_by_site_biomass <- fold_constants$x_by_site_biomass[train_sites, ]
		
		if (zero_inflated) fold_constants$counties_x_occs_pzero_calib_sq <- fold_constants$counties_x_occs_pzero_calib_sq[-test_counties, ]
		fold_constants$w_occs_bias <- fold_constants$w_occs_bias[-test_counties, ]

		# constants: biomass... order matters in assigning these!
		fold_constants$site_index_biomass <- renumSeq(fold_constants$site_index_biomass[fold_constants$site_index_biomass %in% train_sites])
		fold_constants$n_biomass <- fold_constants$n_biomass - length(test_plants)
		fold_constants$n_pheno_sites <- length(train_sites)

		# test data
		x_by_county_occs_test <- constants$x_by_county_occs[test_counties, ]
		x_by_county_biomass_test <- constants$x_by_county_biomass[test_counties, ]

		x_by_site_occs_test <- constants$x_by_site_occs[test_sites, ]
		x_by_site_biomass_test <- constants$x_by_site_biomass[test_sites, ]

		if (zero_inflated) {
			test_x_pzero <- constants$counties_x_occs_pzero_calib_sq[test_counties, ]
		} else {
			test_x_pzero <- NULL
		}

		# inits: occurrences
		fold_inits$lambda_mu_sq <- inits$lambda_mu_sq[-test_counties]
		fold_inits$alpha_occs <- hammer_extract(chains, 'alpha_occs', j = TRUE)
		fold_inits$beta_occs_mu <- hammer_extract(chains, 'beta_occs_mu', j = TRUE)
		
		fold_inits$beta_biomass_mu <- hammer_extract(chains, 'beta_biomass_mu', j = TRUE)

		fold_inits$y_n_ag_sim <- fold_inits$y_n_ag_sim[-test_counties]
		fold_inits$N <- fold_inits$N[-test_counties]

		# inits: biomass
		fold_inits$mu_biomass_site <- hammer_extract(chains, 'mu_biomass_site', j = TRUE)[-test_sites]
		fold_inits$sigma_biomass_within_sites_log <- log(hammer_extract(chains, 'sigma_biomass_within_sites'))
		fold_inits$y_biomass_sim <- fold_inits$y_biomass_sim[-test_plants]

		# inits: integration
		fold_inits$Phi_county <- inits$Phi_county[-test_counties, ]
		fold_inits$Phi_site <- inits$Phi_site[train_sites, ]

		fold_model <- nimbleModel(
			code = model_code,
			constants = fold_constants,
			data = fold_data,
			inits = fold_inits,
			check = FALSE,
			calculate = FALSE,
			buildDerivs = TRUE
		)

		fold_monitors <- c('beta_occs_mu', 'U', 'beta_biomass_mu', 'mu_biomass_site', 'sigma_biomass_within_sites')
		if (!homoscedastic)	fold_monitors <- c(fold_monitors, 'beta_occs_sigma', 'beta_biomass_sigma')
		if (zero_inflated) fold_monitors <- c(fold_monitors, 'beta_pzero')

		fold_conf <- configureMCMC(
			fold_model,
			monitors = fold_monitors,
			print = FALSE,
			enableWAIC = FALSE
		)

		# add no U-turn sampler (Hamiltonian Monte Carlo)
		vars <- 'eta'
		fold_conf$addSampler(target = vars, type = 'NUTS')
		say('NUTS sampler added to ', paste(vars, collapse = ' & '), '.')

		vars <- 'sigma_biomass_within_sites_log'
		fold_conf$addSampler(target = vars, type = 'NUTS')
		say('NUTS sampler added to ', paste(vars, collapse = ' & '), '.')

		vars <- 'alpha_occs'
		fold_conf$addSampler(target = vars, type = 'NUTS')
		say('NUTS sampler added to ', paste(vars, collapse = ' & '), '.')

		var <- 'beta_occs_mu'
		fold_conf$removeSamplers(var)
		fold_conf$addSampler(target = var, type = 'AF_slice')
		say('AF_slice sampler added to ', paste(var, collapse = ' & '), '.')

		var <- 'beta_biomass_mu'	
		fold_conf$removeSamplers(var)
		fold_conf$addSampler(target = var, type = 'AF_slice')
		say('AF_slice sampler added to ', paste(var, collapse = ' & '), '.')

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

		### OCCURRENCE: evaluate predictions
		say('evaluate occurrence predictions')
		preds <- predict_occs_biomass(chains = fold_chains, x_occs = x_by_county_occs_test, x_biomass = x_by_county_biomass_test, homoscedastic = homoscedastic, zero_inflated = zero_inflated, type = 'mu', x_pzero = test_x_pzero)

		preds <- preds$preds_occs

		# if (zero_inflated) preds_pzero <- predict_occs(chains = fold_chains, x = test_x_pzero, homoscedastic = homoscedastic, zero_inflated = zero_inflated, type = 'pzero', x_pzero = NULL)
		
		test_n <- data$y_n_ag[test_counties]
		test_binary <- as.numeric(test_n > 0)

		pres_test_counties <- which(test_binary == 1)
		bg_test_counties <- which(test_binary == 0)

		preds_pres <- preds[ , pres_test_counties]
		preds_bg <- preds[ , bg_test_counties]

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

		cv_occs <- rbind(
			cv_occs,
			data.table(
				k = k,
				n_train = fold_constants$n_counties_occs_calib,
				n_test = length(test_counties),

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

		### BIOMASS: evaluate predictions
		say('evaluate biomass predictions')
		preds <- predict_occs_biomass(chains = fold_chains, x_occs = x_by_site_occs_test, x_biomass = x_by_site_biomass_test, homoscedastic = homoscedastic, zero_inflated = zero_inflated, type = 'mu', x_pzero = test_x_pzero)
		
		preds <- preds$preds_biomass

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

		cv_biomass <- rbind(
			cv_biomass,
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

	# OCCURRENCE: summary
	cv_occs <- rbind(
		cv_occs,
		data.table(
				k = 'means',
				n_train = mean(cv_occs$n_train),
				n_test = mean(cv_occs$n_test),

				n_test_pres = mean(cv_occs$n_test_pres),
				n_test_bg = mean(cv_occs$n_test_bg),

				point_biserial_lower = mean(cv_occs$point_biserial_lower),
				point_biserial_mean = mean(cv_occs$point_biserial_mean),
				point_biserial_upper = mean(cv_occs$point_biserial_upper),

				auc_lower = mean(cv_occs$auc_lower),
				auc_mean = mean(cv_occs$auc_mean),
				auc_upper = mean(cv_occs$auc_upper),

				cbi_lower = mean(cv_occs$cbi_lower),
				cbi_mean = mean(cv_occs$cbi_mean),
				cbi_upper = mean(cv_occs$cbi_upper),

				spearman_lower = mean(cv_occs$spearman_lower),
				spearman_mean = mean(cv_occs$spearman_mean),
				spearman_upper = mean(cv_occs$spearman_upper),

				mae_lower = mean(cv_occs$mae_lower),
				mae_mean = mean(cv_occs$mae_mean),
				mae_upper = mean(cv_occs$mae_upper),

				pzero_acc_lower = mean(cv_occs$pzero_acc_lower),
				pzero_acc_mean = mean(cv_occs$pzero_acc_mean),
				pzero_acc_upper = mean(cv_occs$pzero_acc_upper)
		)
	)

	# BIOMASS: summary
	cv_biomass <- rbind(
		cv_biomass,
		data.table(
				k = 'means',

				n_train_sites = mean(cv_biomass$n_train_sites),
				n_test_sites = mean(cv_biomass$n_test_sites),

				n_train_plants = mean(cv_biomass$n_train_plants),
				n_test_plants = mean(cv_biomass$n_test_plants),
				
				correl_lower = mean(cv_biomass$correl_lower),
				correl_mean = mean(cv_biomass$correl_mean),
				correl_upper = mean(cv_biomass$correl_upper),

				mae_lower = mean(cv_biomass$mae_lower),
				mae_mean = mean(cv_biomass$mae_mean),
				mae_upper = mean(cv_biomass$mae_upper),

				mape_lower = mean(cv_biomass$mape_lower),
				mape_mean = mean(cv_biomass$mape_mean),
				mape_upper = mean(cv_biomass$mape_upper),

				rmse_lower = mean(cv_biomass$rmse_lower),
				rmse_mean = mean(cv_biomass$rmse_mean),
				rmse_upper = mean(cv_biomass$rmse_upper),

				obs_quants_min = mean(cv_biomass$obs_quants_min),
				obs_quants_mean = mean(cv_biomass$obs_quants_mean),
				obs_quants_max = mean(cv_biomass$obs_quants_max),
				obs_quants_prop_in_inner_90 = mean(cv_biomass$obs_quants_prop_in_inner_90)

		)
	)

	# report
	meta_crossvalidation <- list(
		facet = 'occurrence',
		date = date(),
		homoscedastic = homoscedastic,
		zero_inflated = zero_inflated,
		formulae = list(
			formula_occs = formula_occs,
			formula_occs_bias = formula_occs_bias,
			formula_biomass_mu = formula_biomass_mu,
			formula_biomass_sigma = formula_biomass_sigma,
			formula_pzero = formula_pzero
		),
		cv_occs = cv_occs,
		cv_biomass = cv_biomass
	)

	saveRDS(meta_crossvalidation, paste0(out_dir, '/!meta_crossvalidation.rds'))
	sink(paste0(out_dir, '/!meta_crossvalidation.txt'), split = TRUE)
		print(meta_crossvalidation)
	sink()

}

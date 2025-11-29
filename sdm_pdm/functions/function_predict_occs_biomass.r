#' Predict occurrence from chains
#'
#' @param chains
#' @param x_occs Model matrix
#' @param x_biomass Model matrix
#' @param homoscedastic_occs `TRUE` or `FALSE`
#' @param homoscedastic_biomass `TRUE` or `FALSE`
#' @param zero_inflated `TRUE` or `FALSE`
#' @param type 'mu' or 'sigma' or 'pzero'
#' @param x_pzero Model matrix for probability of zero abundance or `NULL`
#'
#' @returns A matrix of predictions. Rows are iterations and columns are sample IDs.
predict_occs_biomass <- function(chains, x_occs, x_biomass, homoscedastic_occs, homoscedastic_biomass, zero_inflated, type = 'mu', x_pzero = NULL) {

	vars <- paste0('beta_occs_', type)
	betas_occs <- hammer_subset(chains, vars, j = TRUE)

	vars <- paste0('beta_biomass_', type)
	betas_biomass <- hammer_subset(chains, vars, j = TRUE)

	sigma_biomass_within_sites <- hammer_subset(chains, 'sigma_biomass_within_sites')

	if (zero_inflated) {
		vars <- paste0('beta_pzero')
		betas_pzero <- hammer_subset(chains, vars, j = TRUE)
	}

	Us <- hammer_subset(chains, 'U', j = TRUE, k = TRUE)

	if (!homoscedastic_occs | !homoscedastic_biomass) {
		stop('not homoscedastic... need code for this')
		# vars <- paste0('beta_biomass_sigma')
		# betas_sigma <- hammer_subset(chains, vars, j = TRUE)
	}

	n_samples <- nrow(x_occs)
	nchains <- hammer_n_chains(chains)
	chain_samples <- hammer_samples(chains)
	iters_per_chain <- nrow(chain_samples[[1]])
	total_iters <- nchains * iters_per_chain

	preds_occs <- preds_biomass <- matrix(NA, nrow = total_iters, ncol = n_samples)

	k <- 1
	for (chain in 1:nchains) {

		for (iter in 1:iters_per_chain) {

			this_beta_occs <- betas_occs$samples[[chain]][iter, ]
			this_beta_occs <- cbind(this_beta_occs)

			this_beta_biomass <- betas_biomass$samples[[chain]][iter, ]
			this_beta_biomass <- cbind(this_beta_biomass)

			this_sigma_biomass_within_sites <- sigma_biomass_within_sites$samples[[chain]][iter, 'sigma_biomass_within_sites']

			# if (zero_inflated) {
			# 	this_beta_pzero <- betas_pzero$samples[[chain]][iter, ]
			# 	this_beta_pzero <- cbind(this_beta_pzero)
			# }

			# predicting mean, homoscedastic, not zero-inflated
			if (type == 'mu' & homoscedastic_occs & homoscedastic_biomass & !zero_inflated) {

				# occurrences
				phi_lambda_mu <- x_occs %*% this_beta_occs
				phi_lambda_mu <- phi_lambda_mu[ , 1]

				# biomass
				phi_biomass_mu <- x_biomass %*% this_beta_biomass
				phi_biomass_mu <- phi_biomass_mu[ , 1]

				# VCV
				U <- Us$samples[[chain]][iter, ]
				U <- matrix(U, ncol = 2)

				# multivariate simulation
				Phis <- matrix(NA, ncol = 2, nrow = n_samples)
				for (i in seq_len(n_samples)) {
					phis <- c(phi_lambda_mu[i], phi_biomass_mu[i])
					Phis[i, ] <- rmnorm_chol(1, mean = phis, cholesky = U, prec_param = FALSE)
				}

				# occurrence
				lambda <- exp(Phis[ , 1])
				pred_occs <- rpois(n_samples, lambda)

				# biomass
				mu_biomass <- exp(Phis[ , 2])
				shape_biomass <- mu_biomass^2 / this_sigma_biomass_within_sites^2
				rate_biomass <- mu_biomass / this_sigma_biomass_within_sites^2

				pred_biomass <- rep(NA_real_, length(mu_biomass))
				for (count in seq_along(pred_biomass)) pred_biomass[count] <- rgamma(1, shape = shape_biomass[count], rate = rate_biomass[count])

			# predicting mean, homoscedastic, zero-inflated
			} else if (type == 'mu' & homoscedastic_occs & homoscedastic_biomass & zero_inflated) {

				# zero inflation
				this_beta_pzero <- betas_pzero$samples[[chain]][iter, ]
				this_beta_pzero <- cbind(this_beta_pzero)
				pzero <- x_pzero %*% this_beta_pzero
				pzero <- pzero[  , 1]
				pzero <- expit(pzero)

				# occurrences
				phi_lambda_mu <- x_occs %*% this_beta_occs
				phi_lambda_mu <- phi_lambda_mu[ , 1]

				# biomass
				phi_biomass_mu <- x_biomass %*% this_beta_biomass
				phi_biomass_mu <- phi_biomass_mu[ , 1]

				this_sigma_biomass_within_sites <- sigma_biomass_within_sites$samples[[chain]][iter, 'sigma_biomass_within_sites']

				# VCV
				U <- Us$samples[[chain]][iter, ]
				U <- matrix(U, ncol = 2)

				# multivariate simulation
				Phis <- matrix(NA, ncol = 2, nrow = n_samples)
				for (i in seq_len(n_samples)) {
					phis <- c(phi_lambda_mu[i], phi_biomass_mu[i])
					Phis[i, ] <- rmnorm_chol(1, mean = phis, cholesky = U, prec_param = FALSE)
				}

				# occurrence
				lambda <- exp(Phis[ , 1])
				pred_occs <- rep(NA_real_, length(lambda))
				for (count in seq_along(lambda)) pred_occs[count] <- rzip(1, lambda = lambda[count], pzero = pzero[count])

				# biomass
				mu_biomass <- exp(Phis[ , 2])
				shape_biomass <- mu_biomass^2 / this_sigma_biomass_within_sites^2
				rate_biomass <- mu_biomass / this_sigma_biomass_within_sites^2

				pred_biomass <- rep(NA_real_, length(exp_Phi_biomass))
				# for (count in seq_along(exp_Phi_biomass)) pred_biomass[count] <- rzigamma(1, shape = shape_biomass, rate = rate_biomass, pzero = pzero[count]) * (pred_occs[count] > 0)
				for (count in seq_along(exp_Phi_biomass)) pred_biomass[count] <- rzigamma(1, shape = shape_biomass, rate = rate_biomass, pzero = pzero[count])

			# predicting mean, not homoscedastic, not zero-inflated
			} else if (type == 'mu' & !homoscedastic_occs & !homoscedastic_biomass & !zero_inflated) {
				stop('Need to code prediction of heteroscedastic occs and biomass.')
				# this_beta_sigma <- betas_sigma$samples[[chain]][iter, ]
				# this_beta_sigma <- cbind(this_beta_sigma)

				# pred_untrans_sigma <- x %*% this_beta_sigma
				# pred_untrans_sigma <- pred_untrans_sigma[ , 1]

				# this_sigma_biomass_among_sites <- exp(pred_untrans_sigma)
				# pred <- rnorm(n_samples, mean = pred_untrans_sigma, sd = this_sigma_biomass_among_sites)
				# pred <- exp(pred)
			
			# predicting pzero
			} else if (type == 'pzero') {

				stop('Need to code prediction of pzero')
				# preds_untrans <- x %*% this_beta
				# pred <- expit(preds_untrans)

			# predicting sigma
			} else if (type == 'sigma') {
				stop('Need to code prediction of sigma')
				# pred <- exp(pred_untrans)			
			}
			preds_occs[k, ] <- pred_occs
			preds_biomass[k, ] <- pred_biomass

			k <- k + 1

		}
		
	}
	list(preds_occs = preds_occs, preds_biomass = preds_biomass)

}


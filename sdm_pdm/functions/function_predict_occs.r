#' Predict occurrence from chains
#'
#' @param chains
#' @param x Model matrix
#' @param homoscedastic `TRUE` or `FALSE`
#' @param zero_inflated `TRUE` or `FALSE`
#' @param type 'mu' or 'sigma'
#' @param x_pzero Model matrix for probability of zero abundance or `NULL`
#'
#' @returns A matrix of predictions. Rows are iterations and columns are sample IDs.
predict_occs <- function(chains, x, homoscedastic, zero_inflated, type = 'mu', x_pzero = NULL) {

	if (type %in% c('mu', 'sigma')) {
		vars <- paste0('beta_occs_', type)
	} else {
		vars <- 'beta_pzero'
	}
	betas <- hammer_subset(chains, vars, j = TRUE)

	if (homoscedastic) {
		lambda_sigma <- hammer_subset(chains, 'lambda_sigma')
	} else {
		# vars <- paste0('beta_biomass_sigma')
		# betas_sigma <- hammer_subset(chains, vars, j = TRUE)
	}
	if (zero_inflated) betas_pzero <- hammer_subset(chains, 'beta_pzero', j = TRUE)

	n_samples <- nrow(x)
	nchains <- hammer_n_chains(chains)
	chain_samples <- hammer_samples(chains)
	iters_per_chain <- nrow(chain_samples[[1]])
	total_iters <- nchains * iters_per_chain

	preds <- matrix(NA, nrow = total_iters, ncol = n_samples)

	k <- 1
	for (chain in 1:nchains) {

		for (iter in 1:iters_per_chain) {

			this_beta <- betas$samples[[chain]][iter, ]
			this_beta <- cbind(this_beta)

			if (zero_inflated) {
				this_beta_pzero <- betas_pzero$samples[[chain]][iter, ]
				this_beta_pzero <- cbind(this_beta_pzero)
			}

			# predicting mean, homoscedastic, not zero-inflated
			if (type == 'mu' & homoscedastic & !zero_inflated) {

				this_lambda_sigma <- lambda_sigma$samples[[chain]][iter, ]

				phi_lambda_mu <- x %*% this_beta
				phi_lambda_mu <- phi_lambda_mu[ , 1]

				log_lambda <- rnorm(n_samples, phi_lambda_mu, sd = this_lambda_sigma)
				lambda <- exp(log_lambda)
				pred <- rpois(n_samples, lambda)

			# predicting mean, homoscedastic, zero-inflated
			} else if (type == 'mu' & homoscedastic & zero_inflated) {

				this_lambda_sigma <- lambda_sigma$samples[[chain]][iter, ]

				phi_lambda_mu <- x %*% this_beta
				phi_lambda_mu <- phi_lambda_mu[ , 1]

				log_lambda <- rnorm(n_samples, phi_lambda_mu, sd = this_lambda_sigma)
				lambda <- exp(log_lambda)

				pzero <- x_pzero %*% this_beta_pzero
				pzero <- pzero[ , 1]
				pzero <- expit(pzero)

				pred <- rep(NA_real_, n_samples)
				for (n in 1:n_samples) pred[n] <- rzip(1, lambda[n], pzero = pzero[n])

			# predicting mean, not homoscedastic, not zero-inflated
			} else if (type == 'mu' & !homoscedastic & !zero_inflated) {
				
				# this_beta_sigma <- betas_sigma$samples[[chain]][iter, ]
				# this_beta_sigma <- cbind(this_beta_sigma)

				# pred_untrans_sigma <- x %*% this_beta_sigma
				# pred_untrans_sigma <- pred_untrans_sigma[ , 1]

				# this_sigma_biomass_among_sites <- exp(pred_untrans_sigma)
				# pred <- rnorm(n_samples, mean = pred_untrans_sigma, sd = this_sigma_biomass_among_sites)
				# pred <- exp(pred)
			
			# predicting pzero
			} else if (type == 'pzero') {

				preds_untrans <- x %*% this_beta
				preds_untrans <- preds_untrans[ , 1]
				pred <- expit(preds_untrans)

			# predicting sigma
			} else if (type == 'sigma') {
				# pred <- exp(pred_untrans)			
			}
			preds[k, ] <- pred

			k <- k + 1

		}
		
	}
	preds

}


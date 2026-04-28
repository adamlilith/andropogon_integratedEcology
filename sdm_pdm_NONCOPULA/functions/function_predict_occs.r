#' Predict occurrence from chains
#'
#' source('C:/Kaji/R/andropogon_integratedEcology/sdm_pdm/functions/function_predict_occs.r')
#' 
#' chains
#' x 				Model matrix
#' zero_inflated 	`TRUE` or `FALSE`
#' x_psi 			Model matrix for probability of zero abundance or `NULL`
#' times			Number of times to make predictions. Since result relies on RNG, may want multiple predictions per datum.
#'
#' Returns a matrix of predictions. Rows are iterations and columns are sample IDs. Will have iterations * times rows.
predict_occs <- function(chains, x, zero_inflated, x_psi = NULL, times = 1) {

	betas <- mc_subset(chains, 'beta_occs', j = TRUE)
	lambda_sigma <- mc_subset(chains, 'lambda_sigma')
	if (zero_inflated) betas_psi <- mc_subset(chains, 'beta_psi', j = TRUE)

	n_samples <- nrow(x)
	nchains <- mc_n_chains(chains)
	chain_samples <- mc_samples(chains)
	iters_per_chain <- nrow(chain_samples[[1]])
	total_iters <- nchains * iters_per_chain

	preds <- matrix(NA, nrow = times * total_iters, ncol = n_samples)

	k <- 1

	for (time in 1:times) {
			
		for (chain in 1:nchains) {

			for (iter in 1:iters_per_chain) {

				this_beta <- betas$samples[[chain]][iter, ]
				this_beta <- cbind(this_beta)

				if (zero_inflated) {
					this_beta_psi <- betas_psi$samples[[chain]][iter, ]
					this_beta_psi <- cbind(this_beta_psi)
				}

				this_lambda_sigma <- lambda_sigma$samples[[chain]][iter, ]
				log_lambda <- x %*% this_beta
				log_lambda <- log_lambda[ , 1]

				lambda <- exp(log_lambda)

				# predicting mean, not zero-inflated
				if (!zero_inflated) {

					pred <- rpois(n_samples, lambda)

				# predicting mean, homoscedastic, zero-inflated
				} else if (zero_inflated) {

					psi <- x_psi %*% this_beta_psi
					psi <- psi[ , 1]
					psi <- expit(psi)

					z <- as.numeric(runif(n_samples) < psi)

					pred <- rep(NA_real_, n_samples)
					for (n in 1:n_samples) pred[n] <- rTruncPseudoZIP(1, lambda[n], z = z[n])

				}

				preds[k, ] <- pred
				k <- k + 1

			}
			
		}

	}
	preds

}


#' Predict occurrence from chains
#'
#' source('C:/Kaji/R/andropogon_integratedEcology/sdm_pdm/functions/function_predict_occs.r')
#' 
#' chains
#' x 				Model matrix
#' x_psi 			Model matrix for probability of zero abundance or `NULL`
#'
#' Returns a matrix of predictions. Rows are iterations and columns are sample IDs.
predict_occs <- function(chains, x, x_psi = NULL) {

	betas <- mc_subset(chains, 'beta_occs', j = TRUE)
	# lambda_sigma <- mc_subset(chains, 'lambda_sigma')
	if (!is.null(x_psi)) betas_psi <- mc_subset(chains, 'beta_psi', j = TRUE)

	n_samples <- nrow(x)
	nchains <- mc_n_chains(chains)
	chain_samples <- mc_samples(chains)
	iters_per_chain <- nrow(chain_samples[[1]])
	total_iters <- nchains * iters_per_chain

	preds <- matrix(NA, nrow = total_iters, ncol = n_samples)

	k <- 1

	for (chain in 1:nchains) {

		for (iter in 1:iters_per_chain) {

			this_beta <- betas$samples[[chain]][iter, ]
			this_beta <- cbind(this_beta)

			# this_lambda_sigma <- lambda_sigma$samples[[chain]][iter, ]
			# phi_lambda_mu <- x %*% this_beta
			# phi_lambda_mu <- phi_lambda_mu[ , 1]
			log_lambda <- x %*% this_beta
			log_lambda <- log_lambda[ , 1]
			# log_lambda <- rnorm(n_samples, phi_lambda_mu, sd = this_lambda_sigma)
			lambda <- exp(log_lambda)

			# predicting zero-UNinflated
			if (is.null(x_psi)) {

				pred <- rpois(n_samples, lambda)

			# predicting zero-inflated
			} else if (!is.null(x_psi)) {

				# zero-inflation
				this_beta_psi <- betas_psi$samples[[chain]][iter, ]
				this_beta_psi <- cbind(this_beta_psi)

				psi <- x_psi %*% this_beta_psi
				psi <- psi[ , 1]
				psi <- expit(psi)

				# z <- as.numeric(runif(n_samples) < psi)

				# predict abundance
				pred <- rep(NA_integer_, n_samples)
				for (i in 1:n_samples) pred[i] <- rHurdlePoisson(1, lambda[i], psi = psi[i])

			}

			preds[k, ] <- pred
			k <- k + 1

		}
		
	}
	preds

}

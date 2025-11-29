#' Predict probability of zero inflation from chains
#'
#' @param chains
#' @param x Model matrix
#'
#' @returns A matrix of predictions. Rows are iterations and columns are sample IDs.
predict_pzero <- function(chains, x) {

	betas <- hammer_subset(chains, 'beta_pzero', j = TRUE)

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

			preds_untrans <- x %*% this_beta
			preds_untrans <- preds_untrans[ , 1]
			pred <- expit(preds_untrans)

			preds[k, ] <- pred

			k <- k + 1

		}
		
	}
	preds

}


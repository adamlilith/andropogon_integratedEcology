#' Predict a biomass model from chains
#'
#' @param chains
#' @param x Model matrix
#' @param zero_inflated Logical indicating if the model is zero-inflated
#' @param type 'mu' or 'sigma' or 'pzero'
#' @param pre_averaged If `TRUE`, then the returned values is the average across each sample unit. If `FALSE`, the returned value is a (chain * iteration)-by-sample matrix. To get the mean, sd, etc. calculate across columns.
#'
#' @returns A matrix of predictions. Rows are iterations and columns are sample IDs.
predict_biomass <- function(chains, x, zero_inflated, type = 'mu', pre_averaged = FALSE) {

	vars <- paste0('beta_', type)
	betas <- hammer_subset(chains, vars, j = TRUE)

	if (type == 'mu' & zero_inflated) {
		vars <- paste0('beta_pzero')
		betas_pzero <- hammer_subset(chains, vars, j = TRUE)
	}

	sigmas <- hammer_subset(chains, 'sigma')

	n_samples <- nrow(x)
	nchains <- length(chains$samples)
	iters_per_chain <- nrow(chains$samples[[1]])
	total_iters <- nchains * iters_per_chain

	if (pre_averaged) {
		preds <- rep(0, n_samples)
	} else {
		preds <- matrix(NA, nrow = total_iters, ncol = n_samples)
	}

	k <- 1
	for (chain in 1:nchains) {
		for (iter in 1:iters_per_chain) {
			
			this_beta <- betas$samples[[chain]][iter, ]
			this_beta <- cbind(this_beta)

			pred_untrans <- x %*% this_beta
			pred_untrans <- pred_untrans[, 1]

			this_sigma <- sigmas$samples[[chain]][iter, 'sigma']

			if (type == 'mu' & !zero_inflated) {

				this_mu <- exp(pred_untrans)

				shape <- this_mu^2 / this_sigma^2
				rate <- this_mu / this_sigma^2

				pred <- rep(NA_real_, nrow(x))
				for (count in seq_along(pred)) {
					pred[count] <- rgamma(1, shape = shape[count], rate = rate[count])
				}

			} else if (type == 'mu' & zero_inflated) {

				# mu
				this_mu <- exp(pred_untrans)

				shape <- this_mu^2 / this_sigma^2
				rate <- this_mu / this_sigma^2

				this_beta_pzero <- betas_pzero$samples[[chain]][iter, ]
				this_beta_pzero <- cbind(this_beta_pzero)

				pred_pzero_untrans <- x %*% this_beta_pzero
				pred_pzero_untrans <- pred_pzero_untrans[, 1]
				pred_pzero <- expit(pred_pzero_untrans)

				pred <- rep(NA_real_, nrow(x))
				for (count in seq_along(pred)) {
					pred[count] <- rzigamma(
						1,
						shape = shape[count],
						rate = rate[count],
						pzero = pred_pzero[count]
					)
				}
				
			} else if (type == 'pzero') {
				pred <- expit(pred_untrans)
			}
			
			if (pre_averaged) {
				preds <- preds + pred
			} else {
				preds[ , k] <- pred
			}
			k <- k + 1
		}
	}
	
	if (pre_averaged) preds <- preds / (k - 1)
	preds

}

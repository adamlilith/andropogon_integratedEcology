#' Predict a biomass model from chains
#'
#' @param chains
#' @param x Model matrix
#' @param homoscedastic `TRUE` or `FALSE`
#' @param type 'mu' or 'sigma' or 'pzero'
#'
#' @returns A matrix of predictions. Rows are iterations and columns are sample IDs.
predict_biomass <- function(chains, x, homoscedastic, type = 'mu') {


	vars <- paste0('beta_biomass_', type)
	betas <- hammer_subset(chains, vars, j = TRUE)

	if (homoscedastic) {
		sigma_biomass_among_sites <- hammer_subset(chains, 'sigma_biomass_among_sites')
	} else {
		vars <- paste0('beta_biomass_sigma')
		betas_sigma <- hammer_subset(chains, vars, j = TRUE)
	}
	sigma_biomass_within_sites <- hammer_subset(chains, 'sigma_biomass_within_sites')

	n_samples <- nrow(x)
	nchains <- length(chains$samples)
	iters_per_chain <- nrow(chains$samples[[1]])
	total_iters <- nchains * iters_per_chain

	preds <- matrix(NA, nrow = total_iters, ncol = n_samples)

	k <- 1
	for (chain in 1:nchains) {

		for (iter in 1:iters_per_chain) {

			this_beta <- betas$samples[[chain]][iter, ]
			this_beta <- cbind(this_beta)

			pred_untrans <- x %*% this_beta
			pred_untrans <- pred_untrans[ , 1]

			this_sigma_biomass_within_sites <- sigma_biomass_within_sites$samples[[chain]][iter, 'sigma_biomass_within_sites']

			if (type == 'mu' & homoscedastic) {
				
				this_sigma_biomass_among_sites <- sigma_biomass_among_sites$samples[[chain]][iter]
				log_mu_biomass <- rnorm(n_samples, mean = pred_untrans, sd = this_sigma_biomass_among_sites)
				mu_biomass <- exp(log_mu_biomass)

				shape_biomass <- mu_biomass^2 / this_sigma_biomass_within_sites^2
				rate_biomass <- mu_biomass / this_sigma_biomass_within_sites^2

				pred <- rep(NA_real_, length(mu_biomass))
				for (count in seq_along(pred)) pred[count] <- rgamma(1, shape = shape_biomass[count], rate = rate_biomass[count])

			} else if (type == 'mu' & !homoscedastic) {
				
				# this_beta_sigma <- betas_sigma$samples[[chain]][iter, ]
				# this_beta_sigma <- cbind(this_beta_sigma)

				# pred_untrans_sigma <- x %*% this_beta_sigma
				# pred_untrans_sigma <- pred_untrans_sigma[ , 1]

				# this_sigma_biomass_among_sites <- exp(pred_untrans_sigma)
				# pred <- rnorm(n_samples, mean = pred_untrans_sigma, sd = this_sigma_biomass_among_sites)
				# pred <- exp(pred)
			
			} else if (type == 'pzero') {
				pred <- expit(pred_untrans)
			} else if (type == 'sigma') {
				pred <- exp(pred_untrans)			
			}
			preds[k, ] <- pred

			k <- k + 1

		}
		
	}
	preds

}

#' Predict a non-biomass model from chains
#'
#' chains			MCMC chains list
#' x 				Model matrix
#' x_psi			Model matrix for psi (asumed to have same form as x). Leave as NULL to use x.
#' resp_distrib 	Named vector of response distribution. This can be 'gamma', 'hGamma' (zero-inflated gamma), 'lognormal', or 'hurdleLN' (zero-inflated lognormal)
#' transform		Named vector of transformations to translate MVN to mean occurrence intensity or biomass: 'identity', 'softplus' or 'exponential'.
#' force_presence	If FALSE, presence/absence dynamically determined using psi. If TRUE, psi set to == 1. Ignored if not zero inflated.
#'
#' @returns A matrix of predictions. Rows are MCMC iterations and columns are samples.
predict_nonbiomass_single_trait <- function(chains, x, x_psi = NULL, resp_distrib, transform, force_presence = FALSE) {

	vars <- paste0('beta_facet')
	betas <- mc_subset(chains, vars, j = TRUE)

	# sigmas_facet_among_sites <- mc_subset(chains, 'sigma_facet_among_sites')
	sigmas_facet_within_sites <- mc_subset(chains, 'sigma_facet_within_sites')

	zero_inflated <- resp_distrib %in% c('hGamma', 'hurdleLN')
	if (zero_inflated) {
		betas_psi <- mc_subset(chains, 'beta_psi', j = TRUE)
		if (is.null(x_psi)) x_psi <- x
	}

	n_samples <- nrow(x)
	if (is.list(chains$samples)) {
		nchains <- mc_n_chains(chains)
		iters_per_chain <- nrow(chains$samples[[1]])
	} else {
		nchains <- 1
		iters_per_chain <- nrow(chains$samples)
	}
	total_iters <- nchains * iters_per_chain

	preds <- matrix(NA, nrow = total_iters, ncol = n_samples)

	k <- 1
	for (chain in 1:nchains) {

		for (iter in 1:iters_per_chain) {

			this_beta <- betas$samples[[chain]][iter, ]
			this_beta <- cbind(this_beta)

			pred_untrans <- x %*% this_beta
			pred_untrans <- pred_untrans[ , 1]

			this_sigma_facet_within_sites <- sigmas_facet_within_sites$samples[[chain]][iter, 'sigma_facet_within_sites']

			# this_sigma_facet_among_sites <- sigmas_facet_among_sites$samples[[chain]][iter, 'sigma_facet_among_sites']
			# log_mu_facet <- rnorm(n_samples, mean = pred_untrans, sd = this_sigma_facet_among_sites)

			# zero-inflation
			if (zero_inflated & !force_presence) {

				beta_psi <- betas_psi$samples[[chain]][iter, ]
				beta_psi <- cbind(beta_psi)
				psi <- x_psi %*% beta_psi
				psi <- psi[ , 1]
				psi <- expit(psi)

			} else if (zero_inflated & force_presence) {
				psi <- rep(1, nrow(x_psi))
			}
			
			# transform mean response
			if (transform == 'identity') {
				# mu_facet <- log_mu_facet
				mu_facet <- pred_untrans
			} else if (transform == 'exponential') {
				# mu_facet <- exp(log_mu_facet)
				mu_facet <- exp(pred_untrans)
			} else if (transform == 'softplus') {
				# mu_facet <- log(1 + exp(log_mu_facet))
				mu_facet <- log(1 + exp(pred_untrans))
			} else {
				stop('Incorrect transform.')
			}

			# prediction
			pred <- rep(NA_real_, nrow(x))
			if (resp_distrib == 'gamma') {

				shape_facet <- mu_facet^2 / this_sigma_facet_within_sites^2
				rate_facet <- mu_facet / this_sigma_facet_within_sites^2

				for (count in seq_along(pred)) pred[count] <- rgamma(1, shape = shape_facet[count], rate = rate_facet[count])

			} else if (resp_distrib == 'hGamma') {

				shape_facet <- mu_facet^2 / this_sigma_facet_within_sites^2
				rate_facet <- mu_facet / this_sigma_facet_within_sites^2

				for (count in seq_along(pred)) pred[count] <- rHurdleGamma(1, shape = shape_facet[count], rate = rate_facet[count], psi = psi[count])

			} else if (resp_distrib == 'lognormal') {
			
				for (count in seq_along(pred)) pred[count] <- rlnorm(1, meanlog = mu_facet[count], sdlog = this_sigma_facet_within_sites)

			} else if (resp_distrib == 'hurdleLN') {
			
				for (count in seq_along(pred)) pred[count] <- rHLN(1, meanlog = mu_facet[count], sdlog = this_sigma_facet_within_sites, psi = psi[count])

			} else {
				stop('Bad `resp_distrib`.')
			}

			preds[k, ] <- pred

			k <- k + 1

		}
		
	}
	preds

}

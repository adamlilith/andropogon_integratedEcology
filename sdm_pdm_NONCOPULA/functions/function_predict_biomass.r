#' Predict a biomass model from chains
#'
#' chains			MCMC chains list
#' x 				Model matrix
#' resp_distrib 	Named vector of response distribution. This can be 'gamma', 'ZIG' (zero-inflated gamma), 'lognormal', or 'ZILN' (zero-inflated lognormal)
#' transform		Named vector of transformations to translate MVN to mean occurrence intensity or biomass: 'identity', 'softplus' or 'exponential'.
#'
#' @returns A matrix of predictions. Rows are MCMC iterations and columns are samples.
predict_biomass <- function(chains, x, resp_distrib, transform = NULL) {

	vars <- paste0('beta_biomass')
	betas <- mc_subset(chains, vars, j = TRUE)
	sigmas_biomass <- mc_subset(chains, 'sigma_biomass')
	betas_psi <- mc_subset(chains, 'beta_psi', j = TRUE)

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

			log_mu_biomass <- x %*% this_beta
			log_mu_biomass <- log_mu_biomass[ , 1]

			this_sigma_biomass <- sigmas_biomass$samples[[chain]][iter, 'sigma_biomass']

			# zero-inflation
			beta_psi <- betas_psi$samples[[chain]][iter, ]
			beta_psi <- cbind(beta_psi)
			psi <- x %*% beta_psi
			psi <- psi[ , 1]
			psi <- expit(psi)
			z <- rbinom(nrow(x), size = 1, prob = psi)
			
			# transform mean response
			if (transform == 'identity') {
				mu_biomass <- log_mu_biomass
			} else if (transform == 'exponential') {
				mu_biomass <- exp(log_mu_biomass)
			} else if (transform == 'softplus') {
				mu_biomass <- log(1 + exp(log_mu_biomass))
			} else {
				stop('Incorrect transform.')
			}

			# prediction
			pred <- rep(NA_real_, nrow(x))
			if (resp_distrib == 'gamma') {

				shape_biomass <- mu_biomass^2 / this_sigma_biomass^2
				rate_biomass <- mu_biomass / this_sigma_biomass^2

				for (count in seq_along(pred)) pred[count] <- rgamma(1, shape = shape_biomass[count], rate = rate_biomass[count])

			} else if (resp_distrib == 'ZIG') {

				shape_biomass <- mu_biomass^2 / this_sigma_biomass^2
				rate_biomass <- mu_biomass / this_sigma_biomass^2

				for (count in seq_along(pred)) pred[count] <- rZIG(1, shape = shape_biomass[count], rate = rate_biomass[count], z = z[count])

			} else if (resp_distrib == 'lognormal') {
			
				for (count in seq_along(pred)) pred[count] <- rlnorm(1, meanlog = log_mu_biomass[count], sdlog = this_sigma_biomass)

			} else if (resp_distrib == 'ZILN') {
			
				for (count in seq_along(pred)) pred[count] <- rZILN(1, meanlog = log_mu_biomass[count], sdlog = this_sigma_biomass, z = z[count])

			} else {
				stop('Bad `resp_distrib`.')
			}

			preds[k, ] <- pred
			k <- k + 1

		}
		
	}
	preds

}

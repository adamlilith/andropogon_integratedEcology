#' Predict occurrence from chains
#'
#' chains
#' x_occs			Model matrix
#' x_biomass		Model matrix
#' x_psi			Model matrix for probability of presence or `NULL`
#' resp_distrib 	Named vector of response distributions. For occurrence, this can be 'Poisson' or 'ZIP'. For biomass this can be 'gamma', 'ZIG' (zero-inflated gamma), 'lognormal', or 'ZILN' (zero-inflated lognormal)
#' transform		Named vector of transformations to translate MVN to mean occurrence intensity or biomass: 'identity', 'softplus' or 'exponential'.
#'
#' Returns a matrix of predictions. Rows are iterations and columns are sample IDs.
predict_occs_biomass <- function(chains, resp_distrib, transform, x_occs, x_biomass, x_psi = NULL) {

	betas_occs <- mc_subset(chains, 'beta_occs', j = TRUE)
	betas_biomass <- mc_subset(chains, 'beta_biomass', j = TRUE)
	sigma_biomass_within_sites <- mc_subset(chains, 'sigma_biomass_within_sites')
	Us <- mc_subset(chains, 'U', j = TRUE, k = TRUE)

	if (!is.null(x_psi)) betas_psi <- mc_subset(chains, 'beta_psi', j = TRUE)

	n_samples <- nrow(x_occs)
	nchains <- mc_n_chains(chains)
	chain_samples <- mc_samples(chains)
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

			# zero-inflation
			if (!is.null(x_psi)) {

				beta_psi <- betas_psi$samples[[chain]][iter, ]
				beta_psi <- cbind(beta_psi)
				psi <- x_psi %*% beta_psi
				psi <- psi[ , 1]
				psi <- expit(psi)
				z <- rbinom(nrow(x_psi), size = 1, prob = psi)

			}

			# occurrence
			lambda <- exp(Phis[ , 1])
			pred_occs <- rpois(n_samples, lambda)

			### BIOMASS
			###########

			# transform mean biomass response
			if (transform['biomass'] == 'identity') {
				mu_biomass <- Phis[ , 2]
			} else if (transform['biomass'] == 'exponential') {
				mu_biomass <- exp(Phis[ , 2])
			} else if (transform['biomass'] == 'softplus') {
				mu_biomass <- log(1 + exp(Phis[ , 2]))
			} else {
				stop('Incorrect transform.')
			}

			pred_biomass <- rep(NA_real_, n_samples)
			if (resp_distrib['biomass'] == 'gamma') {

				shape_biomass <- mu_biomass^2 / this_sigma_biomass_within_sites^2
				rate_biomass <- mu_biomass / this_sigma_biomass_within_sites^2

				for (count in seq_along(pred_biomass)) pred_biomass[count] <- rgamma(1, shape = shape_biomass[count], rate = rate_biomass[count])

			} else if (resp_distrib['biomass'] == 'ZIG') {

				shape_biomass <- mu_biomass^2 / this_sigma_biomass_within_sites^2
				rate_biomass <- mu_biomass / this_sigma_biomass_within_sites^2

				for (count in seq_along(pred_biomass)) pred_biomass[count] <- rZIG(1, shape = shape_biomass[count], rate = rate_biomass[count], z = z[count])

			} else if (resp_distrib['biomass'] == 'lognormal') {
			
				for (count in seq_along(pred_biomass)) pred_biomass[count] <- rlnorm(1, meanlog = mu_biomass[count], sdlog = this_sigma_biomass_within_sites)

			} else if (resp_distrib['biomass'] == 'ZILN') {

				for (count in seq_along(pred_biomass)) pred_biomass[count] <- rZILN(1, meanlog = mu_biomass[count], sdlog = this_sigma_biomass_within_sites, z = z[count])
			
			}

			preds_occs[k, ] <- pred_occs
			preds_biomass[k, ] <- pred_biomass

			k <- k + 1

		}
		
	}
	list(preds_occs = preds_occs, preds_biomass = preds_biomass)

}


#' Predict occurrence from chains
#'
#' source('C:/Kaji/R/andropogon_integratedEcology/sdm_pdm/functions/function_predict_occs_biomass_nonbiomass.r')
#' 
#' chains
#' nonbiomass_facets List of metadata on non-biomass facets
#' x_occs			Model matrix
#' x_biomass		Model matrix
#' x_psi			Model matrix for probability of presence or `NULL`
#' x_nonbiomass     List of model matrices for each non-biomass facet.
#' resp_distrib_biomass 	For biomass this can be 'gamma', 'hGamma' (zero-inflated gamma), 'lognormal', or 'hurdleLN' (zero-inflated lognormal)
#' transform_biomass	Named vector of transformations to translate MVN to mean biomass: 'identity', 'softplus' or 'exponential'.
#' w_occs			Model matrix for bias covariates in occurrence model, or `NULL` if no bias covariates. If not `NULL`, the function will also predict bias-corrected occurrence probabilities.
#' sampled			If `TRUE`, predictions are samples from estimated distributions (ie, as per a posterior predictive node). If `FALSE`, they are the expected (mean) value of the prediction given the environment at the sample. These values will still vary across chains/iterations, but not subject to sampling. In the `sampled = FALSE` case, all probabilities of presence will be forced to 1 if psi >=0.5 and 0 if <0.5.
#' type				'sampled' ==> calculate predictions across each chain/iteration and return all predictions. 'means' ==> calculate predictions using the summary mean across all chains/iterations.
#' force_presence		If `TRUE`, force all predictions assuming presence (N > 0). If FALSE, allow presences and absences to be predicted. This will override the effect of "type" (so if force_presence  TRUE, we will still always force presence at a prediction).
#' 
#' Returns a matrix of predictions. Rows are iterations and columns are sample IDs.
predict_fully_integrated <- function(
	chains,
	nonbiomass_facets,
	resp_distrib_biomass,
	transform_biomass,
	x_occs,
	x_psi,
	x_biomass,
	x_nonbiomass,
	w_occs = NULL,
	sampled = TRUE,
	type = 'samples',
	force_presence = FALSE
) {

	if (type == 'means') {

		# say('Predicting to mean coefficient values across all chains/iterations.')
		summary_all_chains <- chains$summary$all.chains[ , 'Mean']
		chains <- list(samples = list(matrix(summary_all_chains, nrow = 1, dimnames = list(NULL, names(summary_all_chains)))))
		chains$samples[[1]] <- as.mcmc(chains$samples[[1]])
		chains$samples <- as.mcmc.list(chains$samples)
	
	} else if (type == 'samples') {
		# say('Predicting to samples across all chains/iterations.')
	} else {
		stop('Invalid `type` argument. Must be either "samples" or "means".')
	}

	# coefficients: occurrence
	betas_occs <- mc_subset(chains, 'beta_occs', j = TRUE)
	# if (!is.null(w_occs)) alphas_occs <- mc_subset(chains, 'alpha_occs', j = TRUE)

	# coefficients: zero-inflation
	betas_psi <- mc_subset(chains, 'beta_psi', j = TRUE)
	
	# coefficients: biomass
	betas_biomass <- mc_subset(chains, 'beta_biomass', j = TRUE)
	sigma_biomass_within_sites <- mc_subset(chains, 'sigma_biomass_within_sites')

	# coefficients: non-biomass facets
	n_nonbiomass_facets <- length(nonbiomass_facets)
	betas_nonbiomass <- sigmas_facet_within_sites <- list()
	for (f in seq_along(nonbiomass_facets)) {
		betas_nonbiomass[[f]] <- mc_subset(chains, paste0('beta_facet_', f), j = TRUE)
		sigmas_facet_within_sites[[f]] <- mc_subset(chains, paste0('sigma_facet_within_sites_', f))
	}
	names(betas_nonbiomass) <- names(sigmas_facet_within_sites) <- names(nonbiomass_facets)

	# coefficients: integration
	U_stars <- mc_subset(chains, 'U_star', j = TRUE, k = TRUE)
	sigmas <- mc_subset(chains, 'sigmas', j = TRUE)

	n_samples <- nrow(x_occs)
	nchains <- mc_n_chains(chains)
	chain_samples <- mc_samples(chains)
	iters_per_chain <- nrow(chain_samples[[1]])
	total_iters <- nchains * iters_per_chain

	preds_occs <- preds_psi <- preds_bias <- preds_biomass <- matrix(NA_real_, nrow = total_iters, ncol = n_samples)
	preds_nonbiomass <- list()
	for (f in seq_along(nonbiomass_facets)) {
		preds_nonbiomass[[f]] <- matrix(NA_real_, nrow = total_iters, ncol = n_samples)
	}
	names(preds_nonbiomass) <- names(nonbiomass_facets)

	k <- 1
	for (chain in 1:nchains) {

		for (iter in 1:iters_per_chain) {

			# coefficients for this chain/iteration: occurrence
			this_beta_occs <- betas_occs$samples[[chain]][iter, ]
			this_beta_occs <- cbind(this_beta_occs)

			# if (!is.null(w_occs)) {

			# 	this_alpha_occs <- alphas_occs$samples[[chain]][iter, ]
			# 	this_alpha_occs <- cbind(this_alpha_occs)
			
			# }

			# coefficients for this chain/iteration: biomass
			this_beta_biomass <- betas_biomass$samples[[chain]][iter, ]
			this_beta_biomass <- cbind(this_beta_biomass)

			this_sigma_biomass_within_sites <- sigma_biomass_within_sites$samples[[chain]][iter, 'sigma_biomass_within_sites']

			# coefficients for this chain/iteration: non-biomass facets
			this_sigma_facet_within_sites <- rep(NA_real_, length(nonbiomass_facets))
			this_beta_nonbiomass <- list()
			for (f in seq_along(nonbiomass_facets)) {
			
				n_terms <- length(attr(terms(nonbiomass_facets[[f]]$formula), 'term.labels')) + 1
				this_beta_nonbiomass[[f]] <- betas_nonbiomass[[f]]$samples[[chain]][iter, paste0('beta_facet_', f, '[', 1:n_terms,']')]
				this_sigma_facet_within_sites[f] <- sigmas_facet_within_sites[[f]]$samples[[chain]][iter, paste0('sigma_facet_within_sites_', f)]
			
			}

			# predict latent abundance
			phi_lambda_mu <- x_occs %*% this_beta_occs
			phi_lambda_mu <- phi_lambda_mu[ , 1]

			# # predict latent observation probability
			# if (!is.null(w_occs)) {

			# 	p <- w_occs %*% this_alpha_occs
			# 	p <- p[ , 1]
			# 	p <- expit(p)
			
			# }

			# predict latent BIOMASS
			phi_biomass_mu <- x_biomass %*% this_beta_biomass
			phi_biomass_mu <- phi_biomass_mu[ , 1]

			# predict latent NON-BIOMASS facets
			phi_facets_mu <- matrix(NA_real_, nrow = length(nonbiomass_facets), ncol = n_samples)
			for (f in seq_along(nonbiomass_facets)) {
				facet <- names(nonbiomass_facets)[f]
				phi_facets_mu[f, ] <- x_nonbiomass[[facet]] %*% cbind(this_beta_nonbiomass[[f]])
			}

			# Cholesky
			U_star <- U_stars$samples[[chain]][iter, ]
			U_star <- matrix(U_star, ncol = 2 + n_nonbiomass_facets)

			# multivariate simulation
			Phi <- matrix(NA_real_, ncol = 2 + n_nonbiomass_facets, nrow = n_samples)
			n_phis <- 2 + length(nonbiomass_facets)
			for (i in seq_len(n_samples)) {

				z <- if (type == 'samples') {
					rnorm(2 + n_nonbiomass_facets, 0, 1)
				} else {
					rep(0, 2 + n_nonbiomass_facets)
				}

				# abundance
				this_sigmas <- sigmas$samples[[chain]][iter, 'sigmas[1]']
				Phi[i, 1] <- phi_lambda_mu[i] + this_sigmas * inprod(U_star[ , 1], z)
				
				this_sigmas <- sigmas$samples[[chain]][iter, 'sigmas[2]']
				Phi[i, 2] <- phi_biomass_mu[i] + this_sigmas * inprod(U_star[ , 2], z)

				for (f in seq_along(nonbiomass_facets)) {
				
					this_sigmas <- sigmas$samples[[chain]][iter, paste0('sigmas[', f + 2, ']')]
					Phi[i, f + 2] <- phi_facets_mu[f, i] + this_sigmas * inprod(U_star[ , f + 2], z)
				
				}
					
			}

			# zero-inflation
			beta_psi <- betas_psi$samples[[chain]][iter, ]
			beta_psi <- cbind(beta_psi)
			psi <- x_psi %*% beta_psi
			psi <- psi[ , 1]
			psi <- expit(psi)
			if (force_presence) {
				z <- rep(1, n_samples)
			} else if (sampled) {
				z <- rbinom(n_samples, size = 1, prob = psi)
			} else {
				z <- as.numeric(psi >= 0.5)
			}

			preds_psi[k, ] <- psi

			### predict empirical occurrence
			lambda <- exp(Phi[ , 1])
			pred_occs <- rep(NA_real_, n_samples)
			if (sampled) {
				for (i in 1:n_samples) pred_occs[i] <- rHurdlePoisson(1, lambda[i], psi = z[i])
			} else {
				pred_occs <- lambda
			}
			preds_occs[k, ] <- pred_occs

			# ### predict observation probability/number of AG observed
			# if (!is.null(w_occs)) {
			# 	if (sampled) {
			# 		pred_bias <- rbinom(n_samples, size = pred_occs, prob = p)
			# 	} else {
			# 		pred_bias <- p
			# 	}
			# 	preds_bias[k, ] <- pred_bias
			# }

			### predict empirical biomass
			if (transform_biomass == 'identity') {
				mu_biomass <- Phi[ , 2]
			} else if (transform_biomass == 'exponential') {
				mu_biomass <- exp(Phi[ , 2])
			} else if (transform_biomass == 'softplus') {
				mu_biomass <- log(1 + exp(Phi[ , 2]))
			} else {
				stop('Incorrect transform.')
			}

			pred_biomass <- rep(NA_real_, n_samples)
			if (resp_distrib_biomass == 'gamma') {

				shape_biomass <- mu_biomass^2 / this_sigma_biomass_within_sites^2
				rate_biomass <- mu_biomass / this_sigma_biomass_within_sites^2

				if (sampled) {
					for (count in seq_along(pred_biomass)) pred_biomass[count] <- rgamma(1, shape = shape_biomass[count], rate = rate_biomass[count])
				} else {
					pred_biomass <- mu_biomass
				}

			} else if (resp_distrib_biomass == 'hGamma') {

				shape_biomass <- mu_biomass^2 / this_sigma_biomass_within_sites^2
				rate_biomass <- mu_biomass / this_sigma_biomass_within_sites^2

				if (sampled) {
					for (count in seq_along(pred_biomass)) pred_biomass[count] <- rZIG(1, shape = shape_biomass[count], rate = rate_biomass[count], psi = z[count])
				} else {
					pred_biomass <- mu_biomass
				}

			} else if (resp_distrib_biomass == 'lognormal') {
			
				if (sampled) {
					for (count in seq_along(pred_biomass)) pred_biomass[count] <- rlnorm(1, meanlog = mu_biomass[count], sdlog = this_sigma_biomass_within_sites)
				} else {
					pred_biomass <- exp(mu_biomass)
				}

			} else if (resp_distrib_biomass == 'hurdleLN') {

				if (sampled) {
					for (count in seq_along(pred_biomass)) pred_biomass[count] <- rHLN(1, meanlog = mu_biomass[count], sdlog = this_sigma_biomass_within_sites, psi = z[count])
				} else {
					pred_biomass <- exp(mu_biomass)
				}
			
			}
			preds_biomass[k, ] <- pred_biomass

			### predict empirical non-biomass facets
			for (f in seq_along(nonbiomass_facets)) {

				facet <- names(nonbiomass_facets)[f]
				transform <- nonbiomass_facets[[facet]]$transform

				if (transform == 'identity') {
					mu_facet <- Phi[ , 2 + f]
				} else if (transform == 'exponential') {
					mu_facet <- exp(Phi[ , 2 + f])
				} else if (transform == 'softplus') {
					mu_facet <- log(1 + exp(Phi[ , 2 + f]))
				} else {
					stop('Incorrect transform.')
				}

				pred_nonbiomass <- rep(NA_real_, n_samples)
				if (nonbiomass_facets[[facet]]$resp_distrib == 'gamma') {

					shape_facet <- mu_facet^2 / this_sigma_facet_within_sites[f]^2
					rate_facet <- mu_facet / this_sigma_facet_within_sites[f]^2

					if (sampled) {
						for (count in seq_along(mu_facet)) pred_nonbiomass[count] <- rgamma(1, shape = shape_facet[count], rate = rate_facet[count])
					} else {
						pred_nonbiomass <- mu_facet
					}

				} else if (nonbiomass_facets[[facet]]$resp_distrib == 'lognormal') {

					if (sampled) {
						for (count in seq_along(mu_facet)) pred_nonbiomass[count] <- rlnorm(1, meanlog = mu_facet[count], sdlog = this_sigma_facet_within_sites[f])
					} else {
						pred_nonbiomass <- mu_facet
					}
				
				} else if (nonbiomass_facets[[facet]]$resp_distrib == 'hGamma') {

					shape_facet <- mu_facet^2 / this_sigma_facet_within_sites[f]^2
					rate_facet <- mu_facet / this_sigma_facet_within_sites[f]^2

					if (sampled) {
						for (count in seq_along(mu_facet)) pred_nonbiomass[count] <- rZIG(1, shape = shape_facet[count], rate = rate_facet[count], psi = z[count])
					} else {
						pred_nonbiomass <- mu_facet
					}

				} else if (nonbiomass_facets[[facet]]$resp_distrib == 'hurdleLN') {
				
					if (sampled) {
						for (count in seq_along(mu_facet)) pred_nonbiomass[count] <- rHLN(1, meanlog = mu_facet[count], sdlog = this_sigma_facet_within_sites[f], psi = z[count])
					} else {
						pred_nonbiomass <- mu_facet
					}
				
				} else {
					stop('Response distribution for non-biomass facet ', facet, ' not supported.')
				}

				preds_nonbiomass[[facet]][k, ] <- pred_nonbiomass

			}

			k <- k + 1

		}
		
	}

	if (!is.null(w_occs)) {
		list(preds_occs = preds_occs, preds_psi = preds_psi, preds_bias = preds_bias, preds_biomass = preds_biomass, preds_nonbiomass = preds_nonbiomass)
	} else {
		list(preds_occs = preds_occs, preds_psi = preds_psi, preds_biomass = preds_biomass, preds_nonbiomass = preds_nonbiomass)
	}

}


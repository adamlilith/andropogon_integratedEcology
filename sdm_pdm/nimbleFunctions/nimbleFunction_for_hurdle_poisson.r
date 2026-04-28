#' Hurdle Poisson distribution
#'
#' Distribution and random number generator function for a "true" hurdle (zero-inflated) Poisson distribution. The distribution function includes as small "cheat" factor to increase `psi` if it is numerically indistinguishable from 0, and decrease it if it is indistinguishable from 1, to avoid infinite likelihoods. The RNG function does not include this cheat factor, so it can handles cases where `psi` is exactly 0 or 1.
#'
#' x 		Numeric vector.
#' p 		Numeric vector of probabilities.
#' q 		Numeric vector of quantiles.
#' psi 	    Numeric in [0, 1] indicating probability of non-zero value. NOTE this is different from the usual hurdle case where it would represent the probability of a zero value.
#' lambda 	Numeric >= 0. Expected density.
#' log 		Logical (if `TRUE` return the log value), or 0 (`FALSE`)/1 (`TRUE`).
#' n 		Number of values. Can only be 1 when used in a nimble model.
dHurdlePoisson <- nimbleFunction(
    run = function(
        x = double(),
        lambda = double(),
        psi = double(),
        log = integer(default = 0)
    ) {

    returnType(double())

	# if (x == 0) {
	# 	out <- 1 - psi
	# } else {
	# 	out <- psi * (dpois(x, lambda = lambda)) / (1 - exp(-lambda))
	# }

	# if (log) out <- log(out)
	# return(out)

    # clamp psi for safety
    psi <- min(1 - 1e-8, psi)
    psi <- max(1e-8, psi)

    # calculate in log-space because safer if psi ~0 or ~1
    if (x == 0) {
        loglik <- log1p(-psi) # more stable near 0 than log(1 - psi)
    } else {
        lambda <- max(lambda, 1e-8) # clamp lambda for safety

        # normalization by log1p(-exp(-lambda)) to account for fact that any zeros are accounted for by the hurdle component
        loglik <- log(psi) + dpois(x, lambda = lambda, log = 1) - log1p(-exp(-lambda))
    }

    if (log) return(loglik)
    return(exp(loglik))
        
    },
    buildDerivs = 'run'
)

# RNG for truncated, pseudo-zero-inflated Poisson
rHurdlePoisson <- nimbleFunction(
    run = function(
        n = integer(),
        lambda = double(),
        psi = double()
    ) {
        
    returnType(double())

	zero <- runif(1, min = 0, max = 1)
	if (zero < (1 - psi)) {
		return(0)
	} else {
		unif <- runif(1)
		exp_lambda <- exp(-lambda)
		out <- qpois(unif * (1 - exp_lambda) + exp_lambda, lambda)
	}

	return(out)

    } # EOF
)

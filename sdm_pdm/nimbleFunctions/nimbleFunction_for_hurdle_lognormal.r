#' Hurdle (zero-inflated) lognormal distribution
#'
#' Density and random number generator functions for a hurdle (zero-inflated) lognormal distribution. The distribution function includes a "cheat" to obviate issues with taking the log density when `psi` is 0 or 1.
#'
#' x        Numeric vector.
#' p        Numeric vector of probabilities.
#' q        Numeric vector of quantiles.
#' psi      Numeric in range [0, 1] indicating probability of a non-zero value.
#' meanlog,sdlog Numeric. Parameters of the lognormal distribution. `sdlog` must be >=0.
#' log      Logical (if `TRUE` return the log value), or 0 (`FALSE`)/1 (`TRUE`).
#' n        Number of values. Can only be 1.
dHLN <- nimbleFunction(
    run = function(
        x = double(0),
        meanlog = double(0),
        sdlog = double(0),
        psi = double(0),
        log = integer(0, default = 0)
    ) {
    
    returnType(double(0))

    if (log) {
        # clamp psi for safety
        psi <- min(1 - 1e-8, psi)
        psi <- max(1e-8, psi)
        
        if (x == 0) {
            out <- log1p(-psi)
        } else {
            out <- log(psi) + dlnorm(x, meanlog = meanlog, sdlog = sdlog, log = 1)
        }
    } else {
        if (x == 0) {
            out <- 1 - psi
        } else {
            out <- psi * dlnorm(x, meanlog = meanlog, sdlog = sdlog, log = 0)
        }
    }
    return(out)

    },
    buildDerivs = 'run'
)

# RNG for zero-/non-zero forced lognormal
rHLN <- nimbleFunction(
	run = function(
		n = integer(0),
		meanlog = double(0),
		sdlog = double(0),
		psi = double(0)
	) {
		
	returnType(double(0))
	# if (sdlog < 0) {
	# 	print('sdlog must be >= 0.')
	# 	return(NaN)
	# }
	
	# if (n != 1) print('rHLN() only allows n = 1.')
	
    unif <- runif(1)
    if (unif < (1 - psi)) {
        return(0.0)
    } else {
        out <- rlnorm(1, meanlog = meanlog, sdlog = sdlog)
    }

	return(out)
	
	} # EOF
)

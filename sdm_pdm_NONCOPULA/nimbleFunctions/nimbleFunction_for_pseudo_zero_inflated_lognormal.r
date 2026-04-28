#' Hurdle (zero-inflated) lognormal distribution
#'
#' Density and random number generator functions for a hurdle zero-inflated lognormal distribution. If the parameter `z` is 0, then the density of any `x` > 0 is 0, whereas if `x` is 0, the density is 1. If `z` is 1, then the density is 0 if `x` is 0, and the density of `x` > 0 is the density of the typical lognormal density function. In other words, the user forces whether or not the values of interest should have zero or non-zero values.
#'
#' @param x Numeric vector.
#' @param p Numeric vector of probabilities.
#' @param q Numeric vector of quantiles.
#' @param z Numeric, either 0 or 1. Zero indicates "absence" (i.e., x = 0), and 1 "presence" (i.e., x > 0).
#' @param meanlog,sdlog Numeric. Parameters of the lognormal distribution. `sdlog` must be >=0.
#' @param log Logical (if `TRUE` return the log value), or 0 (`FALSE`)/1 (`TRUE`).
#' @param n Number of values. Can only be 1.
dZILN <- nimbleFunction(
    run = function(
        x = double(0),
        meanlog = double(0),
        sdlog = double(0),
        z = double(0),
        log = integer(0, default = 0)
    ) {
    
    returnType(double(0))

	if (sdlog < 0) return(-Inf)

    if (z == 0) {
        if (x == 0) {
            dens <- 1.0
        } else {
            dens <- 0.0
        }
    } else {
        if (x == 0) {
            dens <- 0.0
        } else {
            dens <- dlnorm(x, meanlog = meanlog, sdlog = sdlog)
        }
    }
     
    if (log) dens <- log(dens)
    return(dens)

    },
    buildDerivs = 'run'
)
# RNG for zero-/non-zero forced lognormal
rZILN <- nimbleFunction(
	run = function(
		n = integer(0),
		meanlog = double(0),
		sdlog = double(0),
		z = double(0)
	) {
		
	returnType(double(0))
	if (sdlog < 0) {
		print('sdlog must be >= 0.')
		return(NaN)
	}
	
	if (n != 1) print('dZILN() and rZILN() only allows n = 1.')
	
	if (z == 0) {
		out <- 0.0
	} else {
		out <- rlnorm(1, meanlog = meanlog, sdlog = sdlog)
	}

	return(out)
	
	} # EOF
)

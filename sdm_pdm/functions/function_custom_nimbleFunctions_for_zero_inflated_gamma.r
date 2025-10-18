#' Zero-inflated gamma distribution
#'
#' @param x Numeric vector.
#' @param p Numeric vector of probabilities.
#' @param q Numeric vector of quantiles.
#' @param pzero Numeric >=0 and <= 1. Probability of 0's.
#' @param shape,rate Numeric >0. Parameters of the gamma distribution.
#' @param log Logical (if `TRUE` return the log value), or 0 (`FALSE`)/1 (`TRUE`).
#' @param n Number of values.
dzigamma <- nimbleFunction(
    run = function(
        x = double(0),
        shape = double(0),
        rate = double(0),
        pzero = double(0),
        log = integer(0, default = 0)
    ) {
    
    returnType(double(0))

    if (x == 0) {
        dens <- pzero
	} else if (x > 0) {
        dens <- (1 - pzero) * dgamma(x, shape = shape, rate = rate)
	} else {
		dens <- -Inf  # only defined for x >= 0
	}

    if (log) dens <- log(dens)
    return(dens)

    },
    buildDerivs = 'run'
)

# RNG for zero-inflated gamma
rzigamma <- nimbleFunction(
    run = function(
        n = integer(0),
        shape = double(0),
        rate = double(0),
        pzero = double(0)
    ) {
        
    returnType(double(0))
    
    if (n != 1) print("rzigamma() only allows n = 1")
    
    if (runif(1, 0, 1) < pzero) {
        out <- 0.0
    } else {
        out <- rgamma(1, shape = shape, rate = rate)
    }
    return(out)
    
    } # EOF
)

pzigamma <- nimbleFunction(
    run = function(
        q = double(0),
        shape = double(0),
        rate = double(0),
        pzero = double(0),
        lower.tail = integer(0, default = 1),
        log.p = integer(0, default = 0)
    ) {
    
    returnType(double(0))
        
    if (q < 0) {
        p <- 0.0
    } else if (q == 0) {
        p <- pzero
    } else {
        p <- pzero + (1 - pzero) * pgamma(q, shape = shape, rate = rate, lower.tail = 1, log.p = 0)
    }

    if (!lower.tail) p <- 1 - p
    if (log.p) p <- log(p)
    return(p)

    } # EOF
)

qzigamma <- nimbleFunction(
    run = function(
        p = double(0),
        shape = double(0),
        rate = double(0),
        pzero = double(0),
        lower.tail = integer(0, default = 1),
        log.p = integer(0, default = 0)
    ) {

    returnType(double(0))

    if (log.p) p <- exp(p)
    if (!lower.tail) p <- 1 - p

    if (p <= pzero) {
        out <- 0.0
    } else {
        adj_p <- (p - pzero) / (1 - pzero)
        out <- qgamma(adj_p, shape = shape, rate = rate)
    }
    return(out)

    } # EOF
)

#' Joint Poisson-gamma distribution
#'
#' The distribution and random number generator functions of a joint Poisson and gamma distribution. The distribution function expects the Poisson and gamma values to be both 0 or both non-zero, and the random number function returns numbers that are both zero or both non-zero.
#'
#' @param x Numeric vector with two values. The first is the value for the Poisson and the second for the gamma.
#' @param p Numeric vector of probabilities.
#' @param q Numeric vector of quantiles.
#' @param shape,rate Numeric >0. Parameter of gamma distribution: Parameters of the gamma distribution.
#' @param lambda Numeric >= 0. Parameter of Poisson distribution: Expected density.
#' @param log Logical (if `TRUE` return the log value), or 0 (`FALSE`)/1 (`TRUE`).
#' @param n Number of values.
dpoisgamma1 <- nimbleFunction(
    run = function(
        x = double(1),
		lambda = double(0),
        shape = double(0),
        rate = double(0),
        log = integer(0, default = 0)
    ) {
    
    returnType(double(0))

	x1 <- x[1] # Poisson
	x2 <- x[2] # gamma

    if (x1 == 0 & x2 == 0) {
        # Poisson only since gamma(0, shape, rate, log = 0) is 1
        dens <- dpois(0, lambda, log = 0)
    } else if (x1 > 0 & x2 > 0) {
		dens <- dpois(x1, lambda, log = 0) * dgamma(x2, shape = shape, rate = rate, log = 0)
    } else { # x1 = 0 and x2 > 0 or vice versa
		dens <- 0
	}
    if (log & dens > 0) dens <- log(dens)
    if (log & dens == 0) dens <- -Inf
	return(dens)

    },
    buildDerivs = 'run'
)

# RNG for joint, zero-inflated Poisson and gamma 
rpoisgamma1 <- nimbleFunction(
    run = function(
        n = integer(0),
		lambda = double(0),
        shape = double(0),
        rate = double(0)
    ) {
        
    returnType(double(0))
    
    if (n != 1) print("rpoisgamma1() only allows n = 1")
    
	pois <- rpois(1, lambda)
	g1 <- rgamma(1, shape = shape, rate = rate)
	if (pois == 0) g1 <- 0.0
	out <- c(pois, g1)
    return(out)
    
    } # EOF
)

# pzipoisgamma1 <- nimbleFunction(
#     run = function(
#         q = double(0),
# 		lambda = double(0),
#         shape = double(0),
#         rate = double(0),
#         psi = double(0),
#         lower.tail = integer(0, default = 1),
#         log.p = integer(0, default = 0)
#     ) {
    
#     returnType(double(0))
        
#     if (q < 0) {
#         p <- 0.0
#     } else if (q == 0) {
#         p <- psi
#     } else {
#         p <- psi + (1 - psi) * ppois(q, lambda, log.p = 0) * pgamma(q, shape = shape, rate = rate, lower.tail = 1, log.p = 0)
#     }

#     if (!lower.tail) p <- 1 - p
#     if (log.p) p <- log(p)
#     return(p)

#     } # EOF
# )

# qzigamma <- nimbleFunction(
#     run = function(
#         p = double(0),
#         shape = double(0),
#         rate = double(0),
#         psi = double(0),
#         lower.tail = integer(0, default = 1),
#         log.p = integer(0, default = 0)
#     ) {

#     returnType(double(0))

#     if (log.p) p <- exp(p)
#     if (!lower.tail) p <- 1 - p

#     if (p <= psi) {
#         out <- 0.0
#     } else {
#         adj_p <- (p - psi) / (1 - psi)
#         out <- qgamma(adj_p, shape = shape, rate = rate)
#     }
#     return(out)

#     } # EOF
# )

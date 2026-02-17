#' Hurdle (zero-inflated) inflated gamma distribution
#'
#' Density and random number generator functions for a hurdle zero-inflated gamma distribution. If the parameter `z` is 0, then the density of any `x` > 0 is 0, whereas if `x` is 0, the density is 1. If `z` is 1, then the density of `x` that is 0 is 0, and the density of `x` > 0 is the density of the typical gamma density function. In other words, the user forces whether or not the values of interest should have zero or non-zero values.
#'
#' @param x Numeric vector.
#' @param p Numeric vector of probabilities.
#' @param q Numeric vector of quantiles.
#' @param z Numeric, either 0 or 1. Zero indicates "absence" (i.e., x = 0), and 1 "presence" (i.e., x > 0).
#' @param shape,rate Numeric >0. Parameters of the gamma distribution.
#' @param log Logical (if `TRUE` return the log value), or 0 (`FALSE`)/1 (`TRUE`).
#' @param n Number of values.
dZIG <- nimbleFunction(
    run = function(
        x = double(0),
        shape = double(0),
        rate = double(0),
        z = double(0),
        log = integer(0, default = 0)
    ) {
    
    returnType(double(0))

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
            dens <- dgamma(x, shape = shape, rate = rate)
        }
    }
     
    if (log) dens <- log(dens)
    return(dens)

    },
    buildDerivs = 'run'
)

# RNG for zero-/non-zero forced gamma
rZIG <- nimbleFunction(
    run = function(
        n = integer(0),
        shape = double(0),
        rate = double(0),
        z = double(0)
    ) {
        
    returnType(double(0))
    
    if (n != 1) print('dZIG() only allows n = 1.')
    
    if (z == 0) {
        out <- 0.0
    } else {
        out <- rgamma(1, shape = shape, rate = rate)
    }

    return(out)
    
    } # EOF
)

# # # pzigamma <- nimbleFunction(
# # #     run = function(
# # #         q = double(0),
# # #         shape = double(0),
# # #         rate = double(0),
# # #         psi = double(0),
# # #         lower.tail = integer(0, default = 1),
# # #         log.p = integer(0, default = 0)
# # #     ) {
    
# # #     returnType(double(0))
        
# # #     if (q < 0) {
# # #         p <- 0.0
# # #     } else if (q == 0) {
# # #         p <- psi
# # #     } else {
# # #         p <- psi + (1 - psi) * pgamma(q, shape = shape, rate = rate, lower.tail = 1, log.p = 0)
# # #     }

# # #     if (!lower.tail) p <- 1 - p
# # #     if (log.p) p <- log(p)
# # #     return(p)

# # #     } # EOF
# # # )

# # # qzigamma <- nimbleFunction(
# # #     run = function(
# # #         p = double(0),
# # #         shape = double(0),
# # #         rate = double(0),
# # #         psi = double(0),
# # #         lower.tail = integer(0, default = 1),
# # #         log.p = integer(0, default = 0)
# # #     ) {

# # #     returnType(double(0))

# # #     if (log.p) p <- exp(p)
# # #     if (!lower.tail) p <- 1 - p

# # #     if (p <= psi) {
# # #         out <- 0.0
# # #     } else {
# # #         adj_p <- (p - psi) / (1 - psi)
# # #         out <- qgamma(adj_p, shape = shape, rate = rate)
# # #     }
# # #     return(out)

# # #     } # EOF
# # # )

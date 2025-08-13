#' Zero-inflated Poisson distribution
#'
#' @description See https://r-nimble.org/nimbleExamples/zero_inflated_poisson.html.
#'
#' @param x Numeric vector.
#' @param p Numeric vector of probabilities.
#' @param q Numeric vector of quantiles.
#' @param pzero Numeric >=0 and <= 1. Probability of 0's.
#' @param lambda Numeric >= 0. Expected density.
#' @param log Logical (if `TRUE` return the log value), or 0 (`FALSE`)/1 (`TRUE`).
#' @param n Number of values.
dzip <- nimbleFunction(
    run = function(
        x = double(),
        lambda = double(),
        pzero = double(),
        log = integer(default = 0)
    ) {

    returnType(double())

    # better formulation for Laplacian(?)
    if (x == 0) {
        p <- pzero + (1 - pzero) * dpois(0, lambda)
    } else {
        p <- (1 - pzero) * dpois(x, lambda)
    }
    # p <- pzero * dbinom(x, size = 1, prob = 0) + (1 - pzero) * dpois(x, lambda, log = 0)
    if (log) return(log(p))
    return(p)

    },
    buildDerivs = 'run'
)

# RNG for zero-inflated Poisson
rzip <- nimbleFunction(
    run = function(
        n = integer(),
        lambda = double(),
        pzero = double()
    ) {
        
    returnType(double())

    # is_struct_zero <- rbinom(1, prob = pzero, size = 1)
    # if (is_struct_zero) return(0)
    # return(rpois(1, lambda = lambda))

    if (runif(1, 0, 1) < pzero) {
        return(0)
    } else {
        return(rpois(1, lambda = lambda))
    }
    
    } # EOF
)

# PMF function for a zero-inflated Poisson
pzip <- nimbleFunction(
    run = function(
        q = double(),
        lambda = double(),
        pzero = double(),
        lower.tail = integer(default = 1),
        log.p = integer(default = 0)
    ) {
    
    returnType(double())

    if (q < 0) {
        p <- 0.0
    } else if (lower.tail) {
    
        if (q == 0) {
            p <- pzero + (1 - pzero) * exp(-lambda)
        } else {
            p <- (1 - pzero) * ppois(q, lambda)
        }
    
    } else {
    
        if (q == 0) {
            p <- 1 - (pzero + (1 - pzero) * exp(-lambda))
        } else {
            p <- 1 - (1 - pzero) * ppois(q, lambda)
        }
    
    }

    if (log.p) {
        return(log(p))
    } else {
        return(p)
    }

    } # EOF
)

qzip <- nimbleFunction(
    run = function(
        p = double(),
        lambda = double(),
        pzero = double(),
        lower.tail = integer(default = 1),
        log.p = integer(default = 0)
    ) {

    returnType(integer())

    if (log.p) p[i] <- exp(p[i])

    # invalid
    if (p < 0 | p > 1) {
        out <- NA_integer_
    } else {

        if (!lower.tail) {
            p <- 1 - p
        }

        # mass at zero
        p0 <- pzero + (1 - pzero) * exp(-lambda)

        if (p <= p0) {
            out <- 0L
        } else {
            adj_p <- (p - pzero) / (1 - pzero)
            out <- qpois(p = adj_p, lambda = lambda)
        }
    }
    return(out)

    } # EOF
)

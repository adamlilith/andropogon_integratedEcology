#' Hurdle (zero-inflated) inflated gamma distribution
#'
#' Density and random number generator functions for a hurdle (zero-inflated) gamma distribution.
#'
#' x        Numeric vector.
#' p        Numeric vector of probabilities.
#' q        Numeric vector of quantiles.
#' psi      Numeric in range [0, 1]. Probability of a non-zero value.
#' shape,rate Numeric >0. Parameters of the gamma distribution.
#' log      Logical (if `TRUE` return the log value), or 0 (`FALSE`)/1 (`TRUE`).
#' n        Number of values. Must be 1.
dHurdleGamma <- nimbleFunction(
    run = function(
        x = double(0),
        shape = double(0),
        rate = double(0),
        psi = double(0),
        log = integer(0, default = 0)
    ) {
    
    returnType(double(0))

    # clamp psi for safety
    if (log) {
        psi <- min(1 - 1e-8, psi)
        psi <- max(1e-8, psi)
    }

    if (x == 0) {
        out <- 1 - psi
    } else {
        out <- psi * dgamma(x, shape = shape, rate = rate)
    }

    # if (psi == 0) {
    #     if (x == 0) {
    #         out <- 1.0
    #     } else {
    #         out <- 0.0
    #     }
    # } else {
    #     if (x == 0) {
    #         out <- 0.0
    #     } else {
    #         out <- dgamma(x, shape = shape, rate = rate)
    #     }
    # }
     
    if (log) out <- log(out)
    return(out)

    },
    buildDerivs = 'run'
)

# RNG for zero-/non-zero forced gamma
rHurdleGamma <- nimbleFunction(
    run = function(
        n = integer(0),
        shape = double(0),
        rate = double(0),
        psi = double(0)
    ) {
        
    returnType(double(0))
    
    if (n != 1) print('rHurdleGamma() only allows n = 1.')
    
    unif <- runif(1)
    if (unif < (1 - psi)) {
        return(0)
    } else {
        out <- rgamma(1, shape = shape, rate = rate)
    }

    return(out)
    
    } # EOF
)

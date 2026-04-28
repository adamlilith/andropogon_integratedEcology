#' Distribution and RNG function for joint Poisson/binomial distribution
#'
#' Useful for speeding calculation of likelihood of latent abundance and observed abundance
#' 
#' x        Observed value
#' n        Always = 1
#' lambda   Poisson mean (expected abundance)
#' psi      Hurdle probability (probability of being in the Poisson part of the model)
#' p        Binomial probability (probability of detection/sampling)
#' log      1 (TRUE) or 0 (FALSE)
dHurdlePoissonBinomial <- nimbleFunction(
    run = function(
        x = double(0),
        lambda = double(0),
        psi = double(0),
        p = double(0),
        log = integer(0, default=0)
    ) {
    
    returnType(double())

    # avoid infinite likelihoods through truncation
    psi <- max(1e-8, min(1 - 1e-8, psi))
    lambda <- max(1e-8, lambda)
    p <- max(1e-8, min(1 - 1e-8, p))

    # zero can come from hurdle or from Poisson --> binomial
    if (x == 0) {
        # out <- (1 - psi) + psi * exp(-lambda * p)
        out <- (1 - psi) + psi * (exp(-lambda * p) - exp(-lambda)) / (1 - exp(-lambda))
        return(out)
    } else { # counts > 0
        # calculate in log-space to prevent underflow crashes
        if (log) {
            out <- log(psi) + dpois(x, lambda * p, log = 1) - log1p(-exp(-lambda))
            return(out)
        } else {
            out <- psi * dpois(x, lambda * p, log = 0) / (1 - exp(-lambda))
            return(out)
        }
    }
    # buildDerivs = 'run'
    } # EOF
)

# RNG
rHurdlePoissonBinomial <- nimbleFunction(
    run = function(
        n = integer(0),
        lambda = double(0),
        psi = double(0),
        p = double(0)
    ) {
    
    returnType(double())

    # avoid under/overflows
    psi <- max(1e-8, min(1 - 1e-8, psi))
    lambda <- max(1e-8, lambda)
    p <- max(1e-8, min(1 - 1e-8, p))
    
    # probability of zero
    p_zero <- (1 - psi) + psi * (exp(-lambda * p) - exp(-lambda)) / (1 - exp(-lambda))
    
    u1 <- runif(1)
    if (u1 < p_zero) {
        return(0.0)
    } else {

        # simulate N ~ zero-truncated Poisson by truncating effective lambda at value >0 then applying quantile function
        
        # draw between exp(-lambda) and 1 to force qpois >= 1
        exp_neg_lambda <- exp(-lambda)
        u_trunc <- runif(1)
        N <- qpois(u_trunc * (1 - exp_neg_lambda) + exp_neg_lambda, lambda = lambda)

        # simulate observation from binomial
        return(rbinom(1, size = N, prob = p))

    }
            
    } # EOF
)

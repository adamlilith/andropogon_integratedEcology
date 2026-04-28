#' Custom nimbleFunction()s for truncated Poisson
#'
#' These are distribution and simulation functions for a truncated Poisson distribution where values are prohibited from being less than a user-supplied value. They are useful, for example, when modeling latent abundance as a draw from a Poisson where the mean of the Poisson is in turn a draw from the exponentiation of a normal distribution. In these cases, using an untruncated (default) Poisson would not prohibit the total number of latent individuals to be less than the observed number.
#'
#' @param lambda Expected number.
#' @param y_min Minimum allows value from Poisson.
#' @param log 1 (return log likelihood) or 0 (return arithmetic likelihood)
#' @param n Number of values to simulate. Can only accept 1 for now.
dTruncPoisson <- nimbleFunction(
    # run = function(x = integer(0), lambda = double(0), y_min = integer(0), log = integer(0)) {
    run = function(x = integer(0), lambda = double(0), y_min = integer(0), log = 1) {
        returnType(double(0))
        if (x < y_min) {
            if (log == 1) {
                return(-Inf)
            } else {
                return(0.0)
            }
        } else {

            # log-scale for stability
            log_prob <- dpois(x, lambda, log = 1)
            log_denom <- log(1 - ppois(y_min - 1, lambda = lambda))
            
            min_log <- -36.736800569677 # smallest non-infinite log value
            if (log_denom < min_log) {
                log_denom <- min_log
            }

            log_out <- log_prob - log_denom

            if (log == 1) {
                return(log_out)
            } else {
                return(exp(log_out))
            }
        }
    }
)

# # sampling function using rejection sampling... can be very slow
# rTruncPoisson <- nimbleFunction(
#   run = function(n = integer(0), lambda = double(0), y_min = integer(0)) {
#     returnType(integer(0))
#     value <- 0
#     while (1) {
#       value <- rpois(1, lambda)
#       if (value >= y_min) {
#         return(value)
#       }
#     }
#   }
# )

rTruncPoisson <- nimbleFunction(
    run = function(n = integer(0), lambda = double(0), y_min = integer(0)) {
      
        returnType(integer(0))

        p_upper <- 0.9999
        # p_upper <- 0.999

        # CDF at (y_min - 1)
        if (y_min < 1) {
            p_lower <- 0
        } else {
            p_lower <- ppois(y_min - 1, lambda)
        }

        # need to re-define p_lower in cases when lambda is too small to yield y_min
        # (ppois() can return ~1 if lambda is much smaller than y_min)
        # if (p_lower > p_upper) {
        #     p_lower <- 0.99 * p_upper
        # }

        # need deterministic fallback bc ppois(y_min - 1, lambda) can return ~1 if lambda is much smaller than y_min
        if (p_lower >= p_upper) {
            return(y_min)
        }

        # draw u from uniform across p_lower and p_upper
        u <- runif(1, min = p_lower, max = p_upper)

        # inverse CDF
        x <- qpois(u, lambda)

        if (x < y_min) {
            x <- y_min
        }

        return(x)

    }
)

registerDistributions(list(
  dTruncPoisson = list(
    BUGSdist = 'dTruncPoisson(lambda, y_min)',
    Rdist = 'dTruncPoisson(lambda, y_min)',
    types = c('value = integer(0)', 'lambda = double(0)', 'y_min = integer(0)')
  )
))


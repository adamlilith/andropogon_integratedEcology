#' False hurdle, zero-inflated Poisson distribution
#'
#' Distribution and random number generator function for a hurdle, zero-inflated Poisson distribution. The zero-inflation is reflected in supplying the functions a value `z`, which corresponds to a zero value of the Poisson when zero and a non-zero value when one. The Poisson is truncated because zeros are only possible when `z` is zero (i.e., so exogenous to the distribution function).
#'
#' x 		Numeric vector.
#' p 		Numeric vector of probabilities.
#' q 		Numeric vector of quantiles.
#' z 		Numeric equal to 0 or 1, indicating absence or presence.
#' lambda 	Numeric >= 0. Expected density.
#' log 		Logical (if `TRUE` return the log value), or 0 (`FALSE`)/1 (`TRUE`).
#' n 		Number of values.
dFalseHurdleZIP <- nimbleFunction(
    run = function(
        x = double(),
        lambda = double(),
        z = double(),
        log = integer(default = 0)
    ) {

    returnType(double())

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
			# truncated Poisson (x > 0)
			prob_zero <- exp(-lambda)
			dens <- dpois(x, lambda) / (1 - prob_zero)
		}
	}
	if (log) dens <- log(dens)
	return(dens)

    },
    buildDerivs = 'run'
)

# # # # RNG for truncated, pseudo-zero-inflated Poisson
# # # rtruncpseudozip <- nimbleFunction(
# # #     run = function(
# # #         n = integer(),
# # #         lambda = double(),
# # #         z = double()
# # #     ) {
        
# # #     returnType(double())

# # # 	if (z == 0) {
# # # 		out <- 0.0
# # # 	} else {
# # # 		out <- rpois(1, lambda = lambda)
# # # 		while (out == 0) {
# # # 			out <- rpois(1, lambda = lambda)
# # # 		}
# # # 	}
# # # 	return(out)

# # #     } # EOF
# # # )

# RNG for truncated, pseudo-zero-inflated Poisson
rFalseHurdleZIP <- nimbleFunction(
    run = function(
        n = integer(),
        lambda = double(),
        z = double()
    ) {
        
    returnType(double())

	if (z == 0) {
		out <- 0.0
	} else {

		# CDF across [1, ~Inf)
		p_lower <- ppois(0, lambda)
        # p_upper <- 0.9999
        p_upper <- 0.999999

        if (p_lower >= p_upper) {
            return(1)
        }

        # inverse CDF
		out <- 0.0
		while (out < 1) {

			# draw u from uniform across p_lower and p_upper
			u <- runif(1, min = p_lower, max = p_upper)
	        out <- qpois(u, lambda)
			
		}
	}
	return(out)

    } # EOF
)

# # # # PMF function for a zero-inflated Poisson
# # # pzip <- nimbleFunction(
# # #     run = function(
# # #         q = double(),
# # #         lambda = double(),
# # #         z = double(),
# # #         lower.tail = integer(default = 1),
# # #         log.p = integer(default = 0)
# # #     ) {
    
# # #     returnType(double())

# # #     if (q < 0) {
# # #         p <- 0.0
# # #     } else if (lower.tail) {
    
# # #         if (q == 0) {
# # #             p <- z + (1 - z) * exp(-lambda)
# # #         } else {
# # #             p <- (1 - z) * ppois(q, lambda)
# # #         }
    
# # #     } else {
    
# # #         if (q == 0) {
# # #             p <- 1 - (z + (1 - z) * exp(-lambda))
# # #         } else {
# # #             p <- 1 - (1 - z) * ppois(q, lambda)
# # #         }
    
# # #     }

# # #     if (log.p) {
# # #         return(log(p))
# # #     } else {
# # #         return(p)
# # #     }

# # #     } # EOF
# # # )

# # # qzip <- nimbleFunction(
# # #     run = function(
# # #         p = double(),
# # #         lambda = double(),
# # #         z = double(),
# # #         lower.tail = integer(default = 1),
# # #         log.p = integer(default = 0)
# # #     ) {

# # #     returnType(integer())

# # #     if (log.p) p[i] <- exp(p[i])

# # #     # invalid
# # #     if (p < 0 | p > 1) {
# # #         out <- NA_integer_
# # #     } else {

# # #         if (!lower.tail) {
# # #             p <- 1 - p
# # #         }

# # #         # mass at zero
# # #         p0 <- z + (1 - z) * exp(-lambda)

# # #         if (p <= p0) {
# # #             out <- 0L
# # #         } else {
# # #             adj_p <- (p - z) / (1 - z)
# # #             out <- qpois(p = adj_p, lambda = lambda)
# # #         }
# # #     }
# # #     return(out)

# # #     } # EOF
# # # )

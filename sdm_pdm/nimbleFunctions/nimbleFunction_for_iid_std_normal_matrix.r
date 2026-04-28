#' IID standard normal matrix distribution
#'
#' Density and RNG for a matrix whose elements are independent and identically distributed standard normal values (i.e., y ~ N(0, 1)). Input and outputs are matrices. Used in place of the much slower loop over y[i, j] ~ dnorm() and somewhat slower y[i, 1:n_facets] ~ dmnorm() when the covariance structure is diagonal ones and off-diagonal zeros. dIIDStandardNorm() does the sampling across the entire matrix, whereas dIIDStandardNorm_row does it for a single row.
#' x            Numeric matrix with dimensions n_row x n_col. Values to evaluate.
#' n_row        Integer. Number of rows in x.
#' n_col        Integer. Number of columns in x.
#' log          Integer (0/1) or logical. If 1/TRUE, return log-density.
#' n            Must be 1.
#'
#' Use like: z_county[1:n_counties, 1:n_facets] ~ dIIDStandardNorm(n_row = n_counties, n_col = n_facets)
dIIDStandardNorm <- nimbleFunction(
    run = function(
        x = double(2),
        n_row = integer(0),
        n_col = integer(0),
        log = integer(0, default = 0)
    ) {

    returnType(double(0))

    n_vals <- n_row * n_col
    log_dens <- -0.5 * n_vals * log(2.0 * pi) - 0.5 * sum(x[1:n_row, 1:n_col]^2)

    if (log) {
        return(log_dens)
    } else {
        return(exp(log_dens))
    }

    },
    buildDerivs = 'run'
)

# RNG for dIIDStandardNorm. Use like: z_county[1:n_counties, 1:n_facets] ~ dIIDStandardNorm(n_row = n_counties, n_col = n_facets)
rIIDStandardNorm <- nimbleFunction(
    run = function(
        n = integer(0), # must be 1
        n_row = integer(0),
        n_col = integer(0)
    ) {

    returnType(double(2))

    n_vals <- n_row * n_col
    vals <- rnorm(n_vals, mean = 0.0, sd = 1.0)
    out <- matrix(vals, nrow = n_row, ncol = n_col)
    return(out)

    }
)

dIIDStandardNorm_row <- nimbleFunction(
    run = function(
        x = double(1),
        n_col = integer(0),
        log = integer(0, default = 0)
    ) {

    returnType(double(0))
    log_dens <- -0.5 * n_col * log(2.0 * pi) - 0.5 * sum(x[1:n_col]^2)
    if (log) return(log_dens)
    return(exp(log_dens))

    }
)

rIIDStandardNorm_row <- nimbleFunction(
    run = function(
        n = integer(0), # must be 1
        n_col = integer(0)
    ) {

    returnType(double(1))
    return(rnorm(n_col, mean = 0.0, sd = 1.0))

    }
)

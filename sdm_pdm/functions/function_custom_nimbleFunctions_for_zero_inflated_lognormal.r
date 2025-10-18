#' Zero-inflated lognormal distribution
#'
#' @param x Numeric vector.
#' @param p Numeric vector of probabilities.
#' @param q Numeric vector of quantiles.
#' @param pzero Numeric >=0 and <= 1. Probability of 0's.
#' @param meanlog Mean of of the log of the variable.
#' @param sdlog Standard deviation the log of the variable.
#' @param log Logical (if `TRUE` return the log value), or 0 (`FALSE`)/1 (`TRUE`).
#' @param n Number of values.
dziln <- nimbleFunction(
	run = function(
		x = double(0),
		meanlog = double(0),
		sdlog = double(0),
		pzero = double(0),
		log = integer(0, default = 0)
	) {

	returnType(double(0))

	if (sdlog <= 0) {
		if (log) {
			return(-Inf)
		} else {
			return(0.0)
		}
	}

	if (x == 0) {
		dens <- pzero
	} else if (x > 0) {
		dens <- (1 - pzero) * dlnorm(x, meanlog = meanlog, sdlog = sdlog)
	} else {
		dens <- -Inf  # only defined for x >= 0
	}
    if (log) dens <- log(dens)
    return(dens)

	},
    buildDerivs = 'run'
)

rziln <- nimbleFunction(
	run = function(
		n = integer(0),
		meanlog = double(0),
		sdlog = double(0),
		pzero = double(0)
	) {
	
	returnType(double(0))

	if (sdlog <= 0) {
		if (log) return(-Inf) else return(0.0)
	}

	if (n != 1) print('rziln(0) only allows n = 1.')
	if (runif(1, 0, 1) < pzero) {
		out <- 0.0
	} else {
		out <- rlnorm(1, meanlog = meanlog, sdlog = sdlog)
	}
	return(out)
	
	} # EOF
)

pziln <- nimbleFunction(
	run = function(
		q = double(0), 
		meanlog = double(0), 
		sdlog = double(0), 
		pzero = double(0), 
		lower.tail = integer(0, default = 1), 
		log.p = integer(0, default = 0)
	) {
	
	returnType(double(0))

	if (sdlog <= 0) {
		if (log) return(-Inf) else return(0.0)
	}

	p <- 0.0
	if (q < 0) {
		p <- 0.0
	} else if (q == 0) {
		p <- pzero
	} else {
		p <- pzero + (1 - pzero) * plnorm(q, meanlog = meanlog, sdlog = sdlog, lower.tail = 1, log.p = 0)
	}

	if (!lower.tail) p <- 1 - p
	if (log.p) p <- log(p)

	return(p)
	
	}
)

qziln <- nimbleFunction(
	run = function(
		p = double(0), 
		meanlog = double(0), 
		sdlog = double(0), 
		pzero = double(0), 
		lower.tail = integer(0, default = 1), 
		log.p = integer(0, default = 0)
	) {
	
	returnType(double(0))

	if (sdlog <= 0) {
		if (log) return(-Inf) else return(0.0)
	}

	if (log.p) p <- exp(p)
	if (!lower.tail) p <- 1 - p

	x <- 0.0
	if (p <= pzero) {
		x <- 0.0
	} else {
		adj_p <- (p - pzero) / (1 - pzero)
		x <- qlnorm(adj_p, meanlog = meanlog, sdlog = sdlog, lower.tail = 1, log.p = 0)
	}
	return(x)

	} # EOF
)

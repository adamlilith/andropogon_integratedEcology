#' Calculate logged versions of variables
#' 
#' This function takes a data frame and calculates logged versions of specified variables, adding them as new columns to the data frame. The logged values are calculated using the formula log10(x + 1) to handle zero values appropriately.
#'
#' x		data.frame, data.table, or SpatVector
#' vars		Character vector of variable names to log-transform (default: paste0('bio', c(12:14, 16:19)))
calculate_logged_vars <- function(x, vars = paste0('bio', c(12:14, 16:19))) {

	for (var in vars) {
		y <- x[[var]]
		if (inherits(y, c('data.frame', 'data.table', 'SpatVector'))) y <- unlist(y)
		y <- log10(y + 1)
		x$DUMMY <- y
		names(x)[ncol(x)] <- paste0(var, '_log10p1')
	}
	x

}

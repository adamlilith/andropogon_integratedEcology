#' Delineate range core
#'
#' This function takes as its main argument a SpatVector with predicted values and returns a subset that represents the range core.`
#' pred_vect 		A SpatVector with predicted values, including a column with the number of estimated AG.
#' column			Name of column used to define core.
#' core_quant	A numeric value between 0 and 1 representing the quantile threshold for core delineation.
#' rule '>=' or '<='
#' @returns A SpatVector.
delineate_range_core <- function(pred_vect, column, core_quant = 0.95, rule = '>=') {
	
	x <- pred_vect[[column]]
	x <- x[ , 1, drop = TRUE]
	# if (any(names(x) == 'focal_region')) x <- x[pred_vect$focal_region]
	threshold <- quantile(x, core_quant, na.rm = TRUE)
	range_core <- if (rule == '>=') {
		pred_vect[pred_vect[[column]] >= threshold]
	} else if (rule == '<=') {
		pred_vect[pred_vect[[column]] <= threshold]
	} else if (rule == '<') {
		pred_vect[pred_vect[[column]] < threshold]
	} else {
		stop('Bad `rule` for delineating range core.')
	}
	range_core <- aggregate(range_core)
	range_core <- buffer(range_core, 1000)
	range_core <- buffer(range_core, -1000)
	range_core <- aggregate(range_core)
	range_core

}

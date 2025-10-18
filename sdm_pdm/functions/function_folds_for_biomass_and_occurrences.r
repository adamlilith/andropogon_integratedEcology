#' Define folds for biomass for use in nimble::runCrossValidate().
#'
#' @param i Number in 1:4: fold.
folds_for_biomass_and_occurrences <- function(i) {

	folds_biomass <- folds_for_biomass(i)
	folds_occs <- folds_for_occurences(i)
	c(folds_occs, folds_biomass)

}


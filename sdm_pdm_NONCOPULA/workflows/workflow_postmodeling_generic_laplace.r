#' Post-modeling workflow for any model.
#'
#' @param facet Name of the facet being modeled (e.g., occurrence, biomass, plant height, etc.)
#' @param laplace `list` from `runLaplace()`.
#' @param formulae `List` of model formula.
#' @param descrip Textual description of the model
#' @param homoscedastic `TRUE` if model is homoscedastic.
#' @param zero_inflated `TRUE` if model is zero-inflated.
#' @param out_dir Folder in which to save laplace.
workflow_postmodeling_generic_laplace <- function(facet, laplace, formulae, descrip, homoscedastic, zero_inflated, out_dir) {

	meta <- list(
		facets = facet,
		descrip = descrip,
		homoscedastic = homoscedastic,
		zero_inflated = zero_inflated,
		formulae = formulae,
		log_lik = laplace$summary$logLik,
		df = laplace$summary$df,
		convergence = laplace$MLE$convergence,
		message = laplace$MLE$message,
		coefficients = laplace$summary$params
	)

	saveRDS(meta, paste0(out_dir, '/!meta_generic.rds'))
	sink(paste0(out_dir, '/!meta_generic.txt'), split = TRUE)
		print(meta)
	sink()

}

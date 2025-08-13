#' Creates an array of model matrices for a response curve. If there is >1 covariate, then the output is a 3D array, with one "page" per covariate. If there is just one covariate, then the output is a single model matrix. Within each model matrix, the focal covariate increases from low to high across the rows with their range given by the range across present-day and future environments, while the other covariates are held constant at their median value across counties with at least one observed AG.
#'
#' @param formula: RHS formula with intercept
#' @param  centers: named numeric vector of medians
#' @param  scales: named numeric vector of standard deviations
#' @param  ag_vect_sq SpatVector of AG presences with covariates
#' @param  vects: a `list` of other `SpatVector`s with the same names as the covariates in the formula. These are the covariates that will be used to create the range of the response curve (along with `ag_vect_sq`). Ignored if `NULL`.
#' @param n_response_curve_values Number of rows (values along the focal variable) in the response curve matrix. Default is 200.
create_response_curve_array <- function(
	formula,
	centers,
	scales,
	ag_vect_sq,
	vects = NULL,
	n_response_curve_values = 200
) {

	terms <- terms(formula)
	terms <- attr(terms, 'term.labels')
	linear_terms <- terms
	linear_terms <- linear_terms[!grepl(linear_terms, pattern = '\\^2')]
	linear_terms <- linear_terms[!grepl(linear_terms, pattern = '\\:')]
	linear_terms <- linear_terms[!grepl(linear_terms, pattern = '\\*')]

	n_linear_terms <- length(linear_terms)

	if (n_linear_terms == 0) {

		response_curve_x <- matrix(1, nrow = n_response_curve_values, ncol = 1, dimnames = list(1:n_response_curve_values, '(Intercept)'))
		unscaled <- NA

	} else {

		ag_vect_sq_just_occs <- ag_vect_sq
		ag_vect_sq_just_occs <- as.data.frame(ag_vect_sq_just_occs)
		ag_vect_sq_just_occs <- ag_vect_sq_just_occs[ , linear_terms, drop = FALSE]

		# get range of each variable
		mins <- rep(Inf, n_linear_terms)
		maxs <- rep(-Inf, n_linear_terms)
		names(mins) <- names(maxs) <- linear_terms

		for (linear_term in linear_terms) {

			mins[linear_term] <- min(ag_vect_sq_just_occs[[linear_term]], na.rm = TRUE)
			maxs[linear_term] <- max(ag_vect_sq_just_occs[[linear_term]], na.rm = TRUE)

		}

		if (!is.null(vects)) {
			for (i in seq_along(vects)) {
				for (linear_term in linear_terms) {
					vals <- unlist(as.vector(vects[[i]][[linear_term]]))
					vals <- vals[!is.infinite(vals)]
					if (linear_term == 'aridity') vals <- vals[vals <= quantile(vals, 0.99)]
					mins[linear_term] <- min(mins[linear_term], min(vals, na.rm = TRUE))
					maxs[linear_term] <- max(maxs[linear_term], max(vals, na.rm = TRUE))
				}
			}
		}

		# create 3D array for cases where >1 predictor
		if (n_linear_terms > 1) {
		
			env_array <- array(
				NA,
				dim = c(n_response_curve_values, n_linear_terms, n_linear_terms),
				dimnames = list(
					1:n_response_curve_values,
					linear_terms,
					linear_terms
				)
			)

			unscaled <- matrix(NA_real_, nrow = n_response_curve_values, ncol = n_linear_terms, dimnames = list(1:n_response_curve_values, linear_terms))

			# calculate median of each variable across counties with AG presences
			medians <- apply(ag_vect_sq_just_occs, 2, median, na.rm = TRUE)

			for (i in seq_along(linear_terms)) {
				for (j in seq_along(linear_terms)) {
					env_array[ , i, j] <- medians[i]
				}
			}

			for (i in seq_along(linear_terms)) {
				unscaled[ , i] <- env_array[ , i, i] <- seq(mins[i], maxs[i], length.out = n_response_curve_values)
			}

			# scale
			for (i in seq_along(linear_terms)) {

				linear_term <- linear_terms[i]
				env_array[ , , i] <- sweep(env_array[ , , i], 2, centers[linear_terms], '-')
				env_array[ , , i] <- sweep(env_array[ , , i], 2, scales[linear_terms], '/')
				
			}

			# make into model matrix
			mm_array <- list()
			for (i in seq_along(linear_terms)) {
				mm_array[[i]] <- model.matrix(formula, as.data.frame(env_array[ , , i, drop = TRUE]))
			}

			nc <- ncol(mm_array[[1]])
			response_curve_x <- array(
				as.numeric(unlist(mm_array)),
				dim = c(n_response_curve_values, nc, n_linear_terms),
				dimnames = list(1:n_response_curve_values, colnames(mm_array[[1]]), linear_terms)
			)

		} else {

			## 2D matrix for just one predictor variable
			unscaled <- env_array <- matrix(
				seq(mins, maxs, length.out = n_response_curve_values),
				nrow = n_response_curve_values,
				ncol = n_linear_terms,
				dimnames = list(1:n_response_curve_values, linear_terms)
			)

			# scale
			env_array <- scale(env_array, center = centers[linear_terms], scale = scales[linear_terms])
			colnames(env_array) <- linear_terms

			# make into model matrix
			response_curve_x <- model.matrix(formula, as.data.frame(env_array))

		}

	}


	list(
		resp_curves_x_unscaled = unscaled,			# unscaled model matrix for response curve
		response_curves_x_scaled = response_curve_x	# scaled model matrix for response curve
	)

}

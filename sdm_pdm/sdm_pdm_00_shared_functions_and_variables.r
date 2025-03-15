### MODELING ANDROPOGON GERARDI DISTRIBUTION, PHENOTYPE, PHYSIOLOGY, GENOTYPE, and ASSOCIATED MICROBIAL COMMUNITIES
### Erica Newman | Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2023-12
###
### This script contains shared functions and variables for the integrated SDM/PDM scripts.
###
### source('C:/Ecology/R/andropogon_integratedEcology/sdm_pdm/sdm_pdm_00_data_alignment.r')
### source('C:/Subarashi/R/andropogon_integratedEcology/sdm_pdm/sdm_pdm_00_data_alignment.r')
###
### CONTENTS ###
### setup ###
### create_response_curve_array() ###

#############
### setup ###
#############

	library(bayesplot) # MCMC
	library(coda) # MCMC
	library(cowplot) # combining ggplots
	library(data.table) # fast data tables
	library(enmSdmX) # GIS & SDMing
	library(enmSdmX) # GIS & SDMing
	library(ggplot2) # plotting
	library(ggnewscale) # colors
	library(ggspatial) # spatial plotting
	library(nimbleHMC) # Bayesian modeling
	library(omnibus) # utilities
	library(predicts) # GIS & SDMing
	library(readxl) # Excel
	library(terra) # spatial objects


#####################################
### create_response_curve_array() ###
#####################################

	# create_response_curve_array() creates an array of model matrices for a response curve. It is multi-dimensional, with dimensions for each variable in the model. Each slice of the array is a model matrix for a different variable. All variables are held at their mean values, except for the variable of interest, which is varied across a range of values set by the range of environments within counties with at least one observed AG.

	# NB This function assumes that there is more than one variable in the model. If there is only one variable, the function will not work because it will produce a 3D array, whereas a single variable requires a matrix.

	# form: RHS formula
	# centers: named numeric vector of means
	# scales: named numeric vector of standard deviations
	# ag_pres: data frame of counties
	create_response_curve_array <- function(form, centers, scales, ag_pres) {

		terms <- terms(form)
		terms <- attr(terms, 'term.labels')
		linear_terms <- terms
		linear_terms <- linear_terms[!grepl(linear_terms, pattern = '\\^2')]
		linear_terms <- linear_terms[!grepl(linear_terms, pattern = '\\:')]
		linear_terms <- linear_terms[!grepl(linear_terms, pattern = '\\*')]
	
		n_linear_terms <- length(linear_terms)

		env_array <- array(
			NA,
			dim = c(n_response_curve_values, n_linear_terms, n_linear_terms),
			dimnames = list(
				1:n_response_curve_values,
				linear_terms,
				linear_terms
			)
		)

		# calculate means of each variable across counties with AG presences
		ag_pres <- ag_vect_sq[ag_vect_sq$any_ag_quality_1_to_3 > 0]
		ag_pres <- as.data.frame(ag_pres)[ , linear_terms, drop = FALSE]
		means <- colMeans(ag_pres)

		for (i in seq_along(linear_terms)) {
			for (j in seq_along(linear_terms)) {
				env_array[ , i, j] <- means[i]
			}
		}

		mins <- apply(ag_pres[ , linear_terms, drop = FALSE], 2, min)
		maxs <- apply(ag_pres[ , linear_terms, drop = FALSE], 2, max)

		for (i in seq_along(linear_terms)) {
			env_array[ , i, i] <- seq(mins[i], maxs[i], length.out = n_response_curve_values)
		}

		# scale
		env_array_scaled <- env_array
		for (i in seq_along(linear_terms)) {
			env_array_scaled[ , , i] <- sweep(env_array_scaled[ , , i], 2, centers[linear_terms], '-')
			env_array_scaled[ , , i] <- sweep(env_array_scaled[ , , i], 2, scales[linear_terms], '/')
		}

		# make into model matrix
		mm_array <- list()
		for (i in seq_along(linear_terms)) {
			mm_array[[i]] <- model.matrix(form, as.data.frame(env_array_scaled[ , , i, drop = TRUE]))
		}

		nc <- ncol(mm_array[[1]])
		response_curve_x <- array(
			as.numeric(unlist(mm_array)),
			dim = c(n_response_curve_values, nc, n_linear_terms),
			dimnames = list(1:n_response_curve_values, colnames(mm_array[[1]]), linear_terms)
		)

		response_curve_x

	}


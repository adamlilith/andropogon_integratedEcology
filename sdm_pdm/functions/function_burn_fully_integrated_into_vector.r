#' Burn predictions from an occurrence/biomass/non-biomass model into a SpatVector of North America or Dust Bowl map
#'
#' source('C:/Kaji/R/andropogon_integratedEcology/sdm_pdm/functions/function_burn_fully_integrated_into_vector.r')
#'
#' demesne				'nam' (North America) or '1930s' (Dust Bowl area)
#' chains				Chains from NIMBLE
#' formula_occs			Formula for occurrences
#' formula_bias	Formula for bias in sampling occurrences
#' 
#' formula_psi			Formula for presence/absences (assumes that we use the log/unlogged version of precipitation as per formula_occs)
#' 
#' formula_biomass		Formula for biomass
#' resp_distrib_biomass	Distribution for biomass (eg, 'gamma', 'hGamma' (zero-inflated gamma), 'lognormal', or 'hurdleLN' (zero-inflated lognormal))
#' transform_biomass	Named vector of transformations to translate MVN to mean biomass: 'identity', 'softplus' or 'exponential'. The names of the vector should be the same as the response variable for biomass in formula_biomass.
#' log_precip_biomass	TRUE/FALSE: log precipitation predictors for biomass submodel
#' 
#' nonbiomass_facets	Named list of non-biomass facets
#' force_presence		If `TRUE`, force all predictions assuming presence (N > 0). If FALSE, allow presences and absences to be predicted.
burn_fully_integrated_into_vector <- function(
	demesne,
	chains,
	formula_occs,
	formula_bias,
	formula_psi,
	formula_biomass,
	resp_distrib_biomass,
	transform_biomass,
	log_precip_biomass,
	nonbiomass_facets,
	force_presence = FALSE
) {

	### occurrence data
	###################

	pred_vect <- if (demesne == 'nam') {
		pred_vect <- data_occs_counties$ag_vect_sq
		pred_vect <- simplifyGeom(pred_vect, tolerance = 1000)
	} else if (demesne == '1930s') {
		pred_vect <- vect(paste0('./data_from_adam_and_loretta/andropogon_gerardi_occurrences_with_environment_1931_1940_prism.gpkg'))
	}

	### predict present
	###################
	say('   predicting present occurrence, biomass, and non-biomass... (force_presence = ', force_presence, ')')

	if (demesne == 'nam') {
		x_occs <- data_occs_counties$counties_x_sq
		x_psi <- data_occs_counties$counties_x_sq
		x_biomass <- data_biomass_counties$counties_x_sq
	} else if (demesne == '1930s') {
		x_occs <- data_occs_counties$counties_x_thirties
		x_psi <- data_occs_counties$counties_x_thirties
		x_biomass <- data_biomass_counties$counties_x_thirties
	}

	x_nonbiomass <- list()
	for (f in seq_along(nonbiomass_facets)) {
		
		facet <- names(nonbiomass_facets)[f]
		x_nonbiomass[[f]] <- if (demesne == 'nam') {
			data_nonbiomass_counties[[facet]]$counties_x_sq
		} else if (demesne == '1930s') {
			data_nonbiomass_counties[[facet]]$counties_x_thirties
		}

	}
	names(x_nonbiomass) <- names(nonbiomass_facets)

	preds <- predict_fully_integrated(
		chains = chains,
		nonbiomass_facets = nonbiomass_facets,
		resp_distrib_biomass = resp_distrib_biomass,
		transform_biomass = transform_biomass,
		x_occs = x_occs,
		x_biomass = x_biomass,
		x_psi = x_psi,
		x_nonbiomass = x_nonbiomass,
		w_occs = if (demesne == 'nam') data_occs_counties$counties_w_sq else data_occs_counties$counties_w_thirties,
		force_presence = force_presence
	)

	preds_occs <- preds$preds_occs
	preds_biomass <- preds$preds_biomass
	preds_psi <- preds$preds_psi

	pred_vect$DUMMY1 <- colMeans(preds_psi)
	pred_vect$DUMMY2 <- apply(preds_psi, 2, inner_quant)

	pred_vect$DUMMY3 <- colMeans(preds_occs)
	pred_vect$DUMMY4 <- apply(preds_occs, 2, median)
	pred_vect$DUMMY5 <- apply(preds_occs, 2, inner_quant)

	pred_vect$DUMMY6 <- colMeans(preds_biomass)
	pred_vect$DUMMY7 <- apply(preds_biomass, 2, median)
	pred_vect$DUMMY8 <- apply(preds_biomass, 2, inner_quant)

	index <- (ncol(pred_vect) - 7):ncol(pred_vect)
	if (demesne == 'nam') {
		names(pred_vect)[index] <- c(
			'psi_sq', 'psi_inner_quant_sq',
			'N_ag_mean_sq', 'N_ag_median_sq', 'N_ag_inner_quant_sq',
			'biomass_mean_sq', 'biomass_median_sq', 'biomass_inner_quant_sq'
		)
	} else if (demesne == '1930s') {
		names(pred_vect)[index] <- c(
			'psi_1930s', 'psi_inner_quant_1930s',
			'N_ag_mean_1930s', 'N_ag_median_1930s', 'N_ag_inner_quant_1930s',
			'biomass_mean_1930s', 'biomass_median_1930s', 'biomass_inner_quant_1930s'
		)
	}

	for (f in seq_along(nonbiomass_facets)) {

		facet <- names(nonbiomass_facets)[f]
		pred_nonbiomass <- preds$preds_nonbiomass[[facet]]

		pred_vect$DUMMY1 <- colMeans(pred_nonbiomass)
		pred_vect$DUMMY2 <- apply(pred_nonbiomass, 2, median)
		pred_vect$DUMMY3 <- apply(pred_nonbiomass, 2, inner_quant)

		index <- (ncol(pred_vect) - 2):ncol(pred_vect)
		if (demesne == 'nam') {
			names(pred_vect)[index] <- paste0(facet, c('_mean_sq', '_median_sq', '_inner_quant_sq'))
		} else if (demesne == '1930s') {
			names(pred_vect)[index] <- paste0(facet, c('_mean_1930s', '_median_1930s', '_inner_quant_1930s'))
		}

	}

	### future occurrences--only for entire North American domain
	##############################################################
	if (demesne == 'nam') {

		### predict mean to each future
		###############################
		for (fut in futs) {

			say('   predicting ', fut, ' occurrence, biomass, and non-biomass...')
			
			x_occs <- data_occs_counties[paste0('counties_x_', fut)]
			x_occs <- x_occs[[1]]
		
			x_biomass <- data_biomass_counties[paste0('counties_x_', fut)]
			x_biomass <- x_biomass[[1]]
		
			x_psi <- data_psi_counties[paste0('counties_x_', fut)]
			x_psi <- x_psi[[1]]
		
			x_nonbiomass <- list()
			for (f in seq_along(nonbiomass_facets)) {
				
				facet <- names(nonbiomass_facets)[f]
				x_nonbiomass[[f]] <- data_nonbiomass_counties[[facet]][[paste0('counties_x_', fut)]]

			}
			names(x_nonbiomass) <- names(nonbiomass_facets)

			preds <- predict_fully_integrated(
				chains = chains,
				nonbiomass_facets = nonbiomass_facets,
				resp_distrib_biomass = resp_distrib_biomass,
				transform_biomass = transform_biomass,
				x_occs = x_occs,
				x_biomass = x_biomass,
				x_psi = x_psi,
				x_nonbiomass = x_nonbiomass,
				force_presence = force_presence
			)

			# zero-inflation
			preds_psi <- preds$preds_psi
			pred_vect$DUMMY1 <- colMeans(preds_psi)
			pred_vect$DUMMY2 <- apply(preds_psi, 2, inner_quant)

			index <- (ncol(pred_vect) - 1):ncol(pred_vect)
			names(pred_vect)[index] <- paste0(c('psi_', 'psi_inner_quant_'), fut)

			# abundance
			preds_occs <- preds$preds_occs
			pred_vect$DUMMY1 <- colMeans(preds_occs)
			pred_vect$DUMMY2 <- apply(preds_occs, 2, median)
			pred_vect$DUMMY3 <- apply(preds_occs, 2, inner_quant)

			index <- (ncol(pred_vect) - 2):ncol(pred_vect)
			names(pred_vect)[index] <- paste0(c('N_ag_mean_', 'N_ag_median_', 'N_ag_inner_quant_'), fut)

			# biomass
			preds_biomass <- preds$preds_biomass
			pred_vect$DUMMY4 <- colMeans(preds_biomass)
			pred_vect$DUMMY5 <- apply(preds_biomass, 2, median)
			pred_vect$DUMMY6 <- apply(preds_biomass, 2, inner_quant)

			index <- (ncol(pred_vect) - 2):ncol(pred_vect)
			names(pred_vect)[index] <- paste0(c('biomass_mean_', 'biomass_median_', 'biomass_inner_quant_'), fut)

			for (f in seq_along(nonbiomass_facets)) {

				facet <- names(nonbiomass_facets)[f]
				pred_nonbiomass <- preds$preds_nonbiomass[[facet]]

				pred_vect$DUMMY1 <- colMeans(pred_nonbiomass)
				pred_vect$DUMMY2 <- apply(pred_nonbiomass, 2, median)
				pred_vect$DUMMY3 <- apply(pred_nonbiomass, 2, inner_quant)

				index <- (ncol(pred_vect) - 2):ncol(pred_vect)
				names(pred_vect)[index] <- paste0(facet, c('_mean_', '_median_', '_inner_quant_'), fut)

			}

		} 

	} # if demesne is North America

	pred_vect

}

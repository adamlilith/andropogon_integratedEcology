### MODELING ANDROPOGON GERARDI DISTRIBUTION, PHENOTYPE, PHYSIOLOGY, GENOTYPE, and ASSOCIATED MICROBIAL COMMUNITIES
### Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2026-03
###
### This script conducts external validation of Andropogon gerardi models.
###
### source('C:/Kaji/R/andropogon_integratedEcology/sdm_pdm/sdm_pdm_80_external_validation.r')
###
### CONTENTS ###
### setup ###
### compare models against AG height from Figure 3 in McMillan (1964) ###
### compare predictions to data from Griffan-Nolan & Sandel (2023) and Griffan-Nolan et al. (2025) ###

#############
### setup ###
#############

	rm(list = ls())

	setwd('C:/Kaji/Research/Andropogon/Andropogon')
	source('C:/Kaji/R/andropogon_integratedEcology/sdm_pdm/sdm_pdm_00_shared_functions_and_variables.r')

	set.seed(1)

say('#########################################################################')
say('### compare models against AG height from Figure 3 in McMillan (1964) ###')
say('#########################################################################')

	library(wCorr) # weighted rank correlation

	# Compare model predictions at sites from where McMillan collected AG seed to AG height data in Figure 3 of McMillan (1964). We will use the slope of a linear mixed model (obs ~ predicted | site) and the polyserial correlation (rank variable vs continuous variable) between observed and predicted. Observed values were binned into ordinal classes by McMillan, so the regression assumes a Poisson response.

	test_data <- vect('./data_from_mcmillan/mcmillan_1964_fig_3_ag_prism_1952_1961_soilgrids_top_5_cm.gpkg')
	test_data <- calculate_logged_vars(test_data)
	obs <- test_data[ , c('id', paste0('height', 1:3))]
	obs <- as.data.table(obs)

	centers_scales <- fread('./outputs_loretta/integrated_sdm_pdm/centers_and_scales_for_covariates.csv')

	### unimodal models
	###################
		
		model_folders <- listFiles('./outputs_loretta/integrated_sdm_pdm/models_height', pattern = '\\[height\\~')

		summary <- data.table()
		for (i in seq_along(model_folders)) {

			model_folder <- model_folders[i]
			say(model_folder)

			# posterior samples
			chains <- readRDS(paste0(model_folder, '/chains.rds'))

			# metadata
			meta_facet <- readRDS(paste0(model_folder, '/!meta_height.rds'))
			formula_facet <- meta_facet$formulae$formula_facet
			resp_distrib <- meta_facet$resp_distrib
			transform <- meta_facet$transform

			terms <- terms(formula_facet)
			terms <- attr(terms, 'term.labels')

			covariates <- terms
			covariates <- covariates[!grepl(covariates, pattern = '\\^2')]
			covariates <- covariates[!grepl(covariates, pattern = '\\:')]
			covariates <- covariates[!grepl(covariates, pattern = '\\*')]

			mm <- test_data[ , covariates, drop = FALSE]
			mm <- as.data.table(mm)

			for (covariate in covariates) {

				center <- centers_scales$center_sites[centers_scales$variable == covariate]
				scale <- centers_scales$scale_sites[centers_scales$variable == covariate]

				mm[ , covariate] <- (mm[[covariate]] - center) / scale

			}

			x <- model.matrix(formula_facet, data = mm)

			# predict facet model
			preds <- predict_nonbiomass_single_trait(chains, x = x, resp_distrib = resp_distrib, transform = transform)
			predictions <- colMeans(preds, na.rm = TRUE)

			# calculate rank correlation between McMillan's height classes and site-level height predictions using weighted correlation to reflect different sample sizes at McMillan sites

			# classify predictions by scheme used in McMillan (1964)
			obs_pred <- melt(
				as.data.table(test_data[ , c('id', 'height1', 'height2', 'height3')]),
				id.vars = 'id',
				measure.vars = c('height1', 'height2', 'height3'),
				value.name = 'observed_class'
			)[ , .(id, observed_class)]
			
			obs_pred[ , id := factor(id)]

			obs_pred[ , predicted := predictions[match(obs_pred$id, test_data$id)]]
			obs_pred <- obs_pred[complete.cases(obs_pred)]
			obs_pred[ , site_n := .N, by = id]
			obs_pred[ , w := 1 / site_n]
			obs_pred[ , observed_class := factor(observed_class, level = 1:10, ordered = TRUE)]

			rho <- weightedCorr(obs_pred$predicted, obs_pred$observed_class, weights = obs_pred$w, method = 'Polyserial')

			formula_facet_char <- paste(as.character(formula_facet), collapse = ' ')

			summary <- rbind(
				summary,
				data.table(
					facet = 'height',
					model_type = 'unimodal',
					model_folder = basename(model_folder),
					formula_facet = formula_facet_char, 
					resp_distrib = resp_distrib,
					transform = transform,
					rho = rho
				)
			)

		} # next model

		fwrite(summary, './outputs_loretta/integrated_sdm_pdm/models_height/validation_vs_mcmillan_1964_height_unimodal_models.csv')

	### integrated models
	#####################

		# only using morphological model since that's the only one that included height

		model_dir <- './outputs_loretta/integrated_sdm_pdm/models_integrated/non_centered_MVN([occs~hurdlepoisson_offset]_[biomass~hurdleln]_[bla_can_hei])_eta=2'

		say('Predicting to:\n', model_dir)

		formulae <- readRDS(paste0(model_dir, '/formulae.rds'))

		this_nonbiomass_facets <- nonbiomass_facets[c('blade_width', 'canopy_diameter', 'height')]

		chains <- readRDS(paste0(model_dir, '/chains.rds'))

		### prepare data at sites from which McMillan collected clones
		##############################################################

		# prepare predictors for occurrence
		terms <- attr(terms(formulae$formula_occs), 'term.labels')
		terms_in_df <- terms[terms %in% names(test_data)]
		centers <- centers_scales$center_occs[match(terms_in_df, centers_scales$variable)]
		scales <- centers_scales$scale_occs[match(terms_in_df, centers_scales$variable)]
		names(centers) <- names(scales) <- terms_in_df

		x_occs <- as.data.table(test_data)
		x_occs <- x_occs[ , ..terms_in_df]
		x_occs <- scale(x_occs, center = centers, scale = scales)
		x_occs <- as.data.frame(x_occs)
		x_occs <- model.matrix(formulae$formula_occs, x_occs)

		# prepare predictors for biomass
		terms <- attr(terms(formulae$formula_biomass), 'term.labels')
		terms_in_df <- terms[terms %in% names(test_data)]
		centers <- centers_scales$center_sites[match(terms_in_df, centers_scales$variable)]
		scales <- centers_scales$scale_sites[match(terms_in_df, centers_scales$variable)]
		names(centers) <- names(scales) <- terms_in_df

		x_biomass <- as.data.table(test_data)
		x_biomass <- x_biomass[ , ..terms_in_df]
		x_biomass <- scale(x_biomass, center = centers, scale = scales)
		x_biomass <- as.data.frame(x_biomass)
		x_biomass <- model.matrix(formulae$formula_biomass, x_biomass)

		# prepare predictors for non-biomass facets
		x_nonbiomass <- list()
		for (f in seq_along(this_nonbiomass_facets)) {

			facet <- names(this_nonbiomass_facets)[f]
			formula_facet <- this_nonbiomass_facets[[facet]]$formula

			terms <- attr(terms(formula_facet), 'term.labels')
			terms_in_df <- terms[terms %in% names(test_data)]
			centers <- centers_scales$center_sites[match(terms_in_df, centers_scales$variable)]
			scales <- centers_scales$scale_sites[match(terms_in_df, centers_scales$variable)]
			names(centers) <- names(scales) <- terms_in_df

			x_facet <- as.data.table(test_data)
			x_facet <- x_facet[ , ..terms_in_df]
			x_facet <- scale(x_facet, center = centers, scale = scales)
			x_facet <- as.data.frame(x_facet)
			x_facet <- model.matrix(formula_facet, x_facet)

			x_nonbiomass[[facet]] <- x_facet

		}

		preds <- predict_fully_integrated(
			chains = chains,
			nonbiomass_facets = this_nonbiomass_facets,
			resp_distrib_biomass = formulae$meta_biomass$resp_distrib_biomass,
			transform_biomass = formulae$meta_biomass$transform_biomass,
			x_occs = x_occs,
			x_biomass = x_biomass,
			x_psi = x_occs,
			x_nonbiomass = x_nonbiomass,
			w_occs = NULL,
			force_presence = TRUE
		)

		preds_height <- preds$preds_nonbiomass$height
		predictions <- colMeans(preds_height, na.rm = TRUE)

		# compare facet model predictions to observed height classes
		obs_pred <- melt(
			as.data.table(test_data[ , c('id', 'height1', 'height2', 'height3')]),
			id.vars = 'id',
			measure.vars = c('height1', 'height2', 'height3'),
			value.name = 'observed_class'
		)[ , .(id, observed_class)]
		
		obs_pred[ , id := factor(id)]

		obs_pred[ , predicted := predictions[match(obs_pred$id, test_data$id)]]
		obs_pred <- obs_pred[complete.cases(obs_pred)]
		obs_pred[ , site_n := .N, by = id]
		obs_pred[ , w := 1 / site_n]
		obs_pred[ , observed_class := factor(observed_class, level = 1:10, ordered = TRUE)]

		rho <- weightedCorr(obs_pred$predicted, obs_pred$observed_class, weights = obs_pred$w, method = 'Polyserial')

		formula_facet <- paste('~', as.character(formulae$nonbiomass_facets[['height']]$formula))[2]
		resp_distrib <- nonbiomass_facets[['height']]$resp_distrib
		transform <- nonbiomass_facets[['height']]$transform

		summary <- data.table(
			facet = 'height',
			model_type = 'integrated',
			model_folder = basename(model_dir),
			formula_facet = formula_facet, 
			resp_distrib = resp_distrib,
			transform = transform,
			rho = rho
		)
		
	fwrite(summary, paste0(model_dir, '/validation_vs_mcmillan_1964_height_integrated_model.csv'))

# say('######################################################################################################')
# say('### compare predictions to data from Griffan-Nolan & Sandel (2023) and Griffan-Nolan et al. (2025) ###')
# say('######################################################################################################')

# 	# Compare model predictions to data on AG height from Griffan-Nolan, R.J. and Sandel, B. 2023. Global intraspecific trait-climate relationships for grasses are linked to a species' typical form and function. Ecography 2023:e06586. Acquired directly from Brody Sandel. For records that have the same longitude/latitude, we will assume they are from the same site and collected at the same time, so will use the reciprocal of the number of samples at a "site" (same long/lat) as a weight in the model so that sites with more samples count the same as those with fewer samples.

	# # for cluster-robust standard errors from linear models
	# library(sandwich)
	# library(lmtest)

# 	# data
# 	test_data <- fread('./data_from_brody/data/Height.csv')
# 	test_data <- test_data[complete.cases(test_data[ , c('Longitude', 'Latitude')])]
# 	test_data[ , Height := 100 * Height] # convert to cm
# 	if (any(test_data$Longitude > 0)) test_data$Longitude[test_data$Longitude > 0] <- -1 * test_data$Longitude[test_data$Longitude > 0] # coercing one pos long to neg
# 	test_data[ , site := .GRP, by = .(Longitude, Latitude)]
# 	test_data[ , site_n := .N, by = site]
	
# 	# extract ClimateNA climate values at these sites
# 	test_vect <- vect(test_data, geom = c('Longitude', 'Latitude'), keepgeom = TRUE, crs = getCRS('WGS84'))
# 	prec_dms <- coordImprecision(test_vect, dms = TRUE)
# 	prec_decimal <- coordImprecision(test_vect, dms = FALSE)
# 	prec_m <- pmax(prec_dms, prec_decimal, rep(1000, length(prec_dms))) # minimum buffer of 1 km
# 	test_vect <- buffer(test_vect, prec_m)

# 	# test_vect <- test_vect[test_vect$id != 1] # removing this record bc it encompasses another with better precision
# 	cna <- listFiles('C:/Kaji/Research Data/ClimateNA/v 7.3 AdaptWest/1991-2020', pattern = '.tif')
# 	cna <- rast(cna)
# 	names(cna)[grepl(sources(cna), pattern = 'PPT')] <- paste0('ppt', prefix(1:12, 2))
# 	names(cna)[grepl(sources(cna), pattern = 'Tave')] <- paste0('tmean', prefix(1:12, 2))
# 	names(cna)[grepl(sources(cna), pattern = 'Tmax')] <- paste0('tmax', prefix(1:12, 2))
# 	names(cna)[grepl(sources(cna), pattern = 'Tmin')] <- paste0('tmin', prefix(1:12, 2))
	
# 	test_vect <- project(test_vect, cna)
# 	test_clim <- extract(cna, test_vect, exact = TRUE, fun = mean, ID = FALSE)
# 	names(test_clim) <- names(cna)

# 	# BIOCLIMS
# 	ppt <- test_clim[ , paste0('ppt', prefix(1:12, 2))]
# 	tmin <- test_clim[ , paste0('tmin', prefix(1:12, 2))]
# 	tmax <- test_clim[ , paste0('tmax', prefix(1:12, 2))]

# 	ppt <- as.matrix(ppt)
# 	tmin <- as.matrix(tmin)
# 	tmax <- as.matrix(tmax)

# 	bioclims <- predicts::bcvars(ppt, tmin, tmax)
# 	bioclims <- as.data.frame(bioclims)
# 	names(bioclims) <- paste0('bio', 1:19)
# 	bioclims$aridity <- (bioclims$bio1 + 10) / ((bioclims$bio12 + 1) / 10)
# 	bioclims <- calculate_logged_vars(bioclims)

# 	# soil from SoilGrids 2.0
# 	sand <- rast('C:/Kaji/Research Data/SoilGrids/SoilGrids 2.0/sand_0-5cm_mean_northAmerica.tif')
# 	nitrogen <- rast('C:/Kaji/Research Data/SoilGrids/SoilGrids 2.0/nitrogen_0-5cm_mean_northAmerica.tif')
# 	ph <- rast('C:/Kaji/Research Data/SoilGrids/SoilGrids 2.0/phh2o_0-5cm_mean_northAmerica.tif')
# 	soil <- c(sand, nitrogen, ph)
# 	names(soil) <- c('sand', 'nitrogen', 'ph')

# 	test_vect <- project(test_vect, soil)
# 	soil <- crop(soil, test_vect)
# 	test_soil <- extract(soil, test_vect, exact = TRUE, fun = mean, ID = FALSE, na.rm = TRUE)
# 	test_soil$sand <- test_soil$sand / 1000
# 	test_soil$nitrogen <- test_soil$nitrogen / 1000
# 	test_soil$ph <- test_soil$ph / 10

# 	# data frame with raw (uncentered) climate and soil data for each measurement
# 	x_raw <- cbind(bioclims, test_soil)
# 	centers_scales <- fread('./outputs_loretta/integrated_sdm_pdm/centers_and_scales_for_covariates.csv')

# 	### unimodal models
# 	###################
		
# 		model_folders <- listFiles('./outputs_loretta/integrated_sdm_pdm/models_height', pattern = '\\[height\\~')

# 		summary <- data.table()
# 		for (i in seq_along(model_folders)) {

# 			model_folder <- model_folders[i]
# 			say(model_folder)

# 			# posterior samples
# 			chains <- readRDS(paste0(model_folder, '/chains.rds'))

# 			# metadata
# 			meta_facet <- readRDS(paste0(model_folder, '/!meta_height.rds'))
# 			formula_facet <- meta_facet$formulae$formula_facet
# 			resp_distrib <- meta_facet$resp_distrib
# 			transform <- meta_facet$transform

# 			terms <- terms(formula_facet)
# 			terms <- attr(terms, 'term.labels')

# 			covariates <- terms
# 			covariates <- covariates[!grepl(covariates, pattern = '\\^2')]
# 			covariates <- covariates[!grepl(covariates, pattern = '\\:')]
# 			covariates <- covariates[!grepl(covariates, pattern = '\\*')]

# 			mm <- x_raw[ , covariates, drop = FALSE]
# 			mm <- as.data.table(mm)

# 			for (covariate in covariates) {

# 				center <- centers_scales$center_sites[centers_scales$variable == covariate]
# 				scale <- centers_scales$scale_sites[centers_scales$variable == covariate]

# 				mm[ , covariate] <- (mm[[covariate]] - center) / scale

# 			}

# 			x <- model.matrix(formula_facet, data = mm)

# 			# predict facet model
# 			preds <- predict_nonbiomass_single_trait(chains, x = x, resp_distrib = resp_distrib, transform = transform)
# 			predictions <- colMeans(preds, na.rm = TRUE)
			
# 			# compare facet model predictions to observed height classes
# 			obs_pred <- data.table(
# 				site = factor(test_data$site),
# 				observed = test_data$Height,
# 				predicted = predictions,
# 				w = 1 / test_data$site_n
# 			)
# 			obs_pred <- obs_pred[is.finite(predicted)]

# 			# calibration model
# 			calib <- lm(observed ~ predicted, data = obs_pred, weights = w)

# 			r_squared <- summary(calib)$r.squared
# 			coefs <- coefficients(calib)
# 			intercept <- coefs['(Intercept)']
# 			slope <- coefs['predicted']

# 			# cluster-robust SEs			
# 			vcv <- vcovCL(calib, cluster = ~ site)
# 			robust_se <- coeftest(calib, vcov = vcv, type = 'HC2')
# 			se_slope <- robust_se['predicted', 'Std. Error']
# 			se_intercept <- robust_se['(Intercept)', 'Std. Error']

# 			intercept_lower <- intercept - 1.96 * se_intercept
# 			intercept_upper <- intercept + 1.96 * se_intercept
# 			slope_lower <- slope - 1.96 * se_slope
# 			slope_upper <- slope + 1.96 * se_slope

# 			slope_sig <- ifelse(slope_lower > 1 | slope_upper < 1, '*', '-')
# 			intercept_sig <- ifelse(intercept_lower > 1 | intercept_upper < 1, '*', '-')

# 			formula_facet_char <- paste(as.character(formula_facet), collapse = ' ')

# 			summary <- rbind(
# 				summary,
# 				data.table(
# 					facet = 'height',
# 					model_type = 'unimodal',
# 					model_folder = basename(model_folder),
# 					formula_facet = formula_facet_char, 
# 					resp_distrib = resp_distrib,
# 					transform = transform,
# 					intercept = intercept,
# 					intercept_lower = intercept_lower,
# 					intercept_upper = intercept_upper,
# 					intercept_sig = intercept_sig,
# 					slope = slope,
# 					slope_lower = slope_lower,
# 					slope_upper = slope_upper,
# 					slope_sig = slope_sig,
# 					r_squared = r_squared
# 				)
# 			)

# 		} # next model

# 		fwrite(summary, './outputs_loretta/integrated_sdm_pdm/models_height/validation_vs_sandel_and_griffan_nolan_height_unimodal_models.csv')

# 	### integrated model
# 	####################

# 		# only using morphological model since that's the only one that included height

# 		model_dir <- './outputs_loretta/integrated_sdm_pdm/models_integrated/non_centered_MVN([occs~hurdlepoisson_offset]_[biomass~hurdleln]_[bla_can_hei])_eta=2'
# 		say('Predicting to:\n', model_dir)

# 		formulae <- readRDS(paste0(model_dir, '/formulae.rds'))
# 		this_nonbiomass_facets <- nonbiomass_facets[c('blade_width', 'canopy_diameter', 'height')]

# 		chains <- readRDS(paste0(model_dir, '/chains.rds'))

# 		### prepare data at sites from which McMillan collected seeds
# 		#############################################################

# 		# prepare predictors for occurrence
# 		terms <- attr(terms(formulae$formula_occs), 'term.labels')
# 		terms_in_df <- terms[terms %in% names(x_raw)]
# 		centers <- centers_scales$center_occs[match(terms_in_df, centers_scales$variable)]
# 		scales <- centers_scales$scale_occs[match(terms_in_df, centers_scales$variable)]
# 		names(centers) <- names(scales) <- terms_in_df

# 		x_occs <- as.data.table(x_raw)
# 		x_occs <- x_occs[ , ..terms_in_df]
# 		x_occs <- scale(x_occs, center = centers, scale = scales)
# 		x_occs <- as.data.frame(x_occs)
# 		x_occs <- model.matrix(formulae$formula_occs, x_occs)

# 		# prepare predictors for biomass
# 		terms <- attr(terms(formulae$formula_biomass), 'term.labels')
# 		terms_in_df <- terms[terms %in% names(x_raw)]
# 		centers <- centers_scales$center_sites[match(terms_in_df, centers_scales$variable)]
# 		scales <- centers_scales$scale_sites[match(terms_in_df, centers_scales$variable)]
# 		names(centers) <- names(scales) <- terms_in_df

# 		x_biomass <- as.data.table(x_raw)
# 		x_biomass <- x_biomass[ , ..terms_in_df]
# 		x_biomass <- scale(x_biomass, center = centers, scale = scales)
# 		x_biomass <- as.data.frame(x_biomass)
# 		x_biomass <- model.matrix(formulae$formula_biomass, x_biomass)

# 		# prepare predictors for non-biomass facets
# 		x_nonbiomass <- list()
# 		for (f in seq_along(this_nonbiomass_facets)) {

# 			facet <- names(this_nonbiomass_facets)[f]
# 			formula_facet <- this_nonbiomass_facets[[facet]]$formula

# 			terms <- attr(terms(formula_facet), 'term.labels')
# 			terms_in_df <- terms[terms %in% names(x_raw)]
# 			centers <- centers_scales$center_sites[match(terms_in_df, centers_scales$variable)]
# 			scales <- centers_scales$scale_sites[match(terms_in_df, centers_scales$variable)]
# 			names(centers) <- names(scales) <- terms_in_df

# 			x_facet <- as.data.table(x_raw)
# 			x_facet <- x_facet[ , ..terms_in_df]
# 			x_facet <- scale(x_facet, center = centers, scale = scales)
# 			x_facet <- as.data.frame(x_facet)
# 			x_facet <- model.matrix(formula_facet, x_facet)

# 			x_nonbiomass[[facet]] <- x_facet

# 		}

# 		preds <- predict_fully_integrated(
# 			chains = chains,
# 			nonbiomass_facets = this_nonbiomass_facets,
# 			resp_distrib_biomass = formulae$meta_biomass$resp_distrib_biomass,
# 			transform_biomass = formulae$meta_biomass$transform_biomass,
# 			x_occs = x_occs,
# 			x_biomass = x_biomass,
# 			x_psi = x_occs,
# 			x_nonbiomass = x_nonbiomass,
# 			w_occs = NULL,
# 			force_presence = TRUE
# 		)

# 		preds_height <- preds$preds_nonbiomass$height
# 		predictions <- colMeans(preds_height, na.rm = TRUE)

# 		# compare facet model predictions to observed height classes
# 		obs_pred <- data.table(
# 			site = factor(test_data$site),
# 			observed = test_data$Height,
# 			predicted = predictions,
# 			w = 1 / test_data$site_n
# 		)
# 		obs_pred <- obs_pred[is.finite(predicted)]

# 		# calibration model
# 		calib <- lm(observed ~ predicted, data = obs_pred, weights = w)

# 		r_squared <- summary(calib)$r.squared
# 		coefs <- coefficients(calib)
# 		intercept <- coefs['(Intercept)']
# 		slope <- coefs['predicted']

# 		# cluster-robust SEs			
# 		vcv <- vcovCL(calib, cluster = ~ site)
# 		robust_se <- coeftest(calib, vcov = vcv, type = 'HC2')
# 		se_slope <- robust_se['predicted', 'Std. Error']
# 		se_intercept <- robust_se['(Intercept)', 'Std. Error']

# 		intercept_lower <- intercept - 1.96 * se_intercept
# 		intercept_upper <- intercept + 1.96 * se_intercept
# 		slope_lower <- slope - 1.96 * se_slope
# 		slope_upper <- slope + 1.96 * se_slope

# 		slope_sig <- ifelse(slope_lower > 1 | slope_upper < 1, '*', '-')
# 		intercept_sig <- ifelse(intercept_lower > 1 | intercept_upper < 1, '*', '-')

# 		formula_facet <- paste('~', as.character(formulae$nonbiomass_facets[['height']]$formula))[2]
# 		resp_distrib <- nonbiomass_facets[['height']]$resp_distrib
# 		transform <- nonbiomass_facets[['height']]$transform

# 		summary <- data.table(
# 			facet = 'height',
# 			model_type = 'integrated',
# 			model_folder = basename(model_dir),
# 			formula_facet = formula_facet, 
# 			resp_distrib = resp_distrib,
# 			transform = transform,
# 			intercept = intercept,
# 			intercept_lower = intercept_lower,
# 			intercept_upper = intercept_upper,
# 			intercept_sig = intercept_sig,
# 			slope = slope,
# 			slope_lower = slope_lower,
# 			slope_upper = slope_upper,
# 			slope_sig = slope_sig,
# 			r_squared = r_squared
# 		)
		
# 	fwrite(summary, paste0(model_dir, '/validation_vs_sandel_and_griffan_nolan_height_integrated_model.csv'))


say('DONE', level = 1, deco = '@')

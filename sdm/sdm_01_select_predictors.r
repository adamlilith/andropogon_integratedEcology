### MODELING ANDROPOGON GERARDI DISTRIBUTION, MORPHOLOGY, PHYSIOLOGY, and ASSOCIATED MICROBIAL COMMUNITIES
### Erica Newman | Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2023-11
###
### This script analyzes candidate predictors for collinearity.
###
### source('C:/Ecology/R/andropogon_integratedEcology/sdm/sdm_01_select_predictors.r')
### source('C:/Subarashi/R/andropogon_integratedEcology/sdm/sdm_01_select_predictors.r')
###
### CONTENTS ###
### setup ###
### cluster analysis on candidate predictors based on correlation matrix ###

#############
### setup ###
#############

	rm(list = ls())

	drive <- 'C:/Ecology/'
	# drive <- 'C:/Subarashi/'

	setwd(paste0(drive, '/Research/Andropogon/Andropogon'))

	library(cluster)
	library(omnibus)
	library(terra)

# say('############################################################################')
# say('### cluster analysis on candidate predictors based on correlation matrix ###')
# say('############################################################################')

	ag_vect <- vect('./data_from_adam_and_loretta/andropogon_gerardi_occurrences_with_environment_1961_2020.gpkg')

	preds <- c(paste0('bio', c(1:7, 10:12, 15, 18:19)), 'aridity', 'ph', 'soc', 'sand', 'silt', 'clay')
	vars <- as.data.frame(ag_vect)
	vars <- vars[ , preds]
	vars <- vars[complete.cases(vars), ]
	cor <- cor(vars, method = 'spearman')
	dists <- 1 - abs(cor)
	dists <- as.dist(dists)

	clust <- agnes(dists, method = 'average')

	png('./outputs_loretta/Correlations between SDM Predictors - Cluster Diagram.png', width = 1200, height = 800)
	par(cex = 2)
	plot(clust, which.plot = 2, ylab = '1 - abs(correlation)', xlab = '', main = 'All counties in North America')
	abline(h = 0.3, col = 'red')
	dev.off()

say('#########################################################')
say('### preliminary analysis of candidate predictor terms ###')
say('#########################################################')

	ag_vect <- vect('./data_from_adam_and_loretta/andropogon_gerardi_occurrences_with_environment_1961_2020.gpkg')

	# preds <- c(paste0('bio', c(1:7, 10:12, 15, 18:19)), 'aridity', 'ph', 'soc', 'sand', 'silt', 'clay')
	preds <- c(paste0('bio', c(12, 1, 18, 5, 6, 7, 15, 2, 10)), 'aridity', 'ph', 'soc', 'sand', 'silt', 'clay')

	numer <- ag_vect$any_ag_quality_1_to_3
	denom <- ag_vect$num_poaceae_records

	resp <- numer / denom

	models <- list()
	results <- data.frame()

	for (pred in preds) {

		pred_data <- scale(ag_vect[[pred]])
		pred_data_sq <- pred_data^2
		model <- glm(resp ~ pred_data + pred_data_sq, family = binomial)
		models[[pred]] <- summary(model)

		results <- rbind(
			results,
			data.frame(
				pred = pred,
				converged = model$converged,
				model$deviance,
				model$null.deviance,
				aic = AIC(model),
				coeff_linear = model$coefficients[2],
				coeff_quad = model$coefficients[3]
			)
		)
	
	}

	rownames(results) <- NULL
	results <- results[order(results$aic), ]

	write.csv(results, './outputs_loretta/simple_quadratic_models_on_ratio_of_ag_to_poaceae.csv', row.names = FALSE)



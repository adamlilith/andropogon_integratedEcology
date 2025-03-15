### MODELING ANDROPOGON GERARDI DISTRIBUTION, PHENOTYPE, PHYSIOLOGY, GENOTYPE, and ASSOCIATED MICROBIAL COMMUNITIES
### Erica Newman | Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2024-11
###
### This script constructs a network diagram of microbial taxa that co-occur.
###
### source('C:/Ecology/R/andropogon_integratedEcology/sdm_admixture/admixture_01_sdm_species_level_priors.r')
### source('C:/Subarashi/R/andropogon_integratedEcology/sdm_admixture/admixture_01_sdm_species_level_priors.r')
###
### CONTENTS ###
### setup ###
### integrated SDM-ADMIXTURE model ###

say('#######################', pre = 1)
say('### response curves ###')
say('#######################')

	# chains <- readRDS(paste0(out_dir, '/chains.rds'))
	mcmc <- chains$samples

	# 1st is species
	colors <- c('black', 'red', 'green', 'cyan', 'orange', 'yellow', 'violet', 'floralwhite')

	# make a plot for how AG SDM and each ancestral population responds to each predictor
	responses <- list()
	for (pred in 1:n_predictors) {

		# unscaled predictor value
		x <- env_array[ , pred, pred]

		# create data frame with SDM response
		pars <- paste0('response_curves_sdm[', 1:n_response_curve_rows, ', ', pred, ']')
		sdm_response_mean <- chains$summary$all.chains[pars, 'Mean']
		sdm_response_lower <- chains$summary$all.chains[pars, '95%CI_low']
		sdm_response_upper <- chains$summary$all.chains[pars, '95%CI_upp']

		# data frames to hold predictions in long format
		df_mean <- data.frame(
			x = x,
			response = sdm_response_mean / max(sdm_response_mean),
			type = 'Species-level'
		)

		df_lower <- data.frame(
			x = x,
			response = sdm_response_lower,
			type = 'Species-level'
		)

		df_upper <- data.frame(
			x = x,
			response = sdm_response_upper,
			type = 'Species-level'
		)

		# add to data frames ADMIXTURE responses
		admix_response_mean <- admix_response_lower <- admix_response_upper <- list()
		for (k in 1:K) {

			# response of this ancestral population to this predictor
			pars <- paste0('response_curves_admix[', 1:n_response_curve_rows, ', ', k, ', ', pred, ']')
			admix_response_mean <- chains$summary$all.chains[pars, 'Mean']
			admix_response_lower <- chains$summary$all.chains[pars, '95%CI_low']
			admix_response_upper <- chains$summary$all.chains[pars, '95%CI_upp']

			# data frames to hold predictions in long format
			df_mean <- rbind(
				df_mean,
				data.frame(
					x = x,
					response = admix_response_mean,
					type = paste0('K', k)
				)
			)

			df_lower <- rbind(
				df_lower,
				data.frame(
					x = x,
					response = admix_response_lower,
					type = paste0('K', k)
				)
			)

			df_upper <- rbind(
				df_upper,
				data.frame(
					x = x,
					response = admix_response_upper,
					type = paste0('K', k)
				)
			)

		}


		responses[[pred]] <- ggplot()
		types <- unique(df_mean$type)

		# plot CI first
		for (i in seq_along(types)) {

			type <- types[i]

			this_df_upper <- df_upper[df_mean$type == type, ]
			this_df_lower <- df_lower[df_mean$type == type, ]

			max_score <- max(this_df_upper$response)

			this_df_upper$response <- this_df_upper$response / max_score
			this_df_lower$response <- this_df_lower$response / max_score

			this_df_ci <- data.frame(
				x = c(x, rev(x)),
				y = c(this_df_upper$response, rev(this_df_lower$response)),
				type = type
			)

			responses[[pred]] <- responses[[pred]] +
				geom_polygon(
					data = this_df_ci,
					mapping = aes(x = x, y = y),
					color = NA,
					fill = alpha(colors[i], 0.1)
				)

		}

		# plot trendlines second
		for (i in seq_along(types)) {

			type <- types[i]

			this_df_mean <- df_mean[df_mean$type == type, ]
			this_df_upper <- df_upper[df_mean$type == type, ]
			this_df_lower <- df_lower[df_mean$type == type, ]

			max_score <- max(this_df_mean$response, this_df_upper$response)

			this_df_mean$response <- this_df_mean$response / max_score

			responses[[pred]] <- responses[[pred]] +
				geom_line(
					data = this_df_mean,
					mapping = aes(x = x, y = response),
					color = colors[i]
				)

		}

		if (predictor_names[pred] == 'aridity') {
			nice_title <- 'Aridity (temperature / precipitation)'
			nice_axis <- 'Aridity (°C / mm)'
		} else if (predictor_names[pred] == 'bio7') {
			nice_title <- 'Temperature annual range (BIO 07)'
			nice_axis <- 'Temperature annual range (°C)'
		} else if (predictor_names[pred] == 'bio12') {
			nice_title <- 'Total annual precipitation (BIO 12)'
			nice_axis <- 'Total annual precipitation (mm)'
		} else if (predictor_names[pred] == 'bio15') {
			nice_title <- 'Precipitation seasonality (BIO 15)'
			nice_axis <- 'Precipitation seasonality'
		} else if (predictor_names[pred] == 'bio15') {
			nice_title <- 'Precipitation seasonality (BIO 15)'
			nice_axis <- 'Precipitation seasonality'
		} else if (predictor_names[pred] == 'ph') {
			nice_title <- 'Soil pH'
			nice_axis <- 'pH'
		} else if (predictor_names[pred] == 'sand') {
			nice_title <- 'Soil proportion sand'
			nice_axis <- 'Proportion sand'
		} else if (predictor_names[pred] == 'silt') {
			nice_title <- 'Soil proportion silt'
			nice_axis <- 'Proportion silt'
		} else if (predictor_names[pred] == 'soc') {
			nice_title <- 'Soil organic matter'
			nice_axis <- 'Soil organic matter'
		}

		responses[[pred]] <- responses[[pred]] +
			ylim(0, 1) +
			scale_color_discrete(name = 'Entity') +
			scale_fill_discrete(name = 'Entity') +
			xlab(nice_axis) +
			ylab('Scaled ancestry proportion or\nSDM-predicted abundance') +
			ggtitle(nice_title) +
			theme(
				plot.title = element_text(size = 10)
			)

	} # next predictor

	if (length(predictor_names) <= 4) {
		nrow <- 1
		width <- 12
		height <- 4
	} else {
		nrow <- 2
		width <- 12
		height <- 10	
	}

	responses <- plot_grid(plotlist = responses, nrow = nrow)

	ggsave(plot = responses, filename = paste0(out_dir, '/response_curves.png'), width = width, height = height, dpi = 600)

### MODELING ANDROPOGON GERARDI DISTRIBUTION, MORPHOLOGY, PHYSIOLOGY, and ASSOCIATED MICROBIAL COMMUNITIES
### Erica Newman | Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2023-11
###
### This script has shared functions and variables for modeling microbes associated with Andropogon gerardi.
###
### source('C:/Ecology/Research/Andropogon/Andropogon/andropogon_integratedEcology/hmsc_vs_microbes/microbes_00_shared_variables_and_functions.r')
### source('E:/Adam/Research/Andropogon/Andropogon/andropogon_integratedEcology/hmsc_vs_microbes/microbes_00_shared_variables_and_functions.r')
###
### CONTENTS ###

	# number of years across which "long-term" variables are calculated
	MASTER_long_term_years <- 10

	# start/end month-day of seasons
	MASTER_summer_month_day_start <- '06-01'
	MASTER_summer_month_day_end <- '08-31'

	MASTER_winter_month_day_start <- '12-01'
	MASTER_winter_month_day_end <- '02-28'
	
	MASTER_climate_years <- (2022 - MASTER_long_term_years + 1):2022 # MASTER_climate_years across which to calculate long-term variables

	sites <- read_excel('./data_from_sonny_and_loretta/AGER_site_date_sampling_03JUL2024.xlsx', sheet = 'Sheet1')
	sampling_dates <- as.Date(sites$SAMPLING_DATE)
	MASTER_mean_sampling_date <- mean(sampling_dates)

	MASTER_first_sampling_date <- min(sites$SAMPLING_DATE)
	MASTER_mean_sampling_date <- mean(sites$SAMPLING_DATE)


####################################################
### function to compile selected climate rasters ###
####################################################

	# pr_dir 			folder with PRISM data, should end with "PRISM/working/an81"
	# sampling_date		
	compile_climate_rasters <- function(pr_dir, sampling_date) {

		clim_rasts <- list()

		env_combined <- rbind(env_rhizo, env_bulk)
		env_combined_spatial_nad83 <- vect(env_combined, geom = c('longitude', 'latitude'), crs = getCRS('NAD83'))

		say('   summer mean temperature')

			x <- list()
			for (year in MASTER_climate_years) {
			
				start_date <- paste0(year, '-', MASTER_summer_month_day_start)
				end_date <- paste0(year, '-', MASTER_summer_month_day_end)
			
				dates <- c(start_date, end_date)
				
				x[[length(x) + 1]] <- prStack(
					pr_dir,
					vars = 'tmean',
					dates = dates,
					by = 'month',
					span = TRUE,
					res = 800,
					rastSuffix = 'tif'
				)
				
			}
			
			x <- do.call(c, x)
			x <- mean(x)

			names(x) <- 'summer_tmean_c'
			clim_rasts <- c(clim_rasts, x)
		
		say('total annual precipitation')

			first_sampling_year <- year(MASTER_first_sampling_date)
			first_sampling_month <- month(MASTER_first_sampling_date)

			mean_sampling_month <- month(MASTER_mean_sampling_date)

			start_date <- paste0(first_sampling_year - MASTER_long_term_years, '-', prefix(first_sampling_month, 2), '-01')
			end_date <- paste0(max(MASTER_climate_years), '-', prefix(mean_sampling_month, 2), '-28')
			
			dates <- c(start_date, end_date)
			
			x <- prStack(
				pr_dir,
				vars = 'ppt',
				dates = dates,
				by = 'month',
				span = TRUE,
				res = 800,
				rastSuffix = 'tif'
			)
			
			months <- nlyr(x)
			x <- sum(x)
			x <- 12 * x / months
			
			names(x) <- 'annual_ppt_mm'
			clim_rasts <- c(clim_rasts, x)

		say('daily precipitation CV', level = 2)

			first_sampling_year <- year(MASTER_first_sampling_date)
			first_sampling_month <- month(MASTER_first_sampling_date)

			mean_sampling_month <- month(MASTER_mean_sampling_date)

			start_date <- paste0(first_sampling_year - MASTER_long_term_years, '-', prefix(first_sampling_month, 2), '-01')
			end_date <- paste0(max(years), '-', prefix(mean_sampling_month, 2), '-28')
			
			dates <- c(start_date, end_date)
			
			x <- prStack(
				pr_dir,
				vars = 'ppt',
				dates = dates,
				by = 'day',
				span = TRUE,
				res = 800,
				rastSuffix = 'tif'
			)
			
			avg <- mean(x)
			sd <- stdev(x, pop = FALSE)
			x <- sd / avg

			names(x) <- 'ppt_cv'
			clim_rasts <- c(clim_rasts, x)

		say('precipitation before sampling', level = 2)

		start_date <- MASTER_mean_sampling_date - 6

		x <- prStack(
			pr_dir,
			vars = 'ppt',
			dates = c(start_date, MASTER_mean_sampling_date),
			by = 'day',
			span = TRUE,
			res = 800,
			rastSuffix = 'tif'
		)

		x <- sum(x)

		names(x) <- 'sampling_ppt_mm'
		clim_rasts <- c(clim_rasts, x)

	say('temperature before sampling', level = 2)

		start_date <- MASTER_mean_sampling_date - 29

		x <- prStack(
			pr_dir,
			vars = 'tmean',
			dates = c(start_date, MASTER_mean_sampling_date),
			by = 'day',
			span = TRUE,
			res = 800,
			rastSuffix = 'tif'
		)

		x <- mean(x)

		names(x) <- 'sampling_tmean_C'
		clim_rasts <- c(clim_rasts, x)
		
		clim_rasts

	}


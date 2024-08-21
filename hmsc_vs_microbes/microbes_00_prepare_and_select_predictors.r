### MODELING ANDROPOGON GERARDI DISTRIBUTION, MORPHOLOGY, PHYSIOLOGY, and ASSOCIATED MICROBIAL COMMUNITIES
### Erica Newman | Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2023-11
###
### This script compiles environmental data for modeling of microbial communities associated with Andropogon gerardi.
###
### source('C:/Ecology/Research/Andropogon/Andropogon/andropogon_integratedEcology/hmsc_vs_microbes/microbes_00_prepare_and_select_predictors.r')
### source('E:/Adam/Research/Andropogon/Andropogon/andropogon_integratedEcology/hmsc_vs_microbes/microbes_00_prepare_and_select_predictors.r')
###
### CONTENTS ###
### setup ###
### settings ###
### prepare candidate predictors ###
### correlation between predictors ###
### extract values from raster cells for spatial predictions ###

#############
### setup ###
#############

rm(list = ls())

drive <- 'C:/Ecology/'
# drive <- 'E:/Adam/'

setwd(paste0(drive, '/Research/Andropogon/Andropogon'))

# library(airUpThere)
devtools::load_all(paste0(drive, '/R/airUpThere'))
library(cluster)
library(data.table)
# library(enmSdmX)
library(omnibus)
library(readxl)
library(terra)

################
### settings ###
################

	# number of years across which "long-term" variables are calculated
	long_term_years <- 10

	# start/end month-day of seasons
	summer_month_day_start <- '06-01'
	summer_month_day_end <- '08-31'

	winter_month_day_start <- '12-01'
	winter_month_day_end <- '02-28'
	
	# pr_dir <- 'D:/ecology/PRISM/working/an81'
	pr_dir <- 'E:/ecology/PRISM/working/an81'
	# pr_dir <- 'F:/ecology/PRISM/working/an81'
	
	years <- (2022 - long_term_years + 1):2022 # years across which to calculate long-term variables

# say('####################################')
# say('### prepare candidate predictors ###')
# say('####################################')

	# sink('./data/data_from_sonny_and_loretta/AGER_site_date_sampling_03JUL2024_extracted_environment_readme.txt', split = TRUE)

	# say('CANDIDATE PREDICTORS ATTACHED TO THE SITE SAMPLES FROM WHICH WE HAVE MICROBIAL DATA')
	# say('Data in file "AGER_site_date_sampling_03JUL2024_extracted_environment.csv".')
	# say(date(), post = 1)
	# say('"Long-term" climate variables:', pre = 1)
	# say('Candidate "long-term" predictors are calculated from monthly or daily rasters and include:')
	# say('* summer_tmean_C: Summer mean temperature (JJA) from mean monthly temperature averaged across years from PRISM AN81m across the 10 yr prior to and including sampling (C)')
	# say('* winter_tmean_C: Winter mean temperature (DJF) from mean monthly temperature averaged across years from PRISM AN81m across the 10 yr prior to and including sampling (C)')
	# say('* tmean_sd_C: Temperature standard deviation from mean monthly temperature across all 12 months in all years from PRISM AN81m across the 10 yr prior to and including sampling (C)')
	# say('* summer_ppt_mm: Summer precipitation (JJA) from total monthly precipitation, averaged across years from PRISM AN81m across the 10 yr prior to and including sampling (mm')
	# say('* annual_ppt_mm: Mean annual precipitation from total monthly precipitation across all years from PRISM AN81m across the 10 yr prior to and including sampling (C)')
	# say('* ppt_cv: Precipitation CV from daily precipitation across all days across all years from PRISM AN81d across the 10 yr prior to and including sampling (unitless')
	# say('* conseq_days_sans_precip: Mean length of consecutive runs of days with no precipitation from daily precipitation across all days in all years from PRISM AN81d across the 10 yr prior to and including sampling (days)')
	# say('  e.g., 1 0 0 1 1 0 1 0 1 0 0 0 would have 2 runs (of 2 and 3 d each)')
	
	# say('"Short-term" weather variables:', pre = 1)
	# say('* sampling_ppt_mm: Total precipitation in the 7 days prior to and including sampling from PRISM AN81d (mm).')
	# say('* sampling_tmean_C: Mean of mean daily temperature in the 30 days prior to and including sampling from PRISM AN81d (C).')

	# say('Soil variables:', pre = 1)
	# say('* ph_soilgrids: pH from SoilGrids 2.0 (10*)')
	# say('* nitrogen_soilgrids: N (nitrogen) from SoilGrids 2.0 (g / kg)')
	# say('* soc_soilgrids: SOC (soil organic carbon) from SoilGrids 2.0 (g / kg)')
	# say('* cec_soilgrids: CEC (cation exchange capacity) from SoilGrids 2.0 (cmol(c) / kg)')
	# say('* sand_soilgrids, silt_soilgrids, clay_soilgrids: Sand/silt/clay percentages from SoilGrids 2.0 (proportion)')
	# say('* aluminum_usgs: Aluminum from USGS Smith et al. 2019 (weight %)')
	# say('* calcium_usgs: Calcium from USGS Smith et al. 2019 (weight %)')
	# say('* sodium_usgs: Sodium from USGS Smith et al. 2019 (weight %)')
	# say('* phosphorus_usgs: Phosphorus from USGS Smith et al. 2019 (mg / kg)')

	# sink()

	# sites <- read_excel('./data/data_from_sonny_and_loretta/AGER_site_date_sampling_03JUL2024.xlsx', sheet = 'Sheet1')
	# sites <- vect(sites, geom = c('LONGITUDE', 'LATITUDE'), crs = enmSdmX::getCRS('WGS84'), keepgeom = TRUE)

	# ### dates for calculating variables
	# # using mean sampling date
	# first_sampling_date <- min(sites$SAMPLING_DATE)
	# first_sampling_year <- year(first_sampling_date)
	# first_sampling_month <- month(first_sampling_date)

	# mean_sampling_date <- mean(sites$SAMPLING_DATE)
	# mean_sampling_month <- month(mean_sampling_date)

	# ### summer mean temperature
	# ###########################
	# say('summer mean temperature', level = 2)

		# stack <- list()
		# for (year in years) {
		
			# start_date <- paste0(year, '-', summer_month_day_start)
			# end_date <- paste0(year, '-', summer_month_day_end)
		
			# dates <- c(start_date, end_date)
			
			# stack[[length(stack) + 1]] <- prStack(
				# pr_dir,
				# vars = 'tmean',
				# dates = dates,
				# by = 'month',
				# span = TRUE,
				# res = 800,
				# rastSuffix = 'tif'
			# )
			
		# }
		
		# stack <- do.call(c, stack)
		# x <- extract(stack, sites, ID = FALSE)
		# x <- rowMeans(x)
		
		# sites$summer_tmean_C <- x
		
	# ### winter mean temperature
	# ###########################
	# say('winter mean temperature', level = 2)

		# stack <- list()
		# for (year in years) {
		
			# start_date <- paste0(year - 1, '-', winter_month_day_start)
			# end_date <- paste0(year, '-', winter_month_day_end)
		
			# dates <- c(start_date, end_date)
			
			# stack[[length(stack) + 1]] <- prStack(
				# pr_dir,
				# vars = 'tmean',
				# dates = dates,
				# by = 'month',
				# span = TRUE,
				# res = 800,
				# rastSuffix = 'tif'
			# )
			
		# }
		
		# stack <- do.call(c, stack)
		# x <- extract(stack, sites, ID = FALSE)
		# x <- rowMeans(x)
		
		# sites$winter_tmean_C <- x
		
	# ### temperature SD across entire year
	# #####################################
	# say('temperature SD across entire year', level = 2)

		# start_date <- paste0(first_sampling_year - long_term_years, '-', prefix(first_sampling_month, 2), '-01')
		# end_date <- paste0(max(years), '-', prefix(mean_sampling_month, 2), '-28')
		
		# dates <- c(start_date, end_date)
		
		# stack <- prStack(
			# pr_dir,
			# vars = 'tmean',
			# dates = dates,
			# by = 'month',
			# span = TRUE,
			# res = 800,
			# rastSuffix = 'tif'
		# )
		
		# x <- extract(stack, sites, ID = FALSE)
		# x <- apply(x, 1, sd)
		
		# sites$tmean_sd_C <- x
		
	# ### summer total precipitation
	# ##############################
	# say('summer total precipitation', level = 2)

		# stack <- list()
		# for (year in years) {
		
			# start_date <- paste0(year, '-', summer_month_day_start)
			# end_date <- paste0(year, '-', summer_month_day_end)
		
			# dates <- c(start_date, end_date)
			
			# stack[[length(stack) + 1]] <- prStack(
				# pr_dir,
				# vars = 'ppt',
				# dates = dates,
				# by = 'month',
				# span = TRUE,
				# res = 800,
				# rastSuffix = 'tif'
			# )
			
			# stack[[length(stack)]] <- sum(stack[[length(stack)]])
			
		# }
		
		# stack <- do.call(c, stack)
		# x <- extract(stack, sites, ID = FALSE)
		# x <- rowSums(x)
		
		# sites$summer_ppt_mm <- x
		
	# ### total annual precipitation
	# ##############################
	# say('total annual precipitation', level = 2)

		# first_sampling_date <- min(sites$SAMPLING_DATE)
		# first_sampling_year <- year(first_sampling_date)
		# first_sampling_month <- month(first_sampling_date)

		# mean_sampling_date <- mean(sites$SAMPLING_DATE)
		# mean_sampling_month <- month(mean_sampling_date)

		# start_date <- paste0(first_sampling_year - long_term_years, '-', prefix(first_sampling_month, 2), '-01')
		# end_date <- paste0(max(years), '-', prefix(mean_sampling_month, 2), '-28')
		
		# dates <- c(start_date, end_date)
		
		# stack <- prStack(
			# pr_dir,
			# vars = 'ppt',
			# dates = dates,
			# by = 'month',
			# span = TRUE,
			# res = 800,
			# rastSuffix = 'tif'
		# )
		
		# months <- nlyr(stack)
		# x <- extract(stack, sites, ID = FALSE)
		# x <- rowSums(x)
		# x <- 12 * x / months
		
		# sites$annual_ppt_mm <- x
		
	# ### daily precipitation CV
	# ##########################
	# say('daily precipitation CV', level = 2)

		# first_sampling_date <- min(sites$SAMPLING_DATE)
		# first_sampling_year <- year(first_sampling_date)
		# first_sampling_month <- month(first_sampling_date)

		# mean_sampling_date <- mean(sites$SAMPLING_DATE)
		# mean_sampling_month <- month(mean_sampling_date)

		# start_date <- paste0(first_sampling_year - long_term_years, '-', prefix(first_sampling_month, 2), '-01')
		# end_date <- paste0(max(years), '-', prefix(mean_sampling_month, 2), '-28')
		
		# dates <- c(start_date, end_date)
		
		# stack <- prStack(
			# pr_dir,
			# vars = 'ppt',
			# dates = dates,
			# by = 'day',
			# span = TRUE,
			# res = 800,
			# rastSuffix = 'tif'
		# )
		
		# extr <- extract(stack, sites, ID = FALSE)
		# avg <- rowMeans(extr)
		# sd <- apply(extr, 1, sd)
		# x <- sd / avg
		
		# sites$ppt_cv <- x
		
	# ### median consecutive runs with no precipitation
	# #################################################
	# say('median consecutive runs with no precipitation', level = 2)

		# first_sampling_date <- min(sites$SAMPLING_DATE)
		# first_sampling_year <- year(first_sampling_date)
		# first_sampling_month <- month(first_sampling_date)

		# mean_sampling_date <- mean(sites$SAMPLING_DATE)
		# mean_sampling_month <- month(mean_sampling_date)

		# start_date <- paste0(first_sampling_year - long_term_years, '-', prefix(first_sampling_month, 2), '-01')
		# end_date <- paste0(max(years), '-', prefix(mean_sampling_month, 2), '-28')
		
		# dates <- c(start_date, end_date)
		
		# stack <- prStack(
			# pr_dir,
			# vars = 'ppt',
			# dates = dates,
			# by = 'day',
			# span = TRUE,
			# res = 800,
			# rastSuffix = 'tif'
		# )

		# # consecutive days
		# x <- extract(stack, sites, ID = FALSE)
		# for (i in 1:nrow(x)) {
		
			# for (j in 1:(ncol(x) - 1)) {
			
				# if (x[i, j] == 0 & x[i, j + 1] == 0) {
					# x[i, j] <- 1
				# } else {
					# x[i, j] <- 0
				# }
			
			# }
		
		# }
		# x[ , ncol(x)] <- 0
		
		# x <- apply(x, 1, mean)
		# sites$conseq_days_sans_precip <- x
		
	# ### precipitation before sampling
	# #################################
	# say('precipitation before sampling', level = 2)

	# x <- prExtractRelativeDaily(
		# pr_dir,
		# x = sites,
		# vars = 'ppt',
		# date = 'SAMPLING_DATE',
		# res = 800,
		# rastSuffix = 'tif',
		# windowYears = 0,
		# windowDays = 6,
		# verbose = TRUE
	# )
	
	# x <- rowSums(x)
	
	# sites$sampling_ppt_mm <- x
	
	# ### temperature before sampling
	# ###############################
	# say('temperature before sampling', level = 2)
	
	# x <- prExtractRelativeDaily(
		# pr_dir,
		# x = sites,
		# vars = 'tmean',
		# date = 'SAMPLING_DATE',
		# res = 800,
		# rastSuffix = 'tif',
		# windowYears = 0,
		# windowDays = 29,
		# verbose = TRUE
	# )
	
	# x <- rowMeans(x)
	
	# sites$sampling_tmean_C <- x
	
	# ### soil pH
	# ###########
	# say('pH', level = 2)
	
	# x <- rast(paste0(drive, '/Research Data/SoilGrids/SoilGrids 2.0/phh2o_0-5cm_mean_northAmerica.tif'))
	# x <- extract(x, sites, ID = FALSE)
	# sites$ph_soilgrids <- x[ , 1] / 10

	# ### soil N
	# ###########
	# say('N', level = 2)
	
	# x <- rast(paste0(drive, '/Research Data/SoilGrids/SoilGrids 2.0/nitrogen_0-5cm_mean_northAmerica.tif'))
	# x <- extract(x, sites, ID = FALSE)
	# sites$nitrogen_soilgrids <- x[ , 1]

	# ### SOC
	# #######
	# say('SOC', level = 2)
	
	# x <- rast(paste0(drive, '/Research Data/SoilGrids/SoilGrids 2.0/soc_0-5cm_mean.tif'))
	# x <- extract(x, sites, ID = FALSE)
	# sites$soc_soilgrids <- x[ , 1]

	# ### CEC
	# #######
	# say('CEC', level = 2)
	
	# x <- rast(paste0(drive, '/Research Data/SoilGrids/SoilGrids 2.0/cec_0-5cm_mean_northAmerica.tif'))
	# x <- extract(x, sites, ID = FALSE)
	# sites$cec_soilgrids <- x[ , 1]

	# ### sand
	# ########
	# say('sand', level = 2)
	
	# x <- rast(paste0(drive, '/Research Data/SoilGrids/SoilGrids 2.0/sand_0-5cm_mean_northAmerica.tif'))
	# x <- extract(x, sites, ID = FALSE)
	# sites$sand_soilgrids <- x[ , 1] / 1000

	# ### silt
	# ########
	# say('silt', level = 2)
	
	# x <- rast(paste0(drive, '/Research Data/SoilGrids/SoilGrids 2.0/silt_0-5cm_mean_northAmerica.tif'))
	# x <- extract(x, sites, ID = FALSE)
	# sites$silt_soilgrids <- x[ , 1] / 1000

	# ### clay
	# ########
	# say('clay', level = 2)
	
	# x <- rast(paste0(drive, '/Research Data/SoilGrids/SoilGrids 2.0/clay_0-5cm_mean_northAmerica.tif'))
	# x <- extract(x, sites, ID = FALSE)
	# sites$clay_soilgrids <- x[ , 1] / 1000

	# ### Al
	# ######
	# say('aluminum', level = 2)
	
	# x <- rast(paste0(drive, '/Research Data/Soil Elements of the USA Smith et al 2019 Geochemical and Mineralogical Maps/Elements/Al_tif/Top5_Al.tif'))
	# x <- extract(x, sites, ID = FALSE)
	# sites$aluminum_usgs <- x[ , 1]

	# ### Ca
	# ######
	# say('calcium', level = 2)
	
	# x <- rast(paste0(drive, '/Research Data/Soil Elements of the USA Smith et al 2019 Geochemical and Mineralogical Maps/Elements/Ca_tif/Top5_Ca.tif'))
	# x <- extract(x, sites, ID = FALSE)
	# sites$calcium_usgs <- x[ , 1]

	# ### Na
	# ######
	# say('sodium', level = 2)
	
	# x <- rast(paste0(drive, '/Research Data/Soil Elements of the USA Smith et al 2019 Geochemical and Mineralogical Maps/Elements/Na_tif/Top5_Na.tif'))
	# x <- extract(x, sites, ID = FALSE)
	# sites$sodium_usgs <- x[ , 1]

	# ### P
	# ######
	# say('phosphorus', level = 2)
	
	# x <- rast(paste0(drive, '/Research Data/Soil Elements of the USA Smith et al 2019 Geochemical and Mineralogical Maps/Elements/P_tif/Top5_P.tif'))
	# x <- extract(x, sites, ID = FALSE)
	# sites$phosphorus_usgs <- x[ , 1]

	# ### save
	# sites <- as.data.frame(sites)
	# write.csv(sites, './data/data_from_sonny_and_loretta/AGER_site_date_sampling_03JUL2024_extracted_environment.csv', row.names = FALSE)

# say('######################################')
# say('### correlation between predictors ###')
# say('######################################')

# 	sites <- read.csv('./data/data_from_sonny_and_loretta/AGER_site_date_sampling_03JUL2024_extracted_environment.csv')
	
# 	vars <- c('summer_tmean_C', 'winter_tmean_C', 'tmean_sd_C', 'summer_ppt_mm', 'annual_ppt_mm', 'ppt_cv', 'conseq_days_sans_precip', 'sampling_ppt_mm', 'sampling_tmean_C', 'ph_soilgrids', 'nitrogen_soilgrids', 'soc_soilgrids', 'cec_soilgrids', 'sand_soilgrids', 'silt_soilgrids', 'clay_soilgrids', 'aluminum_usgs', 'calcium_usgs', 'sodium_usgs', 'phosphorus_usgs')
	
# 	site_vars <- sites[ , vars]
# 	cors <- cor(site_vars, method = 'spearman')
	
# 	cor_dist <- 1 - abs(cors)
# 	dists <- as.dist(cor_dist)
	
# 	clust <- agnes(dists, method = 'average')

# 	png('./outputs_sonny/Correlations between Microbe Predictors - Cluster Diagram.png', width = 1200, height = 800)
# 	par(cex = 2)
# 	plot(clust, which.plot = 2, ylab = '1 - abs(correlation)', xlab = '')
# 	abline(h = 0.25, col = 'red')
# 	dev.off()

say('################################################################')
say('### extract values from raster cells for spatial predictions ###')
say('################################################################')

	# selected predictors
	# 'sampling_ppt_mm', 'summer_tmean_C', 'annual_ppt_mm', 'ppt_cv', 'pH', 'soil_c', 'soil_n', 'clay', 'sand'

	# holds predictor rasters
	clim_rasts <- list()
	soil_rasts <- list()

	sites <- read_excel('./data/data_from_sonny_and_loretta/AGER_site_date_sampling_03JUL2024.xlsx', sheet = 'Sheet1')
	sampling_dates <- as.Date(sites$SAMPLING_DATE)
	mean_sampling_date <- mean(sampling_dates)

	### summer mean temperature
	###########################
	say('summer mean temperature', level = 2)

		x <- list()
		for (year in years) {
		
			start_date <- paste0(year, '-', summer_month_day_start)
			end_date <- paste0(year, '-', summer_month_day_end)
		
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

		names(x) <- 'summer_tmean_C'
		clim_rasts <- c(clim_rasts, x)
		
	# ### winter mean temperature
	# ###########################
	# say('winter mean temperature', level = 2)

	# 	x <- list()
	# 	for (year in years) {
		
	# 		start_date <- paste0(year - 1, '-', winter_month_day_start)
	# 		end_date <- paste0(year, '-', winter_month_day_end)
		
	# 		dates <- c(start_date, end_date)
			
	# 		x[[length(x) + 1]] <- prStack(
	# 			pr_dir,
	# 			vars = 'tmean',
	# 			dates = dates,
	# 			by = 'month',
	# 			span = TRUE,
	# 			res = 800,
	# 			rastSuffix = 'tif'
	# 		)
			
	# 	}
		
	# 	x <- do.call(c, x)
	# 	x <- mean(x)

	# 	names(x) <- 'winter_tmean_C'
	# 	clim_rasts <- c(clim_rasts, x)
		
	# ### temperature SD across entire year
	# #####################################
	# say('temperature SD across entire year', level = 2)

	# 	first_sampling_date <- min(sites$SAMPLING_DATE)
	# 	first_sampling_year <- year(first_sampling_date)
	# 	first_sampling_month <- month(first_sampling_date)

	# 	mean_sampling_date <- mean(sites$SAMPLING_DATE)
	# 	mean_sampling_month <- month(mean_sampling_date)

	# 	start_date <- paste0(first_sampling_year - long_term_years, '-', prefix(first_sampling_month, 2), '-01')
	# 	end_date <- paste0(max(years), '-', prefix(mean_sampling_month, 2), '-28')
		
	# 	dates <- c(start_date, end_date)
		
	# 	x <- prStack(
	# 		pr_dir,
	# 		vars = 'tmean',
	# 		dates = dates,
	# 		by = 'month',
	# 		span = TRUE,
	# 		res = 800,
	# 		rastSuffix = 'tif'
	# 	)
		
	# 	x <- do.call(c, x)
	# 	x <- stdev(x, pop = FALSE)

	# 	names(x) <- 'tmean_sd_C'
	# 	clim_rasts <- c(clim_rasts, x)

	# ### summer total precipitation
	# ##############################
	# say('summer total precipitation', level = 2)

	# 	x <- list()
	# 	for (year in years) {
		
	# 		start_date <- paste0(year, '-', summer_month_day_start)
	# 		end_date <- paste0(year, '-', summer_month_day_end)
		
	# 		dates <- c(start_date, end_date)
			
	# 		x[[length(x) + 1]] <- prStack(
	# 			pr_dir,
	# 			vars = 'ppt',
	# 			dates = dates,
	# 			by = 'month',
	# 			span = TRUE,
	# 			res = 800,
	# 			rastSuffix = 'tif'
	# 		)
			
	# 		x[[length(x)]] <- sum(x[[length(x)]])
			
	# 	}

	# 	x <- do.call(c, stack)
	# 	x <- sum(stack)
	# 	x <- 12 * x / months

	# 	names(x) <- 'summer_ppt_mm'
	# 	clim_rasts <- c(clim_rasts, x)
		
	### total annual precipitation
	##############################
	say('total annual precipitation', level = 2)

		first_sampling_date <- min(sites$SAMPLING_DATE)
		first_sampling_year <- year(first_sampling_date)
		first_sampling_month <- month(first_sampling_date)

		mean_sampling_date <- mean(sites$SAMPLING_DATE)
		mean_sampling_month <- month(mean_sampling_date)

		start_date <- paste0(first_sampling_year - long_term_years, '-', prefix(first_sampling_month, 2), '-01')
		end_date <- paste0(max(years), '-', prefix(mean_sampling_month, 2), '-28')
		
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

	### daily precipitation CV
	##########################
	say('daily precipitation CV', level = 2)

		first_sampling_date <- min(sites$SAMPLING_DATE)
		first_sampling_year <- year(first_sampling_date)
		first_sampling_month <- month(first_sampling_date)

		mean_sampling_date <- mean(sites$SAMPLING_DATE)
		mean_sampling_month <- month(mean_sampling_date)

		start_date <- paste0(first_sampling_year - long_term_years, '-', prefix(first_sampling_month, 2), '-01')
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
		
	# # ### median consecutive runs with no precipitation
	# # #################################################
	# # say('median consecutive runs with no precipitation', level = 2)

	# # 	first_sampling_date <- min(sites$SAMPLING_DATE)
	# # 	first_sampling_year <- year(first_sampling_date)
	# # 	first_sampling_month <- month(first_sampling_date)

	# # 	mean_sampling_date <- mean(sites$SAMPLING_DATE)
	# # 	mean_sampling_month <- month(mean_sampling_date)

	# # 	start_date <- paste0(first_sampling_year - long_term_years, '-', prefix(first_sampling_month, 2), '-01')
	# # 	end_date <- paste0(max(years), '-', prefix(mean_sampling_month, 2), '-28')
		
	# # 	dates <- c(start_date, end_date)
		
	# # 	stack <- prStack(
	# # 		pr_dir,
	# # 		vars = 'ppt',
	# # 		dates = dates,
	# # 		by = 'day',
	# # 		span = TRUE,
	# # 		res = 800,
	# # 		rastSuffix = 'tif'
	# # 	)

	# # 	# consecutive days
	# # 	x <- extract(stack, sites, ID = FALSE)
	# # 	for (i in 1:nrow(x)) {
		
	# # 		for (j in 1:(ncol(x) - 1)) {
			
	# # 			if (x[i, j] == 0 & x[i, j + 1] == 0) {
	# # 				x[i, j] <- 1
	# # 			} else {
	# # 				x[i, j] <- 0
	# # 			}
			
	# # 		}
		
	# # 	}
	# # 	x[ , ncol(x)] <- 0
		
	# # 	x <- apply(x, 1, mean)
	# # 	sites$conseq_days_sans_precip <- x
		
	### precipitation before sampling
	#################################
	say('precipitation before sampling', level = 2)

		start_date <- mean_sampling_date - 6

		x <- prStack(
			pr_dir,
			vars = 'ppt',
			dates = c(start_date, mean_sampling_date),
			by = 'day',
			span = TRUE,
			res = 800,
			rastSuffix = 'tif'
		)

		x <- sum(x)

		names(x) <- 'sampling_ppt_mm'
		clim_rasts <- c(clim_rasts, x)
		
	### temperature before sampling
	###############################
	say('temperature before sampling', level = 2)

		start_date <- mean_sampling_date - 29

		x <- prStack(
			pr_dir,
			vars = 'tmean',
			dates = c(start_date, mean_sampling_date),
			by = 'day',
			span = TRUE,
			res = 800,
			rastSuffix = 'tif'
		)

		x <- mean(x)

		names(x) <- 'sampling_tmean_C'
		clim_rasts <- c(clim_rasts, x)
	
	### soil pH
	###########
	say('pH', level = 2)
		
		x <- rast(paste0(drive, '/Research Data/SoilGrids/SoilGrids 2.0/phh2o_0-5cm_mean_northAmerica.tif'))
		x <- x / 10
		names(x) <- 'pH'
		soil_rasts <- c(soil_rasts, x)

	### soil N
	###########
	say('N', level = 2)
		
		x <- rast(paste0(drive, '/Research Data/SoilGrids/SoilGrids 2.0/nitrogen_0-5cm_mean_northAmerica.tif'))
		x <- crop(x, soil_rasts[[1]])
		names(x) <- 'soil_n'
		soil_rasts <- c(soil_rasts, x)


	### SOC
	#######
	say('SOC', level = 2)
	
		x <- rast(paste0(drive, '/Research Data/SoilGrids/SoilGrids 2.0/soc_0-5cm_mean.tif'))
		names(x) <- 'soc'
		x <- crop(x, soil_rasts[[1]])
		soil_rasts <- c(soil_rasts, x)

	# ### CEC
	# #######
	# # say('CEC', level = 2)
	
		# x <- rast(paste0(drive, '/Research Data/SoilGrids/SoilGrids 2.0/cec_0-5cm_mean_northAmerica.tif'))
		# 	names(x) <- 'soil_cec'
		# 	soil_rasts <- c(soil_rasts, x)

	## sand
	#######
	say('sand', level = 2)
	
			x <- rast(paste0(drive, '/Research Data/SoilGrids/SoilGrids 2.0/sand_0-5cm_mean_northAmerica.tif'))
			names(x) <- 'sand'
			soil_rasts <- c(soil_rasts, x)

	# ### silt
	# ########
	# say('silt', level = 2)
	
	# 	x <- rast(paste0(drive, '/Research Data/SoilGrids/SoilGrids 2.0/silt_0-5cm_mean_northAmerica.tif'))
	# 	names(x) <- 'silt'
	# 	soil_rasts <- c(soil_rasts, x)

	### clay
	########
	say('clay', level = 2)
		
		x <- rast(paste0(drive, '/Research Data/SoilGrids/SoilGrids 2.0/clay_0-5cm_mean_northAmerica.tif'))
		names(x) <- 'clay'
		soil_rasts <- c(soil_rasts, x)

	### amalgamate
	clim_rasts <- do.call(c, clim_rasts)
	soil_rasts <- do.call(c, soil_rasts)

	soil_rasts <- project(soil_rasts, clim_rasts, threads = 2)
	rasts <- c(clim_rasts, soil_rasts)

	rasts <- aggregate(rasts, 16, na.rm = TRUE) # aggregate to larger cells to save runtime

	conus_env <- as.data.table(rasts, xy = TRUE, cells = TRUE, na.rm = TRUE)
	saveRDS(conus_env, './outputs_sonny/conus_environment_as_per_prism_aggregated_by_16.rds')

say('DONE!', level = 1, deco = '!')

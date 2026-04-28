#' Get short title and axis names for predictors
#'
#' @param pred Name of predictor in R-friendly format.
#' @returns A title string, and axis title string, and a unit string
get_nice_predictor <- function(pred) {

	# get traits
	if (pred == 'aridity') {
		short <- 'Aridity'
		long <- 'Aridity ((temp. + 10) / ((1 + precip.) / 1000))'
		units <- ''
		range_fx <- standard_range_fx
	} else if (pred == 'bio1') {
		short <- 'Mean Annual Temperature (BIO01)'
		long <- 'Mean annual temperature (°C)'
		units <- '°C'
		range_fx <- standard_range_fx
	} else if (pred == 'bio5') {
		short <- 'Temperature of the Hottest Month (BIO05)'
		long <- 'Temperature of the hottest month (°C)'
		units <- '°C'
		range_fx <- standard_range_fx
	} else if (pred == 'bio6') {
		short <- 'Temperature of Coldest Month (BIO06)'
		long <- 'Temperature of the coldest month (°C)'
		units <- '°C'
		range_fx <- standard_range_fx
	} else if (pred == 'bio7') {
		short <- 'Temperature Annual Range (BIO07)'
		long <- 'Temperature annual range (°C)'
		units <- '°C'
		range_fx <- standard_range_fx
	} else if (pred %in% c('bio12', 'bio12_log10p1')) {
		short <- 'Total Annual Precipitation (BIO12)'
		long <- 'Total annual precipitation (mm)'
		units <- 'mm'
		range_fx <- bottom0_range_fx
	} else if (pred == 'bio15') {
		short <- 'Precipitation Seasonality (BIO15)'
		long <- 'Precipitation seasonality'
		units <- ''
		range_fx <- bottom0_range_fx
	} else if (pred %in% c('bio12', 'bio18_log10p1')) {
		short <- 'Precipitation of Warmest Quarter (BIO18)'
		long <- 'Precipitation of warmest quarter (mm)'
		units <- 'mm'
		range_fx <- bottom0_range_fx
	} else if (pred == 'ph') {
		short <- 'Soil pH'
		long <- 'pH'
		units <- ''
		range_fx <- standard_range_fx
	} else if (pred == 'sand') {
		short <- 'Soil Percent Sand'
		long <- 'Percent sand'
		units <- '%'
		range_fx <- percent_range_fx
	} else if (pred == 'silt') {
		short <- 'Soil Percent Silt'
		long <- 'Percent silt'
		units <- '%'
		range_fx <- percent_range_fx
	} else if (pred == 'clay') {
		short <- 'Soil Percent Clay'
		long <- 'Percent clay'
		units <- 'mm'
		range_fx <- percent_range_fx
	} else if (pred == 'soc') {
		short <- 'Soil Organic Matter'
		long <- 'Soil organic matter'
		units <- '???'
		range_fx <- percent_range_fx
	} else if (pred %in% c('site_nitrogen', 'nitrogen')) {
		short <- 'Soil Nitrogen'
		long <- 'Soil nitrogen'
		units <- '%'
		range_fx <- percent_range_fx
	} else if (pred %in% c('insolation_2000_growing_season_kWh_per_m2', 'insolation_2023_growing_season_kWh_per_m2')) {
		short <- 'Solar Insolation'
		long <- 'Growing season solar insolation'
		units <- 'kWh∙m¯²'
		range_fx <- standard_range_fx
	}

	list(short = short, long = long, units = units, range_fx = range_fx)

}

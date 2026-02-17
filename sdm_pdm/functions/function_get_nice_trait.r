#' Get R-friendly trait name from raw data column name
#'
#' @param raw Name of trait in raw data tables.
#' @returns The R-friendly version of `raw`.
get_rfriendly_trait_name <- function(raw) {

	# get traits
	raws <- c('Biomass', 'Delta13C', 'N_conc', 'CN_ratio', 'Height', 'BladeWidth', 'LeafThick', 'SPAD', 'CanopyDiam', 'WatPot', 'PhotoRate', 'StomCond', 'IntCO2', 'TranspRate')
	
	rfriendlies <- c('biomass', 'delta13c', 'n_concentration', 'cn_ratio', 'height', 'blade_width', 'leaf_thickness', 'spad', 'canopy_diameter', 'water_potential', 'photosynthetic_rate', 'stomatal_conductance', 'internal_co2', 'transpiration_rate')

	rfriendlies[match(raw, raws)]

}

#' Get "raw" "trait name from "R-friendly" name
#'
#' @param rfriendly Name of trait in R-friendly format.
#' @returns The "raw" version of `rfriendly`.
get_raw_trait_name_from_rfriendly <- function(rfriendly) {

	raws <- c('Biomass', 'Delta13C', 'N_conc', 'CN_ratio', 'Height', 'BladeWidth', 'LeafThick', 'SPAD', 'CanopyDiam', 'WatPot', 'PhotoRate', 'StomCond', 'IntCO2', 'TranspRate')

	rfriendlies <- c('biomass', 'delta13c', 'n_concentration', 'cn_ratio', 'height', 'blade_width', 'leaf_thickness', 'spad', 'canopy_diameter', 'water_potential', 'photosynthetic_rate', 'stomatal_conductance', 'internal_co2', 'transpiration_rate')

	raws[match(rfriendly, rfriendlies)]

}

#' Title, extended title, and units of a trait for plots/reports
#'
#' @param rfriendly Name of trait in R-friendly format.
#' @returns Plot/report/table title, axis title, a nice legend title (with line breaks and units), and units of trait, plus a function useful for returning a nice range of the variable. Each function requires values of the variable, plus a "mult" argument for increasing/decreasing the min/max values by a small amount (default = 0.05, representing an increase/decrease of 5%).
get_nice_trait <- function(rfriendly) {

	if (rfriendly == 'occs') {
		short <- 'Abundance'
		long <- 'Abundance'
		units <- 'Abundance'
		legend_title <- 'Abundance'
		range_fx <- bottom0_range_fx
	} else if (rfriendly == 'biomass') {
		short <- 'Biomass'
		long <- 'Aboveground biomass'
		units <- 'Biomass'
		legend_title <- 'Biomass\n(g)'
		range_fx <- bottom0_range_fx
	} else if (rfriendly == 'delta13c') {
		short <- 'δ¹³C'
		long <- 'δ¹³C'
		units <- '???'
		legend_title <- 'δ¹³C'
		range_fx <- standard_range_fx
	} else if (rfriendly == 'n_concentration') {
		short <- 'N concentration'
		long <- 'Leaf N concentration (%)'
		units <- '% dry weight'
		legend_title <- 'N conc\n(%).'
		range_fx <- percent_range_fx
	} else if (rfriendly == 'cn_ratio') {
		short <- 'C:N ratio'
		long <- 'Leaf C:N ratio'
		units <- ''
		legend_title <- 'C:N'
		range_fx <- standard_range_fx
	} else if (rfriendly == 'height') {
		short <- 'Height'
		long <- 'Height (cm)'
		units <- 'cm'
		legend_title <- 'Height\n(cm)'
		range_fx <- bottom0_range_fx
	} else if (rfriendly == 'blade_width') {
		short <- 'Blade width'
		long <- 'Blade width (cm)'
		units <- 'cm'
		legend_title <- 'Blade\nwidth\n(cm)'
		range_fx <- bottom0_range_fx
	} else if (rfriendly == 'leaf_thickness') {
		short <- 'Leaf thickness'
		long <- 'Leaf thickness (mm)'
		units <- 'mm'
		legend_title <- 'Leaf\nthickness\n(mm)'
		range_fx <- bottom0_range_fx
	} else if (rfriendly == 'spad') {
		short <- 'SPAD'
		long <- 'Soil-plant analyses development (SPAD)'
		units <- ''
		legend_title <- 'SPAD'
		range_fx <- bottom0_range_fx
	} else if (rfriendly == 'canopy_diameter') {
		short <- 'Canopy diameter'
		long <- 'Canopy diameter (cm)'
		units <- 'cm'
		legend_title <- 'Canopy\ndia.\n(cm)'
		range_fx <- bottom0_range_fx
	} else if (rfriendly == 'water_potential') {
		short <- 'Water potential'
		long <- 'Mid-day water potential (bars)'
		units <- '(bars)'
		legend_title <- 'Water\npotential\n(bars)'
		range_fx <- top0_range_fx
	} else if (rfriendly == 'photosynthetic_rate') {
		short <- 'Photosynthetic rate'
		long <- 'Photosynthetic rate (mol CO₂ m⁻² s⁻¹)'
		units <- 'mol CO₂ m⁻² s⁻¹'
		legend_title <- 'Photo.\nrate\'(mol CO₂ m⁻² s⁻¹)'
		range_fx <- bottom0_range_fx
	} else if (rfriendly == 'stomatal_conductance') {
		short <- 'Stomatal conductance'
		long <- 'Stomatal conductance (mol H2O₂ m⁻² s⁻¹)'
		units <- 'mol H2O₂ m⁻² s⁻¹'
		legend_title <- 'Stomatal\nconduct.\n()mol H2O₂ m⁻² s⁻¹)'
		range_fx <- bottom0_range_fx
	} else if (rfriendly == 'internal_co2') {
		short <- 'CO₂ concentration'
		long <- 'Internal CO₂ concentration (μmol CO₂ mol air⁻¹)'
		units <- 'μmol CO₂ mol air⁻¹'
		legend_title <- 'Internal\nCO₂\n(μmol CO₂ mol air⁻¹)'
		range_fx <- bottom0_range_fx
	} else if (rfriendly == 'transpiration_rate') {
		short <- 'Transpiration rate'
		long <- 'Transpiration rate (μmol s⁻¹)'
		units <- 'μmol s⁻¹'
		legend_title <- 'Transpir.\nrate\n(μmol s⁻¹)'
		range_fx <- bottom0_range_fx
	}

	list(short = short, long = long, units = units, legend_title = legend_title, range_fx = range_fx)

}

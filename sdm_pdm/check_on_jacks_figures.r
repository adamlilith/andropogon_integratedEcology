rm(list = ls())

setwd('C:/Kaji/Research/Andropogon/Andropogon')
library(data.table) # fast data tables
library(enmSdmX) # GIS/SDM
library(ggplot2) # graphics
library(patchwork) # graphics
library(terra) # GIS

source('C:/Kaji/R/andropogon_integratedEcology/sdm_pdm/function_prepare_biomass.r')
source('C:/Kaji/R/andropogon_integratedEcology/sdm_pdm/function_create_response_curve_array.r')
# future scenarios	
futs <- c(
	'ssp245_2041_2070',
	'ssp245_2071_2100',
	'ssp370_2041_2070',
	'ssp370_2071_2100'
)
biomass_data <- prepare_biomass_data(formula = ~ 1 + ph)

biomass_data_raw <- biomass_data$raw_data_biomass

graph <- ggplot(biomass_data_raw, aes(x = site_ph, y = Biomass)) +
	geom_point() +
	geom_smooth(method = 'lm', formula = y ~ x, se = TRUE) +
	xlab('Field sampled pH') +
	ylab('Biomass (g)') +
	ggtitle('Biomass') +
	theme_minimal()

library(nlme)
lme <- lme(Biomass ~ site_ph, random = ~ 1|SITE, data = biomass_data_raw)
summary(lme)

library(lme4)
lmer <- lmer(Biomass ~ site_ph + (1|SITE), data = biomass_data_raw)
summary(lmer)
anova(lmer)

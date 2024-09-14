### MODELING ANDROPOGON GERARDI DISTRIBUTION, MORPHOLOGY, PHYSIOLOGY, and ASSOCIATED MICROBIAL COMMUNITIES
### Erica Newman | Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2023-11
###
### This script analyzes candidate predictors for collinearity.
###
### source('C:/Ecology/Drive/Research/Andropogon/Andropogon/andropogon_integratedEcology/sdm_01_select_predictors.r')
### source('E:/Adam/Research/Andropogon/Andropogon/andropogon_integratedEcology/sdm_01_select_predictors.r')
###
### CONTENTS ###
### setup ###
### cluster analysis on candidate predictors based on correlation matrix ###

#############
### setup ###
#############

rm(list = ls())

drive <- 'C:/Ecology/Drive/'
# drive <- 'E:/Adam/'

.libPaths(paste0(drive, '/R/libraries'))
setwd(paste0(drive, '/Research/Andropogon/Andropogon'))

library(cluster)
library(omnibus)
library(terra)

say('############################################################################')
say('### cluster analysis on candidate predictors based on correlation matrix ###')
say('############################################################################')

ag_vect <- vect('./data/occurrence_data/andropogon_gerardi_occurrences_with_environment.gpkg')
ag_vect <- ag_vect[!is.na(ag_vect$any_ag_quality1to3), ]

preds <- c(paste0('bio', c(1:7, 10:12, 15, 18:19)), 'gdd_5_deg', 'climatic_moisture_index', 'pet_warmest_quarter_mm', 'aridity', 'ph', 'cec', 'sand', 'silt', 'clay', 'soc')
vars <- as.data.frame(ag_vect)
vars <- vars[ , preds]
vars <- vars[complete.cases(vars), ]
cor <- cor(vars, method = 'spearman')
dists <- 1 - abs(cor)
dists <- as.dist(dists)

clust <- agnes(dists, method = 'average')

png('./outputs/Correlations between SDM Predictors - Cluster Diagram.png', width = 1200, height = 800)
par(cex = 2)
plot(clust, which.plot = 2, ylab = '1 - abs(correlation)', xlab = '')
abline(h = 0.3, col = 'red')
dev.off()


### MODELING OF MICROBES RELATED TO ANDROPOGON GERARDI: CONSTRUCT MODEL(S)
### Erica Newman & Adam B. Smith | 2024-07
###
### source('C:/Ecology/Research/Andropogon/Andropogon/andropogon_integratedEcology/hmsc_vs_microbes/S2_fit_models.r')
### source('E:/Adam/Research/Andropogon/Andropogon/andropogon_integratedEcology/hmsc_vs_microbes/S2_fit_models.r')
###
###	SCRIPT INPUT: Unfitted model(s) created using Hmsc(), ALL PLACED IN THE './data' FOLDER:
### SCRIPT OUTPUT: Fitted models saved in file './models/fit_model.rds'

### setup
#########

rm(list = ls())
set.seed(1)

# workDir <- 'C:/Ecology/Research/Andropogon/Andropogon/andropogon_integratedEcology/hmsc_vs_microbes'
workDir <- 'E:/Adam/Research/Andropogon/Andropogon/andropogon_integratedEcology/hmsc_vs_microbes'

setwd(workDir)

library(BayesLogit)
library(Hmsc)
library(ape) # we need this to construct a taxonomic tree
library(ggplot2)
library(omnibus)

### settings
############

# number of MCMC iterations
samples <- 1000

# burn-in
transient <- 100

# thinning rate
thin <- 1

nParallel <- NULL # default: nParallel = nChains, set to 1 to disable parallel execution [ABS:: does not make sense given below]
nChains <- 2

if (is.null(nParallel)) nParallel <- nChains

model <- readRDS('./models/unfitted_model.rds')

### MCMC
########

fit <- sampleMcmc(
	model, samples = samples, thin = thin,
	adaptNf = rep(ceiling(0.4 * samples * thin), model$nr), 
	transient = ceiling(0.5 * samples * thin),
	nChains = nChains, nParallel = nParallel,
	verbose = TRUE
)

saveRDS(fit, file = './models/fit_model.rds')

say('DONE', deco = '!', level = 1)

### MODELING ANDROPOGON GERARDI DISTRIBUTION, PHENOTYPE, PHYSIOLOGY, GENOTYPE, and ASSOCIATED MICROBIAL COMMUNITIES
### Erica Newman | Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2023-12
###
### This script constructs a non-integrated model for AG geographic distribution.
###
### source('C:/Ecology/R/andropogon_integratedEcology/sdm_02_nonintegrated_model_n_mixture.r')
### source('E:/Adam/R/andropogon_integratedEcology/sdm_02_nonintegrated_model_n_mixture.r')
###
### CONTENTS ###
### setup ###
### non-integrated SDM ###
### model diagnostics ###
### map of current distribution ###

#############
### setup ###
#############

rm(list = ls())

drive <- 'C:/Ecology/'
# drive <- 'E:/Adam/'

.libPaths(paste0(drive, '/R/libraries'))
setwd(paste0(drive, '/Research/Andropogon/Andropogon'))

library(bayesplot) # graphing
library(coda) # Bayesian diagnostics
library(nimble) # Bayes
# library(nimbleHMC) # Hamiltonian Monte Carlo samplers
library(omnibus)
library(scales) # for plotting transparency
library(terra) # spatial objects

sink('./outputs/nonintegrated_sdm.txt', split = TRUE)
say()
say('##########################')
say('### non-integrated SDM ###')
say('##########################')
say(date(), post = 1)

say('This model is for the spatial distribution of AG. Currently, it assumes distribution is driven only by climate. Occurrences are at the county level, so county area is used as an offset.', breaks = 60, post = 1)

### MCMC settings
niter <- 90000
nburnin <- 10000
thin <- 80
nchains <- 4

say('MCMC settings:', level = 2)
say('niter ......', niter)
say('nburnin ... ', nburnin)
say('thin .......', thin)
say('nchains ... ', nchains, post = 1)

### model predictors and terms
say('Predictors and model formula:', level = 2)

# predictor_names <- c('aridity', 'bio7', 'bio2', 'sand', 'bio12', 'bio15', 'cec', 'soc', 'silt')
predictor_names <- c('aridity', 'bio7', 'ph', 'sand', 'bio12', 'bio15', 'cec', 'soc', 'silt')
say('We are using predictors: ', paste(predictor_names, collapse = ' '), post = 1)

preselected_model_terms <- read.csv('./outputs/sig_coeffs_elastic_net_2024_06_04.csv')

# get vector of linear predictors... we need to ensure all predictors appear at least as linear terms to respect marginality
terms <- preselected_model_terms$term
terms <- terms[terms %notin% c('(Intercept)', 'log_area_km2')]
linear_terms <- gsub(terms, pattern = 'I\\(', replacement = '')
linear_terms <- gsub(linear_terms, pattern = '\\^2)', replacement = '')
linear_terms <- strsplit(linear_terms, ':')
linear_terms <- unlist(linear_terms)
linear_terms <- unique(linear_terms)
linear_terms <- sort(linear_terms)

# get non-linear terms
terms <- terms[terms %notin% linear_terms]

# create model formula
form <- paste0(' ~ 1 + ', paste(c(linear_terms, terms), collapse = ' + '))
say('Model formula:')
say(form, post = 1, breaks = 80)

### load AG data
ag_vect <- vect('./data/occurrence_data/andropogon_gerardi_occurrences_with_environment.gpkg')

fields <- c('area_km2', 'any_ag_quality1to3', 'num_poaceae_records', predictor_names)
ag_vect <- ag_vect[ , fields]

ag <- as.data.frame(ag_vect)
completes <- complete.cases(as.data.frame(ag_vect))
ag_focus <- ag[completes, ]
ag_vect_focus <- ag_vect[completes, ]

### collate data

### county area... used as an offset to make underlying model fit an IPP
area_km2 <- ag_focus$area_km2
log_area_km2 <- log(area_km2)
log_area_km2_scaled <- scale(log_area_km2)
log_area_center <- attributes(log_area_km2_scaled)$`scaled:center`
log_area_scale <- attributes(log_area_km2_scaled)$`scaled:scale`
log_area_km2_scaled <- log_area_km2_scaled[ , 1]

### number of Poaceae records... used to model sampling bias
log_num_poaceae_records <- log1p(ag_focus$num_poaceae_records) # log(x + 1)
log_num_poaceae_records_scaled <- scale(log_num_poaceae_records)
log_num_poaceae_records_scaled <- log_num_poaceae_records_scaled[ , 1]

### subset, scale, and manipulate predictors into model frame
form <- as.formula(form)
x_raw <- as.data.frame(ag_focus[ , c('area_km2', predictor_names)])
log_area_km2 <- log(x_raw$area_km2)
log_area_km2 <- data.frame(log_area_km2 = log_area_km2)
x_raw$area_km2 <- NULL
x_raw <- cbind(log_area_km2, x_raw)
x_raw_scaled <- scale(x_raw)
x_raw_scaled <- as.data.frame(x_raw_scaled)
x <- model.matrix(form, as.data.frame(x_raw_scaled))

### inputs for nimble
say('Inputs:', level = 2)
data <- list(
	y = ag_focus$any_ag_quality1to3 # number of AG records in each county
)

n_counties <- nrow(ag_focus)
n_sdm_terms <- ncol(x)

constants <- list(
  n_counties = n_counties,
  n_sdm_terms = n_sdm_terms,
  x = x,
  log_area_km2_scaled = log_area_km2_scaled,
  log_num_poaceae_records_scaled = log_num_poaceae_records_scaled
)

inits <- list(
  beta = rep(0, n_sdm_terms),
  alpha0_sampling = 0,
  alpha_area = 0,
  alpha_poaceae = 0,
  N = 100 * ag_focus$any_ag_quality1to3,
  lambda = 100 * ag_focus$any_ag_quality1to3
)

say('Data:')
print(str(data))

say('Constants:', pre = 1)
print(str(constants))

say('Initializations:', pre = 1)
print(str(inits))

### define model
say('nimbleCode():', level = 2)
say('We assume an N-mixture model (latent, real abundance ~ Poisson, and observations ~ binomial draws from latent abundance.')
model_code <- nimbleCode({
  
	# priors for relationship to environment
	beta[1] ~ dnorm(0, sd = 10) # intercept... not regularized
	for (j in 2:n_sdm_terms) {
		beta[j] ~ ddexp(0, rate = 1) # regularization toward 0
		# beta[j] ~ dnorm(0, sd = 10) # broad prior
	}

	# priors for sampling bias
	alpha0_sampling ~ dnorm(0, sd = 10)
	alpha_area ~ dnorm(0, sd = 10)
	alpha_poaceae ~ dnorm(0, sd = 10)

	# likelihood
	for (i in 1:n_counties) {    # this specifies estimates be made for each county
		
		### actual abundance (latent--unobserved)
		N[i] ~ dpois(lambda[i])

		# relationship between latent abundance and environment
		log(lambda[i]) <- inprod(beta[1:n_sdm_terms], x[i, 1:n_sdm_terms])

		### observed number of AG
		y[i] ~ dbin(prob = p[i], size = N[i])

		# sampling bias
		logit(p[i]) <- alpha0_sampling + alpha_area * log_area_km2_scaled[i] + alpha_poaceae * log_num_poaceae_records_scaled[i]

	}
  
})

say('nimbleModel():', level = 2)
model_species <- nimbleModel(
	code = model_code, # our model
	constants = constants, # constants
	data = data, # data
	inits = inits, # initialization values
	check = TRUE, # any errors?
	calculate = FALSE
	# buildDerivs = TRUE # need for Hamiltonian Monte Carlo
)

say('$initializeInfo() and $calculate():', level = 2)
model_species$initializeInfo()
model_species$calculate()

say('configureMCMC():', level = 2)
monitors <- c('beta', 'alpha0_sampling', 'alpha_area', 'alpha_poaceae', 'lambda')
conf <- configureMCMC(
  model_species,
  monitors = monitors,
  print = TRUE,
  enableWAIC = FALSE
)

# # add no U-turn sampler (Hamiltonian Monte Carlo)
# conf$addSampler(target = c('beta', 'alpha_area', 'alpha_poaceae'), type = 'NUTS')

### compile/build/run model/save MCMC
build <- buildMCMC(conf)

compiled <- compileNimble(model_species, build, showCompilerOutput = FALSE)

chains <- runMCMC(
  compiled$build,
  niter = niter,
  nburnin = nburnin,
  thin = thin,
  nchains = nchains,
  inits = inits,
  progressBar = TRUE,
  samplesAsCodaMCMC = TRUE,
  summary = TRUE,
  WAIC = FALSE,
  perChainWAIC = FALSE
)

saveRDS(chains, './outputs/nonintegrated_sdm_chains.rds')

say(date())
sink()

say('#########################')
say('### model diagnostics ###')
say('#########################')

chains <- readRDS('./outputs/nonintegrated_sdm_chains.rds')

# graphing trace plots for all betas
png('./outputs/nonintegrated_sdm_beta_trace.png', width = 1800, height = 1200)
pars <- paste0('beta[', 1:n_sdm_terms, ']')
print(mcmc_trace(chains$samples, pars = pars))
dev.off()

# graphing trace plots for all "extra" betas
png('./outputs/nonintegrated_sdm_beta_extra_trace.png', width = 1500, height = 800)
pars <- c('alpha0_sampling', 'alpha_poaceae', 'alpha_area')
print(mcmc_trace(chains$samples, pars = pars))
dev.off()

# graphing density plots for all betas
png('./outputs/nonintegrated_sdm_beta_density.png', width = 1800, height = 1200)
pars <- paste0('beta[', 1:n_sdm_terms, ']')
print(mcmc_dens_overlay(chains$samples, pars = pars))
dev.off()

# graphing trace plots for all "extra" betas
png('./outputs/nonintegrated_sdm_beta_extra_density.png', width = 2000, height = 800)
pars <- c('alpha0_sampling', 'alpha_poaceae', 'alpha_area')
print(mcmc_dens_overlay(chains$samples, pars = pars))
dev.off()

say('###################################')
say('### map of current distribution ###')
say('###################################')

chains <- readRDS('./outputs/nonintegrated_sdm_chains.rds')

# subset chain summary to just the lambda's associated with background sites
summary <- chains$summary$all.chains

which_lambda <- grepl(rownames(summary), pattern = 'lambda')
lambda <- summary[which_lambda, ]

# # get just counties with data
ag_vect_focus$lambda_mean <- lambda[ , 'Mean']
ag_vect_focus$lambda_0.05ci <- lambda[ , '95%CI_low']
ag_vect_focus$lambda_0.95ci <- lambda[ , '95%CI_upp']
ag_vect_focus$lambda_ci <- ag_vect_focus$lambda_0.95ci - ag_vect_focus$lambda_0.05ci

ag_vect_outline <- crop(ag_vect, ext(ag_vect_focus))

quants <- quantile(ag_vect_focus$lambda_mean, c(0.25, 0.5, 0.75, 0.95))
ag_vect_focus$quant_col <- NA
ag_vect_focus$quant_col[ag_vect_focus$lambda_mean < quants[1]] <- 'gray90'
ag_vect_focus$quant_col[ag_vect_focus$lambda_mean >= quants[1] & ag_vect_focus$lambda_mean < quants[2]] <- alpha('forestgreen', 0.25)
ag_vect_focus$quant_col[ag_vect_focus$lambda_mean >= quants[2] & ag_vect_focus$lambda_mean < quants[3]] <- alpha('forestgreen', 0.5)
ag_vect_focus$quant_col[ag_vect_focus$lambda_mean >= quants[3] & ag_vect_focus$lambda_mean < quants[4]] <- alpha('forestgreen', 0.75)
ag_vect_focus$quant_col[ag_vect_focus$lambda_mean >= quants[4]] <- alpha('forestgreen', 1)

nam <- vect(paste0(drive, '/Research Data/GADM/Version 4.1/High Res North America Level 1 sans Great Lakes SpatVector WGS84.gpkg'))
nam <- project(nam, ag_vect)
png('./outputs/nonintegrated_sdm_lambda.png', width = 1200, height = 1000, res = 300)

	par(oma = rep(0, 4), mar = rep(0, 4))

	plot(ag_vect_focus, col = NA, border = NA)
	plot(ag_vect, col = 'gray95', border = NA, add = TRUE)
	col <- ag_vect_focus$quant_col
	plot(ag_vect_focus, col = col, border = NA, lwd = 0.01, add = TRUE)
	plot(nam, border = 'gray65', lwd = 0.01, add = TRUE)

dev.off()

png('./outputs/nonintegrated_sdm_lambda_90ci_vs_mean.png', width = 1200, height = 1000, res = 300)

	par(oma = rep(0, 4), mar = rep(0, 4))

	plot(ag_vect_focus, col = NA, border = NA)
	plot(ag_vect, col = 'gray95', border = NA, add = TRUE)
	col <- ag_vect_focus$lambda_ci / ag_vect_focus$lambda_mean
	col <- col - min(col)
	col <- col / max(col)
	plot(ag_vect_focus, col = alpha('darkred', col), border = NA, lwd = 0.01, add = TRUE)
	plot(nam, border = 'gray65', lwd = 0.01, add = TRUE)

dev.off()

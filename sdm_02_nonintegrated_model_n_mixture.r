### MODELING ANDROPOGON GERARDI DISTRIBUTION, PHENOTYPE, PHYSIOLOGY, GENOTYPE, and ASSOCIATED MICROBIAL COMMUNITIES
### Erica Newman | Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2023-12
###
### This script constructs a non-integrated model for AG geographic distribution.
###
### source('C:/Ecology/Drive/Research/Andropogon/Andropogon/andropogon_integratedEcology/sdm_02_nonintegrated_model_n_mixture.r')
### source('E:/Adam/Research/Andropogon/Andropogon/andropogon_integratedEcology/sdm_02_nonintegrated_model_n_mixture.r')
###
### CONTENTS ###
### setup ###
### user-defined constants ###


#############
### setup ###
#############

rm(list = ls())

drive <- 'C:/Ecology/Drive/'
# drive <- 'E:/Adam/'

.libPaths(paste0(drive, '/R/libraries'))
setwd(paste0(drive, '/Research/Andropogon/Andropogon'))

library(bayesplot) # graphing
library(coda) # Bayesian diagnostics
library(nimble) # Bayes
library(nimbleHMC) # Hamiltonian Monte Carlo samplers
library(omnibus)
library(scales) # for plotting transparency
library(terra) # spatial objects

say('##############################')
say('### user-defined constants ###')
say('##############################')

# predictor terms selected using glmnet
# predictor_names <- c('bio2', 'bio3', 'bio4', 'bio5', 'bio6', 
                     # 'bio7', 'bio10', 'bio12', 'bio15', 'bio18', 'bio19', 
                     # 'pet_warmest_quarter_mm', 'gdd_5_deg',
                     # 'climatic_moisture_index')

predictor_names <- c('bio2', 'bio7', 'bio12', 'bio15', 'aridity')

say('We are using predictors: ', paste(predictor_names, collapse = ' '))

say('##########################')
say('### non-integrated SDM ###')
say('##########################')

say('This model is for the spatial distribution of AG. Currently, it assumes distribution is driven only by climate. We should later consider adding soil. Occurrences are at the county level, so county area is used as an offset.', breaks = 60)

### load data
ag_vect <- vect('./data/occurrence_data/andropogon_gerardi_occurrences_with_environment.gpkg')

fields <- c('area_km2', 'any_ag_quality1to3', 'num_poaceae_records', predictor_names)
ag_vect <- ag_vect[ , fields]

ag <- as.data.frame(ag_vect)
completes <- complete.cases(as.data.frame(ag_vect))
ag_focus <- ag[completes, ]
ag_vect_focus <- ag_vect[completes, ]

### collate data

### county area
area_km2 <- ag_focus$area_km2
log_area_km2 <- log(area_km2)
log_area_km2_scaled <- scale(log_area_km2)
log_area_center <- attributes(log_area_km2_scaled)$`scaled:center`
log_area_scale <- attributes(log_area_km2_scaled)$`scaled:scale`

log_area_km2_scaled <- log_area_km2_scaled[ , 1]

### number of Poaceae records
log_num_poaceae_records <- log1p(ag_focus$num_poaceae_records) # log(x + 1)
log_num_poaceae_records_scaled <- scale(log_num_poaceae_records)
log_num_poaceae_records_scaled <- log_num_poaceae_records_scaled[ , 1]

### PCA on predictors
x_raw <- as.data.frame(ag_focus[ , predictor_names])
x_scaled <- scale(x_raw)
centers <- attributes(x_scaled)$`scaled:center`
scales <- attributes(x_scaled)$`scaled:scale`

### construct model matrix
form <- paste0('~ 1 + ',
	paste(predictor_names, sep = '', collapse = ' + '),
	' + ',
	paste('I(', predictor_names, '^2)', sep = '', collapse = ' + ')
)

for (i in 1:(length(predictor_names) - 1)) {
	for (j in (i + 1):length(predictor_names)) {
		form <- paste0(form, ' + ', predictor_names[i], ':', predictor_names[j])
	}
}
form <- as.formula(form)

x <- model.matrix(form, as.data.frame(x_scaled))

# okay, we've just made the data frame with all the predictors we're going to
# use. Now we're going to run the SDM in Nimble

##########################################

### make nimble lists (data, constants, initializations)

data <- list(
	y = ag_focus$any_ag_quality1to3
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

### define model
model_code <- nimbleCode({
  
  # priors
  beta[1] ~ dnorm(0, sd = 10) # intercept... not regularized
  for (j in 2:n_sdm_terms) {
    beta[j] ~ ddexp(0, rate = 1) # regularization toward 0
    # beta[j] ~ dnorm(0, sd = 10) # broad prior
  }
  
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
	logit(p[i]) <- alpha0_sampling +
		alpha_area * log_area_km2_scaled[i] +
		alpha_poaceae * log_num_poaceae_records_scaled[i]

  }
  
  
})

model_species <- nimbleModel(
	code = model_code, # our model
	constants = constants, # constants
	data = data, # data
	inits = inits, # initialization values
	check = TRUE, # any errors?
	calculate = FALSE
	# buildDerivs = TRUE # need for Hamiltonian Monte Carlo
)

model_species$initializeInfo()
model_species$calculate()

set.seed(2069)

say('Turn on WAIC for final runs!')
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
  niter = 40000,
  nburnin = 20000,
  thin = 20,
  nchains = 4,
  inits = inits,
  progressBar = TRUE,
  samplesAsCodaMCMC = TRUE,
  summary = TRUE,
  WAIC = FALSE,
  perChainWAIC = FALSE
)

saveRDS(chains, './outputs/sdm_chains.rds')

# graphing trace plots for all betas
png('./outputs/sdm_beta_trace.png', width = 1800, height = 1200)
pars <- paste0('beta[', 1:n_sdm_terms, ']')
print(mcmc_trace(chains$samples, pars = pars))
dev.off()

# graphing trace plots for all "extra" betas
png('./outputs/sdm_beta_extra_trace.png', width = 1500, height = 800)
pars <- c('alpha0_sampling', 'alpha_poaceae', 'alpha_area')
print(mcmc_trace(chains$samples, pars = pars))
dev.off()

# graphing density plots for all betas
png('./outputs/sdm_beta_density.png', width = 1800, height = 1200)
pars <- paste0('beta[', 1:n_sdm_terms, ']')
print(mcmc_dens_overlay(chains$samples, pars = pars))
dev.off()

# graphing trace plots for all "extra" betas
png('./outputs/sdm_beta_extra_density.png', width = 2000, height = 800)
pars <- c('alpha0_sampling', 'alpha_poaceae', 'alpha_area')
print(mcmc_dens_overlay(chains$samples, pars = pars))
dev.off()

# # graph the posterior probability of occurrence
# pars <- paste0('beta[', 1:3, ']')
# #mcmc_combo(chains$samples, pars = 'lambda[1]')
# mcmc_combo(chains$samples, pars = pars)
# #dev.off()

################
# now we have a 'chains' object with the lambda's, which are the county-level 
# probability of occurrence

# str(chains)
# corner(chains$summary$all.chains)
# head(chains$summary$all.chains, 40)

### map of current distribution
###############################

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
png('./outputs/sdm_lambda.png', width = 1200, height = 1000, res = 300)

	par(oma = rep(0, 4), mar = rep(0, 4))

	plot(ag_vect_focus, col = NA, border = NA)
	plot(ag_vect, col = 'gray95', border = NA, add = TRUE)
	col <- ag_vect_focus$quant_col
	plot(ag_vect_focus, col = col, border = NA, lwd = 0.01, add = TRUE)
	plot(nam, border = 'gray65', lwd = 0.01, add = TRUE)

dev.off()

png('./outputs/sdm_lambda_90ci_vs_mean.png', width = 1200, height = 1000, res = 300)

	par(oma = rep(0, 4), mar = rep(0, 4))

	plot(ag_vect_focus, col = NA, border = NA)
	plot(ag_vect, col = 'gray95', border = NA, add = TRUE)
	col <- ag_vect_focus$lambda_ci / ag_vect_focus$lambda_mean
	col <- col - min(col)
	col <- col / max(col)
	plot(ag_vect_focus, col = alpha('darkred', col), border = NA, lwd = 0.01, add = TRUE)
	plot(nam, border = 'gray65', lwd = 0.01, add = TRUE)

dev.off()

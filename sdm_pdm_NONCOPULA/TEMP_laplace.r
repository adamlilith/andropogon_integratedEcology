

	rm(list = ls())

	drive <- 'C:/Kaji/'

	setwd(paste0(drive, '/Research/Andropogon/Andropogon'))
	source(paste0(drive, '/R/andropogon_integratedEcology/sdm_pdm/sdm_pdm_00_shared_functions_and_variables.r'))

library(nimbleHMC) # Bayesian modeling

	formula_biomass = ~ 1 + bio12
	data_biomass <- prepare_biomass_data(formula_biomass = formula_biomass, n_response_curve_values = n_response_curve_values, calib = FALSE)


data <- list(y_biomass = data_biomass$y_biomass)
constants <- list(

	n_biomass = data_biomass$n_biomass,
	x_by_site_biomass = data_biomass$x_by_site_biomass,
	n_pheno_sites = data_biomass$n_pheno_sites,
	site_index_biomass = data_biomass$site_index_biomass,

	n_counties = data_biomass$n_counties,

	n_terms_biomass = 2
)

inits <- list(
	beta_biomass = c(-2, 2),
	site_biomass_plant_sigma_log = 1,
	site_biomass_mean_sigma = 1,
	log_site_biomass_mu = rep(1, 26),
	y_biomass_sim = data_biomass$y_biomass

)

code <- nimbleCode({

	# BIOMASS: priors for relationship of site-level mean biomass to environment
	for (i in 1:n_terms_biomass) {
		beta_biomass[i] ~ dnorm(0, sd = 10) # broad prior
	}

	# prior for sd of site-level biomass on lognormal (~ half-Cauchy), ==> vague
	site_biomass_plant_sigma_log ~ dnorm(0, sd = cross_site_biomass_plant_sigma_log_prior_sd)
	site_biomass_plant_sigma_log <- exp(site_biomass_plant_sigma_log)

	# prior for sd of mean of biomass at a site on lognormal (~ half-Cauchy), ==> vague
	sigma_biomass_among_sites ~ dnorm(0, sd = sigma_biomass_among_sites_log_prior_sd)
	site_biomass_mean_sigma <- exp(sigma_biomass_among_sites)

	# BIOMASS: parameters of biomass distribution are latent and functions of environment
	# individual plant biomasses are samples from the site-level distribution defined by the site-level distribution (next chunk after this one)
	for (i in 1:n_pheno_sites) {

		# relationship of biomass to the environment
		log_site_biomass_mu[i] ~ dnorm(log_site_biomass_mu_mean[i], sd = site_biomass_mean_sigma)
		log_site_biomass_mu_mean[i] <- inprod(beta_biomass[1:n_terms_biomass], x_by_site_biomass[i, 1:n_terms_biomass])

		# site-level mean biomass
		mu_biomass_site[i] <- exp(log_site_biomass_mu[i])

		# moment matching to get dgamma() parameters
		shape_biomass[i] <- mu_biomass_site[i]^2 / site_biomass_plant_sigma_log^2
		rate_biomass[i] <- mu_biomass_site[i] / site_biomass_plant_sigma_log^2

	}

	# BIOMASS: likelihood of individual plants
	for (i in 1:n_biomass) {

		# likelihood
		y_biomass[i] ~ dgamma(shape = shape_biomass[site_index_biomass[i]], rate = rate_biomass[site_index_biomass[i]])

		# simulated values for unconditional DHARMa residuals
		y_biomass_sim[i] ~ dgamma(shape = shape_biomass[site_index_biomass[i]], rate = rate_biomass[site_index_biomass[i]])

	}

	# BIOMASS: posterior samplers for predictions of to counties in status quo and future
	for (i in 1:n_counties) {

		# biomass: status quo
		log_biomass_mu_county_mu_mean_sq[i] <-
			inprod(beta_biomass[1:n_terms_biomass], counties_x_biomass_sq[i, 1:n_terms_biomass])
		log_biomass_mu_county_sq[i] ~ dnorm(log_biomass_mu_county_mu_mean_sq[i], sd = site_biomass_mean_sigma)

		# biomass: future
		log_biomass_county_mu_ssp245_2041_2070[i] <-
			inprod(beta_biomass[1:n_terms_biomass], counties_x_biomass_ssp245_2041_2070[i, 1:n_terms_biomass])
		mu_biomass_county_ssp245_2041_2070[i] <- exp(log_biomass_county_mu_ssp245_2041_2070[i])
		
		log_biomass_county_mu_ssp245_2071_2100[i] <-
			inprod(beta_biomass[1:n_terms_biomass], counties_x_biomass_ssp245_2071_2100[i, 1:n_terms_biomass])
		mu_biomass_county_ssp245_2071_2100[i] <- exp(log_biomass_county_mu_ssp245_2071_2100[i])
		
		log_biomass_county_mu_ssp370_2041_2070[i] <-
			inprod(beta_biomass[1:n_terms_biomass], counties_x_biomass_ssp370_2041_2070[i, 1:n_terms_biomass])
		mu_biomass_county_ssp370_2041_2070[i] <- exp(log_biomass_county_mu_ssp370_2041_2070[i])

		log_biomass_county_mu_ssp370_2071_2100[i] <-
			inprod(beta_biomass[1:n_terms_biomass], counties_x_biomass_ssp370_2071_2100[i, 1:n_terms_biomass])
		mu_biomass_county_ssp370_2071_2100[i] <- exp(log_biomass_county_mu_ssp370_2071_2100[i])

	}


})


model <- nimbleModel(
	code,
	data = data,
	constants = constants,
	# inits = inits,
	calculate = FALSE,
    buildDerivs = TRUE
)


model$calculate()

compiled <- compileNimble(model)

laplace <- buildLaplace(model, randomEffectsNodes = c('log_site_biomass_mu', 'log_biomass_mu_county_sq'))
compiled_laplace <- compileNimble(laplace, project = model)

results <- runLaplace(compiled_laplace)
results$summary

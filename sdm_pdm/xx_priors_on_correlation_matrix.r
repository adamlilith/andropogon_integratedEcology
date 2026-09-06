library('nimble')
rm(list = ls())

# 1. Custom distribution for the vector of upper-triangular elements
d_custom_cor_vec <- nimbleFunction(
  run = function(x = double(1),
                 mu_vec = double(1),
                 sd_vec = double(1),
                 n_dims = integer(0),
                 log = integer(0, default = 0)) {
    
    returnType(double(0))
    
    # Bound check for correlations
    for(i in 1:length(x)) {
      if(x[i] < -1 | x[i] > 1) return(-Inf)
    }
    
    # Build a temporary matrix to check positive definiteness
    mat <- identityMatrix(n_dims)
    counter <- 1
    for(i in 1:(n_dims - 1)) {
      for(j in (i + 1):n_dims) {
        mat[i, j] <- x[counter]
        mat[j, i] <- x[counter]
        counter <- counter + 1
      }
    }
    
    # Safe positive definite check using eigenvalues
    eig_vals <- eigen(mat, symmetric = TRUE)$values
    for(i in 1:n_dims) {
      if(eig_vals[i] <= 0) return(-Inf)
    }
    
    # Sum log-probabilities based on the custom priors
    log_prob <- 0
    for(i in 1:length(x)) {
      log_prob <- log_prob + dnorm(x[i], mean = mu_vec[i], sd = sd_vec[i], log = 1)
    }
    
    if(log) return(log_prob)
    else return(exp(log_prob))
  }
)

# Dummy generator required by NIMBLE
r_custom_cor_vec <- nimbleFunction(
  run = function(n = integer(0), mu_vec = double(1), sd_vec = double(1), n_dims = integer(0)) {
    returnType(double(1))
    return(rep(0, length(mu_vec))) # Initialize with 0s
  }
)

# 2. Deterministic function to rebuild the symmetric matrix for the likelihood
build_cor_mat <- nimbleFunction(
  run = function(rho_vec = double(1), n_dims = integer(0)) {
    returnType(double(2))
    mat <- identityMatrix(n_dims)
    counter <- 1
    for(i in 1:(n_dims - 1)) {
      for(j in (i + 1):n_dims) {
        mat[i, j] <- rho_vec[counter]
        mat[j, i] <- rho_vec[counter]
        counter <- counter + 1
      }
    }
    return(mat)
  }
)

# Register distributions and functions
registerDistributions(list(
  d_custom_cor_vec = list(
    BUGSdist = 'd_custom_cor_vec(mu_vec, sd_vec, n_dims)',
    types = c('value = double(1)', 'mu_vec = double(1)', 'sd_vec = double(1)', 'n_dims = integer(0)')
  )
))

# 3. Simulate data with a known 5x5 correlation structure
set.seed(42)
dim_size <- 5
n_obs <- 2000
n_tri <- (dim_size * (dim_size - 1)) / 2

true_cor <- diag(dim_size)
true_cor[1, 2:5] <- 0.3
true_cor[2:5, 1] <- 0.3
true_cor[2:5, 2:5] <- 0.7
diag(true_cor) <- 1

cholesky_true <- chol(true_cor)
sim_data <- matrix(rnorm(n_obs * dim_size), nrow = n_obs, ncol = dim_size) %*% cholesky_true

# 4. Construct prior vectors
# Elements 1:4 are correlations with Axis 1
# Elements 5:10 are correlations between Axes 2:5
mu_prior_vec <- c(rep(0, 4), rep(0.7, 6))
sd_prior_vec <- c(rep(0.5, 4), rep(0.05, 6))

# 5. Define the NIMBLE model
model_code <- nimbleCode({
  for(j in 1:n_dims) {
    mu[j] ~ dnorm(0, sd = 10)
  }
  
  # Sample the 10 upper-triangular elements as a single vector
  rho_vec[1:n_tri] ~ d_custom_cor_vec(mu_vec[1:n_tri], 
                                      sd_vec[1:n_tri], 
                                      n_dims)
  
  # Deterministically rebuild the matrix
  cor_mat[1:n_dims, 1:n_dims] <- build_cor_mat(rho_vec[1:n_tri], n_dims)
  
  # MVN likelihood
  for(i in 1:n_obs) {
    y[i, 1:n_dims] ~ dmnorm(mu[1:n_dims], cov = cor_mat[1:n_dims, 1:n_dims])
  }
})

model_constants <- list(
  n_dims = dim_size,
  n_tri = n_tri,
  n_obs = n_obs,
  mu_vec = mu_prior_vec,
  sd_vec = sd_prior_vec
)

model_data <- list(
  y = sim_data
)

# Initialize rho_vec with 0s, mapping to an identity matrix
model_inits <- list(
  mu = rep(0, dim_size),
  rho_vec = rep(0, n_tri)
)

# 6. Build, configure, and compile
nimble_model <- nimbleModel(code = model_code, 
                            constants = model_constants, 
                            data = model_data, 
                            inits = model_inits)

mcmc_conf <- configureMCMC(nimble_model, monitors = c('mu', 'cor_mat', 'rho_vec'))

# Assign block sampler to the vector to navigate the positive-definite bounds efficiently
mcmc_conf$removeSamplers('rho_vec')
mcmc_conf$addSampler(target = 'rho_vec', type = 'RW_block')

nimble_mcmc <- buildMCMC(mcmc_conf)

compiled_model <- compileNimble(nimble_model)
compiled_mcmc <- compileNimble(nimble_mcmc, project = nimble_model)

# 7. Run MCMC
mcmc_out <- runMCMC(compiled_mcmc, niter = 5000, nburnin = 1000, samplesAsCodaMCMC = TRUE)

# The output for cor_mat and rho_vec will now contain valid posterior estimates
summary(mcmc_out)

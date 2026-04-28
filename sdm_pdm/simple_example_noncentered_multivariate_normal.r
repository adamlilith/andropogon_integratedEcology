library(nimble)
library(mvtnorm)

n <- 100
s <- 2
x <- rnorm(n * s)
X <- matrix(x, n, s)

mu <- c(-2, 1)

rho <- -0.7
Sigma <- matrix(c(1, rho, rho, 1), nrow = 2)

y <- rmvnorm(n, mean = mu, sigma = Sigma)


data <- list(y = y)

constants <- list(n = n, n_resp = 2)

inits <- list(mu = c(0, 0))
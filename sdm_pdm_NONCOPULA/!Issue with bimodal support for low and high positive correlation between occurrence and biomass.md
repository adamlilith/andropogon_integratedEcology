# I constructed a multivariate normal model for two responses using nimble. I found that the correlation between the two responses was both strongly positive and weakly positive, depending on the value of the predictor variable. I discovered this because the trace plots show convergence around either strong positive or weak positive correlation, but not between. Is there a way to incorporate a predictor into the covariance of a multivariate normal model?

Yes, you *can* let a predictor enter the covariance of a multivariate normal – but you need to re-parameterize carefully so the covariance matrix stays positive-definite.

Let me split this into two parts:

1. **What your strange trace plots are probably telling you**
2. **How to build a MVN model with predictor-dependent covariance in `nimble`**

---

## 1. Interpreting the “strong vs weak correlation” traces

If your model is

[
\mathbf{y}_i \mid x_i \sim \mathcal{N}_2\big(\boldsymbol{\mu}(x_i), \Sigma\big),
]

with **one constant** covariance matrix (\Sigma) for all observations, then the correlation parameter (\rho) is *not* allowed to change with (x). So if the trace of (\rho) shows:

* one mode near, say, 0.2 (weak positive),
* another mode near, say, 0.8 (strong positive),
* and the chain doesn’t explore values *between* those modes,

that usually means:

* the data are **not very informative** about (\rho), **and/or**
* the current parameterization gives the sampler a nasty, multimodal posterior surface,
* or the mean structure is misspecified and the model is trying to “explain” patterns with covariance instead.

So the trace alone doesn’t *prove* that correlation truly varies with the predictor, but it’s a plausible scientific hypothesis to test by extending the model.

Which brings us to…

---

## 2. Letting a predictor enter the covariance in a 2-D MVN

For two responses, the covariance matrix at observation (i) is

[
\Sigma_i =
\begin{pmatrix}
\sigma_{1i}^2 & \rho_i \sigma_{1i} \sigma_{2i} \
\rho_i \sigma_{1i} \sigma_{2i} & \sigma_{2i}^2
\end{pmatrix}.
]

The trick is to model the **unconstrained** parameters as linear functions of predictors, then transform them back to valid variances/correlations:

* Variances/SDs must be **> 0**
  (\Rightarrow) model (\log \sigma_{1i}), (\log \sigma_{2i}) linearly.

* Correlation must be in **(–1, 1)**
  (\Rightarrow) model an unconstrained (\eta_{\rho,i}) linearly and transform:
  [
  \rho_i = 2\text{ilogit}(\eta_{\rho,i}) - 1
  ]
  where (\text{ilogit}(z) = 1/(1+\exp(-z))).

### Generic model structure

Let (x_i) be your predictor. Then you can write:

[
\begin{aligned}
\mu_{1i} &= \beta_{10} + \beta_{11} x_i, \
\mu_{2i} &= \beta_{20} + \beta_{21} x_i, \
\log \sigma_{1i} &= \alpha_{10} + \alpha_{11} x_i,\
\log \sigma_{2i} &= \alpha_{20} + \alpha_{21} x_i,\
\eta_{\rho,i} &= \gamma_0 + \gamma_1 x_i,\
\rho_i &= 2\text{ilogit}(\eta_{\rho,i}) - 1.
\end{aligned}
]

Then build (\Sigma_i) as above and use it in a `dmnorm` for each observation.

---

## 3. Sketch in `nimble` (2 responses, covariance depends on x)

Here’s a minimal pattern in `nimbleCode` style; adjust priors as you like. I’ll keep everything in a single block so you can see the structure:

```r
library(nimble)

code <- nimbleCode({
  # Likelihood
  for(i in 1:N) {
    ## Means as functions of x[i]
    mu[i, 1] <- beta10 + beta11 * x[i]
    mu[i, 2] <- beta20 + beta21 * x[i]
    
    ## Standard deviations as functions of x[i]
    log_sigma1[i] <- alpha10 + alpha11 * x[i]
    log_sigma2[i] <- alpha20 + alpha21 * x[i]
    sigma1[i] <- exp(log_sigma1[i])
    sigma2[i] <- exp(log_sigma2[i])
    
    ## Correlation as function of x[i], mapped to (-1, 1)
    eta_rho[i] <- gamma0 + gamma1 * x[i]
    rho[i] <- 2 * ilogit(eta_rho[i]) - 1   # or tanh(eta_rho[i])
    
    ## Build covariance matrix Sigma[i, , ]
    Sigma[i, 1, 1] <- pow(sigma1[i], 2)
    Sigma[i, 2, 2] <- pow(sigma2[i], 2)
    Sigma[i, 1, 2] <- rho[i] * sigma1[i] * sigma2[i]
    Sigma[i, 2, 1] <- Sigma[i, 1, 2]
    
    ## Bivariate normal for y[i, 1:2]
    y[i, 1:2] ~ dmnorm(mean = mu[i, 1:2],
                       cov  = Sigma[i, 1:2, 1:2])
  }
  
  # Priors for mean parameters
  beta10 ~ dnorm(0, sd = 10)
  beta11 ~ dnorm(0, sd = 10)
  beta20 ~ dnorm(0, sd = 10)
  beta21 ~ dnorm(0, sd = 10)
  
  # Priors for sd regressions
  alpha10 ~ dnorm(0, sd = 5)
  alpha11 ~ dnorm(0, sd = 5)
  alpha20 ~ dnorm(0, sd = 5)
  alpha21 ~ dnorm(0, sd = 5)
  
  # Priors for correlation regression
  gamma0 ~ dnorm(0, sd = 2)
  gamma1 ~ dnorm(0, sd = 2)
})
```

Key points:

* `Sigma` is a 3-D array `[N, 2, 2]` of deterministic nodes.
* `dmnorm` is called *per observation* with that observation’s slice `Sigma[i, 1:2, 1:2]`.
* The parameterization **guarantees** positive-definiteness:

  * `sigma1[i], sigma2[i] > 0` because of `exp(...)`,
  * `|rho[i]| < 1` because of `2 * ilogit(...) - 1`.

If this compiles and runs fine but is slow, it’s because you’re doing a Cholesky for each `i` each iteration. For only 2 dimensions it’s usually acceptable, but if `N` is huge there are more efficient formulations (e.g., reparameterizing via conditional normals instead of full `dmnorm`).

---

## 4. Practical suggestions and diagnostics

1. **Start simpler.**
   Before putting the predictor into *all* of `σ1`, `σ2` and `ρ`, you might:

   * first allow only the **variances** to depend on `x`,
   * keep `ρ` constant: `rho[i] <- rho; rho ~ dunif(-1, 1)`.
     This often already helps with weird bimodality.

2. **Check residual correlations by predictor level.**
   Fit a simpler model (constant covariance), compute residuals for each response, and look at sample correlations in low- vs high-`x` groups. That will tell you whether a predictor-dependent correlation is really needed.

3. **If you stick with predictor-dependent ρ**, monitor:

   * `gamma1`’s posterior (does it clearly differ from 0?),
   * `rho[i]` at representative values of `x` (say, min/median/max).

---

If you’re willing to share how you’ve currently specified the bivariate normal in `nimble` (mean structure + priors), I can help you translate it directly into a version with predictor-dependent covariance and also comment on why the traces look bimodal in *that* specific setup.

# Currently, I am using an LJK prior for the multivariate normal. How would I translate what you suggest above into this form?

Right, this is the key subtlety:

An LKJ prior is defined on **a single correlation matrix**.
As soon as you let the **correlation itself depend on a predictor**, you no longer have *one* correlation matrix, so you can’t literally keep a standard LKJ prior on “the” correlation matrix.

So there are really two different modeling choices:

---

## Option A (recommended): keep LKJ, let *variances* depend on the predictor

You can still use the LKJ prior exactly as you do now, but treat the *marginal SDs* as functions of your predictor (x), while the **correlation is global**. This often fixes the “bimodal ρ” problem once you give the model enough flexibility in the mean and variance.

Conceptually:

* Global correlation matrix (R) with LKJ prior.
* Observation-specific SDs (\sigma_{1i}(x_i), \sigma_{2i}(x_i)).
* Observation-specific covariance:
  [
  \Sigma_i = D_i, R, D_i,\quad
  D_i = \mathrm{diag}\big(\sigma_{1i},,\sigma_{2i}\big).
  ]

In `nimble` with LKJ on a 2×2 correlation matrix:

```r
code <- nimbleCode({
  ## LKJ prior on the (global) correlation matrix
  L_corr[1:2, 1:2] ~ dlkj_corr_cholesky(eta, K = 2)
  R[1:2, 1:2] <- chol2cor(L_corr[1:2, 1:2])
  rho <- R[1, 2]   # if you want to monitor it
  
  ## Priors for mean structure
  beta10 ~ dnorm(0, sd = 10)
  beta11 ~ dnorm(0, sd = 10)
  beta20 ~ dnorm(0, sd = 10)
  beta21 ~ dnorm(0, sd = 10)
  
  ## Priors for variance regressions
  alpha10 ~ dnorm(0, sd = 5)
  alpha11 ~ dnorm(0, sd = 5)
  alpha20 ~ dnorm(0, sd = 5)
  alpha21 ~ dnorm(0, sd = 5)
  
  ## Likelihood
  for(i in 1:N) {
    ## Means
    mu[i, 1] <- beta10 + beta11 * x[i]
    mu[i, 2] <- beta20 + beta21 * x[i]
    
    ## SDs as functions of x[i]
    log_sigma1[i] <- alpha10 + alpha11 * x[i]
    log_sigma2[i] <- alpha20 + alpha21 * x[i]
    sigma1[i] <- exp(log_sigma1[i])
    sigma2[i] <- exp(log_sigma2[i])
    
    ## Build covariance using global R and local SDs
    Sigma[i, 1, 1] <- pow(sigma1[i], 2)
    Sigma[i, 2, 2] <- pow(sigma2[i], 2)
    Sigma[i, 1, 2] <- R[1, 2] * sigma1[i] * sigma2[i]
    Sigma[i, 2, 1] <- Sigma[i, 1, 2]
    
    ## Bivariate normal
    y[i, 1:2] ~ dmnorm(mean = mu[i, 1:2],
                       cov  = Sigma[i, 1:2, 1:2])
  }
  
  ## Hyperprior for LKJ shape
  eta ~ dexp(1)   # or fix eta to something like 1 or 2
})
```

* Here you **retain the LKJ prior** on (R).
* The predictor enters the *marginal* variances, not the correlation.
* This is the cleanest “LKJ-compliant” way to make the covariance depend on predictors.

Given your trace behavior, I’d actually do this first; often what looks like “ρ changes with x” is partly heteroscedasticity.

---

## Option B: Let correlation depend on the predictor → LKJ no longer natural

If you want what I wrote earlier:

[
\begin{aligned}
\log \sigma_{1i} &= \alpha_{10} + \alpha_{11} x_i,\
\log \sigma_{2i} &= \alpha_{20} + \alpha_{21} x_i,\
\eta_{\rho,i} &= \gamma_0 + \gamma_1 x_i,\
\rho_i &= 2\operatorname{ilogit}(\eta_{\rho,i}) - 1,
\end{aligned}
]

then each observation has its own correlation ( \rho_i). There is no single correlation matrix to put an LKJ prior on; instead, your *primitive parameters* are the regression coefficients (\gamma_0, \gamma_1). You put priors directly on those, e.g.

```r
gamma0 ~ dnorm(0, sd = 2)
gamma1 ~ dnorm(0, sd = 2)
```

and then:

```r
for(i in 1:N) {
  ...
  eta_rho[i] <- gamma0 + gamma1 * x[i]
  rho[i] <- 2 * ilogit(eta_rho[i]) - 1
  
  Sigma[i,1,1] <- pow(sigma1[i], 2)
  Sigma[i,2,2] <- pow(sigma2[i], 2)
  Sigma[i,1,2] <- rho[i] * sigma1[i] * sigma2[i]
  Sigma[i,2,1] <- Sigma[i,1,2]
  
  y[i,1:2] ~ dmnorm(mu[i,1:2], cov = Sigma[i,1:2,1:2])
}
```

That gives you a **covariate-dependent correlation**, but it’s **not** LKJ anymore; there’s simply no single correlation matrix to which LKJ can be applied.

You *could* try to choose priors on (\gamma_0, \gamma_1) that make the induced marginal prior on (ρ_i) roughly match LKJ at “typical” (x), but that’s approximate and messy. For 2D LKJ(η) the marginal density of ρ is proportional to ((1 - ρ^2)^{η-1}); translating that into priors on a linear predictor is not straightforward.

---

### Summary

* **If you want to keep LKJ**:
  → Let the **marginal SDs** depend on (x), keep a single global correlation matrix (R) with LKJ prior (Option A). This is the direct translation of your current LKJ setup.

* **If you truly want ρ(x)**:
  → Drop LKJ, treat the regression on ρ as the primitive parameterization (Option B). LKJ, by construction, doesn’t extend cleanly to an (x)-varying correlation in the way nimble currently implements it.

If you paste your current LKJ-based `nimbleCode` I can show you exactly where to insert the SD regression (Option A) and what minimal edits would be needed to try Option B side-by-side.

# Nimble's distribution function does not use the "cov" argument. Nimble also does not have a chol2cor() function.

You’re absolutely right on both counts:

* `dmnorm()` in **nimble** uses a *precision* (`prec`) argument, not `cov`.
* There is no built-in `chol2cor()`.

My earlier sketch was more “generic MVN” than “nimble-correct.” Let me translate everything into actual nimble patterns.

---

## 1. Using LKJ in nimble + predictor in the *variances* (keep correlation global)

This is the cleanest thing you can do **while still using `dlkj_corr_cholesky`**.

### Key ideas

* `L_corr[ , ] ~ dlkj_corr_cholesky(eta, K)` gives you a **Cholesky factor of a correlation matrix** (R):
  (R = L_{\text{corr}} L_{\text{corr}}^\top).
* For each observation (i), define SDs as functions of a predictor (x_i):
  (\log \sigma_{1i} = \alpha_{10} + \alpha_{11} x_i), etc.
* The Cholesky factor of the **covariance** matrix for obs (i) is

  [
  L_{\Sigma,i} = D_i, L_{\text{corr}}, \quad
  D_i = \text{diag}(\sigma_{1i},, \sigma_{2i}),
  ]

  and you feed that to `dmnorm_cholesky()` with `prec_param = 0`.

Because we’re in 2D, we can write this analytically (no `%*%`, no `chol`, no `cov`):

Let

```text
L_corr = [ l11   0
           l21  l22 ]
```

Then

```text
L_Sigma[i] = D_i %*% L_corr
           = [ sigma1[i]*l11      0
               sigma2[i]*l21   sigma2[i]*l22 ]
```

Here’s a minimal nimble model:

```r
code <- nimbleCode({

  #### LKJ prior for the 2x2 correlation matrix ####
  L_corr[1:2, 1:2] ~ dlkj_corr_cholesky(eta, K = 2)
  eta ~ dexp(1)            # or fix eta to, say, 1 or 2

  # If you want to monitor rho, reconstruct R and extract R[1,2]:
  R[1:2, 1:2] <- L_corr[1:2, 1:2] %*% t(L_corr[1:2, 1:2])
  rho <- R[1, 2]

  #### Priors for mean structure ####
  beta10 ~ dnorm(0, sd = 10)
  beta11 ~ dnorm(0, sd = 10)
  beta20 ~ dnorm(0, sd = 10)
  beta21 ~ dnorm(0, sd = 10)

  #### Priors for variance regressions ####
  alpha10 ~ dnorm(0, sd = 5)
  alpha11 ~ dnorm(0, sd = 5)
  alpha20 ~ dnorm(0, sd = 5)
  alpha21 ~ dnorm(0, sd = 5)

  #### Likelihood ####
  for(i in 1:N) {

    # Means
    mu[i, 1] <- beta10 + beta11 * x[i]
    mu[i, 2] <- beta20 + beta21 * x[i]

    # SDs as functions of x[i]
    log_sigma1[i] <- alpha10 + alpha11 * x[i]
    log_sigma2[i] <- alpha20 + alpha21 * x[i]
    sigma1[i] <- exp(log_sigma1[i])
    sigma2[i] <- exp(log_sigma2[i])

    # Cholesky of covariance for obs i:
    # L_Sigma = diag(sigma) %*% L_corr, written elementwise
    L_Sigma[i, 1, 1] <- sigma1[i] * L_corr[1, 1]
    L_Sigma[i, 1, 2] <- 0
    L_Sigma[i, 2, 1] <- sigma2[i] * L_corr[2, 1]
    L_Sigma[i, 2, 2] <- sigma2[i] * L_corr[2, 2]

    # Bivariate normal with Cholesky of covariance
    y[i, 1:2] ~ dmnorm_cholesky(
      mean      = mu[i, 1:2],
      chol      = L_Sigma[i, 1:2, 1:2],
      prec_param = 0          # 0 => 'chol' is Cholesky of *covariance*
    )
  }
})
```

This:

* Keeps a **single global correlation structure** with LKJ prior.
* Allows heteroscedasticity via `sigma1[i], sigma2[i]` depending on `x[i]`.
* Never uses a `cov=` argument or `chol2cor()`.

For many real datasets, this is enough to cure the “ρ is weirdly bimodal” issue you saw, because a lot of what looks like “ρ changes across x” is actually changing marginal variances.

---

## 2. If you truly want ρ(x): analytic Cholesky, no LKJ

Once you want the **correlation itself** to depend on a predictor, there is *no single correlation matrix* anymore, so the usual LKJ construction doesn’t apply.

In 2D, the cleanest nimble-friendly parameterization is:

[
\Sigma_i =
\begin{pmatrix}
\sigma_{1i}^2 & \rho_i \sigma_{1i}\sigma_{2i} \
\rho_i \sigma_{1i}\sigma_{2i} & \sigma_{2i}^2
\end{pmatrix}
]

and you feed `dmnorm_cholesky()` an analytic Cholesky of Σ(_i). For a 2×2 covariance matrix, if

[
L_i =
\begin{pmatrix}
a_i & 0 \
b_i & c_i
\end{pmatrix},\quad \Sigma_i = L_i L_i^\top,
]

then

* (a_i^2 = \sigma_{1i}^2 \Rightarrow a_i = \sigma_{1i}),
* (a_i b_i = \rho_i \sigma_{1i}\sigma_{2i} \Rightarrow b_i = \rho_i \sigma_{2i}),
* (b_i^2 + c_i^2 = \sigma_{2i}^2 \Rightarrow c_i = \sigma_{2i} \sqrt{1-\rho_i^2}).

So you never need `cov` or `chol` at all.

You can combine that with a regression for ρ(x):

[
\eta_{\rho,i} = \gamma_0 + \gamma_1 x_i,\qquad
\rho_i = 2 ,\text{ilogit}(\eta_{\rho,i}) - 1
]

so (|\rho_i|<1) automatically.

In nimble:

```r
code <- nimbleCode({

  #### Priors for mean structure ####
  beta10 ~ dnorm(0, sd = 10)
  beta11 ~ dnorm(0, sd = 10)
  beta20 ~ dnorm(0, sd = 10)
  beta21 ~ dnorm(0, sd = 10)

  #### Priors for SD regressions ####
  alpha10 ~ dnorm(0, sd = 5)
  alpha11 ~ dnorm(0, sd = 5)
  alpha20 ~ dnorm(0, sd = 5)
  alpha21 ~ dnorm(0, sd = 5)

  #### Priors for correlation regression ####
  gamma0 ~ dnorm(0, sd = 2)
  gamma1 ~ dnorm(0, sd = 2)

  for(i in 1:N) {

    # Means
    mu[i, 1] <- beta10 + beta11 * x[i]
    mu[i, 2] <- beta20 + beta21 * x[i]

    # SDs
    log_sigma1[i] <- alpha10 + alpha11 * x[i]
    log_sigma2[i] <- alpha20 + alpha21 * x[i]
    sigma1[i] <- exp(log_sigma1[i])
    sigma2[i] <- exp(log_sigma2[i])

    # Correlation as function of x[i]
    eta_rho[i] <- gamma0 + gamma1 * x[i]
    rho[i] <- 2 * ilogit(eta_rho[i]) - 1   # in (-1, 1)

    # Analytic Cholesky of covariance
    L_Sigma[i, 1, 1] <- sigma1[i]
    L_Sigma[i, 1, 2] <- 0
    L_Sigma[i, 2, 1] <- rho[i] * sigma2[i]
    L_Sigma[i, 2, 2] <- sigma2[i] * sqrt(1 - pow(rho[i], 2))

    # Bivariate normal
    y[i, 1:2] ~ dmnorm_cholesky(
      mean      = mu[i, 1:2],
      chol      = L_Sigma[i, 1:2, 1:2],
      prec_param = 0
    )
  }
})
```

This is fully nimble-legal, no `cov=` and no `chol2cor()`. The trade-off:

* **Pros**: allows a *predictor-dependent* correlation; analytically simple in 2D.
* **Cons**: you’ve left the world of LKJ; priors on `gamma0, gamma1` now implicitly define the prior on ρ(x).

---

## 3. What I’d do next, given your situation

Given your bimodal ρ traces:

1. **First try**: LKJ + predictor in the SDs (Section 1).

   * Minimal change from what you’re already doing.
   * Often resolves odd behavior in ρ because you’re no longer forcing a single variance structure across all x.

2. If you still see clear evidence that ρ varies with x (e.g., residual correlations by x-strata), then move to **Section 2** and explicitly model ρ(x) with `dmnorm_cholesky` and analytic Cholesky, dropping the LKJ.

If you want, you can paste your current LKJ nimble code and I can do an in-place edit that (a) keeps LKJ but adds predictor-dependent SDs, and (b) shows the minimal change you’d make if you decide to go to ρ(x).


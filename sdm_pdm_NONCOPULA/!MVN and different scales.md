# I have measured 16 different traits on a set of plants. These can be on different scales so numerically can vary between one another by an order of magnitude or two. I am estimating the means of these traits using a multivariate normal with an LJK prior in nimble. I do have to transform the normal values using exp() to get each trait's mean. Do I need to be concerned that different traits can vary by an order of magnitude or more. For example, one trait can have values in the range of about ~0.1 to 2, but another in the range of 200-300 or so.

Short answer (well… *stat prof short*): **yes, you should be concerned—but not for the reason people usually say “scale your data.”** The issue is not *statistical validity*, it’s **geometry and numerics**, especially given `exp()` and an LKJ prior.

Let me unpack that carefully.

---

## 1. Conceptually: different scales are *fine*

From a modeling point of view, having trait means that differ by orders of magnitude is not inherently a problem.

A multivariate normal on a vector of latent means is perfectly happy if:

* one component implies a mean around 0.2
* another implies a mean around 250

There is no violation of probability theory or Bayesian logic here.

So: **this is not a “you must standardize or else your inference is wrong” situation.**

---

## 2. The real problem: `exp()` + MVN + LKJ = nasty curvature

Where things *do* get dangerous is here:

> “I transform the normal values using `exp()` to get each trait’s mean.”

That means your latent MVN is living on the **log scale**, and the actual mean is

[
\mu_j = \exp(\eta_j)
]

Now think about the geometry:

* A difference of **1 unit** in (\eta) multiplies (\mu) by **e ≈ 2.7**
* A difference of **3 units** multiplies (\mu) by **~20**
* A difference of **5 units** multiplies (\mu) by **~150**

So if one trait lives near (\eta \approx \log(0.5) \approx -0.7)
and another lives near (\eta \approx \log(250) \approx 5.5)

then:

* the posterior curvature in (\eta)-space is wildly different across dimensions
* the Hessian becomes badly scaled
* correlations become hard to estimate
* Laplace approximations become fragile
* MCMC mixing can suffer (though it’s more forgiving)

This is **not about the data scale**; it’s about the **induced scale on the latent MVN**.

---

## 3. Why the LKJ prior makes this more sensitive

The LKJ prior operates on the **correlation matrix**, but correlations are only meaningful relative to **comparable marginal variances**.

When:

* trait A needs (\eta) to vary by ±0.2
* trait B needs (\eta) to vary by ±2 or more

then:

* the implied covariance matrix has very uneven marginal scales
* the LKJ prior interacts with those scales indirectly via the Cholesky
* small numerical errors in large-scale dimensions bleed into correlation estimates

This is especially problematic if:

* some traits are tightly estimated
* others are weakly informed
* and you’re integrating over random effects with Laplace

---

## 4. What *actually* goes wrong in practice

Here’s what I’d expect to see if scale is biting you:

* Variance components for large-scale traits get underestimated
* Correlations involving those traits get shrunk toward zero
* Laplace estimates disagree with MCMC
* Optimizers complain or converge inconsistently
* Posterior correlations look suspiciously “too clean”

None of this means the model is wrong—just **poorly conditioned**.

---

## 5. Best practice: center *on the log scale*

The clean fix is **not** to rescale the observed data directly.

Instead:

### Step 1: Work with *centered log-means*

Define:
[
\eta_j = \log(\mu_j) - c_j
]

where (c_j) is a trait-specific offset, e.g.:

* (\log(\text{typical value}))
* (\log(\text{sample mean}))
* or even a rough prior guess

Then:
[
\mu_j = \exp(\eta_j + c_j)
]

Now:

* all (\eta_j) live near 0
* the MVN operates on comparable scales
* the LKJ prior behaves much better
* Laplace curvature is far more reasonable

This is *exactly* the same model statistically—it’s a reparameterization—but numerically it’s night and day.

---

## 6. Alternative: trait-specific scales in the MVN

Another valid approach (more complex, sometimes necessary):

* Give each trait its own scale parameter
* Model the covariance as:
  [
  \Sigma = D R D
  ]
  with (D = \text{diag}(\sigma_1, \dots, \sigma_{16}))

This is standard, but:

* you still want the latent variables roughly centered
* otherwise optimization remains painful

---

## 7. Do you *have* to fix this?

If you’re:

* only using MCMC
* mixing is good
* effective sample sizes are strong
* no divergences or pathologies

then you *might* get away with it.

But if you’re:

* using Laplace
* comparing Laplace vs MCMC
* fitting hierarchical structure
* estimating correlations you care about

then **yes, you should absolutely be concerned**—and the solution is reparameterization, not data scaling.

---

## Bottom line

* Different trait magnitudes are fine **in principle**
* `exp()` turns “order of magnitude differences” into **severe curvature differences**
* LKJ + Laplace magnifies those numerical issues
* Centering the latent log-means is the right fix

If you want, you can paste a sketch of your nimble model block and I can suggest a numerically safer parameterization that keeps the scientific interpretation identical.

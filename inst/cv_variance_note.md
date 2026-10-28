# The Asymptotic Variance of the Sample Coefficient of Variation: A Note on Equations 9 and 22 of Chattopadhyay and Kelley (2016)

This note derives the asymptotic variance of the sample coefficient of
variation by the delta method, shows how it relates to the expression
printed as Equation 9 (and estimated by Equation 22) of

> Chattopadhyay, B., & Kelley, K. (2016). Estimation of the coefficient of
> variation with minimum risk: A sequential method for minimizing sampling
> error and study cost. *Multivariate Behavioral Research, 51*(5), 627–648.
> [doi:10.1080/00273171.2016.1203279](https://doi.org/10.1080/00273171.2016.1203279)

and explains the `variance` argument of `seq_mr_cv()`. The two expressions
agree for a normal population and for any population with no excess
kurtosis; they differ otherwise. The default of `seq_mr_cv()`,
`variance = "corrected"`, uses the general expression; `variance =
"published"` uses Equation 22 of the article, so that the article's worked
example and tables can be reproduced.

The formulas are written as plain text, so that they read the same in any
viewer: `^` marks a power, `_` a subscript, and `sqrt()` a square root.

## Setup

Let `X_1, ..., X_n` be independent observations from a population with
mean `μ > 0`, variance `σ^2`, and third and fourth central moments `μ_3`
and `μ_4`. The coefficient of variation is `κ = σ / μ`, estimated by
`k_n = s_n / m_n`, where `m_n` is the sample mean and `s_n^2` the sample
variance.

By the central limit theorem for the sample mean and the sample variance,
`sqrt(n) (m_n - μ, s_n^2 - σ^2)` is asymptotically normal with mean zero
and covariance matrix

```
Σ = | σ^2    μ_3         |
    | μ_3    μ_4 - σ^4   |
```

## The Delta Method

Write `k_n = g(m_n, s_n^2)` with `g(a, b) = sqrt(b) / a`. The gradient of
`g` at `(μ, σ^2)` is

```
D = ( -σ / μ^2,  1 / (2 σ μ) )
```

so `sqrt(n) (k_n - κ)` is asymptotically normal with mean zero and
variance

```
ξ^2 = D' Σ D
    = (σ^2 / μ^4) σ^2
      + 2 (-σ / μ^2) (1 / (2 σ μ)) μ_3
      + (1 / (4 σ^2 μ^2)) (μ_4 - σ^4)
```

Collecting terms,

```
ξ^2 = μ_4 / (4 μ^2 σ^2) - σ^2 / (4 μ^2) - μ_3 / μ^3 + σ^4 / μ^4
```

This is Equation 30 of Kelley, Darku, and Chattopadhyay (2018), who note
that it is the asymptotic variance found by Albrecher, Ladoucette, and
Teugels (2010). In terms of the skewness `γ_1 = μ_3 / σ^3` and the
kurtosis `β_2 = μ_4 / σ^4` it reads

```
ξ^2 = κ^4 - γ_1 κ^3 + κ^2 (β_2 - 1) / 4
```

## Equation 9 of the 2016 Article

Equation 9 of Chattopadhyay and Kelley (2016, p. 631) gives

```
ξ_2016^2 = μ_4 / (4 μ^4) + σ^4 / (4 μ^4) - μ_3 / μ^3 + σ^2 / (2 μ^2)
         = κ^4 (β_2 + 1) / 4 - γ_1 κ^3 + κ^2 / 2
```

and Equation 22 (p. 634) is its sample version. Subtracting,

```
ξ^2 - ξ_2016^2 = (κ^2 / 4) (1 - κ^2) (β_2 - 3)
```

The two expressions are therefore equal when `β_2 = 3`, which includes
every normal population (both then equal `κ^4 + κ^2 / 2`), and when
`κ = 1`. For a population with heavier tails than the normal (`β_2 > 3`)
and `κ < 1`, Equation 9 is the smaller of the two, and a sequential rule
built on it stops with fewer observations than the number that minimizes
the risk.

## A Numerical Demonstration

For a gamma population with shape `a` (and scale 1), `μ = σ^2 = a`,
`μ_3 = 2 a`, and `μ_4 = 3 a^2 + 6 a`, so `κ^2 = 1 / a`,
`ξ^2 = κ^2 (1 + κ^2) / 2`, and
`ξ_2016^2 = ξ^2 - (3 / 2) κ^4 (1 - κ^2)`. With `a = 1.5625` (`κ = 0.8`)
these are 0.5248 and 0.3036. The code below computes both, and compares
them with `n` times the simulated mean squared error of `k_n` at
`n = 400`.

```r
library(SMSD)

set.seed(2016)
shape <- 1.5625                       # gamma population, kappa = 0.8
kappa <- 1 / sqrt(shape)
n <- 400

# Population moments of a gamma distribution with this shape (scale 1).
mu <- shape
sigma2 <- shape
mu3 <- 2 * shape
mu4 <- 3 * shape^2 + 6 * shape

# The general delta method variance, and Equation 9 of the 2016 article.
general <- mu4 / (4 * mu^2 * sigma2) - sigma2 / (4 * mu^2) -
  mu3 / mu^3 + sigma2^2 / mu^4
equation_9 <- mu4 / (4 * mu^4) + sigma2^2 / (4 * mu^4) -
  mu3 / mu^3 + sigma2 / (2 * mu^2)
c(general = general, equation_9 = equation_9)
#>    general equation_9
#>   0.524800   0.303616

# n times the mean squared error of the sample coefficient of variation.
k <- replicate(50000, {
  x <- rgamma(n, shape = shape)
  sd(x) / mean(x)
})
squared_error <- n * (k - kappa)^2
c(n_mse = mean(squared_error),
  se = sd(squared_error) / sqrt(length(squared_error)))
#>       n_mse          se
#> 0.521506927 0.003396951

# The two estimators on one large sample.
x <- rgamma(200000, shape = shape)
seq_mr_cv(x, A = 400, c = 1)$value[3]
#> [1] 0.5394619
seq_mr_cv(x, A = 400, c = 1, variance = "published")$value[3]
#> [1] 0.3119073
```

The simulated value, 0.522 with a standard error of 0.003, agrees with the
general expression (0.525) and not with Equation 9 (0.304). On a single
sample of 200,000 observations the two estimators return 0.539 and 0.312:
each is close to the population quantity it estimates.

## Consequences for the Sequential Procedure

The stopping rule (Equation 23 of the article) and the optimal fixed sample
size `n_c = (A / c)^(1/2) ξ` (Equation 15) are unchanged; only the value
plugged in for `ξ^2` differs.

* For a normal population the two choices lead to the same procedure, up to
  sampling error in the estimates.
* For the gamma population above, `n_c` is proportional to `ξ`, so a rule
  built on Equation 9 aims at `sqrt(0.3036 / 0.5248) = 0.76` of the sample
  size that minimizes the risk. In a simulation of the full procedure with
  `A` = 5,000,000 and `c` = 10, the design of the article's Table 1
  (p. 642), with 5,000 replications, the mean final sample size was 375.9
  with Equation 22 (the article's Table 1 reports 374.04) and 480.7 with
  the general estimator (standard errors 0.9 and 1.2), against a risk
  minimizing size of 512.2. The simulation is
  `seq_mr_cv_sensitivity(0.8, A = 5e6, c = 10, distribution = "gamma",
  variance = "published", G = 5000, seed = 2016)`, and the same call with
  `variance = "corrected"`.
* The `n_c` column of the article's Table 1 is computed from Equation 9.
  For its gamma row with `κ = 0.4`, for example, it is 174; the general
  expression gives 215.4.

## What the Package Does

`seq_mr_cv(x, A, c)` uses the general expression by default
(`variance = "corrected"`), through the estimator of Kelley, Darku, and
Chattopadhyay (2018, Equations 31 and 32), which is also the estimator
behind `seq_aipe_cv()`. `seq_mr_cv(x, A, c, variance = "published")` uses
Equation 22 of the 2016 article. With it the package reproduces the
article's worked example (p. 638: risk 5964.345 at `n = 18`, 6224.861 at
`n = 19`, and the rule first met at `n = 35` with risk 4891.284) and
agrees with `mr_cv()` in the 'DMAR' package and `mr.cv()` in the 'MBESS'
package, both of which use the published form. (Two cases, both of which
arise in small samples, are treated differently. When either estimate is
zero, negative, or not finite, `seq_mr_cv()` reports the stopping rule as
not met, where the rule as printed would stop. When Equation 22 is positive
but below `n^-3`, `seq_mr_cv()` uses `n^-3` in its place, as it does for
the general estimator, where `mr_cv()` uses Equation 22 itself; the
reported risk then differs, and so does the decision for those values of
`A / c` that put `n^2` between the two criteria. See `?seq_mr_cv`.) The
choice in use is recorded on the returned object as `attr(out, "variance")`.

## References

Albrecher, H., Ladoucette, S. A., & Teugels, J. L. (2010). Asymptotics of
the sample coefficient of variation and the sample dispersion. *Journal of
Statistical Planning and Inference, 140*, 358–368.

Chattopadhyay, B., & Kelley, K. (2016). Estimation of the coefficient of
variation with minimum risk: A sequential method for minimizing sampling
error and study cost. *Multivariate Behavioral Research, 51*(5), 627–648.
[doi:10.1080/00273171.2016.1203279](https://doi.org/10.1080/00273171.2016.1203279)

Kelley, K., Darku, F. B., & Chattopadhyay, B. (2018). Accuracy in parameter
estimation for a general class of effect sizes: A sequential approach.
*Psychological Methods, 23*(2), 226–243.
[doi:10.1037/met0000127](https://doi.org/10.1037/met0000127)

# Correction to Equation 26 of Kelley, Darku, and Chattopadhyay (2018)

This note reproduces the derivation of the asymptotic variance of the
standardized mean difference under the delta method, correcting a
typographical omission in Equation 26 of

> Kelley, K., Darku, F. B., & Chattopadhyay, B. (2018). Accuracy in
> parameter estimation for a general class of effect sizes: A sequential
> approach. *Psychological Methods, 23*(2), 226–243.
> [doi:10.1037/met0000127](https://doi.org/10.1037/met0000127)

The published Equation 26 (p. 232), and the sample version printed under
Equation 27, omit a factor of 1/2 in the first cross term. The corrected
form is what `seq_aipe_smd()` computes. The authors' typeset derivation is
posted at
<https://kenkelley.org/wp-content/uploads/articles/Kelley_Darku_Chattopadhyay_Equation_26.pdf>
and is installed with the package as `eq26_correction.pdf`
(`system.file("eq26_correction.pdf", package = "SMSD")`); the text below
reproduces it.

The correction matters only when the two means differ and the two groups
have different third central moments, because the corrected term is the
product of the mean difference and the difference of the third central
moments. For symmetric distributions, whenever the two groups are skewed
alike, and whenever the two means are equal, the two forms agree. As a
check, for two populations with unit variance, a standardized mean
difference of 0.5, and skewness of `sqrt(8/3)` in opposite directions (a
standardized chi square with 3 degrees of freedom and its mirror image),
the corrected variance is 1.3710 and the published form gives 0.5545; the
simulated value of `n` times the variance of the estimate at `n = 200` per
group is 1.36 (standard error 0.03, 4,000 replications). The package's
tests hold the corrected form to the matrix product `D' Σ D` below, to that
population value, and to that sampling variance.

The formulas are written as plain text, so that they read the same in any
viewer: `^` marks a power, `_` a subscript, and `sqrt()` a square root.

## Setup

For two independent groups with means `μ_1` and `μ_2`, a common variance
`σ^2`, and third and fourth central moments `μ_k3` and `μ_k4` for group
`k = 1, 2`, let

```
θ_1 = μ_1 - μ_2,    θ_2 = σ^2,    g_1(x) = x,    g_2(x) = sqrt(x)
```

Then `δ = g_1(θ_1) / g_2(θ_2) = (μ_1 - μ_2) / σ` is the population
standardized mean difference. The delta method gives

```
ξ^2 = D' Σ D
```

with

```
D = ( g_1'(θ_1) / g_2(θ_2),  -g_1(θ_1) g_2'(θ_2) / g_2(θ_2)^2 )
  = ( 1 / σ,  -(μ_1 - μ_2) / (2 σ^3) )
```

and

```
Σ = | 2 σ^2                 (μ_13 - μ_23) / 2            |
    | (μ_13 - μ_23) / 2     (μ_14 + μ_24 - 2 σ^4) / 4    |
```

## Expansion

Using `x' A x = a_11 x_1^2 + 2 a_12 x_1 x_2 + a_22 x_2^2`,

```
ξ^2 = D' Σ D
    = 2 σ^2 (1 / σ^2)
      + 2 ((μ_13 - μ_23) / 2) (1 / σ) (-(μ_1 - μ_2) / (2 σ^3))
      + ((μ_14 + μ_24 - 2 σ^4) / 4) ((μ_1 - μ_2)^2 / (4 σ^6))
```

The first term equals 2. The second term simplifies to

```
-(μ_1 - μ_2) (μ_13 - μ_23) / (2 σ^4)
```

The third term equals

```
((μ_1 - μ_2)^2 / (4 σ^6)) ((μ_14 + μ_24) / 4 - σ^4 / 2)
```

Defining

```
t_1 = (1 / 2) (μ_1 - μ_2) (μ_13 - μ_23) / σ^4
t_2 = (μ_1 - μ_2)^2 / (4 σ^6)
t_3 = (μ_14 + μ_24) / 4 - σ^4 / 2
```

the asymptotic variance has the three term decomposition

```
ξ^2 = 2 - t_1 + t_2 t_3
```

The factor of 1/2 in `t_1` is the correction; Equation 26 as published had
this factor missing.

## Sample Analog

For two samples of `n` observations each, `seq_aipe_smd()` substitutes

* the two sample means for `μ_1` and `μ_2`,
* the pooled variance `(s_1^2 + s_2^2) / 2` for `σ^2`,
* the unbiased estimates of the third central moment of each group for
  `μ_13` and `μ_23`, and
* the unbiased estimates of the fourth central moment of each group for
  `μ_14` and `μ_24`,

and uses the larger of the resulting estimate and `n^-3`, the lower bound
of the article's Equation 27. An estimate that is zero, negative, or not
finite is not a usable estimate of the variance: `seq_aipe_smd()` then
reports the stopping rule as not met (see its help page). A positive
estimate is used as the article prescribes, through the bound.

# Minimum Risk Point Estimation of the Standardized Mean Difference

Given the data collected so far from two independent groups, evaluates
the purely sequential minimum risk point estimation (MRPE) stopping rule
of Chattopadhyay and Kelley (2017) for the standardized mean difference
\\\delta = (\mu_1 - \mu_2)/\sigma\\. The procedure is for a study in
which each observation costs money and a less accurate estimate also has
a cost: it stops sampling when the estimated total of the two is as
small as it can be. It targets a risk, not a confidence interval, so
there is no confidence level.

## Usage

``` r
seq_mr_smd(group_1, group_2, A, c, gamma = 0.49, na.rm = TRUE)
```

## Arguments

- group_1, group_2:

  Numeric vectors: the observations collected so far in the first and
  the second group, at least four in each. Four is the smallest pilot
  size the article allows, and
  [`seq_mr_smd_pilot`](https://yelleknek.github.io/SMSD/reference/seq_mr_smd_pilot.md)
  never returns less.

- A:

  The loss constant: the price of one unit of squared estimation error,
  in the same currency as `c`. The article builds it as the amount one
  is willing to pay to estimate the standardized mean difference to
  within \\\epsilon\\, divided by \\\epsilon^2\\ (for example \\2500 /
  0.5^2 = 10000\\). A single positive number.

- c:

  The cost of sampling one observation: a single positive number when an
  observation costs the same in both groups, or two positive numbers,
  `c(c1, c2)`, the cost in group 1 and the cost in group 2, in that
  order. Either group may be the costlier one.

- gamma:

  The tuning constant of the correction term, strictly between 0 and
  1/2; default 0.49. Use the same value as in
  [`seq_mr_smd_pilot`](https://yelleknek.github.io/SMSD/reference/seq_mr_smd_pilot.md).

- na.rm:

  Should missing values be removed before the computation? `TRUE` (the
  default) removes them; `FALSE` stops with an error if any are present.

## Value

An object of class `c("smsd_check", "data.frame")` with columns `term`
and `value` and one row for each of

- `smd`:

  the standardized mean difference, group 1 minus group 2, divided by
  the pooled standard deviation;

- `risk`:

  the approximate risk (Equation 9, or Equation 17 with two costs) at
  the current group sizes, with the estimate in place of \\\delta\\;

- `mean_1`, `mean_2`:

  the two sample means;

- `optimal_n1`, `optimal_n2`:

  the criterion sample size for each group. With one cost both equal the
  right side of Equation 15, rounded up. With two costs, the less costly
  group's value is the right side of its stopping rule, rounded up, and
  the costlier group's value is the size the article pairs with the
  larger of the less costly group's current size and its criterion. The
  rule is met exactly when each current size is at least its criterion;

- `current_n1`, `current_n2`:

  the number of observations supplied in each group;

- `next_n1`, `next_n2`:

  the size each group should have at the next check. When the rule is
  met they equal the current sizes: nothing more is called for.
  Otherwise, with one cost, each group is brought up to one more than
  the smaller current size. With two costs the less costly group gets
  one more observation and the costlier group is brought up to the size
  the article pairs with the less costly group's new size; but if the
  costlier group is short of the size paired with the less costly
  group's current size, only that shortfall is called for. See "Sampling
  with two costs" in Details.

Attributes: `is_satisfied` (`TRUE` or `FALSE`, never `NA`: whether to
stop sampling), `current_n` (the total sample size), `criterion_n` (the
sum of the two criteria), `method`, `loss_A`, `cost_1`, `cost_2`,
`gamma`, `variance_unusable` (always `FALSE` here, because \\V_n^2\\ is
never below 2; carried for uniformity with the procedures whose variance
estimate can be zero, negative, or not finite), and `notes` (a character
vector, printed under the result). `conf_level` and `width` are `NULL`.
The rule is decided group by group, so the two totals are not to be
compared with each other: `current_n` can be at least `criterion_n`
while one group is still short of its own criterion. In that case a note
names the group. If every observation within each group is the same, or
the same to ten significant digits, there is no pooled standard
deviation to divide by and the standardized mean difference is not
estimable: `smd`, `risk`, the two criteria, and `criterion_n` are `NA`,
the rule is reported as not met, and a note says why (`next_n1` and
`next_n2` do not depend on the estimate and are still reported). That
test is relative to the size of the values, so observations that do vary
but share a large common offset (each group's range is at most
\\10^{-10}\\ times its largest absolute value, as it can be for times
recorded from a distant origin) are treated as the same too. The note
then says to subtract the same constant from both groups, which leaves
the standardized mean difference unchanged, and to check again. See
[`smsd_check_methods`](https://yelleknek.github.io/SMSD/reference/smsd_check_methods.md)
for the print, summary, and conversion methods.

## Details

The function does not collect data or keep state. Call
[`seq_mr_smd_pilot`](https://yelleknek.github.io/SMSD/reference/seq_mr_smd_pilot.md)
for the starting sample sizes, collect that many observations, and call
`seq_mr_smd()`. While `attr(out, "is_satisfied")` is `FALSE`, bring the
two groups up to the sizes in the rows `next_n1` and `next_n2` and call
it again.

The estimate is \\d = (\bar{X}\_1 - \bar{X}\_2)/s\\, with \\s\\ the
pooled standard deviation (Equations 2 and 3 of Chattopadhyay and
Kelley, 2017). The article assumes a normal population with a common
variance in each group; the stopping rule rests on the normal theory
mean squared error of \\d\\. Use
[`seq_mr_smd_sensitivity`](https://yelleknek.github.io/SMSD/reference/seq_mr_sensitivity.md)
to see what a skewed population costs. The estimate, and with it the
stopping decision, is the same whatever the common units of the two
groups (seconds or minutes, say), up to rounding error.

**One cost for both groups.** The approximate risk (their Equation 9) is
\$\$R = A\left\\\frac{1}{n_1} + \frac{1}{n_2} + \frac{\delta^2}{2(n_1 +
n_2)}\right\\ + c\\(n_1 + n_2),\$\$ the price of the squared error plus
the cost of sampling. It is smallest when each group has \\n_c =
\sqrt{A/(2c)}\\\xi\\ observations, with \\\xi^2 = 2 + \delta^2/4\\
(their Equations 10 and 11), where it equals \\4 c\\ n_c\\ (their
Equation 12). Because \\\delta\\ is unknown, the purely sequential rule
(their Equation 15) stops at the first \\n \ge m\\ observations per
group with \$\$n \ge \sqrt{\frac{A}{2c}}\left(V_n + n^{-\gamma}\right),
\qquad V_n^2 = 2 + \frac{d_n^2}{4},\$\$ where \\m\\ is the pilot size,
\\V_n^2\\ is their Equation 13, and \\n^{-\gamma}\\ with \\0 \< \gamma
\< 1/2\\ keeps the rule from stopping too early. In the article's
example (\\d\\ = 1.00, \\A\\ = 10,000, \\c\\ = 2.40; p. 106) the rule is
first met at \\n\\ = 75. The article states the rule for equal group
sizes. If the two groups supplied differ in size, this function uses the
smaller size as \\n\\, so both groups must reach the criterion.

**A different cost in each group.** This is the article's "A General
Scenario" (pp. 106 to 107, Equations 17 to 21, with Appendix B,
Equations 41 to 43). In this paragraph, subscript 1 marks the less
costly group, as in the article. The approximate risk (their Equation
17) replaces the sampling cost above with \\c_1 n_1 + c_2 n_2\\. Setting
its derivatives to zero (their Equation 41) gives the optimal sizes:
\\n_2^\* = n_1^\* \\ r(n_1^\*)\\ with \\r(n) = \\1 + (c_2 - c_1) n^2 /
A\\^{-1/2}\\ (their Equations 19 and 43), and \\n_1^\*\\ solves
\$\$n_1^\* = \sqrt{\frac{A}{c_1}} \left\\1 + \frac{\delta^2}{2\\\\1 +
r(n_1^\*)\\^2}\right\\^{1/2}.\$\$ The sequential procedure samples the
less costly group one observation (or batch) at a time, keeps the
costlier group at \\\lceil n_1\\ r(n_1) \rceil\\ observations (their
Equation 21 and the Stage II text), and stops at the first \\n_1\\ with
\$\$n_1 \ge \sqrt{\frac{A}{c_1}}\left\[ \left\\1 + \frac{d^2}{2\\\\1 +
r(n_1)\\^2}\right\\^{1/2} + n_1^{-\gamma}\right\],\$\$ where \\d\\ is
computed from all observations in both groups.

**Where the function departs from the article as printed.**

1.  Equations 18 and 42, and the stopping rule built on them, are
    printed with \\\\1 + r(n_1)\\\\ to the first power. The first order
    conditions in Equation 41 give the second power shown above. Only
    the second power minimizes the risk in Equation 17, and only it
    reduces to Equation 10 when the costs are equal (as Lemma 6 of the
    article states). The function uses the second power. For \\A\\ =
    50,000, costs of 2 and 8, and \\\delta\\ = 1, the printed form gives
    \\n_1^\*\\ = 183.4; direct minimization of Equation 17 gives 175.7,
    which is what the form above gives.

2.  The article's two cost procedure does not reduce to its one cost
    procedure when \\c_1 = c_2\\: its pilot (Equation 20) uses \\A/c_1\\
    where Equation 14 uses \\A/(2c)\\, and its correction term is
    \\\sqrt{A/c_1}\\n^{-\gamma}\\ where Equation 15 has
    \\\sqrt{A/(2c)}\\n^{-\gamma}\\. The function uses the one cost
    procedure whenever the two costs supplied are equal (so
    `c = c(3, 3)` and `c = 3` give identical results) and the two cost
    procedure whenever they differ. A small difference between the costs
    therefore moves the criterion by a visible amount; this is a
    property of the article's two procedures, not of the software. So
    that rounding error cannot make that choice, two costs whose
    relative difference is at most the square root of machine precision
    (about 1.5e-8) count as equal, and the smaller of the two is used:
    `c = c(24 * 0.1, 2.4)` gives the same result as `c = 2.4`.

3.  The article assumes the first group is the less costly one. The
    function applies the rule with whichever group is less costly in
    that role, so the order of the groups does not matter.

4.  The function is given data, not a sampling plan, so it also requires
    the costlier group to have at least \\\lceil n_1\\ r(n_1) \rceil\\
    observations before reporting the rule as met.

Everything else in the two cost procedure is the article's.

**Sampling with two costs.** The rows `next_n1` and `next_n2` give the
size each group should have at the next check, so the article's sampling
order does not have to be coded by hand: while the rule is not met,
collect observations until the groups have those sizes and call
`seq_mr_smd()` again. Starting from the pilot sizes, each such step adds
one observation to the less costly group and brings the costlier group
up to \\\lceil n_1\\ r(n_1) \rceil\\, with \\n_1\\ the less costly
group's new size, which is the article's procedure with one observation
per stage (its Stage II, p. 107, and the steps after it). The rows are
in the order of the arguments, group 1 and group 2, whichever of them is
the less costly one. If the costlier group has fewer observations than
the article pairs with the less costly group's current size, the rows
ask for that shortfall alone (the less costly group's row equals its
current size). That gives a way to check less often than after every
observation: add several observations to the less costly group and call
`seq_mr_smd()`; if the rows then ask only for the costlier group's
shortfall, collect it and call again. With one cost the rows bring each
group up to one more than the smaller of the two current sizes. Sampling
in another order (adding to the costlier group and topping up the less
costly one) oversamples: the function can then report the rule as met,
because both groups are past their criteria, at a sampling cost well
above the minimum.

**What to expect.** The article proves that the expected final sample
size and the risk approach their optimal values as \\A/c\\ grows; at
practical sizes the agreement is approximate, and the correction term
makes the procedure take somewhat more than the optimal number of
observations. With one cost, this function reproduces the article's
Tables 1 and 3 within Monte Carlo error. For \\A\\ = 50,000, \\c\\ = 1,
and \\\delta = 0\\, for example, the mean final sample size per group in
5,000 replications was 235.04 (the article reports 235.04; the optimal
size is 223.6), and the ratio of the realized risk to the minimum risk
was 1.00. With costs of 2 and 8, \\A\\ = 50,000, and \\\delta = 1\\, the
mean final sizes were 189.2 and 82.7 against optimal sizes of 175.7 and
81.0, and the risk ratio was 0.99 (Monte Carlo standard error 0.01). At
smaller sample sizes the excess is larger: for \\A\\ = 1,000,000, \\c\\
= 500, and \\\delta = 0\\ the mean final size was 50.0 against an
optimal 44.7, as in the article's Table 3.

## References

Chattopadhyay, B., & Kelley, K. (2017). Estimating the standardized mean
difference with minimum risk: Maximizing accuracy and minimizing cost
with sequential estimation. *Psychological Methods*, *22*(1), 94–113.
[doi:10.1037/met0000089](https://doi.org/10.1037/met0000089)

## See also

[`seq_mr_smd_pilot`](https://yelleknek.github.io/SMSD/reference/seq_mr_smd_pilot.md)
for the pilot sizes and
[`seq_mr_smd_sensitivity`](https://yelleknek.github.io/SMSD/reference/seq_mr_sensitivity.md)
to study the procedure before collecting data; the accuracy in parameter
estimation companion
[`seq_aipe_smd`](https://yelleknek.github.io/SMSD/reference/seq_aipe_smd.md);
in DMAR,
[`mr_smd`](https://yelleknek.github.io/DMAR/reference/mr_smd.html) (the
one cost procedure, from summary statistics), the estimator
[`smd`](https://yelleknek.github.io/DMAR/reference/smd.html), and the
fixed sample size planner
[`ss_aipe_smd`](https://yelleknek.github.io/DMAR/reference/ss_aipe_smd.html).

Other seq_mr_smd:
[`seq_mr_sensitivity`](https://yelleknek.github.io/SMSD/reference/seq_mr_sensitivity.md),
[`seq_mr_smd_pilot()`](https://yelleknek.github.io/SMSD/reference/seq_mr_smd_pilot.md)

## Author

Bhargab Chattopadhyay <bhargab@iiitvadodara.ac.in>, Ken Kelley
<kkelley@nd.edu>

## Examples

``` r
# The pilot data of Chattopadhyay and Kelley (2017, p. 105): reading
# scores of children with high and low exposure to same language
# subtitling. Assessing one child costs 2.40 dollars; A = 2500 / 0.5^2.
high_sls <- c(11, 7, 22, 13, 6, 9, 11, 16, 12, 17, 14, 8, 16)
low_sls  <- c(3, 6, 10, 8, 14, 5, 12, 10, 6, 8, 13, 5, 9)
seq_mr_smd_pilot(A = 10000, c = 2.40)
#> [1] 13
seq_mr_smd(high_sls, low_sls, A = 10000, c = 2.40)
#> Sequential check (minimum risk point estimation (standardized mean difference))
#> 
#>   Current n       : 26
#>   Criterion n     : 164
#> 
#>   smd               : 1.021484
#>   risk              : 1801.521
#>   mean_1            : 12.46154
#>   mean_2            : 8.384615
#>   optimal_n1        : 82
#>   optimal_n2        : 82
#>   current_n1        : 13
#>   current_n2        : 13
#>   next_n1           : 14
#>   next_n2           : 14
#> 
#>   Stopping rule met: FALSE

# Two costs: a visit to a child in the second group costs twice as much.
# While the rule is not met, bring each group up to the size in its
# next_n row and check again. The rows follow the article's sampling
# order whichever group is the costlier one, so the loop also runs as
# written with cost <- c(6, 3).
set.seed(113)
A <- 2000
cost <- c(3, 6)                  # group 1, group 2
true_mean <- c(0.5, 0)
m <- seq_mr_smd_pilot(A = A, c = cost)
y <- list(rnorm(m[["group_1"]], mean = true_mean[1]),
          rnorm(m[["group_2"]], mean = true_mean[2]))
check <- seq_mr_smd(y[[1]], y[[2]], A = A, c = cost)
while (!attr(check, "is_satisfied")) {
  next_n <- check$value[check$term %in% c("next_n1", "next_n2")]
  for (j in 1:2) {
    y[[j]] <- c(y[[j]],
                rnorm(next_n[j] - length(y[[j]]), mean = true_mean[j]))
  }
  check <- seq_mr_smd(y[[1]], y[[2]], A = A, c = cost)
}
check
#> Sequential check (minimum risk point estimation (standardized mean difference))
#> 
#>   Current n       : 51
#>   Criterion n     : 51
#> 
#>   smd               : 0.2240416
#>   risk              : 378.5003
#>   mean_1            : 0.5786443
#>   mean_2            : 0.3443293
#>   optimal_n1        : 31
#>   optimal_n2        : 20
#>   current_n1        : 31
#>   current_n2        : 20
#>   next_n1           : 31
#>   next_n2           : 20
#> 
#>   Stopping rule met: TRUE
```

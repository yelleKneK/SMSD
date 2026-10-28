# Pilot Sample Sizes for Minimum Risk Point Estimation of the Standardized Mean Difference

Returns the pilot (first stage) sample size for each group in the purely
sequential minimum risk point estimation (MRPE) procedure of
Chattopadhyay and Kelley (2017) for the standardized mean difference. It
needs only the cost constants, so it can be computed before any data are
collected.

## Usage

``` r
seq_mr_smd_pilot(A, c, m_0 = 4, gamma = 0.49)
```

## Arguments

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

- m_0:

  The smallest pilot size allowed per group: a whole number of at least
  4 (the default), as in the article.

- gamma:

  The tuning constant of the stopping rule's correction term, strictly
  between 0 and 1/2. The default, 0.49, is the value Chattopadhyay and
  Kelley (2017) suggest. Use the same value in
  [`seq_mr_smd`](https://yelleknek.github.io/SMSD/reference/seq_mr_smd.md).

## Value

When `c` is a single number, a length one integer: the pilot size for
each group. When `c` has two elements, an integer vector with names
`group_1` and `group_2`.

## Details

With one cost for both groups the pilot size per group is Equation 14 of
Chattopadhyay and Kelley (2017), \$\$m = \max\left\\m_0,\\ \left\lceil
(A/(2c))^{1/(2 + 2\gamma)} \right\rceil\right\\.\$\$ The article's
examples give 13 for \\A\\ = 10,000 and \\c\\ = 2.40 (p. 105), and 30
for \\A\\ = 50,000 and \\c\\ = 1 (p. 102).

With two different costs the pilot sizes are Equations 20 and 21 of the
article ("A General Scenario", p. 106). Writing \\c_1\\ for the smaller
cost and \\c_2\\ for the larger, the less costly group starts with
\$\$m\_{1} = \max\left\\m_0,\\ \left\lceil (A/c_1)^{1/(2 + 2\gamma)}
\right\rceil\right\\\$\$ observations and the costlier group with
\$\$m\_{2} = \max\left\\m_0,\\ \left\lceil m\_{1} \left(1 + \frac{(c_2 -
c_1)\\m\_{1}^2}{A}\right)^{-1/2} \right\rceil\right\\.\$\$ For \\A\\ =
50,000 and costs of 2 and 8 these are 30 and 29. The article assumes the
first group is the less costly one; this function accepts the costs in
either order and returns the sizes in the order of `c`.

Equation 20 uses \\A/c_1\\ where Equation 14 uses \\A/(2c)\\, so the
article's two cost pilot does not reduce to its one cost pilot when the
costs are equal. When the two costs supplied are equal this function
returns the one cost (Equation 14) size for both groups. Two costs count
as equal when their relative difference is at most the square root of
machine precision (about 1.5e-8), so that costs computed by arithmetic
(`24 * 0.1` against `2.4`) are not sent to the two cost procedure by
rounding error; the smaller of the two is then the cost used.

## References

Chattopadhyay, B., & Kelley, K. (2017). Estimating the standardized mean
difference with minimum risk: Maximizing accuracy and minimizing cost
with sequential estimation. *Psychological Methods*, *22*(1), 94–113.
[doi:10.1037/met0000089](https://doi.org/10.1037/met0000089)

## See also

[`seq_mr_smd`](https://yelleknek.github.io/SMSD/reference/seq_mr_smd.md)
for the stopping rule and
[`seq_mr_smd_sensitivity`](https://yelleknek.github.io/SMSD/reference/seq_mr_sensitivity.md)
to study the procedure before collecting data; in DMAR,
[`mr_smd`](https://yelleknek.github.io/DMAR/reference/mr_smd.html) (the
one cost procedure) and the fixed sample size planner
[`ss_aipe_smd`](https://yelleknek.github.io/DMAR/reference/ss_aipe_smd.html).

Other seq_mr_smd:
[`seq_mr_sensitivity`](https://yelleknek.github.io/SMSD/reference/seq_mr_sensitivity.md),
[`seq_mr_smd()`](https://yelleknek.github.io/SMSD/reference/seq_mr_smd.md)

## Author

Bhargab Chattopadhyay <bhargab@iiitvadodara.ac.in>, Ken Kelley
<kkelley@nd.edu>

## Examples

``` r
# The worked example of Chattopadhyay and Kelley (2017, p. 105): assessing
# one student costs 2.40 dollars, and the team would pay 2,500 dollars to
# estimate the standardized mean difference to within 0.5.
seq_mr_smd_pilot(A = 2500 / 0.5^2, c = 2.40)
#> [1] 13

# Two costs: an observation in the second group costs four times as much.
seq_mr_smd_pilot(A = 50000, c = c(2, 8))
#> group_1 group_2 
#>      30      29 
```

# Plot the Stopping Path of a Sequential Procedure

Draws the path returned by
[`seq_path`](https://yelleknek.github.io/SMSD/reference/seq_path.md)
with base graphics. The `"estimate"` panel shows the estimate, with its
confidence interval when the procedure has one, against the sample size,
so the interval can be seen narrowing toward its target. The
`"criterion"` panel shows the current sample size against the criterion
sample size, so the stopping rule can be seen being met. The first
sample size at which the rule is met is marked in both.

## Usage

``` r
# S3 method for class 'smsd_path'
plot(
  x,
  which = c("estimate", "criterion"),
  main = NULL,
  xlab = NULL,
  ylab = NULL,
  legend = TRUE,
  ...
)
```

## Arguments

- x:

  An object of class `"smsd_path"`.

- which:

  Which panels to draw: `"estimate"`, `"criterion"`, or both (the
  default; the `"criterion"` panel is then left out for a procedure that
  reports no criterion sample size).

- main:

  The title, placed over the first panel drawn, above a line that gives
  the target (the target width or, for a procedure that only bounds the
  length of its interval, the largest length allowed) and the sample
  size at which the rule was first met: one string or an expression. The
  default is built from the method. Use `main = ""` for neither.

- xlab:

  The label of the horizontal axis: one string or an expression. Default
  `"Sample size (n)"`.

- ylab:

  The labels of the vertical axes: one string or expression, or one for
  each panel drawn (a character vector, an expression vector, or a
  list). The defaults are the name of the quantity estimated and
  `"Sample size"`.

- legend:

  Logical; draw a legend in each panel? Default `TRUE`.

- ...:

  Further graphical parameters, applied to every panel drawn: `col`,
  `lwd`, and `lty` set the estimate line (and the current sample size
  line); anything else (for example `xlim`, `ylim`, `log`, `las`,
  `cex.axis`) is passed to
  [`plot.default`](https://rdrr.io/r/graphics/plot.default.html).
  Request a single panel with `which` when a parameter such as `ylim`
  should apply to that panel only. Room for the legend is left above the
  data unless `ylim` or a logarithmic vertical axis is asked for; the
  legend can then lie over the data, and `legend = FALSE` removes it.

## Value

`x`, invisibly. Called for its side effect of drawing on the current
graphics device.

## Details

In the `"estimate"` panel the stages after the stopping point are drawn
in gray. They show what the rule would have said had sampling continued
and are not part of the result: the estimate and interval a sequential
design reports are the ones at the marked sample size (see
[`seq_path`](https://yelleknek.github.io/SMSD/reference/seq_path.md)).

The `"criterion"` panel plots the columns `current_n` and `criterion_n`
of the path, which are on the scale the procedure itself uses (for
example the total over both groups for
[`seq_mr_smd`](https://yelleknek.github.io/SMSD/reference/seq_mr_smd.md)).
For most procedures the rule is met when the current sample size reaches
the criterion; where a procedure adds a further condition, the marked
sample size follows the procedure's own decision. In particular, a stage
at which the variance estimate was zero, negative, or not finite is
never the stopping point, whatever the two lines show there (the column
`variance_unusable` of the path). A procedure that reports no criterion
sample size has no `"criterion"` panel.

When both panels are drawn they are stacked in one column, and the
graphical parameters changed to do so are restored on exit. When a
single panel is requested no graphical parameter is changed, so the
panel can be placed in a layout of the user's own.

## References

Kelley, K., Darku, F. B., & Chattopadhyay, B. (2018). Accuracy in
parameter estimation for a general class of effect sizes: A sequential
approach. *Psychological Methods*, *23*(2), 226–243.
[doi:10.1037/met0000127](https://doi.org/10.1037/met0000127)

## See also

[`seq_path`](https://yelleknek.github.io/SMSD/reference/seq_path.md),
[`smsd_path_methods`](https://yelleknek.github.io/SMSD/reference/smsd_path_methods.md)

Other seq_path:
[`seq_path()`](https://yelleknek.github.io/SMSD/reference/seq_path.md),
[`smsd_path_methods`](https://yelleknek.github.io/SMSD/reference/smsd_path_methods.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
set.seed(113)
# Satisfaction (1 to 7 scale) and weekly hours of product use for 150
# customers, in the order surveyed. How does the interval for their
# correlation narrow toward a full width of 0.30?
use <- rnorm(150, mean = 10, sd = 3)
satisfaction <- 4 + 0.15 * (use - 10) + rnorm(150, sd = 0.8)
path <- seq_path(seq_aipe_r, list(use, satisfaction), width = 0.30)
plot(path)


# One panel, with the vertical axis and its label chosen by the user.
plot(path, which = "estimate", ylim = c(0, 1.2), las = 1,
     ylab = expression(hat(rho)))

```

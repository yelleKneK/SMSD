# Estimator of the asymptotic variance xi^2 of a linear contrast of K
# independent group means, psi_n = sum_k c_k xbar_kn, with a common per
# group sample size n: sqrt(n) (psi_n - psi) is asymptotically N(0, xi^2),
# with xi^2 = sum_k c_k^2 sigma_k^2 (the contrast is Example 6 of Kelley,
# Darku, and Chattopadhyay, 2018, Equation 16; with g_2 = 1 their Theorem 1
# reduces to the variance of the weighted sum of means). The estimator
# substitutes the sample variances. Used by seq_aipe_c(). Internal; not
# exported.
#
# groups: a list of numeric vectors, one per group (a matrix or data frame
#   with one column per group is also accepted).
# contrast_weights: one weight per group.
# na.rm: remove missing values from each group first? Default FALSE.
# floor: apply the lower bound n^-3? Default TRUE.
# Returns a single number: max{sum_k c_k^2 s_k^2, n^-3}, with n the smallest
# group size (the common size, for the balanced data seq_aipe_c() requires),
# or, with floor = FALSE, sum_k c_k^2 s_k^2 itself.
#
# The raw value is NaN when a group with a nonzero weight has no variation
# (a group whose values agree to ten significant digits counts as having
# none; see .is_constant()): the variance of that group has not been
# estimated, so neither has the variance of the contrast, and seq_aipe_c()
# reports a variance estimate that is not yet usable. This is the reading
# seq_budget() gives the same data. Summing the other groups' terms instead
# would leave the group out of the variance, and the rule could be met with
# an interval built on the other groups alone (the test "a group with no
# variation makes the variance of the contrast not yet usable" in
# tests/testthat/test-seq_aipe_c.R works an example by hand; the help of
# seq_aipe_c() gives the coverage this costs and the observations it
# takes). A group whose weight is zero takes no part in the contrast and is
# not looked at.
#
# Unlike the other effect sizes in the family a contrast is not scale free,
# so the bound n^-3 is on the scale of the squared data: a raw value that is
# positive but below it is replaced by the bound, as for the other effect
# sizes, and seq_aipe_c() adds a note (see its help).
#
# The sample variances are computed from the groups divided by a common
# power of two of their magnitude and returned to the units of the data
# (see .aipe_scale()), so that the value is right whenever it can be
# represented: the sum of squares inside a variance can overflow when the
# variance itself does not. The scale is applied twice rather than squared
# for the same reason.
#
# Authors: Bhargab Chattopadhyay, Ken Kelley.
#' @noRd
.var_c <- function(groups, contrast_weights, na.rm = FALSE, floor = TRUE) {
  if (is.data.frame(groups)) {
    groups <- as.list(groups)
  } else if (is.matrix(groups)) {
    groups <- lapply(seq_len(ncol(groups)), function(k) groups[, k])
  }
  if (na.rm) {
    groups <- lapply(groups, function(g) g[!is.na(g)])
  }
  n <- min(lengths(groups))
  flat <- vapply(groups, .is_constant, logical(1))
  if (any(flat & contrast_weights != 0)) {
    return(.var_floor(NaN, n, floor))
  }
  scale <- .aipe_scale(unlist(groups, use.names = FALSE))
  s2 <- vapply(seq_along(groups),
               function(k) {
                 if (flat[k]) 0 else stats::var(groups[[k]] / scale)
               },
               numeric(1))
  raw <- sum(contrast_weights^2 * s2) * scale * scale
  .var_floor(raw, n, floor)
}

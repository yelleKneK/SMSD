# Internal machinery for the standardized mean difference minimum risk point
# estimation procedure of Chattopadhyay and Kelley (2017). Not exported.
#
# With two costs the article labels the less costly group "Group 1"
# ("c1 and c2 (>= c1)", p. 106). Here that group is called "lo" and the
# costlier one "hi", so nothing depends on the order in which the user
# supplies the groups.

# The cost argument of the SMD procedures: one positive number (the same
# cost in both groups) or two (one per group). Returns c(c1, c2) as doubles.
# A one dimensional numeric array counts as a vector (.drop_one_dim()); a
# number with a class is refused, as by the scalar checks.
#
# Whether the two costs are equal decides which of the article's two
# procedures runs, and the two do not meet at equal costs. That decision is
# made once, here, with a tolerance: two costs whose relative difference is
# at most sqrt(.Machine$double.eps) (about 1.5e-8) are one cost, the smaller
# of the two. Costs computed by arithmetic (24 * 0.1 against 2.4) then give
# the same answer as the number typed, and every later `c1 == c2` in this
# file compares identical values.
#' @noRd
.smd_mr_cost <- function(c) {
  c <- .drop_one_dim(c)
  ok_shape <- is.numeric(c) && is.null(dim(c)) && length(c) %in% c(1L, 2L)
  if (!ok_shape) {
    .smsd_stop("'c' must be one positive number (the same cost per ",
               "observation in both groups) or two positive numbers, ",
               "c(c1, c2) (the cost per observation in group 1 and in ",
               "group 2); got ", .smsd_describe(c), ".")
  }
  .stop_if_classed(c, "c", "one or two plain numbers")
  if (length(c) == 1L) {
    c1 <- .check_positive(c, "c")
    return(c(c1, c1))
  }
  c1 <- .check_positive(c[1L], "c[1]")
  c2 <- .check_positive(c[2L], "c[2]")
  if (abs(c1 - c2) <= sqrt(.Machine$double.eps) * max(c1, c2)) {
    c1 <- c2 <- min(c1, c2)
  }
  c(c1, c2)
}

# The approximate risk of Chattopadhyay and Kelley (2017), Equation 9 (one
# cost) and Equation 17 (two costs), at group sizes n1 and n2. `delta2` is
# the squared standardized mean difference (or its estimate).
#' @noRd
.smd_mr_risk <- function(n1, n2, delta2, A, c1, c2) {
  A * (1 / n1 + 1 / n2 + delta2 / (2 * (n1 + n2))) + c1 * n1 + c2 * n2
}

# The ratio n_hi / n_lo that the first order conditions of the risk imply
# when the less costly group has n observations (Equations 19 and 43; the
# same expression sets the costlier group's size in Equation 21).
#' @noRd
.smd_mr_ratio <- function(n, A, c_lo, c_hi) {
  1 / sqrt(1 + (c_hi - c_lo) * n^2 / A)
}

# The costlier group's sample size that the article pairs with n
# observations in the less costly group (Equation 21 and the Stage II text
# on p. 107). It never exceeds n, so it is representable whenever n is.
# Returned as a double; callers that need an integer convert it.
#' @noRd
.smd_mr_companion <- function(n, A, c_lo, c_hi) {
  ceiling(n * .smd_mr_ratio(n, A, c_lo, c_hi))
}

# Right side of the fixed point equation for the less costly group's optimal
# size, n = sqrt(A / c_lo) * sqrt(1 + delta^2 / (2 (1 + r(n))^2)), with
# r(n) = .smd_mr_ratio(n). This is Equations 18 and 42 of the article with
# (1 + r) squared: the first order conditions (Equation 41) give the square,
# and only the squared form minimizes the risk in Equation 17 and reduces to
# Equation 10 when the costs are equal. The article prints the first power.
#' @noRd
.smd_mr_rhs <- function(n, delta2, A, c_lo, c_hi) {
  r <- .smd_mr_ratio(n, A, c_lo, c_hi)
  sqrt(A / c_lo) * sqrt(1 + delta2 / (2 * (1 + r)^2))
}

# The fixed sample sizes (not rounded) that minimize the approximate risk
# for a known squared standardized mean difference, and the risk there.
# Returns c(n1, n2, risk) in the order of the cost arguments.
#' @noRd
.smd_mr_optimum <- function(delta2, A, c1, c2) {
  if (c1 == c2) {
    n <- sqrt(A / (2 * c1)) * sqrt(2 + delta2 / 4)
    return(c(n, n, .smd_mr_risk(n, n, delta2, A, c1, c2)))
  }
  c_lo <- min(c1, c2)
  c_hi <- max(c1, c2)
  # (1 + r)^2 lies in (1, 4), which brackets the root. For a squared effect
  # below about 1e-15 the two ends of the bracket are the same double, and
  # the root is that number.
  lower <- sqrt(A / c_lo) * sqrt(1 + delta2 / 8)
  upper <- sqrt(A / c_lo) * sqrt(1 + delta2 / 2)
  n_lo <- if (!(upper > lower)) {
    lower
  } else {
    stats::uniroot(function(n) n - .smd_mr_rhs(n, delta2, A, c_lo, c_hi),
                   lower = lower, upper = upper,
                   tol = lower * .Machine$double.eps^0.75)$root
  }
  n_hi <- n_lo * .smd_mr_ratio(n_lo, A, c_lo, c_hi)
  n <- if (c1 < c2) c(n_lo, n_hi) else c(n_hi, n_lo)
  c(n, .smd_mr_risk(n[1L], n[2L], delta2, A, c1, c2))
}

# The sample means and the standardized mean difference (the mean difference
# over the pooled standard deviation). `smd` is NA when there is no pooled
# standard deviation to divide by: both groups are constant, or constant up
# to rounding error by the package wide test .is_constant() (a range of at
# most 1e-10 times the largest absolute value), or the pooled standard
# deviation is zero in floating point although one group varies (its
# values are so much smaller than the other group's that their squares
# underflow on the common scale). A zero or rounding-error pooled standard
# deviation therefore never yields an infinite or meaningless estimate.
# `within_tolerance` is TRUE when the groups were called constant by the
# tolerance alone (some observation differs from another within a group);
# the caller then attaches a note of its own, because data that vary and
# share a large common offset are in that case too.
#
# The standardized mean difference does not depend on the common units of
# the two groups, which are divided by one power of two when they are on an
# extreme scale (see .mrpe_scale()). The sums of squares can then neither
# overflow (computed on the data as given, they do for data of magnitude
# 1e154 and beyond, which would leave no estimate) nor underflow (at 1e-160
# the estimate would be wrong in its fifth digit with no sign of trouble,
# and below that there would be none). The means are returned in the units
# of the data.
#' @noRd
.smd_mr_estimate <- function(group_1, group_2) {
  n1 <- length(group_1)
  n2 <- length(group_2)
  r1 <- c(min(group_1), max(group_1))
  r2 <- c(min(group_2), max(group_2))
  scale <- .mrpe_scale(c(r1, r2))
  if (scale != 1) {
    group_1 <- group_1 / scale
    group_2 <- group_2 / scale
    r1 <- r1 / scale
    r2 <- r2 / scale
  }
  mean_1 <- mean(group_1)
  mean_2 <- mean(group_2)
  # .is_constant() needs only the range.
  constant <- .is_constant(r1) && .is_constant(r2)
  s_p <- sqrt(((n1 - 1) * stats::var(group_1) +
                 (n2 - 1) * stats::var(group_2)) / (n1 + n2 - 2))
  estimable <- !constant && s_p > 0
  list(smd = if (estimable) (mean_1 - mean_2) / s_p else NA_real_,
       mean_1 = mean_1 * scale, mean_2 = mean_2 * scale,
       within_tolerance = constant &&
         (r1[1L] < r1[2L] || r2[1L] < r2[2L]))
}

# Pilot sizes, c(m1, m2), in the order of the cost arguments.
#' @noRd
.smd_mr_pilot <- function(A, c1, c2, gamma, m_0) {
  if (c1 == c2) {
    # Equation 14.
    m <- max(m_0, .as_sample_size((A / (2 * c1))^(1 / (2 + 2 * gamma)),
                                  "pilot sample size",
                                  advice = "Check 'A' and 'c'."))
    return(c(m, m))
  }
  c_lo <- min(c1, c2)
  c_hi <- max(c1, c2)
  # Equations 20 and 21.
  m_lo <- max(m_0, .as_sample_size((A / c_lo)^(1 / (2 + 2 * gamma)),
                                   "pilot sample size",
                                   advice = "Check 'A' and 'c'."))
  m_hi <- max(m_0, as.integer(.smd_mr_companion(m_lo, A, c_lo, c_hi)))
  if (c1 < c2) c(m_lo, m_hi) else c(m_hi, m_lo)
}

# One evaluation of the stopping rule from the current estimate `d` and the
# current group sizes. Shared by seq_mr_smd() and seq_mr_smd_sensitivity(),
# so the Monte Carlo wrapper runs exactly the rule the check function
# reports. Returns list(criterion = c(group 1, group 2), satisfied). The
# criteria are on the scale of n and are not converted to integers here (the
# less costly group's can be very large), so this function cannot fail;
# seq_mr_smd() rounds them up for its optimal_n rows, with the rounding
# .mrpe_ceiling() applies here to decide whether a criterion is reached.
#' @noRd
.smd_mr_rule <- function(d, n1, n2, A, c1, c2, gamma) {
  if (c1 == c2) {
    # Equation 15, with V_n = sqrt(2 + d^2 / 4) from Equation 13. The article
    # states the rule for n1 = n2 = n; with unequal group sizes the smaller
    # one is used, so both groups must reach the criterion.
    n <- min(n1, n2)
    criterion <- sqrt(A / (2 * c1)) * (sqrt(2 + d^2 / 4) + n^(-gamma))
    return(list(criterion = c(criterion, criterion),
                satisfied = isTRUE(n >= .mrpe_ceiling(criterion))))
  }
  lo <- if (c1 < c2) 1L else 2L
  hi <- 3L - lo
  n_now <- c(n1, n2)
  cost <- c(c1, c2)
  n_lo <- n_now[lo]
  # The stopping rule on pp. 106 to 107, with the squared term (see
  # .smd_mr_rhs) and the article's correction term sqrt(A / c_lo) n^(-gamma).
  criterion_lo <- .smd_mr_rhs(n_lo, d^2, A, cost[lo], cost[hi]) +
    sqrt(A / cost[lo]) * n_lo^(-gamma)
  target_lo <- .mrpe_ceiling(criterion_lo)
  # "N2 computed using the value of N1" (p. 107): the costlier group must
  # have the size the article pairs with the less costly group's current
  # size.
  needed_hi <- .smd_mr_companion(n_lo, A, cost[lo], cost[hi])
  criterion <- numeric(2L)
  criterion[lo] <- criterion_lo
  # The size paired with the larger of the less costly group's current size
  # and its criterion (a criterion that is not finite is passed on as it
  # is).
  criterion[hi] <- if (is.finite(target_lo)) {
    .smd_mr_companion(max(n_lo, target_lo), A, cost[lo], cost[hi])
  } else {
    target_lo
  }
  list(criterion = criterion,
       satisfied = isTRUE(n_lo >= target_lo) && n_now[hi] >= needed_hi)
}

# The group sizes, c(group 1, group 2), that the procedure calls for at the
# next check when the rule is not met at sizes n1 and n2. `step` is the
# number of observations added between checks (1 in the purely sequential
# procedure of the article; the batch size in seq_mr_smd_sensitivity()).
#
# One cost (Equation 15, stated for equal group sizes): each group is
# brought up to `step` more than the smaller current size.
#
# Two costs (Stage II and the steps that follow it, p. 107): the less costly
# group gets `step` more observations and the costlier group is brought up
# to the size the article pairs with the less costly group's new size
# (.smd_mr_companion). That presupposes groups that are paired now, as they
# are from the pilot on when the article's order is followed. If the
# costlier group is short of the size paired with the less costly group's
# current size, the stage in hand is completed first: the costlier group is
# brought up to that size and the less costly group stays as it is.
#
# The sizes do not depend on the estimate, so they exist when the estimate
# does not. Returned as doubles. Shared by seq_mr_smd(), for its next_n
# rows, and by seq_mr_smd_sensitivity(), so a simulated replication samples
# exactly as the rows tell a user to.
#' @noRd
.smd_mr_next <- function(n1, n2, A, c1, c2, step = 1) {
  n_now <- as.numeric(c(n1, n2))
  if (c1 == c2) {
    return(pmax(n_now, min(n_now) + step))
  }
  lo <- if (c1 < c2) 1L else 2L
  hi <- 3L - lo
  cost <- c(c1, c2)
  n_next <- n_now
  paired_now <- .smd_mr_companion(n_now[lo], A, cost[lo], cost[hi])
  if (n_now[hi] < paired_now) {
    n_next[hi] <- paired_now
  } else {
    n_next[lo] <- n_now[lo] + step
    n_next[hi] <- max(n_now[hi],
                      .smd_mr_companion(n_next[lo], A, cost[lo], cost[hi]))
  }
  n_next
}

# The note for a check that is not met although the total sample size is at
# least the total of the two criteria: exactly one group is short then (and
# both criteria are no larger than the data, so their sum is an integer).
#' @noRd
.note_smd_group_short <- function(current, optimal) {
  short <- which(current < optimal)
  paste0("Group ", short, " has ", current[short], " observations and its ",
         "criterion is ", optimal[short], ". The stopping rule is met only ",
         "when each group has reached its own criterion, so it is not met ",
         "here although the total sample size (", sum(current), ") is at ",
         "least the total of the two criteria (", sum(optimal), ").")
}

# The note for groups that are not literally constant but are constant by
# the package wide tolerance. Two kinds of data are in that case: equal
# values carrying rounding noise, for which collecting more data is the
# thing to do, and data that do vary but share an offset so large that the
# variation is beyond the tenth significant digit, for which it is not. The
# function cannot tell them apart, so the note names both.
#' @noRd
.note_smd_within_tolerance <- function() {
  paste0("Within each group the observations agree to ten significant ",
         "digits, which is treated as no variation, so the standardized ",
         "mean difference is not estimated and the stopping rule is reported ",
         "as not met. If the observations do vary and share a large common ",
         "offset (times measured from a distant origin, for example), ",
         "subtract the same constant from both groups and check again: the ",
         "standardized mean difference does not depend on it. Otherwise ",
         "collect more data and check again.")
}


#' Pilot Sample Sizes for Minimum Risk Point Estimation of the Standardized
#' Mean Difference
#'
#' Returns the pilot (first stage) sample size for each group in the purely
#' sequential minimum risk point estimation (MRPE) procedure of Chattopadhyay
#' and Kelley (2017) for the standardized mean difference. It needs only the
#' cost constants, so it can be computed before any data are collected.
#'
#' @details
#' With one cost for both groups the pilot size per group is Equation 14 of
#' Chattopadhyay and Kelley (2017),
#' \deqn{m = \max\left\{m_0,\ \left\lceil (A/(2c))^{1/(2 + 2\gamma)}
#'   \right\rceil\right\}.}
#' The article's examples give 13 for \eqn{A} = 10,000 and \eqn{c} = 2.40
#' (p. 105), and 30 for \eqn{A} = 50,000 and \eqn{c} = 1 (p. 102).
#'
#' With two different costs the pilot sizes are Equations 20 and 21 of the
#' article ("A General Scenario", p. 106). Writing \eqn{c_1} for the smaller
#' cost and \eqn{c_2} for the larger, the less costly group starts with
#' \deqn{m_{1} = \max\left\{m_0,\ \left\lceil (A/c_1)^{1/(2 + 2\gamma)}
#'   \right\rceil\right\}}
#' observations and the costlier group with
#' \deqn{m_{2} = \max\left\{m_0,\ \left\lceil m_{1}
#'   \left(1 + \frac{(c_2 - c_1)\,m_{1}^2}{A}\right)^{-1/2}
#'   \right\rceil\right\}.}
#' For \eqn{A} = 50,000 and costs of 2 and 8 these are 30 and 29. The article
#' assumes the first group is the less costly one; this function accepts the
#' costs in either order and returns the sizes in the order of \code{c}.
#'
#' Equation 20 uses \eqn{A/c_1} where Equation 14 uses \eqn{A/(2c)}, so the
#' article's two cost pilot does not reduce to its one cost pilot when the
#' costs are equal. When the two costs supplied are equal this function
#' returns the one cost (Equation 14) size for both groups. Two costs count
#' as equal when their relative difference is at most the square root of
#' machine precision (about 1.5e-8), so that costs computed by arithmetic
#' (\code{24 * 0.1} against \code{2.4}) are not sent to the two cost
#' procedure by rounding error; the smaller of the two is then the cost
#' used.
#'
#' @param A The loss constant: the price of one unit of squared estimation
#'   error, in the same currency as \code{c}. The article builds it as the
#'   amount one is willing to pay to estimate the standardized mean
#'   difference to within \eqn{\epsilon}, divided by \eqn{\epsilon^2} (for
#'   example \eqn{2500 / 0.5^2 = 10000}). A single positive number.
#' @param c The cost of sampling one observation: a single positive number
#'   when an observation costs the same in both groups, or two positive
#'   numbers, \code{c(c1, c2)}, the cost in group 1 and the cost in group 2,
#'   in that order. Either group may be the costlier one.
#' @param m_0 The smallest pilot size allowed per group: a whole number of
#'   at least 4 (the default), as in the article.
#' @param gamma The tuning constant of the stopping rule's correction term,
#'   strictly between 0 and 1/2. The default, 0.49, is the value
#'   Chattopadhyay and Kelley (2017) suggest. Use the same value in
#'   \code{\link{seq_mr_smd}}.
#'
#' @return When \code{c} is a single number, a length one integer: the pilot
#'   size for each group. When \code{c} has two elements, an integer vector
#'   with names \code{group_1} and \code{group_2}.
#'
#' @references
#' Chattopadhyay, B., & Kelley, K. (2017). Estimating the standardized mean
#' difference with minimum risk: Maximizing accuracy and minimizing cost
#' with sequential estimation. \emph{Psychological Methods}, \emph{22}(1),
#' 94--113. \doi{10.1037/met0000089}
#'
#' @author Bhargab Chattopadhyay \email{bhargab@@iiitvadodara.ac.in},
#'   Ken Kelley \email{kkelley@@nd.edu}
#'
#' @examples
#' # The worked example of Chattopadhyay and Kelley (2017, p. 105): assessing
#' # one student costs 2.40 dollars, and the team would pay 2,500 dollars to
#' # estimate the standardized mean difference to within 0.5.
#' seq_mr_smd_pilot(A = 2500 / 0.5^2, c = 2.40)
#'
#' # Two costs: an observation in the second group costs four times as much.
#' seq_mr_smd_pilot(A = 50000, c = c(2, 8))
#'
#' @seealso \code{\link{seq_mr_smd}} for the stopping rule and
#'   \code{\link{seq_mr_smd_sensitivity}} to study the procedure before
#'   collecting data; in \pkg{DMAR}, \code{\link[DMAR]{mr_smd}} (the one cost
#'   procedure) and the fixed sample size planner
#'   \code{\link[DMAR]{ss_aipe_smd}}.
#' @family seq_mr_smd
#' @export
seq_mr_smd_pilot <- function(A, c, m_0 = 4, gamma = 0.49) {
  A <- .check_positive(A, "A")
  cost <- .smd_mr_cost(c)
  m_0 <- .check_count(m_0, "m_0", min = 4L)
  gamma <- .mrpe_gamma(gamma)
  m <- .smd_mr_pilot(A, cost[1L], cost[2L], gamma, m_0)
  if (length(c) == 1L) m[1L] else c(group_1 = m[1L], group_2 = m[2L])
}


#' Minimum Risk Point Estimation of the Standardized Mean Difference
#'
#' Given the data collected so far from two independent groups, evaluates
#' the purely sequential minimum risk point estimation (MRPE) stopping rule
#' of Chattopadhyay and Kelley (2017) for the standardized mean difference
#' \eqn{\delta = (\mu_1 - \mu_2)/\sigma}. The procedure is for a study in
#' which each observation costs money and a less accurate estimate also has
#' a cost: it stops sampling when the estimated total of the two is as small
#' as it can be. It targets a risk, not a confidence interval, so there is no
#' confidence level.
#'
#' The function does not collect data or keep state. Call
#' \code{\link{seq_mr_smd_pilot}} for the starting sample sizes, collect that
#' many observations, and call \code{seq_mr_smd()}. While
#' \code{attr(out, "is_satisfied")} is \code{FALSE}, bring the two groups up
#' to the sizes in the rows \code{next_n1} and \code{next_n2} and call it
#' again.
#'
#' @details
#' The estimate is \eqn{d = (\bar{X}_1 - \bar{X}_2)/s}, with \eqn{s} the
#' pooled standard deviation (Equations 2 and 3 of Chattopadhyay and Kelley,
#' 2017). The article assumes a normal population with a common variance in
#' each group; the stopping rule rests on the normal theory mean squared
#' error of \eqn{d}. Use \code{\link{seq_mr_smd_sensitivity}} to see what a
#' skewed population costs. The estimate, and with it the stopping decision,
#' is the same whatever the common units of the two groups (seconds or
#' minutes, say), up to rounding error.
#'
#' \strong{One cost for both groups.} The approximate risk (their Equation 9)
#' is
#' \deqn{R = A\left\{\frac{1}{n_1} + \frac{1}{n_2} +
#'   \frac{\delta^2}{2(n_1 + n_2)}\right\} + c\,(n_1 + n_2),}
#' the price of the squared error plus the cost of sampling. It is smallest
#' when each group has \eqn{n_c = \sqrt{A/(2c)}\;\xi} observations, with
#' \eqn{\xi^2 = 2 + \delta^2/4} (their Equations 10 and 11), where it equals
#' \eqn{4 c\, n_c} (their Equation 12). Because \eqn{\delta} is unknown, the
#' purely sequential rule (their Equation 15) stops at the first
#' \eqn{n \ge m} observations per group with
#' \deqn{n \ge \sqrt{\frac{A}{2c}}\left(V_n + n^{-\gamma}\right),
#'   \qquad V_n^2 = 2 + \frac{d_n^2}{4},}
#' where \eqn{m} is the pilot size, \eqn{V_n^2} is their Equation 13, and
#' \eqn{n^{-\gamma}} with \eqn{0 < \gamma < 1/2} keeps the rule from stopping
#' too early. In the article's example (\eqn{d} = 1.00, \eqn{A} = 10,000,
#' \eqn{c} = 2.40; p. 106) the rule is first met at \eqn{n} = 75. The article
#' states the rule for equal group sizes. If the two groups supplied differ
#' in size, this function uses the smaller size as \eqn{n}, so both groups
#' must reach the criterion.
#'
#' \strong{A different cost in each group.} This is the article's "A General
#' Scenario" (pp. 106 to 107, Equations 17 to 21, with Appendix B, Equations
#' 41 to 43). In this paragraph, subscript 1 marks the less costly group, as
#' in the article. The approximate risk (their Equation 17) replaces the
#' sampling cost above with \eqn{c_1 n_1 + c_2 n_2}. Setting its derivatives
#' to zero (their Equation 41) gives the optimal sizes: \eqn{n_2^* = n_1^*
#' \, r(n_1^*)} with \eqn{r(n) = \{1 + (c_2 - c_1) n^2 / A\}^{-1/2}} (their
#' Equations 19 and 43), and \eqn{n_1^*} solves
#' \deqn{n_1^* = \sqrt{\frac{A}{c_1}}
#'   \left\{1 + \frac{\delta^2}{2\,\{1 + r(n_1^*)\}^2}\right\}^{1/2}.}
#' The sequential procedure samples the less costly group one observation
#' (or batch) at a time, keeps the costlier group at \eqn{\lceil n_1\,
#' r(n_1) \rceil} observations (their Equation 21 and the Stage II text), and
#' stops at the first \eqn{n_1} with
#' \deqn{n_1 \ge \sqrt{\frac{A}{c_1}}\left[
#'   \left\{1 + \frac{d^2}{2\,\{1 + r(n_1)\}^2}\right\}^{1/2}
#'   + n_1^{-\gamma}\right],}
#' where \eqn{d} is computed from all observations in both groups.
#'
#' \strong{Where the function departs from the article as printed.}
#' \enumerate{
#'   \item Equations 18 and 42, and the stopping rule built on them, are
#'     printed with \eqn{\{1 + r(n_1)\}} to the first power. The first order
#'     conditions in Equation 41 give the second power shown above. Only the
#'     second power minimizes the risk in Equation 17, and only it reduces
#'     to Equation 10 when the costs are equal (as Lemma 6 of the article
#'     states). The function uses the second power. For \eqn{A} = 50,000,
#'     costs of 2 and 8, and \eqn{\delta} = 1, the printed form gives
#'     \eqn{n_1^*} = 183.4; direct minimization of Equation 17 gives 175.7,
#'     which is what the form above gives.
#'   \item The article's two cost procedure does not reduce to its one cost
#'     procedure when \eqn{c_1 = c_2}: its pilot (Equation 20) uses
#'     \eqn{A/c_1} where Equation 14 uses \eqn{A/(2c)}, and its correction
#'     term is \eqn{\sqrt{A/c_1}\,n^{-\gamma}} where Equation 15 has
#'     \eqn{\sqrt{A/(2c)}\,n^{-\gamma}}. The function uses the one cost
#'     procedure whenever the two costs supplied are equal (so
#'     \code{c = c(3, 3)} and \code{c = 3} give identical results) and the
#'     two cost procedure whenever they differ. A small difference between
#'     the costs therefore moves the criterion by a visible amount; this is
#'     a property of the article's two procedures, not of the software. So
#'     that rounding error cannot make that choice, two costs whose relative
#'     difference is at most the square root of machine precision (about
#'     1.5e-8) count as equal, and the smaller of the two is used:
#'     \code{c = c(24 * 0.1, 2.4)} gives the same result as
#'     \code{c = 2.4}.
#'   \item The article assumes the first group is the less costly one. The
#'     function applies the rule with whichever group is less costly in that
#'     role, so the order of the groups does not matter.
#'   \item The function is given data, not a sampling plan, so it also
#'     requires the costlier group to have at least \eqn{\lceil n_1\,
#'     r(n_1) \rceil} observations before reporting the rule as met.
#' }
#' Everything else in the two cost procedure is the article's.
#'
#' \strong{Sampling with two costs.} The rows \code{next_n1} and
#' \code{next_n2} give the size each group should have at the next check, so
#' the article's sampling order does not have to be coded by hand: while the
#' rule is not met, collect observations until the groups have those sizes
#' and call \code{seq_mr_smd()} again. Starting from the pilot sizes, each
#' such step adds one observation to the less costly group and brings the
#' costlier group up to \eqn{\lceil n_1\, r(n_1) \rceil}, with \eqn{n_1}
#' the less costly group's new size, which is the article's procedure with
#' one observation per stage (its Stage II, p. 107, and the steps after it).
#' The rows are in the order of the arguments, group 1 and group 2,
#' whichever of them is the less costly one.
#' If the costlier group has fewer observations than the article pairs with
#' the less costly group's current size, the rows ask for that shortfall
#' alone (the less costly group's row equals its current size). That gives a
#' way to check less often than after every observation: add several
#' observations to the less costly group and call \code{seq_mr_smd()}; if
#' the rows then ask only for the costlier group's shortfall, collect it and
#' call again. With one cost the rows bring each group up to one more than
#' the smaller of the two current sizes. Sampling in another order (adding
#' to the costlier group and
#' topping up the less costly one) oversamples: the function can then report
#' the rule as met, because both groups are past their criteria, at a
#' sampling cost well above the minimum.
#'
#' \strong{What to expect.} The article proves that the expected final
#' sample size and the risk approach their optimal values as \eqn{A/c}
#' grows; at practical sizes the agreement is approximate, and the
#' correction term makes the procedure take somewhat more than the optimal
#' number of observations. With one cost, this function reproduces the
#' article's Tables 1 and 3 within Monte Carlo error. For \eqn{A} = 50,000,
#' \eqn{c} = 1, and \eqn{\delta = 0}, for example, the mean final sample size
#' per group in 5,000 replications was 235.04 (the article reports 235.04;
#' the optimal size is 223.6), and the ratio of the realized risk to the
#' minimum risk was 1.00. With costs of 2 and 8, \eqn{A} = 50,000, and
#' \eqn{\delta = 1}, the mean final sizes were 189.2 and 82.7 against
#' optimal sizes of 175.7 and 81.0, and the risk ratio was 0.99 (Monte Carlo
#' standard error 0.01). At smaller sample sizes the excess is larger: for
#' \eqn{A} = 1,000,000, \eqn{c} = 500, and \eqn{\delta = 0} the mean final
#' size was 50.0 against an optimal 44.7, as in the article's Table 3.
#'
#' @param group_1,group_2 Numeric vectors: the observations collected so far
#'   in the first and the second group, at least four in each. Four is the
#'   smallest pilot size the article allows, and
#'   \code{\link{seq_mr_smd_pilot}} never returns less.
#' @inheritParams seq_mr_smd_pilot
#' @param gamma The tuning constant of the correction term, strictly between
#'   0 and 1/2; default 0.49. Use the same value as in
#'   \code{\link{seq_mr_smd_pilot}}.
#' @param na.rm Should missing values be removed before the computation?
#'   \code{TRUE} (the default) removes them; \code{FALSE} stops with an
#'   error if any are present.
#'
#' @return An object of class \code{c("smsd_check", "data.frame")} with
#'   columns \code{term} and \code{value} and one row for each of
#'   \describe{
#'     \item{\code{smd}}{the standardized mean difference, group 1 minus
#'       group 2, divided by the pooled standard deviation;}
#'     \item{\code{risk}}{the approximate risk (Equation 9, or Equation 17
#'       with two costs) at the current group sizes, with the estimate in
#'       place of \eqn{\delta};}
#'     \item{\code{mean_1}, \code{mean_2}}{the two sample means;}
#'     \item{\code{optimal_n1}, \code{optimal_n2}}{the criterion sample size
#'       for each group. With one cost both equal the right side of
#'       Equation 15, rounded up. With two costs, the less costly group's
#'       value is the right side of its stopping rule, rounded up, and the
#'       costlier group's value is the size the article pairs with the
#'       larger of the less costly group's current size and its criterion.
#'       The rule is met exactly when each current size is at least its
#'       criterion;}
#'     \item{\code{current_n1}, \code{current_n2}}{the number of
#'       observations supplied in each group;}
#'     \item{\code{next_n1}, \code{next_n2}}{the size each group should have
#'       at the next check. When the rule is met they equal the current
#'       sizes: nothing more is called for. Otherwise, with one cost, each
#'       group is brought up to one more than the smaller current size. With
#'       two costs the less costly group gets one more observation and the
#'       costlier group is brought up to the size the article pairs with the
#'       less costly group's new size; but if the costlier group is short of
#'       the size paired with the less costly group's current size, only
#'       that shortfall is called for. See "Sampling with two costs" in
#'       Details.}
#'   }
#'   Attributes: \code{is_satisfied} (\code{TRUE} or \code{FALSE}, never
#'   \code{NA}: whether to stop sampling), \code{current_n} (the total
#'   sample size), \code{criterion_n} (the sum of the two criteria),
#'   \code{method}, \code{loss_A}, \code{cost_1}, \code{cost_2},
#'   \code{gamma}, \code{variance_unusable} (always \code{FALSE} here,
#'   because \eqn{V_n^2} is never below 2; carried for uniformity with the
#'   procedures whose variance estimate can be zero, negative, or not
#'   finite), and \code{notes} (a character vector, printed under the
#'   result). \code{conf_level} and \code{width} are \code{NULL}.
#'   The rule is decided group by group, so the two totals are not to be
#'   compared with each other: \code{current_n} can be at least
#'   \code{criterion_n} while one group is still short of its own
#'   criterion. In that case a note names the group.
#'   If every observation within each group is the same, or the same to ten
#'   significant digits, there is no pooled standard deviation to divide by
#'   and the standardized mean difference is not estimable: \code{smd},
#'   \code{risk}, the two criteria, and \code{criterion_n} are \code{NA},
#'   the rule is reported as not met, and a note says why (\code{next_n1}
#'   and \code{next_n2} do not depend on the estimate and are still
#'   reported). That test is relative to the size of the values, so
#'   observations that do vary but share a large common offset (each
#'   group's range is at most \eqn{10^{-10}} times its largest absolute
#'   value, as it can be for times recorded from a distant origin) are
#'   treated as the same too. The note then says to subtract the same
#'   constant from both groups, which leaves the standardized mean
#'   difference unchanged, and to check again. See
#'   \code{\link{smsd_check_methods}} for the print, summary, and conversion
#'   methods.
#'
#' @references
#' Chattopadhyay, B., & Kelley, K. (2017). Estimating the standardized mean
#' difference with minimum risk: Maximizing accuracy and minimizing cost
#' with sequential estimation. \emph{Psychological Methods}, \emph{22}(1),
#' 94--113. \doi{10.1037/met0000089}
#'
#' @author Bhargab Chattopadhyay \email{bhargab@@iiitvadodara.ac.in},
#'   Ken Kelley \email{kkelley@@nd.edu}
#'
#' @examples
#' # The pilot data of Chattopadhyay and Kelley (2017, p. 105): reading
#' # scores of children with high and low exposure to same language
#' # subtitling. Assessing one child costs 2.40 dollars; A = 2500 / 0.5^2.
#' high_sls <- c(11, 7, 22, 13, 6, 9, 11, 16, 12, 17, 14, 8, 16)
#' low_sls  <- c(3, 6, 10, 8, 14, 5, 12, 10, 6, 8, 13, 5, 9)
#' seq_mr_smd_pilot(A = 10000, c = 2.40)
#' seq_mr_smd(high_sls, low_sls, A = 10000, c = 2.40)
#'
#' # Two costs: a visit to a child in the second group costs twice as much.
#' # While the rule is not met, bring each group up to the size in its
#' # next_n row and check again. The rows follow the article's sampling
#' # order whichever group is the costlier one, so the loop also runs as
#' # written with cost <- c(6, 3).
#' set.seed(113)
#' A <- 2000
#' cost <- c(3, 6)                  # group 1, group 2
#' true_mean <- c(0.5, 0)
#' m <- seq_mr_smd_pilot(A = A, c = cost)
#' y <- list(rnorm(m[["group_1"]], mean = true_mean[1]),
#'           rnorm(m[["group_2"]], mean = true_mean[2]))
#' check <- seq_mr_smd(y[[1]], y[[2]], A = A, c = cost)
#' while (!attr(check, "is_satisfied")) {
#'   next_n <- check$value[check$term %in% c("next_n1", "next_n2")]
#'   for (j in 1:2) {
#'     y[[j]] <- c(y[[j]],
#'                 rnorm(next_n[j] - length(y[[j]]), mean = true_mean[j]))
#'   }
#'   check <- seq_mr_smd(y[[1]], y[[2]], A = A, c = cost)
#' }
#' check
#'
#' @seealso \code{\link{seq_mr_smd_pilot}} for the pilot sizes and
#'   \code{\link{seq_mr_smd_sensitivity}} to study the procedure before
#'   collecting data; the accuracy in parameter estimation companion
#'   \code{\link{seq_aipe_smd}}; in \pkg{DMAR}, \code{\link[DMAR]{mr_smd}}
#'   (the one cost procedure, from summary statistics), the estimator
#'   \code{\link[DMAR]{smd}}, and the fixed sample size planner
#'   \code{\link[DMAR]{ss_aipe_smd}}.
#' @family seq_mr_smd
#' @export
seq_mr_smd <- function(group_1, group_2, A, c, gamma = 0.49, na.rm = TRUE) {
  A <- .check_positive(A, "A")
  cost <- .smd_mr_cost(c)
  gamma <- .mrpe_gamma(gamma)
  group_1 <- .mr_check_data(group_1, na.rm, min_n = 4L, arg = "group_1")
  group_2 <- .mr_check_data(group_2, na.rm, min_n = 4L, arg = "group_2")
  n1 <- length(group_1)
  n2 <- length(group_2)
  c1 <- cost[1L]
  c2 <- cost[2L]
  est <- .smd_mr_estimate(group_1, group_2)
  smd <- est$smd
  mean_1 <- est$mean_1
  mean_2 <- est$mean_2

  terms <- c("smd", "risk", "mean_1", "mean_2",
             "optimal_n1", "optimal_n2", "current_n1", "current_n2",
             "next_n1", "next_n2")
  method <- "minimum risk point estimation (standardized mean difference)"
  extra <- list(loss_A = A, cost_1 = c1, cost_2 = c2, gamma = gamma,
                variance_unusable = FALSE)

  if (is.na(smd)) {
    # No estimate, so the rule cannot be evaluated and is not met. The sizes
    # for the next check do not depend on the estimate.
    next_n <- .smd_mr_next(n1, n2, A, c1, c2)
    return(new_smsd_check(
      term = terms,
      value = c(NA_real_, NA_real_, mean_1, mean_2, NA_real_, NA_real_,
                n1, n2, next_n[1L], next_n[2L]),
      method = method,
      is_satisfied = FALSE,
      current_n = n1 + n2,
      criterion_n = NA_real_,
      extra = extra,
      notes = if (est$within_tolerance) {
        .note_smd_within_tolerance()
      } else {
        .note_not_estimable()
      }
    ))
  }

  rule <- .smd_mr_rule(smd, n1, n2, A, c1, c2, gamma)
  advice <- paste0("Check 'A', 'c', and the data: a pooled standard ",
                   "deviation close to zero makes the standardized mean ",
                   "difference, and with it the criterion, very large.")
  optimal <- c(
    .as_sample_size(rule$criterion[1L], "optimal sample size",
                    advice = advice),
    .as_sample_size(rule$criterion[2L], "optimal sample size",
                    advice = advice))
  current <- c(n1, n2)
  # When the rule is met nothing more is called for.
  next_n <- if (rule$satisfied) {
    as.numeric(current)
  } else {
    .smd_mr_next(n1, n2, A, c1, c2)
  }
  # Summed as doubles: two criteria can each be representable as an integer
  # while their total is not.
  criterion_total <- sum(as.numeric(optimal))
  notes <- if (!rule$satisfied && sum(current) >= criterion_total) {
    .note_smd_group_short(current, optimal)
  } else {
    character(0)
  }

  new_smsd_check(
    term = terms,
    value = c(smd, .smd_mr_risk(n1, n2, smd^2, A, c1, c2), mean_1, mean_2,
              optimal[1L], optimal[2L], n1, n2, next_n[1L], next_n[2L]),
    method = method,
    is_satisfied = rule$satisfied,
    current_n = n1 + n2,
    criterion_n = criterion_total,
    extra = extra,
    notes = notes
  )
}

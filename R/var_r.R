# Estimators of the asymptotic variance xi^2 of a sample correlation
# coefficient (sqrt(n) (r_n - rho) is asymptotically N(0, xi^2)), from
# Kelley, Darku, and Chattopadhyay (2019). Used by seq_aipe_r(). Internal;
# not exported.
#
# x, y: paired numeric vectors (at least 4 pairs).
# method: "pearson", "kendall", or "spearman".
# na.rm: drop a pair when either member is missing? Default FALSE.
# floor: apply the lower bound n^-3? Default TRUE.
# Returns a single number: max{V_n^2, n^-3} (their Equation 8, with
# gamma = 3), or, with floor = FALSE, V_n^2 itself. For the Pearson
# estimator that can be zero or negative; for the rank based estimators it
# is zero or positive; and it is NaN when either variable has no variation
# (see .is_constant()).
#
# Pearson (their Equation 3 for xi^2 and Equation 70 for the estimator):
#
#   V_n^2 = (r^2 / 4) { mu40_hat / s_X^4 + mu04_hat / s_Y^4
#                       + 2 mu22_hat / (s_X^2 s_Y^2) + 4 mu22_hat / s_XY^2
#                       - 4 mu31_hat / (s_XY s_X^2)
#                       - 4 mu13_hat / (s_XY s_Y^2) },
#
# with .mu4() for mu40_hat and mu04_hat (Equations 71 and 72) and .mu22(),
# .mu31(), .mu13() (Equations 73 and 74). Under bivariate normality
# xi^2 = (1 - rho^2)^2. Being a difference of moment estimates, the raw
# value can be negative, and at small n it often is.
#
# As printed, Equation 70 divides by the sample covariance s_XY, so it is
# 0 / 0 when the sample covariance is exactly zero, which happens with
# discrete data (1.4% of samples of 10 pairs of independent 5 point items).
# The factor r^2 = s_XY^2 / (s_X^2 s_Y^2) cancels every such division, and
# the estimator is computed with the cancellation carried out,
#
#   V_n^2 = (1 / 4) { r^2 [ mu40_hat / s_X^4 + mu04_hat / s_Y^4
#                           + 2 mu22_hat / (s_X^2 s_Y^2) ]
#                     + 4 mu22_hat / (s_X^2 s_Y^2)
#                     - 4 r mu31_hat / (s_X^3 s_Y)
#                     - 4 r mu13_hat / (s_X s_Y^3) },
#
# which is the same number whenever s_XY is not zero and is the limiting
# value mu22_hat / (s_X^2 s_Y^2) when it is.
#
# Kendall's tau (their Equations 27 to 29): 16 times the sample variance of
#
#   W_i = (2 / n) #{k : R_xk <= R_xi and R_yk <= R_yi}
#         - R_xi / (n + 1) - R_yi / (n + 1),
#
# with R the ranks.
#
# Spearman's rho (their Equations 51 to 53, the estimator of Kojadinovic and
# Yan, 2010): 144 times the sample variance of
#
#   Z_i = R_xi R_yi / (n + 1)^2
#         + (1 / n) sum_k 1{R_xi <= R_xk} R_yk / (n + 1)
#         + (1 / n) sum_k 1{R_yi <= R_yk} R_xk / (n + 1).
#
# The article writes its Spearman stopping rule (Equation 54) with the
# estimator of Genest and Favre (2007; Equations 45 to 50) and then offers
# Equation 51 as "a very simple but also consistent estimator"; Equation 51
# is the one implemented, and it reproduces the article's Table 10.
#
# The two rank based estimators assume continuous data (no ties). With tied
# values they are too small; seq_aipe_r() warns.
# Sorting, cumulative sums, and a cumulative count tree evaluate these same
# estimators in O(n log n) time with O(n) working memory, including ties.
#
# Exact arithmetic for the rank based estimators. Both are zero on paper for
# some data that vary: every W_i equals 2 / n - 1 when the pairs are
# perfectly discordant (Kendall's tau is -1), and every Z_i is the same for,
# among others, two dichotomous variables in perfect opposition. Computed
# from ranks divided by n + 1, such a set of equal numbers has a sample
# variance of about 1e-32 rather than zero, positive for most n, and a
# positive raw estimate is used through the bound n^-3 (R/var_floor.R): the
# check would report the stopping rule as met (with limits of -1.02 and
# -0.98 at n = 10 and a width of 0.4) for a sample that carries no
# information about the variance.
# The estimators are therefore computed from n (n + 1) W_i and
# n (n + 1)^2 Z_i, which are whole multiples of 1/2 and 1/4 (ranks are
# multiples of 1/2) and so exact in floating point; the sample variance of
# equal numbers is then exactly zero, and that of unequal ones an ordinary
# positive number. The results for other data differ from the direct form
# only in the last digits (a relative difference of at most 2e-14 in 3,000
# samples).
#
# Pearson, cancellation. The Pearson estimator is a difference of terms and
# is zero on paper for some data too (collinear pairs whose sample kurtosis
# is (n^2 - 3) / (n - 1)^2, which discrete data can produce); what is
# computed then is rounding noise of either sign. .aipe_zero_if_cancelled()
# returns zero for it.
#
# Authors: Bhargab Chattopadhyay, Ken Kelley.
#' @noRd
.var_r <- function(x, y, method = c("pearson", "kendall", "spearman"),
                   na.rm = FALSE, floor = TRUE) {
  method <- .check_choice(method, c("pearson", "kendall", "spearman"),
                          "method")
  if (na.rm) {
    ok <- !is.na(x) & !is.na(y)
    x <- x[ok]
    y <- y[ok]
  }
  n <- length(x)
  if (.is_constant(x) || .is_constant(y)) {
    # No variation in a variable: the correlation and its variance are
    # undefined.
    return(.var_floor(NaN, n, floor))
  }

  if (method == "pearson") {
    # The estimate does not depend on the scale of either variable; see
    # .aipe_scale().
    x <- x / .aipe_scale(x)
    y <- y / .aipe_scale(y)
    r <- stats::cor(x, y)
    vx <- stats::var(x)
    vy <- stats::var(y)
    sx_sy <- sqrt(vx) * sqrt(vy)
    # The moment estimates in standard deviation units.
    m40 <- .mu4(x) / vx^2
    m04 <- .mu4(y) / vy^2
    m22 <- .mu22(x, y) / (vx * vy)
    m31 <- .mu31(x, y) / (vx * sx_sy)
    m13 <- .mu13(x, y) / (vy * sx_sy)
    vn2 <- (r^2 * (m40 + m04 + 2 * m22) + 4 * m22 - 4 * r * (m31 + m13)) / 4
    vn2 <- .aipe_zero_if_cancelled(
      vn2, c(r^2 * m40, r^2 * m04, 2 * r^2 * m22, 4 * m22, 4 * r * m31,
             4 * r * m13) / 4
    )

  } else if (method == "kendall") {
    rank_x <- rank(x)
    rank_y <- rank(y)
    # The bivariate empirical distribution function at each observation, as
    # a count: the number of pairs at or below it on both margins.
    count <- .rank_lower_count(rank_x, rank_y)
    # n (n + 1) W_i; see "Exact arithmetic" above.
    w <- 2 * (n + 1) * count - n * (rank_x + rank_y)
    vn2 <- 16 * stats::var(w) / (n * (n + 1))^2

  } else { # spearman
    rank_x <- rank(x)
    rank_y <- rank(y)
    # n (n + 1)^2 Z_i; see "Exact arithmetic" above.
    v <- n * rank_x * rank_y + (n + 1) *
      (.rank_upper_sum(rank_x, rank_y) + .rank_upper_sum(rank_y, rank_x))
    vn2 <- 144 * stats::var(v) / (n * (n + 1)^2)^2
  }

  .var_floor(vn2, n, floor)
}

# Sum weights at ranks greater than or equal to each observation's rank.
# The first position of a tied rank includes the entire tied group.
# Sorting and cumulative sums replace the quadratic comparison matrices.
#' @noRd
.rank_upper_sum <- function(ranks, weights) {
  ordered <- order(ranks)
  sums <- rev(cumsum(rev(weights[ordered])))
  sums[match(ranks, ranks[ordered])]
}

# Counts at or below each observation on both margins. A cumulative count
# tree stores the y ranks as observations are visited in x order. Insert a
# whole tied x group before querying it, so the inequalities include ties
# on both margins. Time is O(n log n), working memory O(n).
#' @noRd
.rank_lower_count <- function(rank_x, rank_y) {
  n <- length(rank_x)
  ordered <- order(rank_x)
  y_index <- match(rank_y, sort(unique(rank_y)))
  tree <- numeric(max(y_index))
  count <- numeric(n)
  ends <- c(which(diff(rank_x[ordered]) != 0), n)
  first <- 1L
  for (last in ends) {
    group <- ordered[seq.int(first, last)]
    for (observation in group) {
      index <- y_index[observation]
      while (index <= length(tree)) {
        tree[index] <- tree[index] + 1
        index <- index + bitwAnd(index, -index)
      }
    }
    for (observation in group) {
      index <- y_index[observation]
      total <- 0
      while (index > 0L) {
        total <- total + tree[index]
        index <- index - bitwAnd(index, -index)
      }
      count[observation] <- total
    }
    first <- last + 1L
  }
  count
}

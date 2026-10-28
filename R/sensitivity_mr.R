# Monte Carlo sensitivity wrappers for the minimum risk point estimation
# procedures (seq_mr_cv, seq_mr_smd), with their internal engines. The
# engines call the same rule functions as the exported check functions
# (.mr_cv_rule and .smd_mr_rule), so a replication is exactly the procedure
# a user would run by hand.

# The replication count, batch size, and cap shared by the two wrappers.
# `pilot` is the largest pilot group size; the cap may not be below it.
#' @noRd
.mrpe_sim_args <- function(G, batch, max_n, pilot) {
  G <- .check_count(G, "G", min = 2L)
  batch <- .check_count(batch, "batch", min = 1L)
  max_n <- .check_count(max_n, "max_n", min = 1L)
  if (max_n < pilot) {
    .smsd_stop("'max_n' (", max_n, ") is smaller than the pilot sample ",
               "size (", pilot, "). Raise 'max_n'.")
  }
  list(G = G, batch = batch, max_n = max_n)
}

# The asymptotic variance of sqrt(n) (k_n - kappa) for a population with
# coefficient of variation kappa, skewness `skew`, and kurtosis `kurt`
# (Kelley, Darku, & Chattopadhyay, 2018, Equation 30, written in
# standardized moments).
#' @noRd
.cv_asymptotic_variance <- function(kappa, skew, kurt) {
  kappa^2 * (kurt - 1) / 4 - skew * kappa^3 + kappa^4
}

# The data generating mechanisms of seq_mr_cv_sensitivity(): a sampler and
# the population's asymptotic variance, for a given coefficient of variation.
# A parameter of the generating distribution that is computed from the
# coefficient of variation must be a positive finite number, as in the
# accuracy in parameter estimation sibling (seq_aipe_cv_sensitivity()), and
# is refused with the same message when it is not: the gamma shape
# underflows to zero and the lognormal log standard deviation overflows for
# an extreme value, and the latter is zero when 1 + kappa^2 rounds to 1.
#' @noRd
.mr_cv_mechanism <- function(distribution, kappa) {
  if (distribution == "gamma") {
    shape <- .sensitivity_check_parameter(
      1 / kappa^2, "gamma", "shape, 1 / true_value^2", kappa)
    list(draw = function(k) stats::rgamma(k, shape = shape, scale = 1),
         xi2 = .cv_asymptotic_variance(kappa, 2 * kappa, 3 + 6 * kappa^2))
  } else if (distribution == "lognormal") {
    w <- 1 + kappa^2
    sdlog <- .sensitivity_check_parameter(
      sqrt(log(w)), "lognormal",
      "log standard deviation, sqrt(log(1 + true_value^2))", kappa)
    list(draw = function(k) stats::rlnorm(k, meanlog = 0, sdlog = sdlog),
         xi2 = .cv_asymptotic_variance(kappa, (w + 2) * sqrt(w - 1),
                                       w^4 + 2 * w^3 + 3 * w^2 - 3))
  } else {
    list(draw = function(k) stats::rnorm(k, mean = 1, sd = kappa),
         xi2 = .cv_asymptotic_variance(kappa, 0, 3))
  }
}

# The benchmark of a study (the best fixed sample size or sizes, and the
# risk there) must be positive finite numbers: the summary divides by them.
# A population value at the edge of double precision passes every other
# check and gives a benchmark of zero, Inf, or NaN (the fourth power of a
# coefficient of variation beyond about 1e77 overflows, and of one below
# about 1e-81 underflows), and the summary would then hold rows of Inf and
# NaN beside ordinary ones. `benchmark` is the named benchmark; `inputs`
# the named arguments it was computed from, shown in the message.
#' @noRd
.mrpe_check_benchmark <- function(benchmark, inputs) {
  bad <- which(!(is.finite(benchmark) & benchmark > 0))
  if (!length(bad)) {
    return(invisible(TRUE))
  }
  .smsd_stop("'true_value' is too extreme to simulate with this 'A' and ",
             "'c': the benchmark the procedure is compared with must be ",
             "positive finite numbers, and ", names(benchmark)[bad[1L]],
             " is ", format(benchmark[[bad[1L]]]), "; got ",
             paste0(names(inputs), " = ", vapply(inputs, format, character(1)),
                    collapse = ", "), ".")
}

# One replication of the coefficient of variation procedure. Returns
# c(n, estimate, capped). A stage at which the estimate does not exist, or
# at which the criterion is too large to represent (a sample mean barely
# above zero), is a stage at which the rule is not met: sampling goes on.
#' @noRd
.mr_cv_run_one <- function(draw, pilot, batch, max_n, A, c, gamma,
                           variance) {
  x <- draw(pilot)
  repeat {
    rule <- .mr_cv_rule(x, A, c, gamma, variance)
    n <- length(x)
    if (rule$satisfied || n >= max_n) break
    x <- c(x, draw(min(batch, max_n - n)))
  }
  c(n = n, estimate = rule$cv, capped = as.numeric(!rule$satisfied))
}

# One replication of the standardized mean difference procedure, sampled in
# the order the article describes: after each check that is not met, the
# groups are brought up to the sizes .smd_mr_next() gives, which are the
# sizes seq_mr_smd() reports in its next_n1 and next_n2 rows when `batch` is
# 1. With one cost each stage adds `batch` observations to each group. With
# two costs each stage adds `batch` observations to the less costly group
# and brings the costlier group up to the size the article pairs with it.
# The last stage is shortened so that no group exceeds max_n. The less
# costly group is drawn first. Returns c(n1, n2, estimate, capped).
#' @noRd
.smd_mr_run_one <- function(draw_1, draw_2, pilot, batch, max_n, A, c1, c2,
                            gamma) {
  g1 <- draw_1(pilot[1L])
  g2 <- draw_2(pilot[2L])
  repeat {
    n1 <- length(g1)
    n2 <- length(g2)
    d <- .smd_mr_estimate(g1, g2)$smd
    satisfied <- !is.na(d) &&
      .smd_mr_rule(d, n1, n2, A, c1, c2, gamma)$satisfied
    if (satisfied || max(n1, n2) >= max_n) break
    target <- .smd_mr_next(n1, n2, A, c1, c2,
                           step = min(batch, max_n - max(n1, n2)))
    # Every stage must add at least one observation, or this loop would not
    # end. .smd_mr_next() always asks for one (the step is at least 1 here);
    # the check is what guarantees that a defect there shows up as an error
    # and not as a simulation that never returns.
    if (!any(target > c(n1, n2))) {
      .smsd_stop("Internal error in the simulation of the minimum risk ",
                 "procedure: the next stage adds no observation to either ",
                 "group (current sizes ", n1, " and ", n2, "). Please ",
                 "report this to the maintainer.")
    }
    if (c1 <= c2) {
      g1 <- c(g1, draw_1(target[1L] - n1))
      g2 <- c(g2, draw_2(target[2L] - n2))
    } else {
      g2 <- c(g2, draw_2(target[2L] - n2))
      g1 <- c(g1, draw_1(target[1L] - n1))
    }
  }
  c(n1 = n1, n2 = n2, estimate = d, capped = as.numeric(!satisfied))
}

# Summaries of the final estimates and realized losses shared by the two
# wrappers: one value per replication, NA where a replication stopped by the
# cap had data that give no estimate. Such replications are left out.
#' @noRd
.mrpe_risk_rows <- function(estimate, loss, minimum_risk) {
  used <- !is.na(loss)
  k <- sum(used)
  mean_risk <- if (k > 0L) mean(loss[used]) else NA_real_
  se <- if (k > 1L) stats::sd(loss[used]) / sqrt(k) else NA_real_
  c(mean_estimate = if (k > 0L) mean(estimate[used]) else NA_real_,
    mean_risk = mean_risk,
    risk_efficiency = mean_risk / minimum_risk,
    se_risk_efficiency = se / minimum_risk,
    regret = mean_risk - minimum_risk)
}


#' Monte Carlo Sensitivity Analysis for Minimum Risk Point Estimation
#'
#' These functions run a minimum risk point estimation (MRPE) procedure from
#' its pilot sample to its stopping point \code{G} times, on data simulated
#' from a population whose true value is known, and summarize what happened:
#' how many observations the procedure took compared with the best fixed
#' sample size, and how the realized risk compares with the risk of that
#' best fixed design. Use them before collecting data, to see what a choice
#' of the loss constant \code{A} and the sampling cost \code{c} implies, and
#' to see how the procedure behaves when checking in batches or when the
#' population is not the one the method assumes.
#'
#' The function \code{seq_mr_cv_sensitivity()} runs \code{\link{seq_mr_cv}}
#' (Chattopadhyay & Kelley, 2016), and the function
#' \code{seq_mr_smd_sensitivity()} runs \code{\link{seq_mr_smd}}
#' (Chattopadhyay & Kelley, 2017). Each replication
#' applies the same stopping rule as the check function, to the same kind of
#' growing sample a user would collect.
#'
#' @details
#' \strong{What a replication does.} It draws the pilot sample
#' (\code{\link{seq_mr_cv_pilot}} or \code{\link{seq_mr_smd_pilot}}),
#' evaluates the stopping rule, and, while the rule is not met, adds
#' \code{batch} observations and evaluates it again. For the standardized
#' mean difference with one cost, each stage adds \code{batch} observations
#' to each group. With two costs it follows the article's order: each stage
#' adds \code{batch} observations to the less costly group and then brings
#' the costlier group up to the size the article pairs with it. With
#' \code{batch = 1} these are the sizes \code{\link{seq_mr_smd}} reports in
#' its rows \code{next_n1} and \code{next_n2}, so a replication is the loop
#' that function's help describes. A stage at which the check function would
#' report the rule as not met for want of an estimate (a sample mean that is
#' not positive, constant data, a variance estimate that is zero, negative,
#' or not finite), or at which the criterion is too large to represent, is a
#' stage at which the rule is not met: sampling goes on. The realized loss
#' of a replication is the squared error of its final estimate priced at
#' \code{A}, plus what its observations cost.
#'
#' \strong{Run time.} The stopping rule is evaluated on all the data
#' collected so far after every batch, so the time for one replication grows
#' faster than its final sample size, and toward the square of it for large
#' samples. On the machine used to check the package, with
#' \code{batch = 1}, one replication of the coefficient of variation
#' procedure took about 0.004 seconds at a final sample size near 100, 0.07
#' seconds near 1,100, and 0.5 seconds near 4,000 (about nine minutes for
#' the default \code{G}); \code{batch = 50} brought the last figure to 0.01
#' seconds. When the optimal sample size is in the thousands, use a larger
#' \code{batch}, a smaller \code{G}, or a smaller \code{max_n}.
#'
#' \strong{The benchmark.} \code{optimal_n} (or \code{optimal_n1} and
#' \code{optimal_n2}) is the fixed sample size that minimizes the large
#' sample risk for the simulated population, and \code{minimum_risk} is the
#' risk there. They are not rounded. For the coefficient of variation they
#' are \eqn{n_c = (A\,\xi^2/c)^{1/2}} and \eqn{2 c\, n_c} (Equations 15 and
#' 16 of Chattopadhyay & Kelley, 2016), with \eqn{\xi^2} the asymptotic
#' variance of the simulated population (Equation 30 of Kelley, Darku, &
#' Chattopadhyay, 2018) whichever \code{variance} is requested. For the
#' standardized mean difference they minimize the approximate risk of
#' Chattopadhyay and Kelley (2017, Equation 9, or Equation 17 with two
#' costs). That risk assumes normal populations. For \code{distribution =
#' "skewed"} the benchmark instead uses the large sample mean squared error
#' of the estimate for the skewed population, \eqn{1/n_1 + 1/n_2 +
#' 3\delta^2 / \{2(n_1 + n_2)\}}, so that the risk efficiency shows what the
#' normality assumption costs.
#'
#' \strong{Reading the result.} A procedure that works as intended has
#' \code{ratio_n} and \code{risk_efficiency} near 1. Both guarantees are
#' asymptotic (they hold as \code{A / c} grows). At small sample sizes the
#' correction term makes the procedure take more than the optimal number of
#' observations (for the standardized mean difference with \code{A = 500}
#' and \code{c = 2}, about 19 per group against 16), and for a strongly
#' skewed population the coefficient of variation procedure stops early
#' (for gamma data with a coefficient of variation of 1.3, \code{A = 900},
#' and \code{c = 1}, about 26 observations against 45). The benchmark is a
#' large sample approximation, so at small sample sizes
#' \code{risk_efficiency} can fall below 1; that reflects the approximation,
#' not a procedure better than the best fixed design.
#' \code{risk_efficiency} is also an average of very variable losses, so
#' judge it against \code{se_risk_efficiency}. If \code{pct_hit_max_n} is
#' above zero, some replications were stopped by \code{max_n} and not by
#' the rule, and the other rows are censored.
#'
#' @param true_value The population value. For the coefficient of
#'   variation it must be positive; for the standardized mean difference it
#'   must be less than 100 million in absolute value. A value so extreme
#'   that the population cannot be generated, or that the benchmark is not
#'   a positive finite number, is an error.
#' @param A The loss constant, the price of one unit of squared estimation
#'   error. A single positive number.
#' @param c The cost of sampling one observation. A single positive number
#'   or, for the standardized mean difference, two positive numbers,
#'   \code{c(c1, c2)}, the cost in each group.
#' @param distribution The population to simulate from. For the coefficient
#'   of variation: \code{"gamma"} (the default), \code{"lognormal"}, or
#'   \code{"normal"} (mean 1), each with coefficient of variation
#'   \code{true_value}. For the standardized mean difference:
#'   \code{"normal"} (the default; unit variance in both groups, means
#'   \code{true_value} and 0) or \code{"skewed"} (the same means and
#'   variance, with a standardized chi square distribution on 3 degrees of
#'   freedom in place of the normal).
#' @param variance For \code{seq_mr_cv_sensitivity()}, the variance
#'   estimator passed to \code{\link{seq_mr_cv}}: \code{"corrected"} (the
#'   default) or \code{"published"}.
#' @param G The number of Monte Carlo replications: a whole number of at
#'   least 2. The default is 1000. The examples use far fewer so that they
#'   run quickly; a reported analysis deserves at least the default.
#' @param seed A single whole number that \code{set.seed()} accepts, to
#'   make the result reproducible, or \code{NULL} (the default) to use the
#'   current state of the random number generator. When a seed is given, the
#'   generator's state is restored afterward.
#' @param m_0 The smallest pilot size allowed (per group for the
#'   standardized mean difference): a whole number of at least 4, the
#'   default.
#' @param gamma The tuning constant of the correction term, strictly between
#'   0 and 1/2; default 0.49.
#' @param max_n A cap on the sample size (on each group's size for the
#'   standardized mean difference): a replication that reaches it is
#'   stopped there. A whole number no smaller than the pilot size. The cap
#'   is exact: the last batch is shortened if needed. The default is large;
#'   see "Run time" in Details before using it with constants whose optimal
#'   sample size is in the thousands.
#' @param batch The number of observations added between checks of the
#'   stopping rule: a whole number of at least 1. The default, 1, is the
#'   purely sequential procedure of the articles. For the standardized mean
#'   difference it is the number added per group with one cost, and the
#'   number added to the less costly group with two costs.
#'
#' @return An object of class \code{c("smsd_sensitivity", "data.frame")}
#'   (see \code{\link{smsd_sensitivity_methods}}): a data frame with
#'   columns \code{term} (character) and \code{value} (numeric), one row
#'   per summary, that prints in fixed notation. The values are stored
#'   unrounded. For \code{seq_mr_cv_sensitivity()} the rows are
#'   \describe{
#'     \item{\code{true_value}, \code{loss_A}, \code{cost},
#'       \code{replications}}{the inputs;}
#'     \item{\code{optimal_n}, \code{minimum_risk}}{the benchmark: the best
#'       fixed sample size and the risk there;}
#'     \item{\code{mean_n}, \code{median_n}, \code{sd_n}}{the mean, median,
#'       and standard deviation of the final sample size;}
#'     \item{\code{ratio_n}}{\code{mean_n} divided by \code{optimal_n};}
#'     \item{\code{mean_estimate}}{the mean of the final estimates;}
#'     \item{\code{mean_risk}}{the mean realized loss;}
#'     \item{\code{risk_efficiency}}{\code{mean_risk} divided by
#'       \code{minimum_risk};}
#'     \item{\code{se_risk_efficiency}}{the Monte Carlo standard error of
#'       \code{risk_efficiency};}
#'     \item{\code{regret}}{\code{mean_risk} minus \code{minimum_risk}, in
#'       the currency of \code{c};}
#'     \item{\code{pct_hit_max_n}}{the percentage (0 to 100) of replications
#'       stopped by \code{max_n}.}
#'   }
#'   For \code{seq_mr_smd_sensitivity()} the rows are the same, except that
#'   \code{cost} is \code{cost_1} and \code{cost_2}, and each sample size
#'   row is reported for both groups: \code{optimal_n1}, \code{optimal_n2},
#'   \code{mean_n1}, \code{mean_n2}, \code{median_n1}, \code{median_n2},
#'   \code{sd_n1}, \code{sd_n2}, \code{ratio_n1}, and \code{ratio_n2}. A
#'   replication stopped by \code{max_n} whose final data give no estimate
#'   (a sample mean that is not positive, or constant data) is left out of
#'   \code{mean_estimate} and the four risk rows, which are \code{NA} when
#'   every replication is left out (\code{pct_hit_max_n} is then 100).
#'
#' @references
#' Chattopadhyay, B., & Kelley, K. (2016). Estimation of the coefficient
#' of variation with minimum risk: A sequential method for minimizing
#' sampling error and study cost. \emph{Multivariate Behavioral Research},
#' \emph{51}(5), 627--648. \doi{10.1080/00273171.2016.1203279}
#'
#' Chattopadhyay, B., & Kelley, K. (2017). Estimating the standardized mean
#' difference with minimum risk: Maximizing accuracy and minimizing cost
#' with sequential estimation. \emph{Psychological Methods}, \emph{22}(1),
#' 94--113. \doi{10.1037/met0000089}
#'
#' Kelley, K., Darku, F. B., & Chattopadhyay, B. (2018). Accuracy in
#' parameter estimation for a general class of effect sizes: A sequential
#' approach. \emph{Psychological Methods}, \emph{23}(2), 226--243.
#' \doi{10.1037/met0000127}
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @seealso \code{\link{seq_mr_cv}} and \code{\link{seq_mr_smd}}, the
#'   procedures simulated. The accuracy in parameter estimation
#'   counterparts: \code{\link{seq_aipe_sensitivity}}. How the result
#'   prints: \code{\link{smsd_sensitivity_methods}}.
#'
#' @name seq_mr_sensitivity
NULL

#' @rdname seq_mr_sensitivity
#' @examples
#' # A survey team prices squared error at A = 2,000 and pays 1 per
#' # respondent. If household spending has a coefficient of variation of
#' # 0.5, how many respondents will the procedure take?
#' out <- seq_mr_cv_sensitivity(true_value = 0.5, A = 2000, c = 1, G = 100,
#'                              seed = 113)
#' out
#' out$value[out$term == "ratio_n"]
#'
#' @family seq_mr_cv
#' @export
seq_mr_cv_sensitivity <- function(true_value, A, c,
                                  distribution = c("gamma", "lognormal",
                                                   "normal"),
                                  variance = c("corrected", "published"),
                                  G = 1000, seed = NULL, m_0 = 4,
                                  gamma = 0.49, max_n = 1e5L, batch = 1L) {
  kappa <- .check_positive(true_value, "true_value")
  A <- .check_positive(A, "A")
  c <- .check_positive(c, "c")
  distribution <- .check_choice(distribution,
                                c("gamma", "lognormal", "normal"),
                                "distribution")
  variance <- .check_choice(variance, c("corrected", "published"),
                            "variance")
  seed <- .check_seed(seed)
  gamma <- .mrpe_gamma(gamma)
  pilot <- seq_mr_cv_pilot(A = A, c = c, m_0 = m_0, gamma = gamma)
  sim <- .mrpe_sim_args(G, batch, max_n, pilot)

  mechanism <- .mr_cv_mechanism(distribution, kappa)
  optimal_n <- sqrt(A * mechanism$xi2 / c)
  minimum_risk <- 2 * c * optimal_n
  .mrpe_check_benchmark(c(optimal_n = optimal_n, minimum_risk = minimum_risk),
                        list(true_value = kappa, A = A, c = c))

  runs <- .with_seed(seed, function() {
    vapply(seq_len(sim$G), function(g) {
      .mr_cv_run_one(mechanism$draw, pilot, sim$batch, sim$max_n, A, c,
                     gamma, variance)
    }, numeric(3))
  })
  n <- runs["n", ]
  estimate <- runs["estimate", ]
  loss <- A * (estimate - kappa)^2 + c * n
  risk <- .mrpe_risk_rows(estimate, loss, minimum_risk)

  new_smsd_sensitivity(
    term = c("true_value", "loss_A", "cost", "replications",
             "optimal_n", "minimum_risk",
             "mean_n", "median_n", "sd_n", "ratio_n",
             "mean_estimate",
             "mean_risk", "risk_efficiency", "se_risk_efficiency", "regret",
             "pct_hit_max_n"),
    value = c(kappa, A, c, sim$G,
              optimal_n, minimum_risk,
              mean(n), stats::median(n), stats::sd(n), mean(n) / optimal_n,
              risk[["mean_estimate"]],
              risk[["mean_risk"]], risk[["risk_efficiency"]],
              risk[["se_risk_efficiency"]], risk[["regret"]],
              100 * mean(runs["capped", ]))
  )
}

#' @rdname seq_mr_sensitivity
#' @examples
#' # A training study prices squared error at A = 500 and pays 2 per
#' # participant in each group. For a true standardized mean difference of
#' # 0.5, checking after every pair of participants:
#' seq_mr_smd_sensitivity(true_value = 0.5, A = 500, c = 2, G = 100,
#'                        seed = 113)
#'
#' # The control group costs twice as much per participant, and the rule is
#' # checked after every fifth participant in the less costly group.
#' seq_mr_smd_sensitivity(true_value = 0.5, A = 500, c = c(2, 4), G = 100,
#'                        seed = 113, batch = 5)
#'
#' @family seq_mr_smd
#' @export
seq_mr_smd_sensitivity <- function(true_value, A, c,
                                   distribution = c("normal", "skewed"),
                                   G = 1000, seed = NULL, m_0 = 4,
                                   gamma = 0.49, max_n = 1e5L, batch = 1L) {
  delta <- .check_number(true_value, "true_value")
  # The limit of the accuracy in parameter estimation sibling
  # (seq_aipe_smd_sensitivity()), with its message: beside a larger value the
  # unit variance of the generated data is lost to rounding.
  .sensitivity_check_location(delta)
  A <- .check_positive(A, "A")
  cost <- .smd_mr_cost(c)
  c1 <- cost[1L]
  c2 <- cost[2L]
  distribution <- .check_choice(distribution, c("normal", "skewed"),
                                "distribution")
  seed <- .check_seed(seed)
  m_0 <- .check_count(m_0, "m_0", min = 4L)
  gamma <- .mrpe_gamma(gamma)
  pilot <- .smd_mr_pilot(A, c1, c2, gamma, m_0)
  sim <- .mrpe_sim_args(G, batch, max_n, max(pilot))

  if (distribution == "normal") {
    draw_1 <- function(k) stats::rnorm(k, mean = delta, sd = 1)
    draw_2 <- function(k) stats::rnorm(k, mean = 0, sd = 1)
    kurtosis <- 3
  } else {
    draw_1 <- function(k) .skewed_std(k) + delta
    draw_2 <- function(k) .skewed_std(k)
    kurtosis <- 7    # standardized chi square on 3 df: 3 + 12 / 3
  }
  # Large sample mean squared error of d: 1/n1 + 1/n2 +
  # delta^2 (kurtosis - 1) / (4 (n1 + n2)); for normal data this is
  # Equation 7 of Chattopadhyay and Kelley (2017). The benchmark is the
  # article's optimum at the squared effect delta^2 (kurtosis - 1) / 2.
  optimum <- .smd_mr_optimum(delta^2 * (kurtosis - 1) / 2, A, c1, c2)
  minimum_risk <- optimum[3L]
  .mrpe_check_benchmark(c(optimal_n1 = optimum[[1L]],
                          optimal_n2 = optimum[[2L]],
                          minimum_risk = optimum[[3L]]),
                        list(true_value = delta, A = A, c1 = c1, c2 = c2))

  runs <- .with_seed(seed, function() {
    vapply(seq_len(sim$G), function(g) {
      .smd_mr_run_one(draw_1, draw_2, pilot, sim$batch, sim$max_n, A, c1,
                      c2, gamma)
    }, numeric(4))
  })
  n1 <- runs["n1", ]
  n2 <- runs["n2", ]
  estimate <- runs["estimate", ]
  loss <- A * (estimate - delta)^2 + c1 * n1 + c2 * n2
  risk <- .mrpe_risk_rows(estimate, loss, minimum_risk)

  new_smsd_sensitivity(
    term = c("true_value", "loss_A", "cost_1", "cost_2", "replications",
             "optimal_n1", "optimal_n2", "minimum_risk",
             "mean_n1", "mean_n2", "median_n1", "median_n2",
             "sd_n1", "sd_n2", "ratio_n1", "ratio_n2",
             "mean_estimate",
             "mean_risk", "risk_efficiency", "se_risk_efficiency", "regret",
             "pct_hit_max_n"),
    value = c(delta, A, c1, c2, sim$G,
              optimum[1L], optimum[2L], minimum_risk,
              mean(n1), mean(n2), stats::median(n1), stats::median(n2),
              stats::sd(n1), stats::sd(n2),
              mean(n1) / optimum[1L], mean(n2) / optimum[2L],
              risk[["mean_estimate"]],
              risk[["mean_risk"]], risk[["risk_efficiency"]],
              risk[["se_risk_efficiency"]], risk[["regret"]],
              100 * mean(runs["capped", ]))
  )
}

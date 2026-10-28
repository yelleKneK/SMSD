# Monte Carlo sensitivity analysis for the budget constrained sequential
# procedure. The engine runs on the same internal helpers as seq_budget()
# and seq_budget_next() (.budget_state() and .budget_pick() in
# R/seq_budget.R), so a replication is, step for step, the loop a user
# writes with those two functions; tests/testthat/test-seq_budget.R checks
# that it is.

# One observation scheme for both the engine and the tests: k observations
# from a population with the given mean and standard deviation.
#' @noRd
.budget_draw <- function(k, mean, sd, distribution) {
  z <- if (distribution == "normal") stats::rnorm(k) else .skewed_std(k)
  mean + sd * z
}

# The populations must be ones for which the package's tolerance for "no
# variation" plays no part in the simulated studies. Scores are generated as
# mean + sd * Z, and every check in the package treats scores that agree to
# ten significant digits as showing no variation (.is_constant()). Beside a
# mean M standard deviations from zero, that is agreement within
# M * .smsd_constant_tol standard deviations, and the two scores of a pilot
# of 2 fall that close together with probability about 0.56 times that width
# (the density at zero of the difference of two standard normal draws is
# 1 / sqrt(pi)). The limit is the mean at which the width is 1e-4 of a
# standard deviation, that is, 1e6 standard deviations: about 6 pilots of 2
# in 100,000 are then affected, and far fewer with a larger pilot. In 10,000
# simulated studies of three groups with a pilot of 2, each compared with the
# same study at means of zero, 1 took another path at 1e6, 12 at 1e7, 158 at
# 1e8 (the coverage of the 95 percent interval moved from 0.901 to 0.909),
# and 1,618 at 1e9 (0.942); with a pilot of 5 none did up to 1e9. From about
# 1e11 on no group ever counts as varying, so every study is sampled by
# rounds and the summary would not say so. A group with a weight of zero is
# exempt, because its variance is never used.
#' @noRd
.budget_check_location <- function(means, sds, w) {
  limit <- 1e-4 / .smsd_constant_tol
  bad <- which(w != 0 & abs(means) >= limit * sds)
  if (!length(bad)) {
    return(invisible(TRUE))
  }
  i <- bad[1L]
  .smsd_stop("'true_means' must each be less than ", format(limit),
             " standard deviations ('true_sds') from zero in a group with ",
             "a nonzero weight; entry ", i, " is ", format(means[i]),
             " with a standard deviation of ", format(sds[i]),
             if (length(bad) > 1L) {
               paste0(" (", length(bad), " entries in all are too far)")
             },
             ". The package treats scores that agree to ten significant ",
             "digits as showing no variation, and beside a mean that large, ",
             "scores a small fraction of a standard deviation apart agree ",
             "that closely, so the simulated studies would depend on where ",
             "the scores are located. Subtract a common constant from the ",
             "means; a contrast whose weights sum to zero is unchanged by ",
             "it.")
}

# The largest study seq_budget_sensitivity() simulates, as the number of
# observations the budget pays for in the cheapest group (no study can be
# larger). A replication adds one observation at a time and recomputes a
# standard deviation at each step (.budget_sim_one()), so its time grows
# with the square of its size: at this limit one replication takes minutes
# (the Details of ?seq_budget_sensitivity give run times). seq_budget()
# accepts a budget that pays for up to about two billion observations
# (.budget_check_size()), and a simulation of a study of that size would
# not finish. The limit is far beyond any study the procedure is meant for.
#' @noRd
.budget_sens_max_n <- 1e6

# Stop, before anything is simulated, when the budget pays for more
# observations than that.
#' @noRd
.budget_sens_check_size <- function(budget, cost_per_group) {
  most <- floor(budget / min(cost_per_group))
  if (most > .budget_sens_max_n) {
    .smsd_stop("The budget pays for ", format(most, scientific = FALSE),
               " observations in the cheapest group, and ",
               "seq_budget_sensitivity() simulates studies of at most ",
               format(.budget_sens_max_n, scientific = FALSE),
               ": a study is simulated one observation at a time, and the ",
               "time for one replication grows with the square of its ",
               "size. Simulate a smaller budget to see how the procedure ",
               "behaves, and use seq_budget_optimal() for the allocation ",
               "at the full budget.")
  }
  invisible(TRUE)
}

# One replication: pilot, then one observation at a time to the group
# .budget_pick() names, until it names none. Returns the final sizes, the
# total cost, the contrast estimate, its estimated variance, and whether
# every group met its criterion.
#' @noRd
.budget_sim_one <- function(means, sds, w, a, budget, m_0, distribution,
                            first_passage, order) {
  K <- length(w)
  x <- lapply(seq_len(K), function(i) {
    .budget_draw(m_0, means[i], sds[i], distribution)
  })
  n <- rep.int(m_0, K)
  s <- vapply(x, stats::sd, numeric(1))
  varies <- vapply(x, .budget_varies, logical(1))
  estimable <- vapply(seq_len(K), function(i) {
    .budget_estimable(varies[i], s[i])
  }, logical(1))
  closed <- if (first_passage) logical(K) else NULL
  repeat {
    state <- .budget_state(n, s, estimable, w, a, budget, closed)
    if (first_passage) closed <- state$met
    g <- .budget_pick(n, s, state$target, state$met, state$not_estimable,
                      state$affordable, order)
    if (is.na(g)) break
    x[[g]] <- c(x[[g]], .budget_draw(1L, means[g], sds[g], distribution))
    n[g] <- n[g] + 1L
    s[g] <- stats::sd(x[[g]])
    # A group that varies still varies when an observation is added: its
    # range can only grow, and by at least as much as its largest absolute
    # value does. So only a group that has been constant so far is looked
    # at again, which keeps this step from scanning every group's data.
    if (!varies[g]) varies[g] <- .budget_varies(x[[g]])
    estimable[g] <- .budget_estimable(varies[g], s[g])
  }
  c(n,
    state$total_cost,
    sum(w * vapply(x, mean, numeric(1))),
    sum(w^2 * s^2 / n),
    as.numeric(state$rule_met))
}

# The arguments of seq_budget_sensitivity() that the optimal allocation is
# computed from, for the message when it cannot be represented.
#' @noRd
.budget_sens_arguments <-
  "'contrast_weights', 'true_sds', 'cost_per_group', and 'budget'"

# The tidy summary of a set of replications. `runs` has one column per
# replication, as .budget_sim_one() returns it: the K final sample sizes,
# the total cost, the contrast estimate, its estimated variance, and whether
# every group met its criterion (1 or 0). `margin` is NULL when no
# noninferiority test was asked for.
#' @noRd
.budget_summarize <- function(runs, means, sds, w, a, budget, conf_level,
                              margin, alpha_level) {
  K <- length(w)
  final_n <- runs[seq_len(K), , drop = FALSE]
  cost <- runs[K + 1L, ]
  estimate <- runs[K + 2L, ]
  variance <- runs[K + 3L, ]
  rule_met <- runs[K + 4L, ]

  true_contrast <- sum(w * means)
  optimal <- .budget_optimal(w, sds, a, budget, .budget_sens_arguments)
  # The same limits seq_budget() and the decision function compute.
  interval <- .budget_limits(estimate, variance, (1 - conf_level) / 2)
  covered <- interval$lower <= true_contrast &
    true_contrast <= interval$upper
  if (is.null(margin)) {
    margin_out <- NA_real_
    alpha_out <- NA_real_
    rejection_rate <- NA_real_
  } else {
    margin_out <- margin
    alpha_out <- alpha_level
    rejection_rate <- mean(
      .budget_limits(estimate, variance, alpha_level)$lower > -margin
    )
  }

  new_smsd_sensitivity(
    term = c("true_contrast", "budget", "conf_level", "replications",
             "coverage", "mean_total_cost", "sd_total_cost",
             "pct_within_budget", "pct_rule_met",
             paste0("mean_n_", seq_len(K)),
             paste0("optimal_n_", seq_len(K)),
             "mean_estimated_variance", "mean_true_variance",
             "empirical_variance", "optimal_variance",
             "margin", "alpha_level", "noninferiority_rejection_rate"),
    value = c(true_contrast, budget, conf_level, ncol(runs),
              mean(covered), mean(cost), stats::sd(cost),
              100 * mean(cost <= budget * (1 + 1e-9)),
              100 * mean(rule_met),
              rowMeans(final_n),
              optimal$n_per_group,
              mean(variance), mean(colSums(w^2 * sds^2 / final_n)),
              mean((estimate - true_contrast)^2),
              optimal$minimum_variance,
              margin_out, alpha_out, rejection_rate)
  )
}

#' Monte Carlo Sensitivity Analysis for the Budget Constrained Procedure
#'
#' Runs the budget constrained sequential procedure of Chattopadhyay,
#' Bandyopadhyay, Kelley, and Padalunkal (2025) from pilot to stopping
#' \code{G} times, for population means, standard deviations, costs, and a
#' budget that the user supplies, and summarizes what the procedure
#' delivers: how the budget is spent, how close the final sample sizes and
#' the variance of the contrast estimator come to the optimal allocation of
#' \code{\link{seq_budget_optimal}}, the coverage of the confidence
#' interval, and, for a given margin, how often the noninferiority test
#' rejects. Use it before a study to see what a budget buys.
#'
#' @details
#' Each replication draws \code{m_0} observations per population, then adds
#' one observation at a time to the group that
#' \code{\link{seq_budget_next}} names, re-evaluating
#' \code{\link{seq_budget}} after each, until no group that still needs an
#' observation can be paid for. The total cost can therefore never pass the
#' budget; \code{pct_within_budget} reports this and is 100 by
#' construction.
#'
#' Observations are generated as \code{true_means[i] + true_sds[i] * Z},
#' where \eqn{Z} is standard normal (\code{distribution = "normal"}) or a
#' standardized right skewed chi square variable with 3 degrees of freedom
#' (\code{"skewed"}; skewness about 1.63).
#'
#' The interval and the test are asymptotic, and the size of the shortfall
#' for a design of interest is what this function reports. As a guide, in
#' runs of 10,000 replications (standard error about 0.002) over five
#' designs with a pilot of 10 and about 15 to 120 observations per group,
#' coverage of the 95 percent interval was 0.936 to 0.949 for normal
#' populations and 0.925 to 0.941 for skewed ones, and the rejection rate
#' of the noninferiority test at the margin (nominal 0.05) was
#' 0.050 to 0.061 for normal populations and 0.040 to 0.085 for skewed ones.
#' The default pilot of 5 does a little worse than that when the
#' populations are skewed and the budget is small, and a pilot of 2 does
#' much worse; \code{\link{seq_budget_pilot}} has the figures, and
#' \code{m_0} lets the comparison be made for the design at hand.
#'
#' \code{procedure = "first_passage"} is the package's reading of the
#' article's sequence, and differs from it in three details that the help
#' page of \code{\link{seq_budget_next}} sets out; the main one is that
#' every open group is judged after each observation, where the article
#' judges a population only after an observation drawn from it.
#'
#' Each step of a replication recomputes the standard deviation of the
#' group that received the observation, so the time for one replication
#' grows faster than the number of observations the budget buys:
#' on one core of the computer used to check the package, the default
#' \code{G = 1000} took about 6 seconds when a study ended with about 200
#' observations in all, about 2 minutes with about 3,200, and about 11
#' minutes with about 12,800. For a budget that buys thousands of
#' observations, start with a small \code{G}. A budget that pays for more
#' than one million observations in the cheapest group (\code{budget /
#' min(cost_per_group)}) is refused before anything is simulated: at that
#' limit, with costs of 38, 15, and 10, a study ended with about 415,000
#' observations and one replication took more than two minutes.
#'
#' @param true_means Numeric vector of length \emph{K} (at least 2): the
#'   population means. In a group with a nonzero weight the mean must be
#'   less than \code{1e6} standard deviations (\code{true_sds}) from zero.
#'   The package treats scores that agree to ten significant digits as
#'   showing no variation, and beside a larger mean the simulated studies
#'   would begin to depend on where the scores are located. Subtracting a
#'   common constant from the means leaves a contrast whose weights sum to
#'   zero unchanged.
#' @param true_sds Numeric vector of length \emph{K}: the population
#'   standard deviations, each positive.
#' @param contrast_weights Numeric vector of length \emph{K} giving the
#'   contrast coefficients, not all zero.
#' @param cost_per_group Numeric vector of length \emph{K} giving the cost
#'   of one observation from each population, each positive. The four
#'   vector arguments are matched by position; if more than one is named,
#'   the names must be the same and in the same order.
#' @param budget The total budget, pilot included: a single positive
#'   number. The pilot, which costs \code{m_0 * sum(cost_per_group)}, must
#'   not cost more than the budget, and the budget must not pay for more
#'   than one million observations in the cheapest group (see Details).
#' @param conf_level The confidence level of the two-sided interval whose
#'   coverage is reported (default 0.95).
#' @param distribution The shape of the populations: \code{"normal"}
#'   (default) or \code{"skewed"}; an unambiguous abbreviation is accepted.
#' @param margin Optional noninferiority margin: the size of the margin, a
#'   single positive number, as in
#'   \code{\link{decide_seq_budget_noninferiority}}. When supplied, the
#'   rejection rate of that test at this margin is reported. It is the
#'   power of the test when the true contrast is above \code{-margin} and
#'   its Type I error rate when the true contrast is at or below
#'   \code{-margin}.
#' @param alpha_level The one-sided level of the noninferiority test: a
#'   single number strictly between 0 and 0.5 (default 0.05, the default of
#'   the decision functions). It is used only when \code{margin} is
#'   supplied, but is checked either way.
#' @param G Number of Monte Carlo replications: a whole number of at least
#'   2 (default 1000).
#' @param seed Optional whole number that \code{set.seed()} accepts. When
#'   supplied the replications are reproducible and the random number
#'   generator state of the session is restored on exit; \code{NULL}
#'   (default) uses the current state.
#' @param m_0 Pilot sample size per population: a whole number of at least
#'   2 (default 5, as in \code{\link{seq_budget_pilot}}).
#' @param procedure Which reading of the procedure to run (an unambiguous
#'   abbreviation is accepted).
#'   \code{"reallocate"} (default) is the package's: every group is judged
#'   on the current data at each check and the next observation goes to the
#'   group furthest below its estimated optimal size. The other choice,
#'   \code{"first_passage"}, follows the article's rule: a group is closed
#'   for good the first time it meets its criterion, and the open groups
#'   are visited in rounds (\code{closed} in \code{\link{seq_budget}} with
#'   \code{order = "rounds"} in \code{\link{seq_budget_next}}).
#'
#' @return An object of class \code{c("smsd_sensitivity", "data.frame")}:
#'   a data frame with columns \code{term} and \code{value}, one row per
#'   quantity, that prints in fixed notation (see
#'   \code{\link{smsd_sensitivity_methods}}) and can otherwise be used as
#'   any data frame. With \emph{K} populations the rows are, in order:
#'   \describe{
#'     \item{\code{true_contrast}}{the population contrast,
#'       \code{sum(contrast_weights * true_means)}.}
#'     \item{\code{budget}, \code{conf_level}, \code{replications}}{the
#'       inputs.}
#'     \item{\code{coverage}}{the proportion (0 to 1) of replications whose
#'       interval contains the true contrast.}
#'     \item{\code{mean_total_cost}, \code{sd_total_cost}}{the mean and
#'       standard deviation of the final total cost.}
#'     \item{\code{pct_within_budget}}{the percentage (0 to 100) of
#'       replications whose final cost is no more than the budget.}
#'     \item{\code{pct_rule_met}}{the percentage (0 to 100) of replications
#'       that ended with every group having met its criterion; the others
#'       ended because the budget could not pay for another needed
#'       observation.}
#'     \item{\code{mean_n_1}, ..., \code{mean_n_K}}{the mean final sample
#'       size of each group.}
#'     \item{\code{optimal_n_1}, ..., \code{optimal_n_K}}{for comparison,
#'       the optimal sample sizes for the true standard deviations, as
#'       \code{\link{seq_budget_optimal}} returns them.}
#'     \item{\code{mean_estimated_variance}}{the mean over replications of
#'       the estimated variance of the contrast estimator at stopping (the
#'       \code{variance} row of the final check).}
#'     \item{\code{mean_true_variance}}{the mean over replications of the
#'       true variance of the contrast estimator at the final sample sizes,
#'       \code{sum(contrast_weights^2 * true_sds^2 / n)}. It is never below
#'       \code{optimal_variance}; the excess is the price of not knowing
#'       the standard deviations in advance.}
#'     \item{\code{empirical_variance}}{the mean squared deviation of the
#'       contrast estimate from the true contrast over replications: the
#'       variance the procedure actually delivered. It has more Monte Carlo
#'       error than the two rows above.}
#'     \item{\code{optimal_variance}}{the smallest variance the budget
#'       allows with known standard deviations, \eqn{V_o}.}
#'     \item{\code{margin}, \code{alpha_level}}{the inputs of the
#'       noninferiority test (\code{NA} when no margin was supplied).}
#'     \item{\code{noninferiority_rejection_rate}}{the proportion (0 to 1)
#'       of replications in which inferiority was rejected (\code{NA} when
#'       no margin was supplied).}
#'   }
#'
#' @author
#' Ken Kelley \email{kkelley@@nd.edu},
#' Bhargab Chattopadhyay \email{bhargab@@iiitvadodara.ac.in}
#'
#' @references
#' Chattopadhyay, B., Bandyopadhyay, T., Kelley, K., & Padalunkal, J. J.
#' (2025). A sequential approach for noninferiority or equivalence of a
#' linear contrast under cost constraints. \emph{Psychological Methods},
#' \emph{30}(2), 425--439. \doi{10.1037/met0000570}
#'
#' @examples
#' # A new training format (group 1) against the average of two existing
#' # ones, with unequal costs and a budget of 1,500. G = 100 keeps the
#' # example quick; a reported analysis deserves the default G = 1000 or
#' # more.
#' seq_budget_sensitivity(
#'   true_means = c(50, 50, 50),
#'   true_sds = c(9, 6, 6),
#'   contrast_weights = c(1, -0.5, -0.5),
#'   cost_per_group = c(38, 15, 10),
#'   budget = 1500,
#'   margin = 5,
#'   G = 100,
#'   seed = 113,
#'   m_0 = 10
#' )
#'
#' @seealso \code{\link{seq_budget}}, \code{\link{seq_budget_next}},
#'   \code{\link{seq_budget_optimal}},
#'   \code{\link{decide_seq_budget_noninferiority}}. In \pkg{DMAR},
#'   \code{\link[DMAR]{ss_seq_c_sensitivity}} is the corresponding study
#'   for the sequential fixed-width interval for a contrast.
#' @family seq_budget
#' @export
seq_budget_sensitivity <- function(true_means, true_sds, contrast_weights,
                                   cost_per_group, budget,
                                   conf_level = 0.95,
                                   distribution = c("normal", "skewed"),
                                   margin = NULL, alpha_level = 0.05,
                                   G = 1000, seed = NULL, m_0 = 5,
                                   procedure = c("reallocate",
                                                 "first_passage")) {
  distribution <- .check_choice(distribution, c("normal", "skewed"),
                                "distribution")
  procedure <- .check_choice(procedure, c("reallocate", "first_passage"),
                             "procedure")
  w <- .budget_check_weights(contrast_weights)
  K <- length(w)
  .budget_check_K(K, "contrast_weights")
  per_weight <- "one per entry of 'contrast_weights'"
  means <- .check_vector(true_means, "true_means",
                         "one mean per population", len = K,
                         len_what = per_weight)
  sds <- .check_vector(true_sds, "true_sds",
                       "one standard deviation per population", len = K,
                       len_what = per_weight, positive = TRUE)
  a <- .check_vector(cost_per_group, "cost_per_group",
                     "the cost of one observation in each population",
                     len = K, len_what = per_weight, positive = TRUE)
  .budget_check_location(means, sds, w)
  budget <- .check_positive(budget, "budget")
  conf_level <- .check_conf_level(conf_level)
  .budget_check_size(budget, a)
  .budget_sens_check_size(budget, a)
  .budget_check_names(list(true_means = names(true_means),
                           true_sds = names(true_sds),
                           contrast_weights = names(contrast_weights),
                           cost_per_group = names(cost_per_group)))
  if (!is.null(margin)) {
    margin <- .check_positive(margin, "margin")
  }
  alpha_level <- .budget_check_alpha(alpha_level)
  G <- .check_count(G, "G", min = 2L)
  seed <- .check_seed(seed)
  m_0 <- .check_count(m_0, "m_0", min = 2L)
  pilot_cost <- m_0 * sum(a)
  if (pilot_cost > budget * (1 + 1e-9)) {
    .smsd_stop("The pilot alone (m_0 = ", m_0, " per group) costs ",
               format(pilot_cost), ", which is more than the budget of ",
               format(budget), ". Lower 'm_0' or raise 'budget'.")
  }
  # The summary compares the studies with the optimal allocation. If that
  # cannot be represented, say so now and not after G studies have been run.
  .budget_optimal(w, sds, a, budget, .budget_sens_arguments)

  first_passage <- procedure == "first_passage"
  order <- if (first_passage) "rounds" else "shortfall"
  runs <- .with_seed(seed, function() {
    replicate(G, .budget_sim_one(means, sds, w, a, budget, m_0,
                                 distribution, first_passage, order))
  })
  # Scores so large that a sample variance or the contrast overflows: the
  # check a user would run on such data stops with this message, so the
  # summary does not describe them either.
  if (!all(is.finite(runs))) {
    .budget_stop_scale(
      "The contrast estimate or its variance in a simulated study",
      "'true_means', 'true_sds', and 'contrast_weights'")
  }
  .budget_summarize(runs, means, sds, w, a, budget, conf_level, margin,
                    alpha_level)
}

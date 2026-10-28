# The budget constrained sequential procedure for a linear contrast of K
# population means (Chattopadhyay, Bandyopadhyay, Kelley, and Padalunkal,
# 2025): a priori allocation, pilot size, the check, the choice of the next
# observation, and the noninferiority, superiority, inferiority, and
# equivalence decisions. The Monte Carlo sibling is in
# R/sensitivity_budget.R and runs on the internal helpers defined here.

# ---- Internal helpers ---------------------------------------------------

# The method label carried by every check this family returns. The decision
# functions and seq_budget_next() recognize their input by it.
#' @noRd
.budget_method <- "budget constrained K-population contrast"

# Arguments with one entry per group are matched by position. When two of
# them carry names, the names must agree, in order, so that a vector named
# in another order is not applied to the wrong groups without a word.
# `labels` is a named list of the names() of each argument.
#' @noRd
.budget_check_names <- function(labels) {
  usable <- vapply(labels, function(nm) {
    !is.null(nm) && !anyNA(nm) && all(nzchar(nm))
  }, logical(1))
  if (sum(usable) < 2L) {
    return(invisible(TRUE))
  }
  show <- function(nm) {
    if (length(nm) > 8L) nm <- c(nm[1:8], "...")
    paste(nm, collapse = ", ")
  }
  first <- which(usable)[1L]
  for (i in which(usable)[-1L]) {
    if (!identical(labels[[i]], labels[[first]])) {
      .smsd_stop("'", names(labels)[i], "' is named (", show(labels[[i]]),
                 ") but '", names(labels)[first], "' is named (",
                 show(labels[[first]]), "). Entries are matched to groups ",
                 "by position, not by name, so the names must agree in ",
                 "order. Put '", names(labels)[i], "' in the order of '",
                 names(labels)[first], "', or remove its names with ",
                 "unname().")
    }
  }
  invisible(TRUE)
}

# Stop when a quantity cannot be represented because a product overflowed
# or underflowed. Ordinary data never reach this.
#' @noRd
.budget_stop_scale <- function(quantity, arguments) {
  .smsd_stop(quantity, " cannot be computed: the values in ", arguments,
             " are too large or too small in magnitude to combine. Rescale ",
             "them, for example to another unit, and check again.")
}

# Contrast weights: finite, one per group, not all zero. `len_what` says in
# the message what the length must match.
#' @noRd
.budget_check_weights <- function(x, K = NULL, len_what = NULL) {
  x <- .check_vector(x, "contrast_weights", "one contrast weight per group",
                     len = K, len_what = len_what)
  if (all(x == 0)) {
    .smsd_stop("'contrast_weights' must have at least one nonzero entry; ",
               "all ", length(x), " are zero, so there is no contrast to ",
               "estimate.")
  }
  x
}

# The number of groups implied by a vector argument: at least two.
#' @noRd
.budget_check_K <- function(K, name) {
  if (K < 2L) {
    .smsd_stop("The procedure compares at least 2 populations; '", name,
               "' describes ", K, ".")
  }
  invisible(K)
}

# The advice that closes the message when a sample size implied by the budget
# cannot be represented.
#' @noRd
.budget_size_advice <- "Check 'budget' and 'cost_per_group'."

# The budget must buy a representable number of observations. No criterion
# can exceed what the whole budget buys in the cheapest group, so after this
# check every criterion is representable too.
#' @noRd
.budget_check_size <- function(budget, cost_per_group) {
  .as_sample_size(
    budget / min(cost_per_group),
    what = paste("number of observations that the budget pays for in the",
                 "cheapest group"),
    advice = .budget_size_advice
  )
  invisible(TRUE)
}

# The list of per group data vectors, cleaned.
#' @noRd
.budget_check_groups <- function(groups) {
  if (!is.list(groups) || is.data.frame(groups)) {
    .smsd_stop("'groups' must be a list with one numeric vector per group; ",
               "got ", .smsd_describe(groups), ".")
  }
  .budget_check_K(length(groups), "groups")
  lapply(seq_along(groups), function(i) {
    .check_data(groups[[i]], sprintf("groups[[%d]]", i), min_n = 2L,
                na.rm = NULL)
  })
}

# The groups already closed under the first passage rule: NULL, or one TRUE
# or FALSE per group.
#' @noRd
.budget_check_closed <- function(closed, K) {
  if (is.null(closed)) {
    return(NULL)
  }
  if (!is.logical(closed) || !is.null(dim(closed)) || length(closed) != K ||
      anyNA(closed)) {
    .smsd_stop("'closed' must be NULL or a logical vector with one TRUE or ",
               "FALSE per group (K = ", K, ", none missing), for example ",
               "attr(previous_check, \"met_per_group\"); got ",
               .smsd_describe(closed), ".")
  }
  closed
}

# Do the observations of a group vary, beyond rounding error? The opposite of
# .is_constant(), the one tolerance of the whole package.
#' @noRd
.budget_varies <- function(x) {
  !.is_constant(x)
}

# Can the variance of one group be used? Not when the group does not vary
# (whatever variance was computed is then zero or noise), and not when its
# sample variance is zero or not finite (.var_unusable(), the rule of the
# whole package; it is applied to the standard deviation `s`, which is zero
# or not finite exactly when the variance is). Only the zero case reaches a
# result: a sample variance that is not finite means the scores overflowed,
# and then seq_budget() and seq_budget_sensitivity() stop with the message
# of .budget_stop_scale() instead of returning anything.
#' @noRd
.budget_estimable <- function(varies, s) {
  varies && !.var_unusable(s)
}

# Can one more observation from each group be paid for? Costs need not be
# whole numbers, so the comparison carries a relative tolerance: three
# observations at 0.1 each are not "more than" a budget of 0.3.
#' @noRd
.budget_affordable <- function(total_cost, a, budget) {
  total_cost + a <= budget * (1 + 1e-9)
}

# The smallest whole number N with N >= target for each group, by the
# package's rule for a criterion sample size (.as_sample_size()): a target
# that is a whole number up to representation error (111.00000000000001; the
# rule forgives one part in 1e10) is not rounded up to the next one. (A plain
# loop: the shared helper takes one value at a time, and this runs at every
# step of a simulated study.)
#' @noRd
.budget_criterion <- function(target) {
  criterion <- integer(length(target))
  for (i in seq_along(target)) {
    criterion[i] <- .as_sample_size(target[i],
                                    what = paste("criterion for group", i),
                                    advice = .budget_size_advice)
  }
  criterion
}

# The state of the procedure at the current data. Everything the check, the
# choice of the next observation, and the Monte Carlo engine need is computed
# here, once, from the per group sample sizes and standard deviations.
#
#   n, s        current sample sizes and sample standard deviations
#   estimable   is each group's variance usable? (.budget_estimable())
#   w, a        contrast weights and unit costs
#   budget      the total budget
#   closed      NULL, or for the first passage rule the groups from which
#               sampling has already stopped
#' @noRd
.budget_state <- function(n, s, estimable, w, a, budget, closed = NULL) {
  weighted <- w != 0
  not_estimable <- weighted & !estimable
  # |c_i| s_i sqrt(a_i), the terms of the sum in Equations 12 and 18.
  s[!estimable] <- 0
  share <- abs(w) * s * sqrt(a)
  total_share <- sum(share)
  total_cost <- sum(a * n)
  if (!is.finite(total_share) || !is.finite(total_cost) ||
      any(share[weighted & estimable] <= 0)) {
    .budget_stop_scale(
      "The share of the budget that each group should receive",
      "'contrast_weights', 'cost_per_group', and the data")
  }
  # Equation 18's right hand side, A_0 |c_i| s_i / (sqrt(a_i) sum(...)),
  # computed as what the whole budget would buy in group i times group i's
  # share of the spending, so that very large or very small costs cancel
  # before they can overflow.
  target <- if (total_share > 0) {
    (budget / a) * (share / total_share)
  } else {
    numeric(length(n))
  }
  criterion <- .budget_criterion(target)
  # A group whose variance is unusable needs at least one more observation
  # before its share means anything.
  criterion[not_estimable] <- as.integer(n[not_estimable] + 1L)
  met <- n >= criterion
  if (!is.null(closed)) {
    met <- met | (closed & !not_estimable)
  }
  affordable <- .budget_affordable(total_cost, a, budget)
  open <- !met
  list(
    target = target,
    criterion = criterion,
    not_estimable = not_estimable,
    met = met,
    affordable = affordable,
    total_cost = total_cost,
    remaining = budget - total_cost,
    rule_met = !any(open),
    budget_exhausted = any(open) && !any(open & affordable),
    budget_exceeded = total_cost > budget * (1 + 1e-9),
    variance_unusable = any(not_estimable)
  )
}

# The group that receives the next observation, or NA when no group that
# still needs one can be paid for.
#
#   order = "shortfall": the group furthest below its estimated optimal size
#     (smallest n_i / target_i).
#   order = "rounds": the group with the fewest observations, which visits
#     the open groups one at a time, in rounds.
# Ties go to the larger standard deviation, then to the lower position.
#
# A group whose variance is not yet estimable has no target to be measured
# against, so under either order it takes its turn by rounds: it is next
# when no group whose variance can be estimated (target > 0) has fewer
# observations than it has, and otherwise only when no other group can be
# sampled. It is never sampled ahead of the others, so a group that stays
# constant cannot absorb the budget.
#' @noRd
.budget_pick <- function(n, s, target, met, not_estimable, affordable,
                         order) {
  candidates <- which(!met & affordable)
  if (!length(candidates)) {
    return(NA_integer_)
  }
  unknown <- candidates[not_estimable[candidates]]
  known <- candidates[!not_estimable[candidates]]
  if (length(unknown)) {
    fewest <- min(n[unknown])
    # An open group whose variance can be estimated has a positive target,
    # so `known` not empty implies that some target is positive.
    if (!length(known) || fewest <= min(n[target > 0])) {
      return(unknown[n[unknown] == fewest][1L])
    }
  }
  key <- if (order == "shortfall") {
    n[known] / target[known]
  } else {
    n[known]
  }
  # `known` is in increasing position, so after the two filters the first
  # element is the lowest position among those still tied.
  best <- known[key == min(key)]
  if (length(best) > 1L) {
    best <- best[s[best] == max(s[best])]
  }
  best[1L]
}

# "group 2", "groups 1 and 3", "groups 1, 2, and 4".
#' @noRd
.budget_group_label <- function(which_groups) {
  k <- length(which_groups)
  if (k == 1L) {
    return(paste("group", which_groups))
  }
  if (k == 2L) {
    return(paste("groups", which_groups[1L], "and", which_groups[2L]))
  }
  paste0("groups ", paste(which_groups[-k], collapse = ", "), ", and ",
         which_groups[k])
}

# An amount of money for a note. Display only: the amount is rounded to the
# number of decimal places that floating point error in a sum the size of
# the budget leaves meaningful (about thirteen significant digits of the
# budget), so that dust such as -1.1e-16 prints as 0 while a unit cost of
# 1200.4 against a budget of ten million keeps its decimals.
#' @noRd
.budget_money <- function(x, budget) {
  places <- max(0, floor(13 - log10(budget)))
  x <- round(x, places)
  if (x == 0) {
    return("0")
  }
  fixed <- abs(x) >= 1e-4 && abs(x) < 1e15
  format(x, digits = 15L, scientific = if (fixed) FALSE else NA)
}

# The plain-language notes printed under the stopping line.
#' @noRd
.budget_notes <- function(state, a, budget, first_passage) {
  notes <- character(0)
  open <- which(!state$met)
  flat <- which(state$not_estimable)
  if (state$budget_exceeded) {
    notes <- c(notes, paste0(
      "The data already cost ", .budget_money(state$total_cost, budget),
      ", which is more than the budget of ", .budget_money(budget, budget),
      ", so the budget constraint does not hold for this result. Raise the ",
      "budget or use a smaller pilot; once sampling is under way, choosing ",
      "each new observation with seq_budget_next() keeps the cost within ",
      "the budget."))
  } else if (state$budget_exhausted) {
    notes <- c(notes, paste0(
      "Sampling has ended because the budget cannot pay for another ",
      "observation: ", .budget_money(state$remaining, budget), " of ",
      .budget_money(budget, budget),
      " is left, and one more observation from a group still below its ",
      "criterion (", .budget_group_label(open), ") costs at least ",
      .budget_money(min(a[open]), budget), "."))
  } else if (state$rule_met && first_passage) {
    notes <- c(notes, paste0(
      "Every group has met its criterion at this check or an earlier one. ",
      "Under the first passage rule sampling from a group stops for good ",
      "the first time it does, so a current size can be below the ",
      "current criterion."))
  }
  if (state$variance_unusable) {
    flat_txt <- paste0(
      "The observations in ", .budget_group_label(flat), " are all equal, ",
      "or differ too little relative to their size for a variance to be ",
      "computed")
    if (state$budget_exhausted) {
      # The standard note advises collecting more data, which the budget no
      # longer allows; say the same thing with what holds instead.
      notes <- c(notes, paste0(
        flat_txt, ", so the variance estimate is not usable: the stopping ",
        "rule is reported as not met and the interval is too narrow to rely ",
        "on. More data from ",
        if (length(flat) == 1L) "that group" else "those groups",
        " would take a larger budget."))
    } else {
      notes <- c(notes, .note_variance_unusable(), paste0(
        flat_txt, ", so ",
        if (length(flat) == 1L) "its share" else "their shares",
        " of the budget cannot be estimated yet. seq_budget_next() names ",
        "such a group when no group whose variance can be estimated has ",
        "fewer observations."))
    }
  }
  notes
}

# Stop unless `check` came from seq_budget().
#' @noRd
.budget_require_check <- function(check) {
  ok <- inherits(check, "smsd_check") &&
    identical(attr(check, "method"), .budget_method) &&
    all(c("lower_limit", "contrast", "upper_limit", "variance") %in%
          check$term) &&
    !is.null(attr(check, "n_per_group"))
  if (!ok) {
    got <- if (inherits(check, "smsd_check")) {
      sprintf("a check from another procedure (\"%s\")",
              as.character(attr(check, "method"))[1L])
    } else {
      .smsd_describe(check)
    }
    .smsd_stop("'check' must be an object returned by seq_budget(); got ",
               got, ".")
  }
  invisible(check)
}

# A named term of a check.
#' @noRd
.budget_term <- function(check, term) {
  check$value[match(term, check$term)]
}

# The one-sided level of a decision function.
#' @noRd
.budget_check_alpha <- function(alpha_level) {
  .check_number(alpha_level, "alpha_level", lower = 0, upper = 0.5,
                lower_open = TRUE, upper_open = TRUE)
}

# A margin of a decision function: a single finite number on the scale of
# the contrast. `what` says, for the message, what the margin is when it was
# not supplied.
#' @noRd
.budget_check_margin <- function(margin, name, what) {
  if (missing(margin)) {
    .smsd_stop("'", name, "' must be supplied: ", what, ".")
  }
  .check_number(margin, name)
}

# The limits every decision is computed from: estimate -/+ z sqrt(variance)
# with z the upper alpha_level point of the standard normal distribution,
# from the contrast and variance rows of the check. The confidence limits
# stored on the check, which are two-sided at its conf_level, are never
# read. Refuses a check whose variance is unusable and warns when sampling
# has not ended.
#' @noRd
.budget_decision_limits <- function(check, alpha_level) {
  alpha_level <- .budget_check_alpha(alpha_level)
  if (isTRUE(attr(check, "variance_unusable"))) {
    .smsd_stop("No decision can be made from this check: the observations ",
               "in at least one group are all equal, so the variance ",
               "estimate is not usable and the variance of the contrast is ",
               "not estimable. Collect more data from that group if the ",
               "budget allows.")
  }
  if (!isTRUE(attr(check, "is_satisfied"))) {
    warning("Sampling has not ended at this check (the stopping rule is ",
            "not met), so the decision is provisional. The level of the ",
            "test is stated for one decision, made when sampling has ended.",
            call. = FALSE)
  }
  .budget_limits(.budget_term(check, "contrast"),
                 .budget_term(check, "variance"), alpha_level)
}

# The limits estimate -/+ z sqrt(variance), with z the upper `tail` point of
# the standard normal distribution: tail = (1 - conf_level) / 2 for the
# interval on a check, tail = alpha_level for a one-sided test. One function
# serves the check, the decisions, and the Monte Carlo summary, so the three
# cannot disagree in the last digit.
#' @noRd
.budget_limits <- function(estimate, variance, tail) {
  half <- stats::qnorm(tail, lower.tail = FALSE) * sqrt(variance)
  list(lower = estimate - half, upper = estimate + half)
}

# The optimal allocation of Theorem 1 (Equations 12 and 13) for inputs that
# have already been checked. `arguments` names, for the message, the
# arguments of the calling function that the values came from, so that
# seq_budget_optimal() and seq_budget_sensitivity() each speak of their own.
#' @noRd
.budget_optimal <- function(w, sigmas, a, budget, arguments) {
  # |c_i| sigma_i sqrt(a_i), the terms of the sum in Equations 12 and 13.
  share <- abs(w) * sigmas * sqrt(a)
  total_share <- sum(share)
  # Equation 12 as (what the whole budget would buy in group i) times
  # (group i's share of the spending).
  n_per_group <- (budget / a) * (share / total_share)
  minimum_variance <- total_share^2 / budget
  if (!all(is.finite(n_per_group)) || any(share[w != 0] <= 0) ||
      !is.finite(minimum_variance) || minimum_variance <= 0) {
    .budget_stop_scale("The optimal allocation", arguments)
  }
  list(
    n_per_group = n_per_group,
    minimum_variance = minimum_variance
  )
}


# ---- A priori allocation ------------------------------------------------

#' Optimal Sample Sizes for a Budget Constrained Contrast
#'
#' For planning before any data are collected: given assumed population
#' standard deviations, returns the sample sizes that estimate a linear
#' contrast of \emph{K} population means as precisely as a fixed budget
#' allows when an observation costs a different amount in each population.
#' This is Theorem 1 of Chattopadhyay, Bandyopadhyay, Kelley, and Padalunkal
#' (2025).
#'
#' Write the contrast as \eqn{\psi = \sum_{i=1}^{K} c_i \mu_i}, estimated
#' by \eqn{\hat{\psi} = \sum_{i=1}^{K} c_i \bar{X}_i} with variance
#' \eqn{\sum_{i=1}^{K} c_i^2 \sigma_i^2 / n_i} (Equations 3 to 5 of the
#' article). With a cost of \eqn{a_i} per observation from population
#' \eqn{i} and a budget \eqn{A_0}, the sample sizes that minimize that
#' variance subject to \eqn{\sum_{i=1}^{K} a_i n_i = A_0} are
#' (Equation 12)
#' \deqn{n_{i0} = \frac{A_0 \, |c_i| \, \sigma_i}
#'   {\sqrt{a_i} \, \sum_{l=1}^{K} |c_l| \, \sigma_l \, \sqrt{a_l}},
#' \quad i = 1, \ldots, K,}
#' and the minimized variance is (Equation 13)
#' \deqn{V_o = \frac{1}{A_0}
#'   \left( \sum_{i=1}^{K} |c_i| \, \sigma_i \, \sqrt{a_i} \right)^2.}
#' A population receives more observations when its scores vary more, when
#' its weight in the contrast is larger, and when it is cheaper to sample.
#' Minimizing the variance is what maximizes the power of the article's
#' noninferiority, superiority, and equivalence tests for a given budget.
#'
#' The sizes are returned as they come out of the formula, not rounded.
#' Rounding every size up would cost more than the budget: in the article's
#' example (three populations, 30 per observation, a budget of 10,000;
#' p. 430) each \eqn{n_{i0}} is 111.11, which holds when the products
#' \eqn{|c_i| \sigma_i} are equal, and 111 per group costs 9,990 where 112
#' per group would cost 10,080.
#'
#' In practice the standard deviations are not known. \code{\link{seq_budget}}
#' replaces them with the sample standard deviations and updates the
#' allocation as observations arrive, and
#' \code{\link{seq_budget_sensitivity}} shows by simulation how close that
#' procedure comes to the allocation returned here.
#'
#' The three vector arguments are matched by position. If more than one of
#' them is named, the names must be the same and in the same order;
#' otherwise the function stops, so that a vector named in another order is
#' not applied to the wrong populations.
#'
#' @param contrast_weights Numeric vector of length \emph{K} (at least 2)
#'   giving the contrast coefficients \eqn{c_i}, one per population, not all
#'   zero. The article takes them to sum to zero (a contrast); the
#'   allocation holds for any linear combination, so the sum is not checked.
#'   A population with a weight of zero is allocated no observations.
#' @param sigmas Numeric vector of length \emph{K} giving the assumed
#'   population standard deviations \eqn{\sigma_i}, each positive, in the
#'   order of the weights.
#' @param cost_per_group Numeric vector of length \emph{K} giving the cost
#'   \eqn{a_i} of one observation from each population, each positive, in
#'   the units of \code{budget}. The order is that of the weights.
#' @param budget The total budget \eqn{A_0}: a single positive number.
#'
#' @return A list with two components: \code{n_per_group}, a numeric vector
#'   of length \emph{K} with the optimal sample sizes \eqn{n_{i0}} (not
#'   rounded; their total cost equals \code{budget}), and
#'   \code{minimum_variance}, the minimized variance \eqn{V_o} of the
#'   contrast estimator.
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
#' Luh, W. M., & Guo, J. H. (2016). Sample size planning for the
#' noninferiority or equivalence of a linear contrast with cost
#' considerations. \emph{Psychological Methods}, \emph{21}(1), 13--34.
#' \doi{10.1037/met0000039}
#'
#' @examples
#' # Two new training formats (groups 1 and 3) against the current one
#' # (group 2): the weights 1, -2, 1 give twice the difference between the
#' # average of the new formats and the current format. An observation costs
#' # 38, 15, and 10 in the three groups, and the budget is 5,000 (the normal
#' # row of Table 1 in the article).
#' seq_budget_optimal(
#'   contrast_weights = c(1, -2, 1),
#'   sigmas = c(4.5, 0.9, 1),
#'   cost_per_group = c(38, 15, 10),
#'   budget = 5000
#' )
#'
#' # The article's example on p. 430: equal costs, and equal products of
#' # weight and standard deviation.
#' seq_budget_optimal(
#'   contrast_weights = c(-1, 0.5, 0.5),
#'   sigmas = c(1, 2, 2),
#'   cost_per_group = c(30, 30, 30),
#'   budget = 10000
#' )$n_per_group
#'
#' @seealso \code{\link{seq_budget}} for the sequential procedure when the
#'   standard deviations are unknown. In \pkg{DMAR},
#'   \code{\link[DMAR]{ss_seq_c}} reports the same cost optimal allocation
#'   for a sequential fixed-width interval.
#' @family seq_budget
#' @export
seq_budget_optimal <- function(contrast_weights, sigmas, cost_per_group,
                               budget) {
  w <- .budget_check_weights(contrast_weights)
  K <- length(w)
  .budget_check_K(K, "contrast_weights")
  sigmas_in <- sigmas
  sigmas <- .check_vector(sigmas, "sigmas",
                          "one standard deviation per population", len = K,
                          len_what = "one per entry of 'contrast_weights'",
                          positive = TRUE)
  a <- .check_vector(cost_per_group, "cost_per_group",
                     "the cost of one observation in each population",
                     len = K,
                     len_what = "one per entry of 'contrast_weights'",
                     positive = TRUE)
  budget <- .check_positive(budget, "budget")
  .budget_check_names(list(contrast_weights = names(contrast_weights),
                           sigmas = names(sigmas_in),
                           cost_per_group = names(cost_per_group)))

  .budget_optimal(
    w, sigmas, a, budget,
    "'contrast_weights', 'sigmas', 'cost_per_group', and 'budget'"
  )
}


# ---- Pilot --------------------------------------------------------------

#' Pilot Sample Size for the Budget Constrained Sequential Procedure
#'
#' Returns the number of observations to collect from every population
#' before \code{\link{seq_budget}} is evaluated for the first time. The
#' article gives no formula for the pilot size: it requires at least two
#' observations per population, so that a sample variance exists, and
#' leaves the rest to the researcher. This function checks that
#' requirement and returns the chosen size as an integer.
#'
#' Chattopadhyay, Bandyopadhyay, Kelley, and Padalunkal (2025, p. 431)
#' advise against the minimum of two, because variance estimates from two
#' observations are unstable, and against a pilot so large that it
#' oversamples; in their simulations (p. 432) the pilot was 10 per
#' population with a budget of 5,000. The default of 5 is this package's
#' choice. The pilot
#' is paid for from the same budget, so its cost,
#' \code{m_0 * sum(cost_per_group)}, must be well below the budget; a check
#' whose data already cost more than the budget says so.
#'
#' \strong{What the pilot size costs.} A group whose pilot standard
#' deviation happens to come out small meets its criterion at once, and it
#' is sampled again only if the other standard deviations fall. With a
#' very small pilot this is common: the group is left with two or three
#' observations, and the variance of the contrast, estimated from them,
#' comes out far too small. The table gives the coverage of the 95 percent
#' interval and, in parentheses, the rejection rate of the noninferiority
#' test at nominal level 0.05 with the true contrast exactly at the margin,
#' from 10,000 simulated studies per cell (standard error about 0.002).
#' The first design has the standard deviations, weights, and unit costs
#' of the normal row of the article's Table 1 (costs 38, 15, and 10; a
#' budget of 5,000, which buys about 40 to 100 observations per group);
#' the second has those of the example on the
#' \code{\link{seq_budget_sensitivity}} page (the same costs; a budget of
#' 1,500, about 15 to 28 per group). Each was run with normal and with
#' right skewed populations.
#' \tabular{llccc}{
#'   Budget \tab Populations \tab Pilot of 10 \tab Pilot of 5 \tab Pilot of 2 \cr
#'   5,000 \tab normal \tab 0.947 (0.055) \tab 0.949 (0.053) \tab 0.928 (0.063) \cr
#'   5,000 \tab skewed \tab 0.937 (0.040) \tab 0.939 (0.043) \tab 0.903 (0.060) \cr
#'   1,500 \tab normal \tab 0.938 (0.060) \tab 0.931 (0.060) \tab 0.896 (0.081) \cr
#'   1,500 \tab skewed \tab 0.924 (0.054) \tab 0.907 (0.075) \tab 0.852 (0.121)
#' }
#' With a pilot of 2 the estimated variance of the contrast averaged only
#' 0.59 to 0.76 of the true variance of the estimator, and the true
#' variance was 1.30 to 1.51 times the smallest the budget allows; with
#' pilots of 5 and 10 the two figures were 0.86 to 0.99 and 1.01 to 1.06.
#' In the pilot of 2 runs with normal populations, 6 percent (budget
#' 5,000) and 18 percent (budget 1,500) of the studies ended with a group
#' of two or three observations, and the interval covered the contrast in
#' 0.64 and 0.73 of those studies against 0.94 of the rest. A pilot of 5
#' or 10 costs little when the budget buys 40 or more observations per
#' group. When it buys fewer than about 30, or the populations may be
#' skewed, 10 is the safer choice, and a pilot of 2 should not be used.
#' \code{\link{seq_budget_sensitivity}} gives the figures for a design of
#' interest through its \code{m_0} argument.
#'
#' @param m_0 The pilot sample size per population: a whole number of at
#'   least 2. Default 5.
#'
#' @return A length one integer giving the pilot sample size per
#'   population.
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
#' seq_budget_pilot()
#' seq_budget_pilot(m_0 = 10)
#'
#' @family seq_budget
#' @export
seq_budget_pilot <- function(m_0 = 5) {
  .check_count(m_0, "m_0", min = 2L)
}


# ---- The check ----------------------------------------------------------

#' Sequential Check for a Budget Constrained K-Population Contrast
#'
#' Given the data collected so far from each of \emph{K} populations, the
#' cost of an observation in each, and a total budget, reports whether
#' sampling should continue under the sequential procedure of
#' Chattopadhyay, Bandyopadhyay, Kelley, and Padalunkal (2025), and returns
#' the current estimate of the contrast with its confidence interval. Use
#' it when a contrast of means has to be estimated as precisely as a fixed
#' budget allows and the population variances that determine the best
#' allocation are unknown. \code{\link{seq_budget_next}} says which group
#' the next observation should come from.
#'
#' @details
#' \strong{The criterion.} The contrast \eqn{\psi = \sum_{i=1}^{K} c_i
#' \mu_i} is estimated by \eqn{\hat{\psi} = \sum_{i=1}^{K} c_i \bar{X}_i},
#' whose variance is estimated by \eqn{\hat{V} = \sum_{i=1}^{K} c_i^2
#' s_i^2 / N_i} (Equation 5 of the article with sample variances). The
#' variance is smallest, for a budget \eqn{A_0} and unit costs \eqn{a_i},
#' at the allocation of \code{\link{seq_budget_optimal}}, which depends on
#' the unknown standard deviations. The sequential procedure substitutes
#' the current sample standard deviations: group \eqn{i} has met its
#' criterion when (Equation 18)
#' \deqn{N_i \ge \frac{A_0 \, |c_i| \, s_i}
#'   {\sqrt{a_i} \, \sum_{l=1}^{K} |c_l| \, s_l \, \sqrt{a_l}}.}
#' The right hand side is stored unrounded as \code{target_per_group}; the
#' smallest whole number that satisfies the inequality is
#' \code{criterion_per_group}. As for every criterion sample size in the
#' package, a target that exceeds a whole number by no more than one part
#' in 1e10 is taken to be that number, so that representation error in a
#' target that is a whole number does not add an observation. The interval
#' is \eqn{\hat{\psi} \pm z \sqrt{\hat{V}}} with the standard normal
#' quantile \eqn{z} for \code{conf_level}.
#'
#' \strong{How close the coverage is to nominal.} The coverage of the
#' interval is asymptotic. In samples of realistic size it falls short of
#' nominal, because the estimated variance runs below the true one, and the
#' shortfall depends on the pilot size. In simulations with a pilot of 10
#' per group (the article's simulation setting) and a budget that bought
#' about 40 to 100 observations per group, the 95 percent interval covered
#' the contrast in 0.942 to 0.949 of studies for normal populations
#' and 0.936 to 0.941 for skewed ones; with 15 to 30 per group the figures
#' were 0.936 and 0.925. A pilot of 5 gave about the same with the larger
#' budgets and somewhat less with the small one (0.931 and 0.907 in
#' another set of runs), and the minimum pilot of 2 gave markedly less;
#' \code{\link{seq_budget_pilot}} has those figures.
#' \code{\link{seq_budget_sensitivity}} measures the coverage for a design
#' of interest.
#'
#' \strong{When sampling ends.} Sampling ends when (a) every group has met
#' its criterion or (b) the budget cannot pay for another needed
#' observation (the article, p. 430). Three attributes state which:
#' \describe{
#'   \item{\code{rule_met}}{\code{TRUE} when every group with a nonzero
#'     weight has met its criterion.}
#'   \item{\code{budget_exhausted}}{\code{TRUE} when at least one group is
#'     still below its criterion and what is left of the budget cannot pay
#'     for one more observation from any such group. Only the groups that
#'     still need observations count: a cheaper group that has already met
#'     its criterion does not keep the procedure going.}
#'   \item{\code{is_satisfied}}{\code{TRUE} when sampling has ended, by (a)
#'     or by (b), and the variance of every group could be used. It is
#'     \code{FALSE} while sampling should continue. It is also \code{FALSE}
#'     whenever the variance of a group is unusable, which the attribute
#'     \code{variance_unusable} reports (see below).}
#' }
#' The function \code{\link{seq_budget_next}} has no group to name, and
#' returns \code{NA}, exactly when one of the first two attributes is
#' \code{TRUE}; a loop that adds observations should end on that
#' \code{NA}. Unless a variance is unusable, that is also exactly when
#' \code{is_satisfied} is \code{TRUE}.
#'
#' The targets of Equation 18 always cost exactly the budget in total
#' (multiply each by its \eqn{a_i} and add). Judged on the current data,
#' therefore, every group can be at or above its criterion without passing
#' the budget only when each sample size equals its target exactly. So
#' with the default \code{closed = NULL} the procedure nearly always ends
#' by (b), with less of the budget left than the next needed observation
#' costs and one or more groups still below their criteria. The printed
#' object then carries a note that gives the amount left and names those
#' groups, and its stopping line reads \code{Stopping rule met: FALSE
#' (sampling has ended: the budget is exhausted)}, because the rule itself
#' (every group at its criterion) was not met. \code{is_satisfied} is
#' \code{TRUE} all the same: it says that sampling has ended.
#' A fourth attribute, \code{budget_exceeded},
#' is \code{TRUE} when the data supplied already cost more than the budget
#' (for example a pilot the budget cannot pay for); the result is then
#' flagged with a note, because the budget constraint does not hold for it.
#'
#' \strong{Two readings of "has met its criterion".} In the article a
#' population stops being sampled the first time it meets its criterion
#' and is never sampled again (Stage 1 and Stage 2, pp. 429 to 430), and
#' the open populations are visited one at a time, in rounds. A check sees
#' only the current data, so by default it judges every group afresh: a
#' group that met its criterion earlier and has since fallen below it,
#' because the other standard deviations moved, is sampled again. This
#' default is the package's choice, not the article's procedure. In
#' simulations of three designs from
#' Table 1 of the article (budget 5,000, pilot 10, normal and skewed
#' populations) it left about 10 to 12 of the budget unspent where
#' the first passage rule left 19 to 42, and the true variance of the
#' estimator at the final sample sizes was 0.2 to 1.0 percent smaller;
#' coverage differed by at most 0.006, in either direction. For the
#' first passage rule, pass the \code{met_per_group} attribute of
#' the previous check as \code{closed}, and ask
#' \code{\link{seq_budget_next}} for \code{order = "rounds"}; the help
#' page of that function lists the three respects in which the result
#' still differs from the article's procedure. Under the first passage
#' rule \code{rule_met} can become \code{TRUE} with part of the budget
#' unspent.
#'
#' \strong{A group whose variance is unusable.} When the sample variance of
#' a group with a nonzero weight is zero, no variance has been estimated
#' for that group. As everywhere in the package when a variance estimate is
#' unusable, the stopping rule is then reported as not met (\code{rule_met}
#' and \code{is_satisfied} are \code{FALSE}), \code{attr(out,
#' "variance_unusable")} is \code{TRUE}, and a note says so; see
#' \code{\link{smsd_check_methods}}. This happens when every observation in
#' the group is the same, which a small pilot and coarsely scored data can
#' produce, or the same to about ten significant digits, so that whatever
#' variance is computed is rounding error. (Scores so large that a sample
#' variance overflows are not this case: the function stops with a message
#' that says to rescale them.) The group's share of the budget cannot be
#' estimated and the variance of the contrast is understated; its target is
#' reported as zero and its criterion as one more than its current size.
#' The limits are still computed, from the understated variance, so they
#' are too close together (when every group is constant both equal the
#' estimate) and are not an interval to report. \code{summary()},
#' \code{glance()}, and \code{tidy()} return them without the note, so
#' test \code{attr(out, "variance_unusable")} before using the limits.
#' This is a deliberate departure from the stopping rule as printed in the
#' article, under which such a group has met its criterion (the right hand
#' side of Equation 18 is zero) and is never sampled again, so that a pilot
#' of equal scores would fix that group's variance at zero for the rest of
#' the study.
#'
#' \code{\link{seq_budget_next}} names such a group whenever no group
#' whose variance can be used has fewer observations than it has, so
#' it is sampled at the pace of the least sampled of the other groups
#' until its observations differ. When the equal scores were an accident
#' of a small pilot, one or two further observations usually settle the
#' matter. A group that really is constant, or nearly so (a ceiling or a
#' floor on the measure), keeps that pace for as long as its observations
#' stay equal, which can be to the end of the study. It then receives
#' about as many observations as the least sampled of the other groups,
#' far more than its share when it is expensive, and the study ends with a
#' variance that is still unusable. In a simulated example in which the
#' most expensive of three groups scored at its ceiling 98 percent of the
#' time, that group took 49 percent of the budget where its optimal share
#' was 14 percent, the variance of the estimator was 1.6 times the
#' optimum, and 61 percent of the studies ended with the group's variance
#' still unusable. If a group shows no variation after more than a
#' few observations beyond the pilot, stop and reconsider the measure
#' before spending the rest of the budget.
#'
#' Sampling still ends when the budget cannot pay for another observation:
#' \code{budget_exhausted} is \code{TRUE} and \code{seq_budget_next()}
#' returns \code{NA}. \code{is_satisfied} stays \code{FALSE}, the printed
#' notes say why, and the decision functions refuse the check, because the
#' interval rests on a variance of zero for at least one group. (For large
#' scores that differ only far down their digits, subtracting the same
#' constant from every score restores a usable variance; it leaves a
#' contrast whose weights sum to zero unchanged.)
#'
#' \strong{Names.} \code{groups}, \code{contrast_weights}, and
#' \code{cost_per_group} are matched by position. If more than one of them
#' is named, the names must be the same and in the same order; otherwise
#' the function stops, so that weights or costs named in another order
#' (for example against a list made by \code{split()}, which is in
#' alphabetical order) are not applied to the wrong groups.
#'
#' @param groups A list of \emph{K} (at least 2) numeric vectors, one per
#'   population, holding the observations collected so far. Each needs at
#'   least 2 observations, so that a sample variance exists, and no missing
#'   or infinite values. The groups may have different sizes.
#' @param contrast_weights Numeric vector of length \emph{K} giving the
#'   contrast coefficients \eqn{c_i}, in the order of \code{groups}, not all
#'   zero. Code the contrast so that larger values are better if a
#'   noninferiority or superiority decision is planned. The article takes
#'   the weights to sum to zero; the sum is not checked.
#' @param cost_per_group Numeric vector of length \emph{K} giving the cost
#'   \eqn{a_i} of one observation from each population, each positive, in
#'   the same units as \code{budget} and in the order of \code{groups}.
#' @param budget The total budget \eqn{A_0} for all observations, the pilot
#'   included: a single positive number.
#' @param conf_level The confidence level of the two-sided interval
#'   (default 0.95). It plays no part in the decision functions, which
#'   take a one-sided level of their own.
#' @param closed Optional; leave \code{NULL} (the default) to judge every
#'   group on the current data. For the first passage rule, a
#'   logical vector of length \emph{K} that is \code{TRUE} for each group
#'   from which sampling has already stopped: the \code{met_per_group}
#'   attribute of the previous check. A closed group stays closed, unless
#'   its variance is unusable.
#'
#' @return An object of class \code{c("smsd_check", "data.frame")} with
#'   columns \code{term} and \code{value} and six rows, in this order:
#'   \code{lower_limit}, \code{contrast} (the estimate \eqn{\hat{\psi}}),
#'   \code{upper_limit}, \code{variance} (the estimated variance
#'   \eqn{\hat{V}} of the contrast estimator), \code{current_total_n}, and
#'   \code{current_total_cost}. A group whose variance is unusable adds
#'   nothing to the variance row, so the two limits are then too close
#'   together and are not an interval to report (see Details). Its
#'   attributes are:
#'   \describe{
#'     \item{\code{method}}{the label of the procedure,
#'       \code{"budget constrained K-population contrast"}. The next
#'       observation and decision functions recognize the object by it.}
#'     \item{\code{is_satisfied}, \code{rule_met},
#'       \code{budget_exhausted}}{the three stopping attributes described
#'       under Details.}
#'     \item{\code{budget_exceeded}}{\code{TRUE} when the data already cost
#'       more than the budget.}
#'     \item{\code{variance_unusable}}{\code{TRUE} when the sample variance
#'       of at least one group with a nonzero weight is unusable (zero, or
#'       computed from observations that agree to about ten significant
#'       digits); \code{FALSE} otherwise.}
#'     \item{\code{current_n}, \code{criterion_n}}{the totals over groups
#'       of the current sample sizes and of the criteria.}
#'     \item{\code{conf_level}}{the confidence level of the interval.}
#'     \item{\code{notes}}{a character vector with the notes that
#'       \code{print()} shows (empty when there are none).}
#'     \item{\code{n_per_group}, \code{mean_per_group},
#'       \code{sd_per_group}}{the sample size, mean, and standard deviation
#'       of each group, in the order of \code{groups}.}
#'     \item{\code{target_per_group}}{the right hand side of Equation 18
#'       for each group, not rounded.}
#'     \item{\code{criterion_per_group}}{the smallest whole number at or
#'       above each target, a target within one part in 1e10 above a whole
#'       number counting as that number (an integer vector); for a group
#'       whose variance is unusable, its current size plus one.}
#'     \item{\code{met_per_group}}{\code{TRUE} for each group that has met
#'       its criterion or, under the first passage rule, was closed
#'       earlier.}
#'     \item{\code{not_estimable_per_group}}{\code{TRUE} for each group
#'       with a nonzero weight whose variance is unusable.}
#'     \item{\code{cost_per_group}, \code{contrast_weights},
#'       \code{budget}}{the inputs.}
#'     \item{\code{remaining_budget}}{the budget minus the current total
#'       cost; negative when the budget is exceeded.}
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
#' Luh, W. M., & Guo, J. H. (2016). Sample size planning for the
#' noninferiority or equivalence of a linear contrast with cost
#' considerations. \emph{Psychological Methods}, \emph{21}(1), 13--34.
#' \doi{10.1037/met0000039}
#'
#' Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). \emph{Designing
#' experiments and analyzing data: A model comparison perspective}
#' (4th ed.). Routledge.
#'
#' @examples
#' # Scores after a new training format (group 1) and two existing formats.
#' # An observation costs 38, 15, and 10, and the budget is 1,500.
#' set.seed(113)
#' weights <- c(1, -0.5, -0.5)
#' costs <- c(38, 15, 10)
#' scores <- list(
#'   rnorm(10, mean = 52, sd = 9),
#'   rnorm(10, mean = 50, sd = 6),
#'   rnorm(10, mean = 50, sd = 6)
#' )
#' out <- seq_budget(scores, weights, cost_per_group = costs, budget = 1500)
#' out
#' attr(out, "criterion_per_group")
#' attr(out, "remaining_budget")
#'
#' # Add one observation at a time, to the group seq_budget_next() names,
#' # until it has none to name. The cost never passes the budget.
#' true_means <- c(52, 50, 50)
#' true_sds <- c(9, 6, 6)
#' repeat {
#'   out <- seq_budget(scores, weights, cost_per_group = costs, budget = 1500)
#'   g <- seq_budget_next(out)
#'   if (is.na(g)) break
#'   scores[[g]] <- c(scores[[g]], rnorm(1, true_means[g], true_sds[g]))
#' }
#' out
#' attr(out, "n_per_group")
#'
#' @seealso \code{\link{seq_budget_next}} for the next observation,
#'   \code{\link{seq_budget_pilot}} and \code{\link{seq_budget_optimal}}
#'   for planning, \code{\link{decide_seq_budget}} for the tests, and
#'   \code{\link{seq_budget_sensitivity}} for a simulation of a planned
#'   design. The vignette is \code{vignette("seq-budget", package = "SMSD")}.
#'
#'   In \pkg{DMAR}, \code{\link[DMAR]{ss_seq_c}} is the sequential
#'   fixed-width interval for a contrast, and
#'   \code{\link[DMAR]{ss_aipe_c}} and \code{\link[DMAR]{equivalence_c}}
#'   are the fixed sample size cousins.
#' @family seq_budget
#' @export
seq_budget <- function(groups, contrast_weights, cost_per_group, budget,
                       conf_level = 0.95, closed = NULL) {
  group_names <- if (is.list(groups)) names(groups) else NULL
  groups <- .budget_check_groups(groups)
  K <- length(groups)
  w <- .budget_check_weights(contrast_weights, K,
                             len_what = "one per group in 'groups'")
  a <- .check_vector(cost_per_group, "cost_per_group",
                     "the cost of one observation in each group", len = K,
                     len_what = "one per group in 'groups'", positive = TRUE)
  budget <- .check_positive(budget, "budget")
  conf_level <- .check_conf_level(conf_level)
  .budget_check_size(budget, a)
  closed <- .budget_check_closed(closed, K)
  .budget_check_names(list(groups = group_names,
                           contrast_weights = names(contrast_weights),
                           cost_per_group = names(cost_per_group)))

  n <- lengths(groups)
  means <- vapply(groups, mean, numeric(1))
  sds <- vapply(groups, stats::sd, numeric(1))
  estimable <- vapply(seq_len(K), function(i) {
    .budget_estimable(.budget_varies(groups[[i]]), sds[i])
  }, logical(1))

  state <- .budget_state(n, sds, estimable, w, a, budget, closed)

  contrast <- sum(w * means)
  pieces <- w^2 * sds^2 / n
  variance <- sum(pieces)
  # A group with a nonzero weight and a usable standard deviation must add
  # something to the variance; a zero there is an underflow.
  if (!is.finite(contrast) || !is.finite(variance) ||
      any(pieces[w != 0 & estimable] <= 0)) {
    .budget_stop_scale("The contrast estimate or its variance",
                       "'groups' and 'contrast_weights'")
  }
  limits <- .budget_limits(contrast, variance, (1 - conf_level) / 2)

  new_smsd_check(
    term = c("lower_limit", "contrast", "upper_limit", "variance",
             "current_total_n", "current_total_cost"),
    value = c(limits$lower, contrast, limits$upper, variance, sum(n),
              state$total_cost),
    method = .budget_method,
    # The package's rule: an unusable variance is never "satisfied", even
    # when the budget has ended the sampling.
    is_satisfied = (state$rule_met || state$budget_exhausted) &&
      !state$variance_unusable,
    current_n = sum(n),
    criterion_n = sum(as.numeric(state$criterion)),
    conf_level = conf_level,
    width = NULL,
    extra = list(
      n_per_group = n,
      mean_per_group = means,
      sd_per_group = sds,
      target_per_group = state$target,
      criterion_per_group = state$criterion,
      met_per_group = state$met,
      budget = budget,
      remaining_budget = state$remaining,
      cost_per_group = a,
      contrast_weights = w,
      rule_met = state$rule_met,
      budget_exhausted = state$budget_exhausted,
      budget_exceeded = state$budget_exceeded,
      variance_unusable = state$variance_unusable,
      not_estimable_per_group = state$not_estimable
    ),
    notes = .budget_notes(state, a, budget, first_passage = !is.null(closed))
  )
}


# ---- The next observation -----------------------------------------------

#' Which Group to Sample Next Under a Budget Constraint
#'
#' Given a check from \code{\link{seq_budget}}, returns the group that
#' should receive the next observation, or \code{NA} when sampling has
#' ended. Only a group whose next observation the remaining budget can pay
#' for is ever returned, so a study that adds one observation at a time to
#' the group named here cannot spend more than its budget.
#'
#' @details
#' A group is a candidate when it has not met its criterion (it still needs
#' observations) and one more observation from it costs no more than what
#' is left of the budget. If there is no candidate the function returns
#' \code{NA}, and the check's \code{rule_met} and \code{budget_exhausted}
#' attributes say why. A group that needs observations but cannot be paid
#' for is passed over in favor of one that can, so the budget is used as
#' fully as the unit costs allow.
#'
#' Among the candidates:
#' \describe{
#'   \item{\code{order = "shortfall"}}{(the default, and the package's
#'     choice) the next observation goes to the group furthest below its
#'     estimated optimal size, that is, the one with the smallest ratio of
#'     its current size to its target in Equation 18 of Chattopadhyay,
#'     Bandyopadhyay, Kelley, and Padalunkal (2025). This keeps the
#'     allocation as close as possible to the estimated optimal one at
#'     every step.}
#'   \item{\code{order = "rounds"}}{the next observation goes to the
#'     candidate with the fewest observations. With equal pilot sizes this
#'     visits the open groups one at a time, in rounds, in decreasing order
#'     of their standard deviations, which is the order the article
#'     describes (Stage 2, pp. 429 to 430; one observation per group per
#'     stage). Combine it with the \code{closed} argument of
#'     \code{\link{seq_budget}} for the first passage rule.}
#' }
#' Ties go to the group with the larger sample standard deviation and then
#' to the one listed first.
#'
#' \strong{A group whose variance is unusable}, because every observation
#' in it is equal, has no target to be measured against (the check marks
#' it in \code{not_estimable_per_group}). Under either order it takes its
#' turn by rounds: it is named when no group whose variance can be used
#' has fewer observations than it has, and otherwise only when it is the
#' one group that can still be sampled. It is never sampled ahead of the
#' others, so a group that stays constant cannot draw the whole budget to
#' itself; it does keep pace with the least sampled of the other groups,
#' as the help for \code{\link{seq_budget}} explains. This rule is the
#' package's choice: the article does not treat the case.
#'
#' \strong{How the first passage option differs from the article.} With
#' \code{closed} and \code{order = "rounds"} the package follows the
#' article's rule that a population is not sampled again once it has met
#' its criterion, but three details differ. (1) The article orders the
#' populations by their first stage standard deviations; a check does not
#' know those, so the current ones are used. (2) The article's flowchart
#' (Figure 1) ends the procedure at the first observation that cannot be
#' paid for, and its text does not say whether a cheaper open group may
#' still be sampled when a more expensive one can no longer be afforded;
#' the package continues with the cheaper group, which uses the budget more
#' fully and still cannot exceed it. (3) The article judges a population
#' against its criterion only after an observation has been drawn from that
#' population (Equations 16 and 17); a check judges every open group each
#' time, so a group can be closed after an observation drawn from another
#' group, sooner than the article would close it. The differences are
#' small: in simulations of the normal row of the article's Table 1, the
#' option ended with every group closed in 31 percent of studies where the
#' article's own sequence, programmed separately, did in 29 percent; the
#' mean final sample sizes differed by 0.3 or less and the mean cost by
#' about 1 in a budget of 5,000. Even so, the option is the package's
#' reading of the article's procedure, not a step for step copy of it.
#'
#' @param check An object returned by \code{\link{seq_budget}}.
#' @param order How the next group is chosen among those that still need an
#'   observation and can be paid for: \code{"shortfall"} (default) or
#'   \code{"rounds"}; an unambiguous abbreviation is accepted. See Details.
#'
#' @return A length one integer: the position in \code{groups} (as passed
#'   to \code{\link{seq_budget}}) of the group that should receive the next
#'   observation. It is \code{NA} when sampling has ended, which is exactly
#'   when the check's \code{rule_met} or \code{budget_exhausted} attribute
#'   is \code{TRUE}: either every group has met its criterion, or the
#'   remaining budget cannot pay for an observation from any group that has
#'   not. (That is also when \code{is_satisfied} is \code{TRUE}, unless the
#'   variance of a group is unusable.)
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
#' # Scores after a new training format (group 1) and two existing formats;
#' # an observation costs 38, 15, and 10.
#' set.seed(113)
#' scores <- list(
#'   rnorm(10, mean = 52, sd = 9),
#'   rnorm(10, mean = 50, sd = 6),
#'   rnorm(10, mean = 50, sd = 6)
#' )
#' out <- seq_budget(scores, contrast_weights = c(1, -0.5, -0.5),
#'                   cost_per_group = c(38, 15, 10), budget = 1500)
#' seq_budget_next(out)
#'
#' # Near the end of a smaller study, 24 is left of a budget of 680. Group 1
#' # is still below its criterion but costs 38 per observation, so the next
#' # observation goes to group 3, which also needs one and costs 10.
#' set.seed(113)
#' scores <- list(
#'   rnorm(12, mean = 52, sd = 9),
#'   rnorm(8, mean = 50, sd = 6),
#'   rnorm(8, mean = 50, sd = 6)
#' )
#' out <- seq_budget(scores, contrast_weights = c(1, -0.5, -0.5),
#'                   cost_per_group = c(38, 15, 10), budget = 680)
#' attr(out, "remaining_budget")
#' attr(out, "met_per_group")
#' seq_budget_next(out)
#'
#' @seealso \code{\link{seq_budget}}
#' @family seq_budget
#' @export
seq_budget_next <- function(check, order = c("shortfall", "rounds")) {
  .budget_require_check(check)
  order <- .check_choice(order, c("shortfall", "rounds"), "order")
  affordable <- .budget_affordable(
    .budget_term(check, "current_total_cost"),
    attr(check, "cost_per_group"), attr(check, "budget")
  )
  .budget_pick(
    n = attr(check, "n_per_group"),
    s = attr(check, "sd_per_group"),
    target = attr(check, "target_per_group"),
    met = attr(check, "met_per_group"),
    not_estimable = attr(check, "not_estimable_per_group"),
    affordable = affordable,
    order = order
  )
}


# ---- Decisions ----------------------------------------------------------

#' Noninferiority and Equivalence Decisions for a Budget Constrained Contrast
#'
#' Four tests of the contrast estimated by \code{\link{seq_budget}}, for use
#' once sampling has ended. Each is a one-sided test at level
#' \code{alpha_level} (equivalence is two one-sided tests, each at that
#' level), computed from the contrast estimate and its standard error. The
#' noninferiority and equivalence functions are the tests of Chattopadhyay,
#' Bandyopadhyay, Kelley, and Padalunkal (2025). The superiority and
#' inferiority functions are the one-sided tests against a single margin:
#' the article's superiority test and its noninferiority test are special
#' cases of the first, and the two together are the two halves of its
#' equivalence test.
#'
#' @details
#' Write \eqn{\psi} for the population contrast, \eqn{\hat{\psi}} for its
#' estimate (the \code{contrast} row of the check), \eqn{\hat{V}} for the
#' estimated variance of the estimate (the \code{variance} row), and
#' \eqn{z_{\alpha}} for the upper \eqn{100\alpha} percent point of the
#' standard normal distribution, with \eqn{\alpha} the one-sided level
#' \code{alpha_level}.
#'
#' \strong{The level.} Every decision is computed from those two rows of
#' the check and \code{alpha_level}. None is read from the confidence limits
#' printed on the check, which are two-sided at the check's
#' \code{conf_level}, so a decision does not depend on that confidence
#' level. The limits a decision compares with its margin, \eqn{\hat{\psi}
#' \pm z_{\alpha} \sqrt{\hat{V}}}, are those of a \eqn{1 - 2\alpha}
#' interval: a 90 percent interval for the default \code{alpha_level =
#' 0.05}, the level of the article's tables. This is how the sibling
#' functions \code{\link[DMAR]{ss_seq_c}} and
#' \code{\link[DMAR]{equivalence_c}} in \pkg{DMAR} state the level.
#'
#' \strong{The direction.} The article assumes that a larger response is
#' better (p. 428), so code the contrast so that larger values of
#' \eqn{\psi} favor the treatment of interest. If a smaller response is
#' better, reverse the sign of every contrast weight first.
#'
#' \strong{How close the level is to nominal.} All four tests are
#' asymptotic (the article relies on the central limit theorem, Appendix A),
#' and their level is not exact in samples of realistic size. In simulations
#' of the noninferiority test at nominal level 0.05 (10,000 studies per
#' condition, standard error about 0.002; pilot of 10 per group; the true
#' contrast exactly at the margin), the rejection rate was 0.050 to 0.057
#' for normal populations when the budget bought about 30 to 120
#' observations per group. In a smaller design with 15 to 28 per group it
#' was about 0.058 (0.052 to 0.061 over seven such runs). A one-sided test
#' is more sensitive to skewness than a two-sided interval: for skewed
#' populations the rate ranged from 0.040 to 0.085, the largest where the
#' groups with negative weights were small (about 29 and 48 observations).
#' The inferiority test rejects in the other tail. In two of these designs
#' (the normal row of the article's Table 1 with a budget of 5,000, and the
#' smaller design) its rate at the margin was 0.054 and 0.059 for normal
#' populations and 0.066 and 0.076 for right skewed ones; the
#' noninferiority test in the same runs had 0.056 and 0.059, and 0.043 and
#' 0.055. The equivalence test, with the true contrast on one of its
#' margins, rejects at the rate of the corresponding one-sided test when
#' the other margin is far away, and less often when the two margins are
#' close. A smaller pilot raises the rates further;
#' \code{\link{seq_budget_pilot}} has the figures for the noninferiority
#' test. For normal populations the simulated power agreed with
#' the article's Equation C9 (0.80 to 0.81 in five designs against 0.804).
#' \code{\link{seq_budget_sensitivity}} reports the rejection rate for a
#' design of interest.
#'
#' \strong{When the decision is made.} The level of each test is stated
#' for one decision, made when sampling has ended. If the check says that
#' sampling has not ended (\code{is_satisfied} is \code{FALSE}), the
#' decision is still returned, with a warning that it is provisional. If
#' the variance of a group is unusable (the check's
#' \code{variance_unusable} attribute is \code{TRUE}: the observations in
#' some group are all equal), the functions stop, whether or not the budget
#' has ended the sampling, because the decision would rest on a variance of
#' zero for that group.
#'
#' @section Noninferiority:
#' \code{decide_seq_budget_noninferiority()} is the article's test for
#' noninferiority (Equations 6 and 7; in Appendix C, Equations C1 and C3).
#' With \eqn{\delta} the margin, it tests \eqn{H_0: \psi \le -\delta}
#' against \eqn{H_A: \psi > -\delta}. In words, the null hypothesis is
#' inferiority: the treatment of interest falls short of the comparison by
#' the margin or more. The alternative is noninferiority: any shortfall is
#' smaller than the margin. \code{margin} is the size of the margin, a
#' positive number (the article's \eqn{|\psi_0|}, "the least relevant value
#' of the contrast"); the hypotheses are about \code{-margin}. The null
#' hypothesis is rejected, and noninferiority is declared, when (Equation
#' C3)
#' \deqn{\hat{\psi} > -\delta + z_{\alpha} \sqrt{\hat{V}}.}
#' Noninferiority is the special case \code{margin = -abs(delta)} of the
#' superiority test below: with \code{delta} the noninferiority margin,
#' these two calls return the same decision at any \code{alpha_level}.
#' \preformatted{decide_seq_budget_noninferiority(check, margin = delta)
#' decide_seq_budget_superiority(check, margin = -abs(delta))}
#'
#' @section Superiority:
#' \code{decide_seq_budget_superiority()} is the one-sided test that the
#' contrast is above a margin \eqn{m}: \eqn{H_0: \psi \le m} against
#' \eqn{H_A: \psi > m}. In words, the null hypothesis is that the contrast
#' does not exceed the margin, and the alternative is that it does.
#' \code{margin} is the value \eqn{m} itself, with its sign: it may be
#' positive, zero, or negative. The null hypothesis is rejected when
#' \deqn{\hat{\psi} > m + z_{\alpha} \sqrt{\hat{V}}.}
#' With a positive margin, \eqn{m = |\psi_0|}, this is the article's test
#' for superiority (Equations 8 and 9; Equations C2 and C4): the treatment
#' of interest is better than the comparison by more than the margin. With
#' \code{margin = 0} it is the ordinary one-sided test that the contrast is
#' positive. With a negative margin it is the noninferiority test above.
#' With \eqn{m} the lower equivalence margin \eqn{\psi_1}, it is the first
#' of the two one-sided tests that make up the article's test for
#' equivalence (Equations C11 and C13).
#'
#' @section Inferiority:
#' \code{decide_seq_budget_inferiority()} is the one-sided test that the
#' contrast is below a margin \eqn{m}: \eqn{H_0: \psi \ge m} against
#' \eqn{H_A: \psi < m}. In words, the null hypothesis is that the contrast
#' is at or above the margin, and the alternative, the conclusion the data
#' must establish, is that it is below the margin. \code{margin} is the
#' value \eqn{m} itself, with its sign. The null hypothesis is rejected
#' when
#' \deqn{\hat{\psi} < m - z_{\alpha} \sqrt{\hat{V}}.}
#' This is the mirror image of the superiority test: reversing the sign of
#' every contrast weight and of the margin turns one into the other. In the
#' article these hypotheses and this rejection rule appear, with \eqn{m}
#' the upper equivalence margin \eqn{\psi_2}, as the second of the two
#' one-sided tests that make up the test for equivalence (Equations C12
#' and C14). As a test on its own the article states the hypotheses once,
#' in its footnote 1, with \eqn{m = -|\psi_0|}, as the form its superiority
#' test takes when a smaller response is better. When a larger response is
#' better, \code{margin = -delta} declares the treatment of interest
#' inferior by more than \code{delta}. This is not the article's
#' noninferiority test, in which inferiority is the null hypothesis.
#'
#' @section Equivalence:
#' \code{decide_seq_budget_equivalence()} is the article's test for
#' equivalence (Equation 11; in Appendix C, Equations C10 to C15). With
#' \eqn{\psi_1} the lower and \eqn{\psi_2} the upper equivalence margin,
#' it tests \eqn{H_0: \psi \le \psi_1} or \eqn{\psi \ge \psi_2} against
#' \eqn{H_A: \psi_1 < \psi < \psi_2}. In words, the null hypothesis is
#' nonequivalence (the contrast is at or beyond one of the margins) and the
#' alternative is equivalence (the contrast lies between them).
#' \code{margin_lower} and \code{margin_upper} are the margins
#' \eqn{\psi_1} and \eqn{\psi_2} themselves, with their signs: for
#' equivalence within 3 points either way they are \code{-3} and \code{3}.
#' The article takes \eqn{\psi_1 < 0 < \psi_2}; only \code{margin_lower <
#' margin_upper} is required here. The null hypothesis is rejected when
#' (Equation C15)
#' \deqn{\psi_1 + z_{\alpha} \sqrt{\hat{V}} < \hat{\psi} <
#'   \psi_2 - z_{\alpha} \sqrt{\hat{V}}.}
#' These are two one-sided tests, each at level \eqn{\alpha} (Equations
#' C11 to C14): equivalence is declared exactly when the superiority test
#' at \code{margin_lower} and the inferiority test at \code{margin_upper},
#' each at the same \code{alpha_level}, both reject.
#'
#' @param check An object returned by \code{\link{seq_budget}}. A check
#'   from any other procedure is refused.
#' @param margin A single finite number on the scale of the contrast. For
#'   noninferiority it is the size of the noninferiority margin and must be
#'   positive; inferiority, the null hypothesis, is that the contrast is at
#'   or below \code{-margin}. For superiority and inferiority it is the
#'   value the contrast is compared with, given with its sign (negative,
#'   zero, or positive).
#' @param margin_lower,margin_upper For equivalence: the lower and upper
#'   equivalence margins, single finite numbers on the scale of the
#'   contrast given with their signs, with \code{margin_lower <
#'   margin_upper} (for example \code{-3} and \code{3}).
#' @param alpha_level The one-sided level of the test: a single number
#'   strictly between 0 and 0.5 (default 0.05). For equivalence it is the
#'   level of each of the two one-sided tests.
#'
#' @return A length one logical: \code{TRUE} if the null hypothesis is
#'   rejected, so that the alternative (noninferiority, superiority,
#'   inferiority, or equivalence) is supported, and \code{FALSE} otherwise.
#'   \code{FALSE} is not evidence for the null hypothesis.
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
#' Luh, W. M., & Guo, J. H. (2016). Sample size planning for the
#' noninferiority or equivalence of a linear contrast with cost
#' considerations. \emph{Psychological Methods}, \emph{21}(1), 13--34.
#' \doi{10.1037/met0000039}
#'
#' @examples
#' # Is a cheaper online training format (group 1) no worse, by more than
#' # 3 points, than the average of two established formats?
#' set.seed(113)
#' weights <- c(1, -0.5, -0.5)
#' costs <- c(10, 38, 15)
#' true_means <- c(50, 50, 50)
#' scores <- lapply(true_means, function(m) rnorm(10, mean = m, sd = 6))
#' repeat {
#'   out <- seq_budget(scores, weights, cost_per_group = costs, budget = 3000)
#'   g <- seq_budget_next(out)
#'   if (is.na(g)) break
#'   scores[[g]] <- c(scores[[g]], rnorm(1, mean = true_means[g], sd = 6))
#' }
#' out
#'
#' # Noninferiority with a margin of 3 points, at one-sided level 0.05.
#' decide_seq_budget_noninferiority(out, margin = 3)
#'
#' # The same test, as the special case margin = -abs(3) of superiority.
#' decide_seq_budget_superiority(out, margin = -abs(3))
#'
#' # Is the online format better than the average of the other two at all,
#' # or worse than it?
#' decide_seq_budget_superiority(out, margin = 0)
#' decide_seq_budget_inferiority(out, margin = 0)
#'
#' # Equivalence within 3 points either way: two one-sided tests at 0.05,
#' # and then at a much stricter level.
#' decide_seq_budget_equivalence(out, margin_lower = -3, margin_upper = 3)
#' decide_seq_budget_equivalence(out, margin_lower = -3, margin_upper = 3,
#'                               alpha_level = 0.005)
#'
#' @seealso \code{\link{seq_budget}}. In \pkg{DMAR},
#'   \code{\link[DMAR]{equivalence_c}} gives the fixed sample size
#'   equivalence and noninferiority tests for a contrast with the same
#'   one-sided \code{alpha_level}, and \code{\link[DMAR]{ss_seq_c}} plans a
#'   sequential fixed-width interval for those decisions.
#' @family seq_budget
#' @name decide_seq_budget
NULL

#' @rdname decide_seq_budget
#' @export
decide_seq_budget_noninferiority <- function(check, margin,
                                             alpha_level = 0.05) {
  .budget_require_check(check)
  margin <- .budget_check_margin(
    margin, "margin",
    "the noninferiority margin, a positive number on the scale of the contrast"
  )
  if (margin <= 0) {
    .smsd_stop("'margin' must be greater than 0; got ", format(margin),
               ". Give the size of the margin as a positive number: the ",
               "hypothesis of inferiority is that the contrast is at or ",
               "below -margin.")
  }
  limits <- .budget_decision_limits(check, alpha_level)
  isTRUE(limits$lower > -margin)
}

#' @rdname decide_seq_budget
#' @export
decide_seq_budget_superiority <- function(check, margin, alpha_level = 0.05) {
  .budget_require_check(check)
  margin <- .budget_check_margin(
    margin, "margin",
    "the value the contrast must exceed, on the scale of the contrast"
  )
  limits <- .budget_decision_limits(check, alpha_level)
  isTRUE(limits$lower > margin)
}

#' @rdname decide_seq_budget
#' @export
decide_seq_budget_inferiority <- function(check, margin, alpha_level = 0.05) {
  .budget_require_check(check)
  margin <- .budget_check_margin(
    margin, "margin",
    "the value the contrast must fall below, on the scale of the contrast"
  )
  limits <- .budget_decision_limits(check, alpha_level)
  isTRUE(limits$upper < margin)
}

#' @rdname decide_seq_budget
#' @export
decide_seq_budget_equivalence <- function(check, margin_lower, margin_upper,
                                          alpha_level = 0.05) {
  .budget_require_check(check)
  margin_lower <- .budget_check_margin(
    margin_lower, "margin_lower",
    "the lower equivalence margin, on the scale of the contrast"
  )
  margin_upper <- .budget_check_margin(
    margin_upper, "margin_upper",
    "the upper equivalence margin, on the scale of the contrast"
  )
  if (margin_lower >= margin_upper) {
    .smsd_stop("'margin_lower' must be less than 'margin_upper'; got ",
               "margin_lower = ", format(margin_lower),
               " and margin_upper = ", format(margin_upper), ".")
  }
  limits <- .budget_decision_limits(check, alpha_level)
  isTRUE(limits$lower > margin_lower && limits$upper < margin_upper)
}

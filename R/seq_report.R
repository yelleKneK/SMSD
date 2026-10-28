#' Report a Sequential Procedure in Words
#'
#' Turns the result of a sequential procedure into a short paragraph that can
#' be pasted into a methods or results section: the procedure and its
#' source, its target (the confidence interval width and confidence level,
#' or the loss constants), the sample size at which the stopping rule was
#' met (or that it has not yet been met), and the estimate with its
#' confidence interval. The full reference for the procedure is returned
#' with the paragraph.
#'
#' Given a check object, the paragraph describes that one check: the rule
#' was met at the current sample size, or it has not yet been met. Whether
#' the current sample size is the \emph{first} at which the rule was met is
#' something only the user's sampling record can say. Given a path from
#' \code{\link{seq_path}}, the paragraph describes the stopping stage (the
#' first sample size at which the rule was met or, when it was never met,
#' the last stage checked) and adds where checking began and how often the
#' rule was checked.
#'
#' For an AIPE result whose confidence level is asymptotic, the paragraph
#' states that satisfying the width criterion does not establish nominal
#' coverage at the observed sample size. The simulation vignette
#' \code{vignette("coverage-and-design", package = "SMSD")} gives
#' operating characteristics for difficult designs.
#'
#' The paragraph is a starting point. It reports what the object holds and
#' cannot know how the data were collected, so read it against the study
#' before using it. Notes carried by the object are appended, each
#' introduced by "Note:". One such note says that the variance estimate was
#' zero, negative, or not finite at that sample size, so that no variance
#' had yet been estimated and the stopping rule was reported as not met (see
#' \code{\link{smsd_check_methods}}). For an object flagged in that way
#' (its attribute \code{variance_unusable} is \code{TRUE}) the paragraph
#' gives the estimate without the criterion sample size and without a
#' confidence interval, because neither rests on an estimated variance; the
#' object still holds what the check function returned. For a path, the
#' paragraph also says at how many of the stages before the one it reports
#' the variance estimate was unusable, because the published stopping rules
#' do not set such stages aside.
#'
#' The procedure and its reference are looked up from the \code{method}
#' label of the object. The sources credited are Kelley, Darku, and
#' Chattopadhyay (2018) for the sequential accuracy in parameter estimation
#' procedures for the coefficient of variation, the standardized mean, the
#' standardized mean difference, the regression slope, and the contrast;
#' Chattopadhyay and De (2016) for the Gini index; Kelley, Darku, and
#' Chattopadhyay (2019) for the correlation coefficients and the squared
#' multiple correlation coefficient; Chattopadhyay and Kelley (2016, 2017)
#' for the minimum risk procedures for the coefficient of variation and the
#' standardized mean difference; Chattopadhyay, Bandyopadhyay, Kelley, and
#' Padalunkal (2025) for the budget constrained procedure; and Mukhopadhyay
#' and de Silva (2009) for the classical procedures. An object whose label
#' is not known to the package is still reported, without a citation and
#' without a reference. The constants of a risk based procedure are
#' reported when the object carries them, under the attributes
#' \code{loss_A} (the loss constant), \code{cost_c} or \code{cost_1} and
#' \code{cost_2} (the sampling cost), \code{risk_bound}, and, where the
#' procedure has them, \code{k} (the power of the estimation error in the
#' loss) and \code{t} (the power of the sample size in the cost). A width
#' is reported as the target full width of the confidence interval, except
#' for a procedure that bounds the length of its interval and does not aim
#' at a width (the bounded length interval for a median), where it is
#' reported as the largest length allowed.
#'
#' Sampling under the budget constrained procedure ends in one of three
#' ways, and the paragraph says which: every group reached its criterion
#' sample size, the budget could not pay for another observation from a
#' group still below its criterion, or the data in hand already cost more
#' than the budget. The second is reported as what it is, with the amount
#' of the budget left and the total of the criterion sample sizes that was
#' not reached (left out when the variance of a group could not be
#' estimated), and never as a stopping rule that was met; the third with
#' what the data cost and the budget. The paragraph reads these endings
#' from the attributes \code{rule_met}, \code{budget_exhausted}, and
#' \code{budget_exceeded} that \code{\link{seq_budget}} sets, and the
#' amounts from \code{budget} and \code{remaining_budget}; for a budget
#' constrained check without them it says only that sampling ended, which
#' is true of every ending. The note of \code{\link{seq_budget}} that says
#' the budget is spent, or already exceeded, is not appended when the
#' paragraph has given the amounts, since it would say the same again.
#'
#' @param x For \code{seq_report()}: an object of class
#'   \code{"smsd_check"}, as returned by a sequential check function (for
#'   example \code{\link{seq_aipe_smd}} or \code{\link{seq_mr_cv}}), or of
#'   class \code{"smsd_path"}, as returned by \code{\link{seq_path}}. For
#'   \code{format()} and \code{print()}: the report that
#'   \code{seq_report()} returned.
#' @param digits The number of decimal places shown for the estimate and its
#'   confidence limits. A whole number from 0 to 15; default 3. This affects
#'   the text only: nothing stored in \code{x} is rounded. Constants the
#'   user supplied (the target width, the confidence level, the loss
#'   constants) are shown as supplied, to at most 7 significant digits.
#'
#' @return \code{seq_report()} returns an object of class
#'   \code{c("smsd_report", "character")}: a single character string
#'   holding the paragraph, with the attributes \code{reference} (a
#'   character vector with the full reference for the procedure; of length
#'   zero when the method label is not known) and \code{method} (the method
#'   label of \code{x}). Use \code{as.character()} for the bare string.
#'
#'   \code{print()} shows the paragraph, wrapped to \code{width}
#'   characters, followed by the reference, and returns its input invisibly.
#'   \code{format()} returns the bare paragraph as a plain character vector
#'   of the same length as its input (one string), so a report can be held
#'   in a column of a data frame.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @references
#' Kelley, K., Darku, F. B., & Chattopadhyay, B. (2018). Accuracy in
#' parameter estimation for a general class of effect sizes: A sequential
#' approach. \emph{Psychological Methods}, \emph{23}(2), 226--243.
#' \doi{10.1037/met0000127}
#'
#' Kelley, K., Darku, F. B., & Chattopadhyay, B. (2019). Sequential accuracy
#' in parameter estimation for population correlation coefficients.
#' \emph{Psychological Methods}, \emph{24}(4), 492--515.
#' \doi{10.1037/met0000203}
#'
#' Chattopadhyay, B., & Kelley, K. (2016). Estimation of the coefficient of
#' variation with minimum risk: A sequential method for minimizing sampling
#' error and study cost. \emph{Multivariate Behavioral Research},
#' \emph{51}(5), 627--648. \doi{10.1080/00273171.2016.1203279}
#'
#' Chattopadhyay, B., & Kelley, K. (2017). Estimating the standardized mean
#' difference with minimum risk: Maximizing accuracy and minimizing cost
#' with sequential estimation. \emph{Psychological Methods}, \emph{22}(1),
#' 94--113. \doi{10.1037/met0000089}
#'
#' Chattopadhyay, B., Bandyopadhyay, T., Kelley, K., & Padalunkal, J. J.
#' (2025). A sequential approach for noninferiority or equivalence of a
#' linear contrast under cost constraints. \emph{Psychological Methods},
#' \emph{30}(2), 425--439. \doi{10.1037/met0000570}
#'
#' Chattopadhyay, B., & De, S. K. (2016). Estimation of Gini index within
#' pre-specified error bound. \emph{Econometrics}, \emph{4}(3), 30.
#' \doi{10.3390/econometrics4030030}
#'
#' Mukhopadhyay, N., & de Silva, B. M. (2009). \emph{Sequential methods and
#' their applications}. Chapman & Hall/CRC. ISBN:9780367386535
#'
#' @examples
#' set.seed(113)
#' # Engagement after a training program against a control group.
#' trained <- rnorm(100, mean = 4.8, sd = 0.9)
#' control <- rnorm(100, mean = 4.5, sd = 0.9)
#'
#' # One check, on the data collected so far.
#' out <- seq_aipe_smd(trained, control, width = 0.60)
#' seq_report(out)
#'
#' # A path: the report describes the first sample size at which the rule
#' # was met, and says how often the rule was checked.
#' path <- seq_path(seq_aipe_smd, list(trained, control), width = 0.60,
#'                  by = 5)
#' report <- seq_report(path, digits = 2)
#' report
#' attr(report, "reference")
#'
#' @seealso \code{\link{seq_path}}; \code{\link{smsd_check_methods}} for the
#'   other methods of a check object.
#' @export
seq_report <- function(x, digits = 3L) {
  if (missing(x)) {
    .smsd_stop("'x' must be supplied: a sequential check object (class ",
               "'smsd_check') or a path from seq_path() (class ",
               "'smsd_path').")
  }
  digits <- .check_count(digits, "digits", min = 0L)
  if (digits > 15L) {
    .smsd_stop("'digits' must be a single whole number from 0 to 15 (a ",
               "double holds no more than 15 significant decimal digits); ",
               "got ", digits, ".")
  }
  path <- NULL
  if (inherits(x, "smsd_path")) {
    .path_validate(x)
    path <- x
    check <- attr(x, "final_check", exact = TRUE)
  } else if (inherits(x, "smsd_check")) {
    check <- x
  } else {
    .smsd_stop("'x' must be a sequential check object (class 'smsd_check') ",
               "or a path from seq_path() (class 'smsd_path'); got ",
               .smsd_describe(x), ".")
  }

  s <- summary(check)
  entry <- .report_methods()[[s$method]]
  known <- !is.null(entry)
  if (!known) {
    entry <- list(
      procedure = sprintf("a sequential procedure (%s)", s$method),
      parameter = s$estimate_term, source = NULL, stage_unit = "",
      n_scale = "", two_stage = FALSE, budget = FALSE,
      width_is_bound = FALSE
    )
  }
  source <- if (known) .report_sources()[[entry$source]] else NULL
  met <- isTRUE(s$is_satisfied)
  # How sampling ended under the budget constrained procedure (NULL for
  # every other procedure and for a budget check still under way). Sampling
  # can end there without the rule being met, so the wording follows
  # `ended`, not `met`.
  ending <- if (isTRUE(entry$budget)) .report_budget_ending(check, met)
  ended <- met || !is.null(ending)
  unusable <- isTRUE(attr(check, "variance_unusable", exact = TRUE))
  # The paragraph is English text for a manuscript: a period for the
  # decimal mark and a comma for thousands, whatever getOption("OutDec") is
  # (with a decimal comma the two marks would be the same character).
  num <- function(v) {
    formatC(round(v, digits) + 0, digits = digits, format = "f",
            big.mark = ",", decimal.mark = ".")
  }

  design <- paste0(
    if (ended) "Sample size was determined with " else
      "Sample size is being determined with ", entry$procedure,
    if (!is.null(source)) paste0(" (", source$cite, ")"),
    .report_target(check, s, entry), "."
  )
  stopping <- if (is.null(ending)) {
    .report_stopping(s, entry, met, path, unusable)
  } else {
    .report_budget_stopping(s, entry, ending, unusable,
                            .report_budget_amounts(check))
  }

  estimate <- if (!is.finite(s$estimate)) {
    paste0("The estimate of ", entry$parameter, " is not finite at this ",
           "sample size.")
  } else {
    # No interval for an object that flags its variance estimate as
    # unusable: whatever limits it holds (the fixed width limits, or limits
    # from the variances of the other groups) rest on no estimated variance.
    interval <- ""
    if (!unusable && is.finite(s$lower_limit) && is.finite(s$upper_limit)) {
      interval <- paste0(
        ", ",
        if (is.na(s$conf_level)) "interval" else
          paste0(.report_const(100 * s$conf_level), "% CI"),
        " [", num(s$lower_limit), ", ", num(s$upper_limit), "]"
      )
    }
    if (ended) {
      paste0("At that sample size the estimate of ", entry$parameter,
             " was ", num(s$estimate), interval, ".")
    } else {
      paste0("The interim estimate of ", entry$parameter, " is ",
             num(s$estimate), interval, ".")
    }
  }

  # The check functions explain an unusable variance estimate in a note of
  # their own, which is carried over as it is. The standard note is added
  # only for an object that is flagged and carries no note at all.
  notes <- attr(check, "notes", exact = TRUE)
  if (!is.character(notes)) notes <- character(0)
  if (!is.null(ending)) {
    notes <- .report_budget_notes(notes, ending, .report_budget_amounts(check))
  }
  if (unusable && !length(notes)) {
    notes <- .report_note_unusable(ended)
  }
  notes <- if (length(notes)) paste("Note:", notes) else character(0)

  coverage <- if (identical(attr(check, "coverage_guarantee", exact = TRUE),
                            "asymptotic")) {
    paste0("The confidence level is asymptotic; meeting the width target ",
           "does not establish nominal coverage at this sample size.")
  }
  text <- paste(c(design, stopping, estimate, coverage, .report_earlier_stages(path),
                  notes), collapse = " ")
  structure(text,
            reference = if (is.null(source)) character(0) else
              source$reference,
            method = s$method,
            class = c("smsd_report", "character"))
}

#' @rdname seq_report
#' @param ... Not used.
#' @export
format.smsd_report <- function(x, ...) {
  # One string for each element of x, as format() methods must return, so
  # that a report can sit in a data frame column. Wrapping and the
  # reference belong to print().
  as.character(x)
}

#' @rdname seq_report
#' @param width The line width, in characters, at which \code{print()} wraps
#'   the text (not a confidence interval width). Default: nine tenths of
#'   \code{getOption("width")}.
#' @export
print.smsd_report <- function(x, width = 0.9 * getOption("width"), ...) {
  cat(.report_lines(x, width), sep = "\n")
  invisible(x)
}

# The lines print.smsd_report() writes: the paragraph wrapped to `width`
# and, when the procedure is known, its reference.
#' @noRd
.report_lines <- function(x, width) {
  width <- .check_number(width, "width", lower = 10)
  lines <- strwrap(as.character(x), width = width)
  reference <- attr(x, "reference", exact = TRUE)
  if (length(reference)) {
    lines <- c(lines, "",
               if (length(reference) == 1L) "Reference" else "References")
    for (ref in reference) {
      lines <- c(lines, strwrap(ref, width = width, exdent = 4))
    }
  }
  lines
}

# A constant the user supplied (a width, a confidence percentage, a loss
# constant) or a sample size, shown as supplied: no scientific notation, no
# added zeros, and at most 7 significant digits, so that a width given as
# 1 / 3 reads 0.3333333 and not as sixteen digits. A whole number is shown
# in full whatever its size. `places`, when given, rounds to that many
# decimal places first (an amount of money that is a difference of sums of
# costs, for example). Display only.
#' @noRd
.report_const <- function(v, places = NULL) {
  if (!is.null(places)) v <- round(v, places) + 0
  format(v, big.mark = ",", decimal.mark = ".", scientific = FALSE,
         trim = TRUE, digits = 7L)
}

# A sample size on the scale the procedure reports it.
#' @noRd
.report_n <- function(n, n_scale) {
  n_txt <- .report_const(n)
  switch(n_scale,
         per_group = paste0("n = ", n_txt, " per group"),
         total = paste0("a total sample size of N = ", n_txt),
         paste0("n = ", n_txt))
}

# The target of the procedure: the interval width and confidence level, or
# the loss constants of a risk based procedure when the object carries them.
#' @noRd
.report_target <- function(check, s, entry) {
  scalar <- function(name) {
    v <- attr(check, name, exact = TRUE)
    if (is.numeric(v) && length(v) == 1L && is.finite(v)) v else NULL
  }
  has_width <- !is.na(s$width)
  has_conf <- !is.na(s$conf_level)
  if (has_width && has_conf && isTRUE(entry$width_is_bound)) {
    # A bounded length procedure: the interval's length is not set in
    # advance, only bounded.
    return(paste0(", with a ", .report_const(100 * s$conf_level),
                  "% confidence interval no longer than ",
                  .report_const(s$width)))
  }
  if (has_width && has_conf) {
    return(paste0(", with a target full width of ", .report_const(s$width),
                  " for the ", .report_const(100 * s$conf_level),
                  "% confidence interval"))
  }
  if (has_width) {
    return(paste0(", with a target full width of ", .report_const(s$width),
                  " for the confidence interval"))
  }
  if (has_conf) {
    budget <- scalar("budget")
    return(paste0(", with a ", .report_const(100 * s$conf_level),
                  "% confidence interval",
                  if (!is.null(budget)) {
                    paste0(" and a total budget of ", .report_const(budget))
                  }))
  }
  loss <- scalar("loss_A")
  cost <- scalar("cost_c")
  cost_1 <- scalar("cost_1")
  cost_2 <- scalar("cost_2")
  bound <- scalar("risk_bound")
  # The power of the sample size in the cost, c n^t, for the one procedure
  # that has it (seq_mr_nexp_mean()). Only with t = 1 is c the cost of one
  # observation.
  power <- scalar("t")
  # The power of the estimation error in the loss, A |e|^k, for the
  # procedures that let the user choose it (seq_mr_nexp_mean(),
  # seq_br_norm_mean(), ts_br_norm_mean()). The criterion sample size
  # depends on it, so two analyses that differ in k are not the same
  # design. The other risk based procedures use squared error and carry no
  # such attribute.
  error_power <- scalar("k")
  parts <- c(
    if (!is.null(loss) && !is.null(error_power)) {
      paste0("a loss of A |e|^k for an estimation error of e (A = ",
             .report_const(loss), ", k = ", .report_const(error_power), ")")
    } else if (!is.null(loss)) {
      paste0("a loss constant of A = ", .report_const(loss))
    },
    if (!is.null(cost) && !is.null(power) && power != 1) {
      paste0("a sampling cost of c n^t for n observations (c = ",
             .report_const(cost), ", t = ", .report_const(power), ")")
    } else if (!is.null(cost)) {
      paste0("a sampling cost of c = ", .report_const(cost),
             " per observation")
    } else if (!is.null(cost_1) && !is.null(cost_2)) {
      if (cost_1 == cost_2) {
        paste0("a sampling cost of c = ", .report_const(cost_1),
               " per observation")
      } else {
        paste0("sampling costs of ", .report_const(cost_1), " and ",
               .report_const(cost_2), " per observation in the first and ",
               "second groups")
      }
    },
    if (!is.null(bound)) paste0("a risk bound of w = ", .report_const(bound))
  )
  if (!length(parts)) {
    return("")
  }
  joined <- if (length(parts) == 1L) {
    parts
  } else if (length(parts) == 2L) {
    paste(parts, collapse = " and ")
  } else {
    paste0(paste(parts[-length(parts)], collapse = ", "), ", and ",
           parts[length(parts)])
  }
  paste0(", with ", joined)
}

# How sampling ended under the budget constrained procedure, from the
# attributes seq_budget() sets; NULL while sampling is still under way.
#
#   "exceeded"   the data in hand already cost more than the budget
#   "criterion"  every group reached its criterion sample size
#   "budget"     the budget could not pay for another observation from a
#                group still below its criterion
#   "ended"      the check reports that sampling has ended and carries
#                neither attribute (an object not made by seq_budget())
#
# Sampling that ends because the budget is spent has not reached the
# criterion sample sizes and is never reported as if it had. It can also end
# that way with a variance that cannot be estimated, in which case the check
# reports the rule as not met (is_satisfied is FALSE) although no further
# observation can be taken; `met` alone therefore does not say whether
# sampling has ended.
#' @noRd
.report_budget_ending <- function(check, met) {
  flag <- function(name) isTRUE(attr(check, name, exact = TRUE))
  rule_met <- flag("rule_met")
  exhausted <- flag("budget_exhausted")
  if (!(met || rule_met || exhausted)) {
    return(NULL)
  }
  if (flag("budget_exceeded")) {
    "exceeded"
  } else if (rule_met) {
    "criterion"
  } else if (exhausted) {
    "budget"
  } else {
    "ended"
  }
}

# The budget and what was left of it, from the attributes seq_budget() sets:
# list(budget, remaining, cost), each as text, or NULL for an object that
# does not carry both amounts as single finite numbers. An amount is shown
# as .report_const() shows every constant, so that the budget reads the same
# way wherever the paragraph gives it. What is left is a difference of sums
# of costs, so it is first rounded to the decimal places such a sum the size
# of the budget leaves meaningful (about thirteen significant digits), as
# the notes of seq_budget() do: 12 is not written 11.99999999999.
#' @noRd
.report_budget_amounts <- function(check) {
  amount <- function(name) {
    v <- attr(check, name, exact = TRUE)
    if (is.numeric(v) && length(v) == 1L && is.finite(v)) v else NULL
  }
  budget <- amount("budget")
  remaining <- amount("remaining_budget")
  if (is.null(budget) || is.null(remaining) || budget <= 0) {
    return(NULL)
  }
  places <- max(0, floor(13 - log10(budget)))
  list(budget = .report_const(budget),
       remaining = .report_const(remaining, places),
       cost = .report_const(budget - remaining, places))
}

# The sentence about how sampling ended under the budget constrained
# procedure. Every form begins "Sampling ended at", which is true of each
# ending; the reason follows when the object gives one, with the amounts
# (`amounts`, from .report_budget_amounts()) when the object carries them.
# The total of the criterion sample sizes is left out when a variance could
# not be estimated, because the criterion of such a group is a placeholder
# (one more than its current size); .report_stopping() does the same while
# sampling continues.
#' @noRd
.report_budget_stopping <- function(s, entry, ending, unusable,
                                    amounts = NULL) {
  criterion <- if (is.finite(s$criterion_n) && !unusable) {
    .report_const(s$criterion_n)
  }
  reason <- switch(
    ending,
    exceeded = if (is.null(amounts)) {
      ", which cost more than the budget"
    } else {
      paste0(", which cost ", amounts$cost, ", more than the budget of ",
             amounts$budget, ", so the budget constraint does not hold for ",
             "this result")
    },
    criterion = ", when every group had reached its criterion sample size",
    budget = paste0(
      " because the budget could not pay for another observation from a ",
      "group still below its criterion sample size",
      if (!is.null(amounts)) {
        paste0(": ", amounts$remaining, " of the ", amounts$budget,
               " was left")
      },
      if (!is.null(criterion)) {
        paste0(" (the criterion sample sizes at that point came to ",
               criterion, " in total)")
      }),
    ""
  )
  paste0("Sampling ended at ", .report_n(s$current_n, entry$n_scale), reason,
         ".")
}

# The notes of a budget constrained check that the paragraph carries over.
# seq_budget() explains in a note of its own that the budget is spent, or
# that the data already cost more than the budget; the sentence written by
# .report_budget_stopping() says the same with the same amounts, so that
# note is left out and nothing is said twice. It is recognized by its
# opening words (a test holds them to what seq_budget() writes), and it is
# kept when the object does not carry the amounts, because the sentence
# then cannot give them. Every other note is carried over as it is.
#' @noRd
.report_budget_notes <- function(notes, ending, amounts) {
  opening <- switch(ending,
                    budget = "Sampling has ended because the budget cannot ",
                    exceeded = "The data already cost ",
                    NULL)
  if (is.null(opening) || is.null(amounts)) {
    return(notes)
  }
  notes[!startsWith(notes, opening)]
}

# The note added for an object that flags its variance estimate as unusable
# and carries no note: the package's standard note or, when sampling has
# ended (so that its advice to collect more data cannot be followed), the
# same statement without the advice.
#' @noRd
.report_note_unusable <- function(ended) {
  if (!ended) {
    return(.note_variance_unusable())
  }
  paste0("The variance estimate is not usable at this sample size (it is ",
         "zero, negative, or not finite), so the stopping rule is reported ",
         "as not met and no confidence interval is given.")
}

# For a path: at how many of the stages before the one reported the variance
# estimate was unusable (zero, negative, or not finite), so that the rule was
# reported as not met there. The published rules do not set such stages
# aside, so a write-up should say when it happened. Only a stage the check
# flagged counts: a check that reports no flag (NA in the path) says nothing
# about its variance estimate. NULL when there is no such stage, and for a
# check object.
#' @noRd
.report_earlier_stages <- function(path) {
  if (is.null(path)) {
    return(NULL)
  }
  earlier <- seq_len(.path_stop_row(path) - 1L)
  count <- sum(unclass(path)$variance_unusable[earlier] %in% TRUE)
  if (count == 0L) {
    return(NULL)
  }
  paste0(
    "At ",
    if (length(earlier) == 1L) {
      "the one earlier stage"
    } else if (count < length(earlier)) {
      paste0(.report_const(count), " of the ", .report_const(length(earlier)),
             " earlier stages")
    } else if (count == 2L) {
      "both earlier stages"
    } else {
      paste0("all ", .report_const(count), " earlier stages")
    },
    " the variance estimate was zero, negative, or not finite, so the ",
    "stopping rule was reported as not met there."
  )
}

# The sentence about the stopping rule, for every procedure but a budget
# constrained one whose sampling has ended (.report_budget_stopping()). The
# criterion sample size is left out for an object that flags its variance
# estimate as unusable (`unusable`): the number such an object holds comes
# from no estimated variance (one more than the current sample size, a
# number worked out from the lower bound of the variance estimate, or a
# total that includes such a number for one group), and "the criterion
# sample size at this stage is 13" would tell the reader that 13
# observations will do.
#' @noRd
.report_stopping <- function(s, entry, met, path, unusable) {
  n_txt <- .report_n(s$current_n, entry$n_scale)
  criterion <- if (is.finite(s$criterion_n) && !unusable) {
    .report_const(s$criterion_n)
  }
  if (isTRUE(entry$two_stage)) {
    required <- paste0(
      "the total sample size the procedure requires",
      if (!is.null(criterion)) paste0(" (", criterion, ")")
    )
    return(if (met) {
      paste0("With ", n_txt, ", ", required, " had been reached.")
    } else {
      paste0("With ", n_txt, ", ", required, " has not yet been reached, ",
             "so sampling continues.")
    })
  }
  pending <- paste0(
    " the stopping rule has not yet been met",
    if (!is.null(criterion)) {
      paste0(" (the criterion sample size at this stage is ", criterion, ")")
    },
    ", so sampling continues."
  )
  if (is.null(path)) {
    return(if (met) {
      paste0("The stopping rule was met at ", n_txt, ".")
    } else {
      paste0("At ", n_txt, pending)
    })
  }

  per_group <- if (identical(entry$stage_unit, "per_group")) " per group" else
    ""
  by <- attr(path, "by", exact = TRUE)
  stages <- nrow(path)
  first <- paste0("n = ", .report_const(attr(path, "start", exact = TRUE)),
                  per_group)
  if (identical(attr(path, "start_rule", exact = TRUE), "pilot")) {
    first <- paste0("the pilot sample size (", first, ")")
  }
  checked <- if (stages == 1L) {
    paste0("The stopping rule was evaluated at ", first)
  } else {
    paste0("The stopping rule was evaluated from ", first, " onward, after ",
           if (by == 1L) "each additional observation" else
             paste0("every ", .report_const(by), " additional observations"),
           per_group)
  }
  if (met) {
    paste0(checked, ", and was first met at ", n_txt, ".")
  } else {
    paste0(checked, ". At ", n_txt, pending)
  }
}

# The sources the package's procedures come from: the citation used in the
# running text and the full reference, written as the reference lists of
# the help pages show them: an en dash in a page range (the escape \u2013,
# so that this file stays ASCII), the DOI as doi:..., and the ISBN of the
# book. Plain text otherwise.
#' @noRd
.report_sources <- function() {
  list(
    kdc2018 = list(
      cite = "Kelley, Darku, & Chattopadhyay, 2018",
      reference = paste0(
        "Kelley, K., Darku, F. B., & Chattopadhyay, B. (2018). Accuracy in ",
        "parameter estimation for a general class of effect sizes: A ",
        "sequential approach. Psychological Methods, 23(2), 226\u2013243. ",
        "doi:10.1037/met0000127")
    ),
    kdc2019 = list(
      cite = "Kelley, Darku, & Chattopadhyay, 2019",
      reference = paste0(
        "Kelley, K., Darku, F. B., & Chattopadhyay, B. (2019). Sequential ",
        "accuracy in parameter estimation for population correlation ",
        "coefficients. Psychological Methods, 24(4), 492\u2013515. ",
        "doi:10.1037/met0000203")
    ),
    ck2016 = list(
      cite = "Chattopadhyay & Kelley, 2016",
      reference = paste0(
        "Chattopadhyay, B., & Kelley, K. (2016). Estimation of the ",
        "coefficient of variation with minimum risk: A sequential method ",
        "for minimizing sampling error and study cost. Multivariate ",
        "Behavioral Research, 51(5), 627\u2013648. ",
        "doi:10.1080/00273171.2016.1203279")
    ),
    ck2017 = list(
      cite = "Chattopadhyay & Kelley, 2017",
      reference = paste0(
        "Chattopadhyay, B., & Kelley, K. (2017). Estimating the ",
        "standardized mean difference with minimum risk: Maximizing ",
        "accuracy and minimizing cost with sequential estimation. ",
        "Psychological Methods, 22(1), 94\u2013113. ",
        "doi:10.1037/met0000089")
    ),
    cbkp2025 = list(
      cite = "Chattopadhyay, Bandyopadhyay, Kelley, & Padalunkal, 2025",
      reference = paste0(
        "Chattopadhyay, B., Bandyopadhyay, T., Kelley, K., & Padalunkal, ",
        "J. J. (2025). A sequential approach for noninferiority or ",
        "equivalence of a linear contrast under cost constraints. ",
        "Psychological Methods, 30(2), 425\u2013439. ",
        "doi:10.1037/met0000570")
    ),
    cd2016 = list(
      cite = "Chattopadhyay & De, 2016",
      reference = paste0(
        "Chattopadhyay, B., & De, S. K. (2016). Estimation of Gini index ",
        "within pre-specified error bound. Econometrics, 4(3), 30. ",
        "doi:10.3390/econometrics4030030")
    ),
    mds2009 = list(
      cite = "Mukhopadhyay & de Silva, 2009",
      reference = paste0(
        "Mukhopadhyay, N., & de Silva, B. M. (2009). Sequential methods ",
        "and their applications. Chapman & Hall/CRC. ",
        "ISBN:9780367386535")
    )
  )
}

# Whether the `width` of the procedure with this method label is the
# largest interval length allowed and not a width the interval is aimed at
# (the bounded length interval for a median). The report, the print method
# of a path, and the line under the title of its plot all word the number
# from this one flag.
#' @noRd
.report_width_is_bound <- function(method) {
  if (!(is.character(method) && length(method) == 1L && !is.na(method))) {
    return(FALSE)
  }
  isTRUE(.report_methods()[[method]]$width_is_bound)
}

# Every `method` label the package's check functions produce, with the
# wording of the procedure, the quantity estimated, and the source.
#
# stage_unit: "per_group" when a stage of seq_path() adds one observation to
#   each of several groups; "" otherwise.
# n_scale: the scale of the object's current_n ("per_group", "total" over
#   the groups, or "" for a single sample or paired data).
# two_stage: the procedure fixes a total sample size from a pilot sample
#   rather than checking a rule after every observation.
# budget: sampling can end because a budget is spent, with the rule not met
#   (the budget constrained procedure).
# width_is_bound: the object's `width` is the largest interval length
#   allowed, not the length the interval will have (the bounded length
#   interval for a median).
#
# A label that is not here is reported without a citation, never with a
# guessed one: when a check function is added or its label changes, add or
# change the entry and its case in tests/testthat/test-seq_report.R (the
# test there fails until both are done).
#' @noRd
.report_methods <- function() {
  entry <- function(procedure, parameter, source, stage_unit = "",
                    n_scale = "", two_stage = FALSE, budget = FALSE,
                    width_is_bound = FALSE) {
    list(procedure = procedure, parameter = parameter, source = source,
         stage_unit = stage_unit, n_scale = n_scale, two_stage = two_stage,
         budget = budget, width_is_bound = width_is_bound)
  }
  aipe <- function(target) {
    paste("the sequential accuracy in parameter estimation procedure for",
          target)
  }
  mrpe <- function(target) {
    paste("the sequential minimum risk point estimation procedure for",
          target)
  }
  list(
    # Sequential accuracy in parameter estimation.
    "coefficient of variation" = entry(
      aipe("the coefficient of variation"), "the coefficient of variation",
      "kdc2018"),
    "standardized mean" = entry(
      aipe("the standardized mean"), "the standardized mean", "kdc2018"),
    "standardized mean difference" = entry(
      aipe("the standardized mean difference"),
      "the standardized mean difference", "kdc2018",
      stage_unit = "per_group", n_scale = "per_group"),
    "simple linear regression slope" = entry(
      aipe("the slope of a simple linear regression"),
      "the regression slope", "kdc2018"),
    "contrast" = entry(
      aipe("a linear contrast of means"), "the contrast", "kdc2018",
      stage_unit = "per_group", n_scale = "per_group"),
    "Gini index" = entry(
      aipe("the Gini index"), "the Gini index", "cd2016"),
    "Pearson correlation" = entry(
      aipe("the Pearson correlation coefficient"),
      "the Pearson correlation coefficient", "kdc2019"),
    "Kendall correlation" = entry(
      aipe("Kendall's tau"), "Kendall's tau", "kdc2019"),
    "Spearman correlation" = entry(
      aipe("Spearman's rho"), "Spearman's rho", "kdc2019"),
    "squared multiple correlation, Wald interval" = entry(
      aipe("the squared multiple correlation coefficient (Wald interval)"),
      "the squared multiple correlation coefficient", "kdc2019"),
    "squared multiple correlation, log transformed interval" = entry(
      aipe(paste("the squared multiple correlation coefficient (log",
                 "transformed interval)")),
      "the squared multiple correlation coefficient", "kdc2019"),

    # Minimum risk point estimation.
    "minimum risk point estimation (coefficient of variation)" = entry(
      mrpe("the coefficient of variation"), "the coefficient of variation",
      "ck2016"),
    "minimum risk point estimation (standardized mean difference)" = entry(
      mrpe("the standardized mean difference"),
      "the standardized mean difference", "ck2017",
      stage_unit = "per_group", n_scale = "total"),
    "minimum risk point estimation (normal mean)" = entry(
      mrpe("the mean of a normal distribution"), "the mean", "mds2009"),
    "minimum risk point estimation (exponential mean)" = entry(
      mrpe("the mean of an exponential distribution"), "the mean",
      "mds2009"),
    "minimum risk point estimation (negative exponential location)" = entry(
      mrpe("the location of a negative exponential distribution"),
      "the location parameter", "mds2009"),
    "minimum risk point estimation (distribution free mean, arbitrary variance)" = entry(
      paste("the distribution free sequential minimum risk point",
            "estimation procedure for a mean"), "the mean", "mds2009"),

    # Bounded risk, fixed width, and bounded length (purely sequential).
    "bounded risk estimation (normal mean)" = entry(
      paste("the sequential bounded risk point estimation procedure for",
            "the mean of a normal distribution"), "the mean", "mds2009"),
    "fixed width interval (normal mean)" = entry(
      paste("the sequential fixed width confidence interval procedure for",
            "the mean of a normal distribution"), "the mean", "mds2009"),
    "fixed width interval (distribution free mean)" = entry(
      paste("the distribution free sequential fixed width confidence",
            "interval procedure for a mean"), "the mean", "mds2009"),
    "fixed width interval (negative exponential location)" = entry(
      paste("the sequential fixed width confidence interval procedure for",
            "the location of a negative exponential distribution"),
      "the location parameter", "mds2009"),
    "bounded length confidence interval (distribution free)" = entry(
      paste("the distribution free sequential bounded length confidence",
            "interval procedure for a median"), "the median", "mds2009",
      width_is_bound = TRUE),

    # Two stage procedures.
    "two stage fixed width interval (normal mean)" = entry(
      paste("the two stage fixed width confidence interval procedure for",
            "the mean of a normal distribution"), "the mean", "mds2009",
      two_stage = TRUE),
    "two stage fixed width interval (distribution free mean)" = entry(
      paste("the distribution free two stage fixed width confidence",
            "interval procedure for a mean"), "the mean", "mds2009",
      two_stage = TRUE),
    "two stage fixed width interval (negative exponential location)" = entry(
      paste("the two stage fixed width confidence interval procedure for",
            "the location of a negative exponential distribution"),
      "the location parameter", "mds2009", two_stage = TRUE),
    "two stage bounded risk estimation (normal mean)" = entry(
      paste("the two stage bounded risk point estimation procedure for",
            "the mean of a normal distribution"), "the mean", "mds2009",
      two_stage = TRUE),
    "two stage bounded risk estimation (exponential mean)" = entry(
      paste("the two stage bounded risk point estimation procedure for",
            "the mean of an exponential distribution"), "the mean",
      "mds2009", two_stage = TRUE),

    # Budget constrained.
    "budget constrained K-population contrast" = entry(
      paste("the sequential budget constrained procedure for a linear",
            "contrast of population means"), "the contrast", "cbkp2025",
      n_scale = "total", budget = TRUE)
  )
}

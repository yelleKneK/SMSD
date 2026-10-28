#' Stopping Path of a Sequential Procedure for Data in Hand
#'
#' Applies a sequential check function to data already collected, at each
#' successive sample size, and returns the whole path in one object: the
#' estimate, its confidence limits, the criterion sample size, and the
#' stopping decision at every stage, together with the first sample size at
#' which the stopping rule is met. One call replaces the collect and check
#' loop that would otherwise be written by hand, and the \code{plot()} method
#' draws the picture a preregistration or a methods section needs.
#'
#' @section What the Path Shows, and What It Does Not:
#' The path shows what the stopping rule \emph{would have said} at each
#' sample size for these data in this order. A sequential design decides at
#' each stage from the observations collected so far, without looking ahead:
#' it stops at the first sample size at which the rule is met, and the
#' estimate and interval it reports are the ones at that sample size
#' (\code{stop_n}), not at a later, more favorable one. The rows after
#' \code{stop_n} are shown for illustration only. Choosing a later stage
#' because its estimate or interval looks better is not the sequential
#' procedure, and the properties the procedure was derived to have (the
#' width and coverage of the interval, or the risk) do not hold for an
#' estimate chosen that way. For the same reason the observations must be in
#' the order in which they were, or will be, collected.
#'
#' @section Supplying the Data:
#' At each stage \code{data} is cut to its first \eqn{n} observations and
#' passed to \code{check}.
#' \itemize{
#'   \item A vector is passed as the first argument of \code{check} (the one
#'     sample procedures).
#'   \item A list holds one element for each data argument of \code{check}:
#'     either named after those arguments, as in
#'     \code{list(group_1 = a, group_2 = b)} for two independent groups or
#'     \code{list(x = x, y = y)} for a predictor and a response, or unnamed
#'     and in the order \code{check} expects them (the help page of each
#'     check function gives the names of its data arguments). Each element
#'     is a vector, a matrix or data frame whose rows are observations, or a
#'     list of these (for example the \code{groups} of
#'     \code{\link{seq_aipe_c}}), every member of which is cut to its first
#'     \eqn{n} observations.
#'   \item A data frame is read as a named list of its columns, so the column
#'     names must be argument names of \code{check}.
#' }
#' Every element must hold the same number of observations, because stage
#' \eqn{n} uses the first \eqn{n} of each: \eqn{n} is the sample size for one
#' sample procedures, the per group sample size for group designs, and the
#' number of pairs for paired data. The column \code{current_n} gives the
#' sample size on the scale \code{check} itself reports (for example the
#' total over both groups for \code{\link{seq_mr_smd}}), which is the scale
#' of \code{criterion_n}.
#'
#' \emph{All} of the data go in \code{data}. Only \code{data} is cut to its
#' first \eqn{n} observations; everything in \code{...} is passed whole at
#' every stage. A second group or a second variable written after
#' \code{data} (for example the second group in
#' \code{seq_path(seq_mr_smd, group_1, group_2, ...)}) would therefore be
#' seen in full from the first stage on, which is the looking ahead a
#' sequential design rules out. \code{seq_path()} stops with a message when
#' an argument in \code{...} looks like data: it holds two or more
#' observations and \code{check} still runs when it is cut along with
#' \code{data}. In the rare case that such an argument really is not data,
#' wrap it in \code{\link[base:AsIs]{I}()} to pass it whole.
#'
#' Missing values are refused. The path is indexed by the number of
#' observations collected, and several check functions drop missing values
#' by default (the minimum risk, bounded risk, fixed width, and bounded
#' length functions, whose \code{na.rm} is \code{TRUE} unless set
#' otherwise; the \code{seq_aipe_*} functions stop on a missing value), so
#' with missing values the sample size a stage reports would be smaller
#' than the stage and the first check would come before the pilot sample
#' size had been reached. Remove the incomplete observations first
#' (both members of a pair, the whole row of a matrix or data frame), and
#' trim the elements to a common length if they then differ.
#'
#' @section Where the Path Starts:
#' A sequential procedure begins checking at its pilot sample size, not
#' before, and the path does the same. By default (\code{start = NULL}) the
#' first stage is the procedure's pilot sample size whenever \code{check} is
#' one of this package's check functions: its \code{_pilot} companion (for
#' example \code{\link{seq_aipe_cv_pilot}} for \code{\link{seq_aipe_cv}}) is
#' called with those arguments in \code{...} that it shares with
#' \code{check}, and with its own defaults for the rest. For any other
#' function (a wrapper, or a check function of the user's own) the default
#' is the smallest sample size the function accepts, found by trying.
#' \code{start} can also be given as a number, as a pilot function, or as
#' \code{"earliest"} to see the rule from the smallest sample size the check
#' accepts. A pilot function that needs an argument \code{check} does not
#' take (the minimum pilot sample size \code{m_0}, for example, or a number
#' of predictors) cannot be given that argument through \code{...}: call the
#' pilot function yourself and pass the number it returns as \code{start}.
#'
#' Where the path starts can change where it stops. The stopping rules of
#' \code{\link{seq_aipe_cv}}, \code{\link{seq_aipe_sm}},
#' \code{\link{seq_aipe_smd}}, \code{\link{seq_aipe_r}},
#' \code{\link{seq_aipe_reg_coef}}, \code{\link{seq_aipe_c}},
#' \code{\link{seq_aipe_gini}}, \code{\link{seq_mr_cv}}, and
#' \code{\link{seq_mr_smd}} (with one sampling cost) contain a correction
#' term that keeps the rule from being met below the sample size the pilot
#' formula gives. The pilot sample size is the larger of that and a minimum
#' (4 by default), so for these procedures an earlier start can change
#' \code{stop_n} only when the target is so loose that the formula gives
#' less than the minimum, and then only if the check function accepts fewer
#' observations than the minimum. Several of the classical procedures (for
#' example the fixed width, minimum risk, and bounded risk procedures for a
#' normal mean) have no correction term at all: started below the pilot,
#' their rule can be met by chance at a very small sample size, which is
#' not a sample size at which the procedure itself would have stopped. If a
#' study used a pilot sample size other than the default (a larger minimum,
#' for example), supply that number.
#'
#' @section Procedures That Do Not Fit:
#' The two stage procedures (\code{ts_*}) take one pilot sample and one
#' further sample whose size the pilot determines, so they have no stage by
#' stage path: call the \code{ts_*} function once on the pilot sample to
#' learn the total sample size required, and once more when the additional
#' observations are in. The budget constrained procedure
#' (\code{\link{seq_budget}}) allocates observations unequally across groups
#' according to their costs and variances, so its path is not indexed by one
#' common sample size: call \code{seq_budget()} after each new observation.
#' The minimum risk procedure for the standardized mean difference with two
#' costs (\code{\link{seq_mr_smd}} with \code{c = c(c1, c2)}) is refused for the
#' same reason: its pilot returns one size per group. \code{seq_path()} stops
#' with a message saying so when given one of these.
#'
#' @param check A sequential check function: any function that returns a
#'   sequential check object (class \code{"smsd_check"}), such as
#'   \code{\link{seq_aipe_cv}}, \code{\link{seq_aipe_smd}},
#'   \code{\link{seq_aipe_r}}, or \code{\link{seq_mr_cv}}.
#' @param data The data collected so far, in the order collected: a vector,
#'   a list with one element for each data argument of \code{check}, or a
#'   data frame whose columns are named after those arguments. See
#'   \emph{Supplying the Data}.
#' @param ... The other arguments of \code{check}, passed unchanged and
#'   whole at every stage (for example \code{width} and \code{conf_level},
#'   or \code{A} and \code{c}). Data do not belong here: an argument that
#'   looks like data is refused unless it is wrapped in
#'   \code{\link[base:AsIs]{I}()} (see \emph{Supplying the Data}).
#'   \code{seq_path()} matches its own arguments \code{check} and
#'   \code{data} by exact name or by position only, so an argument such as
#'   \code{c = 75} or \code{d = 0.5} reaches the check function and is not
#'   taken, by partial matching, for \code{check} or \code{data}.
#' @param start The first sample size checked. \code{NULL} (the default)
#'   is the procedure's pilot sample size when \code{check} is one of this
#'   package's check functions, and the smallest sample size \code{check}
#'   accepts otherwise. Other choices: a single whole number; a pilot
#'   function (called with those arguments in \code{...} that it shares
#'   with \code{check}); \code{"pilot"} (the package's pilot function for
#'   \code{check}; an error if there is none); or \code{"earliest"} (the
#'   smallest sample size \code{check} accepts). See \emph{Where the Path
#'   Starts}.
#' @param by The number of observations (per data element) added between
#'   checks. The default, 1, is the purely sequential procedure of the cited
#'   articles. A larger value shows what a study that checks after every
#'   batch of \code{by} observations would have seen; observations beyond
#'   the last complete batch are not used.
#'
#' @return An object of class \code{c("smsd_path", "data.frame")} with one
#'   row per stage and the columns
#'   \describe{
#'     \item{\code{n}}{the number of observations taken from each data
#'       element at this stage (integer).}
#'     \item{\code{estimate}}{the point estimate.}
#'     \item{\code{lower_limit}, \code{upper_limit}}{the confidence limits
#'       (\code{NA} for the risk based procedures, which have no interval).}
#'     \item{\code{current_n}}{the sample size as \code{check} reports it.}
#'     \item{\code{criterion_n}}{the criterion sample size at this stage
#'       (\code{NA} for a procedure that has none).}
#'     \item{\code{is_satisfied}}{whether the stopping rule is met at this
#'       stage.}
#'     \item{\code{variance_unusable}}{whether the variance estimate was
#'       zero, negative, or not finite at this stage. No variance has then
#'       been estimated, and the check functions report the stopping rule as
#'       not met, whatever the other columns show (see
#'       \code{\link{smsd_check_methods}}): the criterion sample size and
#'       any confidence limits in that row rest on no estimated variance.
#'       \code{NA} when \code{check} does not report it.}
#'   }
#'   The attributes are \code{stop_n} (the first \code{n} at which the rule
#'   is met; \code{NA} if it is never met), \code{method},
#'   \code{estimate_term} (the name of the quantity estimated),
#'   \code{conf_level} and \code{width} (absent for procedures without
#'   them), \code{start} (the first sample size checked), \code{start_rule}
#'   (how it was chosen: \code{"pilot"}, \code{"earliest"}, or
#'   \code{"supplied"}), \code{by}, \code{n_observations} (the number of
#'   observations in hand per data element), \code{arguments} (the list
#'   passed through \code{...}, as \code{check} received it), and
#'   \code{final_check}: the \code{smsd_check} object at \code{stop_n} or,
#'   when the rule is never met, at the last stage. Stored values are never
#'   rounded. A warning that \code{check} gives at several stages is
#'   reported once, with the number of stages.
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
#' @examples
#' set.seed(113)
#' # Response times (seconds) of 120 participants, in the order tested. At
#' # what sample size does a 95% confidence interval of full width 0.10 for
#' # the coefficient of variation first meet the stopping rule?
#' response_time <- rgamma(120, shape = 16, rate = 20)
#' path <- seq_path(seq_aipe_cv, response_time, width = 0.10)
#' path
#' attr(path, "stop_n")
#'
#' # Engagement after a training program against a control group, checked
#' # after every 10 participants per group. This study's pilot was 20 per
#' # group, larger than the default pilot, so it is supplied.
#' trained <- rnorm(120, mean = 4.8, sd = 0.9)
#' control <- rnorm(120, mean = 4.5, sd = 0.9)
#' seq_path(seq_aipe_smd, list(group_1 = trained, group_2 = control),
#'          width = 0.60, start = 20, by = 10)
#'
#' # A minimum risk procedure: the loss constant A and the unit cost c
#' # replace the width and the confidence level.
#' seq_path(seq_mr_cv, response_time, A = 400000, c = 75)
#'
#' @seealso \code{\link{plot.smsd_path}} for the plot,
#'   \code{\link{smsd_path_methods}} for the other methods, and
#'   \code{\link{seq_report}} for a sentence that reports the result.
#' @family seq_path
#' @export
seq_path <- function(check, data, ..., start = NULL, by = 1L) {
  # R matches a named argument to `check` or `data` by partial matching, and
  # the check functions take arguments named `c` (the sampling cost) and `d`
  # (a half width), which R would take for `check` and `data`. The arguments
  # are therefore matched again here, by exact name only. `start` and `by`
  # follow `...` in the formals, so R itself matches them exactly.
  supplied <- .path_supplied(sys.call(), parent.frame())
  check <- supplied$check
  data <- supplied$data
  extra <- supplied$extra

  if (!is.function(check)) {
    .smsd_stop("'check' must be a sequential check function such as ",
               "seq_aipe_cv; got ", .smsd_describe(check), ".")
  }
  by <- .check_count(by, "by", min = 1L)
  .path_refuse_family(check)

  .path_refuse_repeated(extra)
  prepared <- .path_prepare(data, check)
  elements <- prepared$elements
  n_obs <- prepared$n_obs
  .path_check_collision(elements, extra, check, .path_hint(prepared))
  # An argument wrapped in I() is the user's word that it is not data.
  as_is <- vapply(extra, inherits, logical(1), what = "AsIs")
  extra[as_is] <- lapply(extra[as_is], .path_unwrap)

  stage <- .path_stage_call(check, elements, extra)
  .path_refuse_data_in_dots(check, stage, extra, as_is, n_obs)
  .path_refuse_pilot_argument(check, extra)

  # Everything in hand first: if the check cannot run on all the data, the
  # arguments or the data are wrong, and the path cannot be computed.
  cache <- list()
  probe <- .path_try(stage, n_obs)
  if (!is.null(probe$error)) {
    .smsd_stop("The check function stopped with an error when given ",
               if (n_obs == 1L) {
                 "the one observation in hand: "
               } else {
                 paste0("all ", n_obs, " observations: ")
               }, .path_sentence(probe$error), .path_hint(prepared))
  }
  .path_check_result(probe$value)
  .path_refuse_label(probe$value)
  estimate_term <- summary(probe$value)$estimate_term
  cache[[as.character(n_obs)]] <- probe

  first <- .path_start(start, check, stage, n_obs, extra)
  if (!is.null(first$tried)) {
    cache[[as.character(first$n)]] <- first$tried
  }
  stages <- seq.int(first$n, n_obs, by = by)
  k <- length(stages)

  estimate <- lower_limit <- upper_limit <- rep(NA_real_, k)
  current_n <- criterion_n <- rep(NA_real_, k)
  is_satisfied <- logical(k)
  variance_unusable <- rep(NA, k)
  warned <- list()
  stop_check <- NULL
  last_check <- NULL

  for (i in seq_len(k)) {
    n <- stages[i]
    tried <- cache[[as.character(n)]]
    if (is.null(tried)) {
      tried <- .path_try(stage, n)
    }
    if (!is.null(tried$error)) {
      .smsd_stop("The check function stopped with an error at n = ", n,
                 ": ", .path_sentence(tried$error),
                 if (i == 1L && first$rule != "earliest") {
                   paste0(" The path starts at n = ", n,
                          if (first$rule == "pilot") {
                            ", the pilot sample size"
                          } else {
                            " because 'start' says so"
                          }, "; a larger 'start' may be needed.")
                 })
    }
    res <- tried$value
    .path_check_result(res)
    row <- .path_stage_values(res, estimate_term)
    estimate[i] <- row$estimate
    lower_limit[i] <- row$lower_limit
    upper_limit[i] <- row$upper_limit
    current_n[i] <- row$current_n
    criterion_n[i] <- row$criterion_n
    is_satisfied[i] <- row$is_satisfied
    variance_unusable[i] <- row$variance_unusable
    for (msg in unique(tried$warnings)) {
      if (!nzchar(msg)) msg <- "(a warning with no message)"
      warned[[msg]] <- c(warned[[msg]], n)
    }
    if (is_satisfied[i] && is.null(stop_check)) {
      stop_check <- res
    }
    last_check <- res
  }

  for (msg in names(warned)) {
    at <- warned[[msg]]
    warning("The check function gave this warning at ", length(at), " of ",
            k, if (k == 1L) " stage" else " stages", " (first at n = ",
            at[1L], "): ", msg, call. = FALSE)
  }

  met <- which(is_satisfied)
  final_check <- if (is.null(stop_check)) last_check else stop_check
  final_summary <- summary(final_check)

  out <- data.frame(
    n = stages,
    estimate = estimate,
    lower_limit = lower_limit,
    upper_limit = upper_limit,
    current_n = current_n,
    criterion_n = criterion_n,
    is_satisfied = is_satisfied,
    variance_unusable = variance_unusable,
    stringsAsFactors = FALSE
  )
  attr(out, "stop_n") <- if (length(met)) stages[met[1L]] else NA_integer_
  attr(out, "method") <- final_summary$method
  attr(out, "estimate_term") <- estimate_term
  attr(out, "conf_level") <- attr(final_check, "conf_level", exact = TRUE)
  attr(out, "width") <- attr(final_check, "width", exact = TRUE)
  attr(out, "start") <- stages[1L]
  attr(out, "start_rule") <- first$rule
  attr(out, "by") <- by
  attr(out, "n_observations") <- n_obs
  attr(out, "arguments") <- extra
  attr(out, "final_check") <- final_check
  class(out) <- c("smsd_path", "data.frame")
  out
}

# The columns of an smsd_path, in order.
#' @noRd
.path_columns <- c("n", "estimate", "lower_limit", "upper_limit",
                   "current_n", "criterion_n", "is_satisfied",
                   "variance_unusable")

# The attributes an smsd_path carries beyond those of a data frame.
#' @noRd
.path_attrs <- c("stop_n", "method", "estimate_term", "conf_level", "width",
                 "start", "start_rule", "by", "n_observations", "arguments",
                 "final_check")

# Evaluate the arguments of a seq_path() call as written and match `check`
# and `data` by exact name or by position, never by partial matching.
# Everything else, apart from `start` and `by`, belongs to the check
# function. Evaluating list(<the arguments as written>) in the caller's
# frame expands a `...` forwarded by a wrapper and evaluates each argument
# once; the promises R created for `check`, `data`, and `...` are never
# forced.
#' @noRd
.path_supplied <- function(call, envir) {
  written <- as.list(call)[-1L]
  written_names <- names(written)
  if (is.null(written_names)) written_names <- rep("", length(written))
  # `start` and `by` are R's to evaluate; here they are only set aside
  # (before evaluation when written out, after it when they arrive through
  # a forwarded `...`).
  written <- written[!(written_names %in% c("start", "by"))]
  values <- tryCatch(
    eval(as.call(c(list(base::list), written)), envir),
    error = function(e) .smsd_stop(conditionMessage(e))
  )
  nms <- names(values)
  if (is.null(nms)) nms <- rep("", length(values))
  own <- nms %in% c("start", "by")
  values <- values[!own]
  nms <- nms[!own]

  take <- function(arg, nms) {
    i <- match(arg, nms)
    if (is.na(i)) i <- match("", nms)
    i
  }
  i <- take("check", nms)
  if (is.na(i)) {
    .smsd_stop("'check' must be supplied: a sequential check function such ",
               "as seq_aipe_cv.")
  }
  check <- values[[i]]
  values <- values[-i]
  nms <- nms[-i]

  i <- take("data", nms)
  if (is.na(i)) {
    .smsd_stop("'data' must be supplied: the observations collected so far, ",
               "as a vector or as a list with one element for each data ",
               "argument of the check function.")
  }
  data <- values[[i]]
  values <- values[-i]

  list(check = check, data = data, extra = values)
}

# The two families whose procedures have no path indexed by one common
# sample size. The package's own functions are recognized by what they are
# (whatever their arguments are called), a function of the user's by the
# arguments those families use, and anything else, as a backstop for a
# wrapper, by the label of the object returned.
#' @noRd
.path_stop_two_stage <- function() {
  .smsd_stop("seq_path() does not apply to the two stage procedures (ts_*). ",
             "A two stage procedure takes one pilot sample and then one ",
             "further sample whose size the pilot determines, so it has no ",
             "stage by stage path. Call the ts_* function once on the pilot ",
             "sample to learn the total sample size required, and once more ",
             "when the additional observations are in.")
}

#' @noRd
.path_stop_budget <- function() {
  .smsd_stop("seq_path() does not apply to the budget constrained ",
             "procedure. seq_budget() allocates observations unequally ",
             "across the groups according to their costs and variances, so ",
             "its path is not indexed by one common sample size. Call ",
             "seq_budget() after each new observation instead; see its help ",
             "page.")
}

#' @noRd
.path_refuse_family <- function(check) {
  switch(.path_family(check),
         two_stage = .path_stop_two_stage(),
         budget = .path_stop_budget())
  invisible(TRUE)
}

# The family without a path that a function belongs to: "two_stage",
# "budget", or "" for neither. A function of the package (`ns`) is known by
# its name: a two stage check function is called ts_*, the budget check
# function seq_budget (their pilot functions are not check functions and
# are left to the general message). Any other function is known by the
# argument names of those families.
#' @noRd
.path_family <- function(check, ns = topenv()) {
  own <- function(pattern) {
    for (name in grep("_pilot$", ls(ns, pattern = pattern), value = TRUE,
                      invert = TRUE)) {
      if (identical(get0(name, envir = ns, inherits = FALSE), check)) {
        return(TRUE)
      }
    }
    FALSE
  }
  arg_names <- names(formals(check))
  if (own("^ts_") || any(c("pilot_data", "additional_data") %in% arg_names)) {
    return("two_stage")
  }
  if (own("^seq_budget$") ||
      any(c("budget", "cost_per_group") %in% arg_names)) {
    return("budget")
  }
  ""
}

# A stage result must be an smsd_check.
#' @noRd
.path_check_result <- function(res) {
  if (!inherits(res, "smsd_check")) {
    .smsd_stop("'check' must be a sequential check function (for example ",
               "seq_aipe_cv or seq_mr_smd): a function that returns an ",
               "object of class 'smsd_check'. The function supplied ",
               "returned ", .smsd_describe(res), ". The pilot (_pilot), ",
               "sensitivity (_sensitivity), and planning functions are not ",
               "check functions.")
  }
  invisible(TRUE)
}

#' @noRd
.path_refuse_label <- function(res) {
  label <- attr(res, "method", exact = TRUE)
  if (is.character(label) && length(label) == 1L) {
    if (grepl("^two stage", label)) .path_stop_two_stage()
    if (grepl("^budget constrained", label)) .path_stop_budget()
  }
  invisible(TRUE)
}

# The quantities one stage contributes to the path, read from the term rows
# and the attributes of the check object. They are the quantities
# summary.smsd_check() reports (a test holds the two together); building its
# one row data frame at every stage would take several times as long as the
# check itself. A stage whose rows do not hold the estimate term is read
# through summary() instead.
#' @noRd
.path_stage_values <- function(res, estimate_term) {
  terms <- .subset2(res, "term")
  values <- .subset2(res, "value")
  pick <- function(name) {
    j <- match(name, terms)
    if (is.na(j)) NA_real_ else values[j]
  }
  unusable <- attr(res, "variance_unusable", exact = TRUE)
  if (!(is.logical(unusable) && length(unusable) == 1L)) unusable <- NA
  if (is.na(match(estimate_term, terms))) {
    s <- summary(res)
    return(list(estimate = s$estimate, lower_limit = s$lower_limit,
                upper_limit = s$upper_limit, current_n = s$current_n,
                criterion_n = s$criterion_n,
                is_satisfied = isTRUE(s$is_satisfied),
                variance_unusable = unusable))
  }
  list(estimate = pick(estimate_term), lower_limit = pick("lower_limit"),
       upper_limit = pick("upper_limit"),
       current_n = attr(res, "current_n", exact = TRUE),
       criterion_n = attr(res, "criterion_n", exact = TRUE),
       is_satisfied = isTRUE(attr(res, "is_satisfied", exact = TRUE)),
       variance_unusable = unusable)
}

# The number of observations in one leaf of the data: a vector, or a matrix
# or data frame whose rows are observations.
#' @noRd
.path_leaf_n <- function(el, label) {
  if (is.data.frame(el)) {
    return(nrow(el))
  }
  if (!is.null(el) && is.atomic(el)) {
    d <- dim(el)
    if (is.null(d) || length(d) == 1L) {
      return(length(el))
    }
    if (length(d) == 2L) {
      return(d[1L])
    }
  }
  .smsd_stop(label, " must be a vector, or a matrix or data frame whose ",
             "rows are observations; got ", .smsd_describe(el), ".")
}

# The number of observations in one leaf of the data that have a missing
# value (NA or NaN): elements of a vector, rows of a matrix or data frame.
#' @noRd
.path_leaf_incomplete <- function(el) {
  if (is.data.frame(el)) {
    if (!ncol(el)) {
      return(0L)
    }
    return(sum(Reduce(`|`, lapply(el, function(column) {
      d <- dim(column)
      if (length(d) == 2L) rowSums(is.na(column)) > 0 else is.na(column)
    }))))
  }
  if (length(dim(el)) == 2L) {
    return(sum(rowSums(is.na(el)) > 0))
  }
  sum(is.na(el))
}

# The first n observations of a data element.
#' @noRd
.path_slice <- function(el, n) {
  if (is.data.frame(el)) {
    return(el[seq_len(n), , drop = FALSE])
  }
  if (is.list(el)) {
    return(lapply(el, .path_slice, n = n))
  }
  if (length(dim(el)) == 2L) {
    return(el[seq_len(n), , drop = FALSE])
  }
  el[seq_len(n)]
}

# Turn `data` into the list of data arguments passed to the check function
# (named, or unnamed and positional), and find the common number of
# observations.
#' @noRd
.path_prepare <- function(data, check) {
  from_frame <- is.data.frame(data)
  from_list <- is.list(data)
  elements <- if (from_frame) {
    as.list(data)
  } else if (from_list) {
    as.list(unclass(data))
  } else {
    list(data)
  }
  if (!length(elements)) {
    .smsd_stop("'data' is empty: supply the observations collected so far.")
  }

  nms <- names(elements)
  if (is.null(nms)) nms <- rep("", length(elements))
  nms[is.na(nms)] <- ""
  named <- nzchar(nms)
  if (any(named) && !all(named)) {
    .smsd_stop("Name every element of 'data' after the argument of the ",
               "check function it belongs to, or name none and give them in ",
               "the order the check function expects.")
  }
  if (all(named)) {
    if (anyDuplicated(nms)) {
      .smsd_stop("The elements of 'data' must have distinct names; got ",
                 paste(nms, collapse = ", "), ".")
    }
    arg_names <- names(formals(check))
    unknown <- setdiff(nms, arg_names)
    if (!is.null(arg_names) && !("..." %in% arg_names) && length(unknown)) {
      .smsd_stop("'data' has ",
                 if (from_frame) "columns" else "elements", " named ",
                 paste(unknown, collapse = ", "), ", which ",
                 if (length(unknown) == 1L) {
                   "is not an argument"
                 } else {
                   "are not arguments"
                 },
                 " of the check function (its arguments are ",
                 paste(arg_names, collapse = ", "), "). ",
                 if (from_frame) {
                   paste0("A data frame is read as a named list of its ",
                          "columns; pass list(...) with the columns the ",
                          "check function needs, named after its arguments.")
                 } else {
                   paste0("Name the elements after the data arguments of ",
                          "the check function, or leave them unnamed and in ",
                          "order.")
                 })
    }
  } else {
    names(elements) <- NULL
  }

  labels <- if (all(named)) paste0("'", nms, "'") else if (from_list) {
    paste("element", seq_along(elements), "of 'data'")
  } else {
    "'data'"
  }
  counts <- integer(0)
  incomplete <- integer(0)
  for (j in seq_along(elements)) {
    el <- elements[[j]]
    if (is.list(el) && !is.data.frame(el)) {
      if (!length(el)) {
        .smsd_stop(labels[j], " is an empty list.")
      }
      member_n <- vapply(seq_along(el), function(g) {
        as.integer(.path_leaf_n(el[[g]],
                                paste0("Member ", g, " of ", labels[j])))
      }, integer(1))
      names(member_n) <- paste0(labels[j], "[[", seq_along(el), "]]")
      counts <- c(counts, member_n)
      incomplete <- c(incomplete,
                      vapply(el, .path_leaf_incomplete, integer(1)))
    } else {
      leaf_n <- as.integer(.path_leaf_n(el, labels[j]))
      names(leaf_n) <- labels[j]
      counts <- c(counts, leaf_n)
      incomplete <- c(incomplete, .path_leaf_incomplete(el))
    }
  }
  if (any(incomplete > 0L)) {
    bad <- incomplete > 0L
    .smsd_stop(
      "'data' holds missing values (NA or NaN): ",
      paste0(names(counts)[bad], ", ", incomplete[bad], " of ", counts[bad],
             " observations", collapse = "; "),
      ". A path is indexed by the number of observations collected, and ",
      "stage n is the first n of them, so every observation must be ",
      "usable: a check function that drops missing values would report a ",
      "sample size smaller than the stage, and the first check would come ",
      "before the pilot sample size had been reached. Remove the ",
      "incomplete observations first (both members of a pair, the whole ",
      "row of a matrix or data frame), and trim the elements of 'data' to ",
      "a common length if they then differ."
    )
  }
  if (length(unique(counts)) > 1L) {
    .smsd_stop("Every element of 'data' must hold the same number of ",
               "observations, because stage n applies the check to the ",
               "first n observations of each; got ",
               paste0(names(counts), ": ", counts, collapse = ", "),
               ". Trim them to a common length. (For a design whose groups ",
               "are meant to grow at different rates, call the check ",
               "function in a loop of your own.)")
  }
  if (counts[1L] < 1L) {
    .smsd_stop("'data' holds no observations.")
  }
  list(elements = elements, n_obs = counts[[1L]],
       positional_list = from_list && !all(named) && length(elements) > 1L)
}

# An argument given twice in `...`. R would stop with "formal argument
# matched by multiple actual arguments" only at the first stage, and the
# checks made before that would describe the second copy as something it is
# not (an argument of the pilot function, say).
#' @noRd
.path_refuse_repeated <- function(extra) {
  nms <- names(extra)
  nms <- nms[!is.na(nms) & nzchar(nms)]
  twice <- unique(nms[duplicated(nms)])
  if (length(twice)) {
    .smsd_stop(paste0("'", twice, "'", collapse = ", "),
               if (length(twice) == 1L) " is" else " are",
               " supplied more than once in '...'. Give each argument of ",
               "the check function once.")
  }
  invisible(TRUE)
}

# An argument may reach the check function from `data` or from `...`, not
# from both. Elements of `data` that carry no name go to the check
# function by position, that is, to its first arguments, so a name in `...`
# that R would match to one of those (exactly or partially) collides too:
# R would give the named value to that argument and push the data on to the
# next one, and the error would be about an argument the user never wrote.
# A name that is exactly another argument's goes to that argument, as R
# matches exact names first. `hint` is added to the message (the reminder
# of how an unnamed list is read, when `data` was one).
#' @noRd
.path_check_collision <- function(elements, extra, check, hint = "") {
  both <- intersect(names(elements), names(extra))
  if (is.null(names(elements))) {
    arg_names <- names(formals(check))
    first <- arg_names[seq_len(min(length(elements), length(arg_names)))]
    extra_names <- names(extra)
    if (length(first) && !("..." %in% first) && !is.null(extra_names)) {
      extra_names <- extra_names[!is.na(extra_names) & nzchar(extra_names)]
      extra_names <- extra_names[!(extra_names %in%
                                     setdiff(arg_names, first))]
      hit <- pmatch(extra_names, first, duplicates.ok = TRUE)
      both <- unique(first[hit[!is.na(hit)]])
    }
  }
  both <- both[nzchar(both)]
  if (length(both)) {
    .smsd_stop("The argument", if (length(both) > 1L) "s", " ",
               paste(both, collapse = ", "), " of the check function ",
               if (length(both) > 1L) "are" else "is",
               " supplied both in 'data' and in '...'. Data belong in ",
               "'data'; everything else in '...'.", hint)
  }
  invisible(TRUE)
}

# A message from R or from a check function, closed with a period when it
# has none, so that a sentence can follow it.
#' @noRd
.path_sentence <- function(msg) {
  msg <- trimws(msg, which = "right")
  if (grepl("[.!?][\"')]*$", msg)) msg else paste0(msg, ".")
}

# Remove the "AsIs" class that I() adds, leaving the object as the check
# function would otherwise have received it.
#' @noRd
.path_unwrap <- function(a) {
  rest <- setdiff(oldClass(a), "AsIs")
  oldClass(a) <- if (length(rest)) rest
  a
}

# The number of observations an argument given in `...` would hold if it
# were data: the length of a vector, the rows of a matrix or data frame, or
# the common count of a list of these. NA for anything else (a function, an
# empty object, a list whose members differ).
#' @noRd
.path_extra_count <- function(a) {
  if (is.data.frame(a)) {
    return(nrow(a))
  }
  if (is.list(a)) {
    if (!length(a)) {
      return(NA_integer_)
    }
    counts <- vapply(a, function(member) {
      if (is.list(member) && !is.data.frame(member)) {
        NA_integer_
      } else {
        as.integer(.path_extra_count(member))
      }
    }, integer(1))
    if (anyNA(counts) || length(unique(counts)) != 1L) {
      return(NA_integer_)
    }
    return(counts[[1L]])
  }
  if (is.null(a) || !is.atomic(a)) {
    return(NA_integer_)
  }
  d <- dim(a)
  if (is.null(d) || length(d) == 1L) {
    return(length(a))
  }
  if (length(d) == 2L) d[1L] else NA_integer_
}

# Data given in `...` in place of `data`. Everything in `...` is passed
# whole at every stage, so a second group or a second variable written
# there, as in seq_path(seq_mr_smd, g1, g2, ...), would be seen in full from
# the first stage on: every stage would look ahead, and for a check that
# accepts groups of unequal size no error would say so. An argument is taken
# for data when it holds two or more observations and the check function
# also runs when that argument and the data are both cut to their first m
# observations, m being one less than the shorter of the two. A vector that
# is not data (contrast weights, a pair of costs) fails that test, because
# the check function does not accept it one element short or does not accept
# so few observations. I() overrides the test.
#' @noRd
.path_refuse_data_in_dots <- function(check, stage, extra, as_is, n_obs) {
  for (j in seq_along(extra)) {
    if (as_is[j]) next
    held <- .path_extra_count(extra[[j]])
    if (is.na(held) || held < 2L) next
    m <- min(held, n_obs) - 1L
    if (m < 1L) next
    whole <- extra[[j]]
    symbol <- stage$extra_symbols[j]
    assign(symbol, .path_slice(whole, m), envir = stage$env)
    tried <- .path_try(stage, m)
    assign(symbol, whole, envir = stage$env)
    if (!is.null(tried$error) || !inherits(tried$value, "smsd_check")) next

    given <- names(extra)[j]
    matched <- .path_extra_names(check, stage, extra)[j]
    which_one <- if (length(given) && !is.na(given) && nzchar(given)) {
      paste0("The argument '", given, "'")
    } else {
      paste0("Unnamed argument ", j, " after 'data'",
             if (nzchar(matched)) {
               paste0(" (the check function's '", matched, "')")
             })
    }
    # The list to pass, in the check function's own argument names when R
    # can match them.
    members <- c(.path_matched_names(check, stage, stage$data_symbols,
                                     names(stage$elements)), matched)
    if (!all(nzchar(members)) || anyDuplicated(members)) {
      members <- c("group_1", "group_2")
    }
    .smsd_stop(
      which_one, " holds ", held, " observations",
      if (held != n_obs) paste0(" ('data' holds ", n_obs, ")"),
      ", and the check function also runs when it is cut along with the ",
      "data, so it looks like data. Everything in '...' is passed whole at ",
      "every stage: each stage would then use observations not yet ",
      "collected at that stage. Put all of the data in 'data', here ",
      "seq_path(check, list(", paste(members, collapse = ", "), "), ...)",
      if (held != n_obs) ", trimmed to a common number of observations",
      ". If the argument is not data, wrap it in I() to pass it whole."
    )
  }
  invisible(TRUE)
}

# An argument of the pilot function given in `...`. The check function does
# not take it, so R's own message would be "unused argument"; say what to do
# instead.
#' @noRd
.path_refuse_pilot_argument <- function(check, extra) {
  nms <- names(extra)
  check_args <- names(formals(check))
  if (is.null(nms) || is.null(check_args) || "..." %in% check_args) {
    return(invisible(TRUE))
  }
  companion <- .path_pilot_companion(check)
  if (is.null(companion)) {
    return(invisible(TRUE))
  }
  # R matches a partial name to an argument of the check function first.
  unmatched <- nms[nzchar(nms) & is.na(pmatch(nms, check_args))]
  pilot_only <- intersect(unmatched, names(formals(companion$fn)))
  if (length(pilot_only)) {
    check_name <- sub("_pilot$", "", companion$name)
    .smsd_stop(
      paste0("'", pilot_only, "'", collapse = ", "),
      if (length(pilot_only) == 1L) " is an argument" else " are arguments",
      " of ", companion$name, "(), not of ", check_name, "(), so ",
      if (length(pilot_only) == 1L) "it" else "they",
      " cannot be passed through '...'. Work out the pilot sample size ",
      "yourself and supply it: start = ", companion$name, "(..., ",
      paste0(pilot_only, " = ", collapse = ", "), ")."
    )
  }
  invisible(TRUE)
}

# A hint added to the first error when an unnamed list may have been meant
# as one list of groups.
#' @noRd
.path_hint <- function(prepared) {
  if (!isTRUE(prepared$positional_list)) {
    return("")
  }
  paste0(" (An unnamed list is read as one element for each data argument ",
         "of the check function, in order. For a function that takes a ",
         "single list of groups, such as seq_aipe_c, use ",
         "data = list(groups = list(...)).)")
}

# The call made at every stage, built once. Each argument is passed as a
# symbol bound in a private environment rather than as a value inside the
# call, so that a message from R or from the check function names the
# argument (data_1, width, and so on) instead of printing the data, and so
# that the check function never has a large object in its call.
#' @noRd
.path_stage_call <- function(check, elements, extra) {
  data_names <- names(elements)
  if (is.null(data_names)) data_names <- rep("", length(elements))
  extra_names <- names(extra)
  if (is.null(extra_names)) extra_names <- rep("", length(extra))
  symbols <- make.unique(c(
    "check_function",
    ifelse(nzchar(data_names), data_names,
           paste0("data_", seq_along(elements))),
    ifelse(nzchar(extra_names), extra_names,
           paste0("argument_", seq_along(extra)))
  ))
  data_symbols <- symbols[1L + seq_along(elements)]
  extra_symbols <- symbols[1L + length(elements) + seq_along(extra)]

  env <- new.env(parent = baseenv())
  assign(symbols[1L], check, envir = env)
  for (j in seq_along(extra)) {
    assign(extra_symbols[j], extra[[j]], envir = env)
  }
  args <- lapply(symbols[-1L], as.name)
  names(args) <- c(data_names, extra_names)
  list(call = as.call(c(list(as.name(symbols[1L])), args)), env = env,
       elements = elements, data_symbols = data_symbols,
       extra_symbols = extra_symbols)
}

# The names of the check function's arguments that the symbols of the stage
# call (`symbols`, with the names they were given in `given`) are matched
# to, as R matches them: an argument given by position, or by a partial
# name, gets the full name of the argument it becomes. A symbol R cannot
# match keeps the name it was given ("" for none).
#' @noRd
.path_matched_names <- function(check, stage, symbols, given) {
  if (is.null(given)) given <- rep("", length(symbols))
  matched <- tryCatch(as.list(match.call(check, stage$call))[-1L],
                      error = function(e) NULL)
  matched_names <- names(matched)
  if (is.null(matched) || is.null(matched_names)) {
    return(given)
  }
  matched_symbols <- vapply(matched, function(v) {
    if (is.symbol(v)) as.character(v) else ""
  }, character(1))
  for (j in seq_along(symbols)) {
    hit <- which(matched_symbols == symbols[j])
    if (length(hit) == 1L && nzchar(matched_names[hit])) {
      given[j] <- matched_names[hit]
    }
  }
  given
}

# The same for the elements of `...`. Used to hand a pilot function the
# arguments it shares with the check function.
#' @noRd
.path_extra_names <- function(check, stage, extra) {
  .path_matched_names(check, stage, stage$extra_symbols, names(extra))
}

# One call of the check function on the first n observations. Returns the
# value, or the error message, with any warnings (which are muffled here and
# reported once by seq_path()).
#' @noRd
.path_try <- function(stage, n) {
  for (j in seq_along(stage$elements)) {
    assign(stage$data_symbols[j], .path_slice(stage$elements[[j]], n),
           envir = stage$env)
  }
  log <- new.env(parent = emptyenv())
  log$warnings <- character(0)
  value <- withCallingHandlers(
    tryCatch(eval(stage$call, stage$env), error = function(e) e),
    warning = function(w) {
      log$warnings <- c(log$warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  if (inherits(value, "error")) {
    return(list(value = NULL, error = conditionMessage(value),
                warnings = character(0)))
  }
  list(value = value, error = NULL, warnings = log$warnings)
}

# The pilot function that goes with a check function of this package: the
# function named <check>_pilot, when `check` is the package's own function
# of that name. NULL for anything else (a wrapper, a user's function).
#' @noRd
.path_pilot_companion <- function(check) {
  ns <- topenv()
  for (pilot_name in ls(ns, pattern = "_pilot$")) {
    owner <- get0(sub("_pilot$", "", pilot_name), envir = ns,
                  inherits = FALSE)
    if (is.function(owner) && identical(owner, check)) {
      return(list(name = pilot_name, fn = get(pilot_name, envir = ns)))
    }
  }
  NULL
}

# The first sample size of the path, and how it was chosen ("pilot",
# "earliest", or "supplied").
#' @noRd
.path_start <- function(start, check, stage, n_obs, extra) {
  advice <- paste0("Supply 'start': the pilot sample size as a number, a ",
                   "pilot function, or \"earliest\" for the smallest sample ",
                   "size the check function accepts.")
  pilot_name <- NULL
  if (is.null(start)) {
    companion <- .path_pilot_companion(check)
    start <- if (is.null(companion)) "earliest" else "pilot"
  }
  if (is.character(start) && length(start) == 1L) {
    start <- .path_also(
      .check_choice(start, c("pilot", "earliest"), "start"),
      "It may also be NULL, a single whole number, or a pilot function."
    )
    if (start == "earliest") {
      for (n in seq_len(n_obs - 1L)) {
        tried <- .path_try(stage, n)
        if (is.null(tried$error)) {
          return(list(n = n, tried = tried, rule = "earliest"))
        }
      }
      return(list(n = n_obs, tried = NULL, rule = "earliest"))
    }
    companion <- .path_pilot_companion(check)
    if (is.null(companion)) {
      .smsd_stop("No pilot function is known for this check function (it ",
                 "is not one of the package's check functions that has a ",
                 "_pilot companion). ", advice)
    }
    pilot_name <- companion$name
    start <- companion$fn
  }

  rule <- "supplied"
  if (is.function(start)) {
    rule <- "pilot"
    called <- if (is.null(pilot_name)) {
      "The pilot function given as 'start'"
    } else {
      paste0("The pilot function ", pilot_name, "()")
    }
    extra_names <- .path_extra_names(check, stage, extra)
    pilot_formals <- formals(start)
    usable <- nzchar(extra_names) & extra_names %in% names(pilot_formals)
    pilot_args <- extra[usable]
    names(pilot_args) <- extra_names[usable]
    # An argument the pilot function requires and the check function does
    # not take (a number of predictors, say) cannot come from '...'.
    required <- names(pilot_formals)[vapply(pilot_formals, function(a) {
      is.symbol(a) && !nzchar(as.character(a))
    }, logical(1))]
    not_shared <- setdiff(required,
                          c(names(pilot_args), names(formals(check)), "..."))
    if (length(not_shared)) {
      .smsd_stop(
        called, " needs ", paste0("'", not_shared, "'", collapse = ", "),
        ", which the check function does not take, so the pilot sample ",
        "size cannot be worked out from the arguments in '...'. Call the ",
        "pilot function yourself and supply the number it returns as ",
        "'start'", if (!is.null(pilot_name)) {
          paste0(": start = ", pilot_name, "(..., ",
                 paste0(not_shared, " = ", collapse = ", "), ")")
        }, "."
      )
    }
    start <- tryCatch(
      do.call(start, pilot_args),
      error = function(e) {
        .smsd_stop(called, " stopped with an error: ",
                   .path_sentence(conditionMessage(e)),
                   " It is called with those arguments in '...' that it ",
                   "shares with the check function. ", advice)
      }
    )
    if (is.numeric(start) && length(start) > 1L) {
      .smsd_stop(called, " returned more than one sample size (",
                 paste(vapply(start, .path_arg_text, character(1)),
                       collapse = ", "), "), one for each group, as ",
                 "it does for a design whose groups are sampled at ",
                 "different rates. A path indexed by one common sample size ",
                 "cannot follow such a design. Pass a single number as ",
                 "'start' to see the rule at equal group sizes, or call the ",
                 "check function in a loop of your own.")
    }
    start <- .check_count(start,
                          "start (the value the pilot function returned)",
                          min = 1L)
  } else {
    # Anything else must be a number. The shared check says so; what it
    # cannot know is that a number is only one of the forms 'start' takes.
    start <- .path_also(
      .check_count(start, "start", min = 1L),
      "It may also be NULL, a pilot function, \"pilot\", or \"earliest\"."
    )
  }
  if (start > n_obs) {
    .smsd_stop(
      if (rule == "pilot") "The pilot sample size is " else "'start' is ",
      start, ", but only ", n_obs,
      if (n_obs == 1L) " observation" else " observations",
      if (length(stage$elements) > 1L) " per data element",
      if (n_obs == 1L) " is" else " are",
      " in hand, so the first check has not yet been reached. Collect ",
      "at least ", start,
      if (rule == "pilot") {
        paste0(", or pass start = \"earliest\" (or a number) to see the ",
               "rule at the sample sizes in hand")
      },
      "."
    )
  }
  list(n = start, tried = NULL, rule = rule)
}

# The value of a shared argument check (`checked`, evaluated here), with one
# sentence (`also`) added to the message of a check that fails: what else
# the argument accepts, which the shared check cannot know.
#' @noRd
.path_also <- function(checked, also) {
  tryCatch(checked, error = function(e) {
    .smsd_stop(conditionMessage(e), " ", also)
  })
}

# The plain stage table: no class, no path attributes.
#' @noRd
.path_strip <- function(x) {
  as.data.frame(unclass(x)[.path_columns], stringsAsFactors = FALSE)
}

# What is wrong with an object that has class smsd_path but is no longer a
# path as seq_path() returned it; NULL when nothing is. Subsetting through
# `[` returns a plain data frame for anything less than the whole path, but
# other operations (rbind(), a data frame verb of another package, setting
# an attribute by hand) can keep the class on part of a path, on two paths
# bound together, or on a table whose attributes are gone. The stages must
# be exactly start, start + by, ..., and stop_n the first of them at which
# the rule is met.
#' @noRd
.path_problem <- function(x) {
  cols <- unclass(x)
  if (!is.list(cols) || !all(.path_columns %in% names(cols))) {
    return("columns are missing")
  }
  required <- setdiff(.path_attrs, c("conf_level", "width"))
  absent <- setdiff(required, names(attributes(x)))
  if (length(absent)) {
    return(paste0(if (length(absent) == 1L) "the attribute " else
      "the attributes ", paste(absent, collapse = ", "),
      if (length(absent) == 1L) " is missing" else " are missing"))
  }
  one_count <- function(v) {
    is.numeric(v) && length(v) == 1L && !is.na(v) && v >= 1 && v == round(v)
  }
  start <- attr(x, "start", exact = TRUE)
  by <- attr(x, "by", exact = TRUE)
  n_obs <- attr(x, "n_observations", exact = TRUE)
  stop_n <- attr(x, "stop_n", exact = TRUE)
  final_check <- attr(x, "final_check", exact = TRUE)
  ok <- one_count(start) && one_count(by) && one_count(n_obs) &&
    start <= n_obs && length(stop_n) == 1L &&
    (is.na(stop_n) || one_count(stop_n)) &&
    inherits(final_check, "smsd_check") &&
    is.character(attr(x, "method", exact = TRUE)) &&
    is.character(attr(x, "estimate_term", exact = TRUE)) &&
    is.character(attr(x, "start_rule", exact = TRUE))
  if (!ok) {
    return("its attributes are not those of a path")
  }
  n <- cols$n
  stages <- seq.int(start, n_obs, by = by)
  if (!is.numeric(n) || length(n) != length(stages) || anyNA(n) ||
      any(n != stages)) {
    return(paste0("its rows are not the stages seq_path() checked, n = ",
                  start, " to ", stages[length(stages)],
                  if (by > 1L) paste0(" by ", by)))
  }
  met <- which(cols$is_satisfied %in% TRUE)
  first_met <- if (length(met)) n[met[1L]] else NA_integer_
  if (!identical(as.integer(stop_n), as.integer(first_met))) {
    return("'stop_n' is not the first stage at which the rule is met")
  }
  NULL
}

# Stop unless `x` is a path as seq_path() returned it.
#' @noRd
.path_validate <- function(x) {
  problem <- .path_problem(x)
  if (!is.null(problem)) {
    .smsd_stop("This object has class 'smsd_path' but is no longer a path ",
               "as seq_path() returned it (", problem, "). It was changed ",
               "afterward, for example by binding paths together, by a ",
               "data frame verb that kept the class on part of the path, ",
               "or by changing an attribute. Use as.data.frame() for the ",
               "plain table of stages, or call seq_path() again.")
  }
  invisible(x)
}

# The row of the stopping stage: the first stage at which the rule is met,
# or the last stage when it is never met.
#' @noRd
.path_stop_row <- function(x) {
  stop_n <- attr(x, "stop_n", exact = TRUE)
  n <- unclass(x)$n
  if (is.na(stop_n)) length(n) else match(stop_n, n)
}

# One of the check's other arguments as short text. A number is shown as it
# would be typed, without scientific notation, while that is short (from
# 1e-5 up to 1e15 in absolute value, and zero); outside that range fixed
# notation would run to hundreds of digits that a double does not hold
# (1e150 as a 151 digit integer), so scientific notation is used, whatever
# the session's scipen option. Anything long is described, not written
# out. The target width and the confidence level in the header and the
# plot are written the same way, so that 95% never reads 9.5e+01%.
#' @noRd
.path_arg_text <- function(a) {
  if (is.numeric(a) && length(a) == 1L && is.null(attributes(a)) &&
      is.finite(a)) {
    fixed <- a == 0 || (abs(a) >= 1e-5 && abs(a) < 1e15)
    return(format(a, scientific = !fixed, digits = 7L, trim = TRUE))
  }
  if (is.function(a)) {
    return("a function")
  }
  if (!is.null(dim(a)) || length(a) > 8L) {
    return(.smsd_describe(a))
  }
  one <- if (is.numeric(a) && length(a) > 1L && is.null(attributes(a)) &&
             all(is.finite(a))) {
    # A short vector of numbers, each written as a single number is.
    paste0("c(", paste(vapply(a, .path_arg_text, character(1)),
                       collapse = ", "), ")")
  } else {
    paste(deparse(a, width.cutoff = 60L, nlines = 1L), collapse = " ")
  }
  if (nchar(one) > 40L) one <- paste0(substr(one, 1L, 37L), "...")
  one
}

# The check's other arguments as short text, one piece ("name = value") for
# each, for the header and the plot.
#' @noRd
.path_args_pieces <- function(args, omit = character(0)) {
  nms <- names(args)
  if (is.null(nms)) nms <- rep("", length(args))
  keep <- !(nms %in% omit)
  args <- args[keep]
  nms <- nms[keep]
  if (!length(args)) {
    return(character(0))
  }
  txt <- vapply(args, .path_arg_text, character(1))
  unname(ifelse(nzchar(nms), paste0(nms, " = ", txt), txt))
}

# The same as one string ("" when there are none).
#' @noRd
.path_args_text <- function(args, omit = character(0)) {
  paste(.path_args_pieces(args, omit), collapse = ", ")
}

# The "Arguments" lines of the printed header: the pieces separated by
# commas and broken only between arguments, never inside one ("f = a
# function" stays on one line). `width` is the longest line allowed; a
# single piece longer than that has a line to itself.
#' @noRd
.path_args_lines <- function(pieces, width = 78L) {
  initial <- "  Arguments       : "
  indent <- strrep(" ", nchar(initial))
  k <- length(pieces)
  pieces <- paste0(pieces, c(rep(",", k - 1L), ""))
  lines <- character(0)
  current <- ""
  for (piece in pieces) {
    if (!nzchar(current)) {
      current <- piece
    } else if (nchar(indent) + nchar(current) + 1L + nchar(piece) <= width) {
      current <- paste(current, piece)
    } else {
      lines <- c(lines, current)
      current <- piece
    }
  }
  lines <- c(lines, current)
  paste0(c(initial, rep(indent, length(lines) - 1L)), lines)
}

# The number of significant digits for display: a whole number from 1 to 22
# (the range base::format() accepts).
#' @noRd
.path_digits <- function(digits) {
  digits <- .check_count(digits, "digits", min = 1L)
  if (digits > 22L) {
    .smsd_stop("'digits' must be a single whole number from 1 to 22; got ",
               digits, ".")
  }
  digits
}

#' Methods for Sequential Path Objects
#'
#' Methods for the \code{smsd_path} objects returned by
#' \code{\link{seq_path}}: the stage by stage record of a sequential
#' procedure applied to data in hand.
#'
#' \code{print()} shows a compact header (the method, the stages checked, the
#' target, the first sample size at which the stopping rule was met, the
#' estimate with its interval there, and any note the check attached at that
#' stage), then the first and last few stages, and closes with a reminder of
#' which estimate is the one to report. With \code{options(OutDec = ",")}
#' the two limits are separated by a semicolon. To keep
#' the table narrow it leaves out columns that carry no information for the
#' path at hand: the confidence limits of a procedure that has none,
#' \code{current_n} when it equals \code{n} at every stage, and
#' \code{variance_unusable} when it is never \code{TRUE}. When it is, the
#' header counts the stages and a line under the table names their sample
#' sizes, because they need not be among the rows shown (stages too
#' scattered to list are counted, with the range they fall in), and says
#' what the column means: at those stages the variance estimate was zero,
#' negative, or not finite, so no variance had yet been estimated and the
#' stopping rule was reported as not met (see
#' \code{\link{smsd_check_methods}}). The header
#' gives the estimate at such a stage without its interval. For a procedure
#' whose interval is bounded in length and not aimed at a width (the bounded
#' length interval for a median) the header gives the bound as the largest
#' length allowed.
#' \code{format()} returns the same lines as a character vector.
#' \code{summary()} returns the stopping stage as a one row data frame: the
#' first stage at which the rule was met or, when it was never met, the last
#' stage checked. \code{as.data.frame()} returns the plain table of stages
#' without the class or the path attributes.
#' Subsetting with \code{[} returns a plain data frame (or vector) whenever
#' the result is no longer the complete path, so part of a path is never
#' printed as if it were the whole of one; a subset that keeps every stage
#' and every column is still the path.
#'
#' Other operations can leave the class on something that is no longer a
#' path: \code{rbind()} of two paths, a data frame verb of another package
#' that keeps the class on some of the rows, or an attribute changed by
#' hand. \code{print()}, \code{format()}, \code{summary()}, \code{plot()},
#' and \code{\link{seq_report}} stop with a message that says so, because
#' anything they reported would describe a path that was never computed.
#' \code{as.data.frame()} always works and gives the plain table.
#'
#' @param x,object An object of class \code{"smsd_path"}.
#' @param digits Number of significant digits to show: a whole number from
#'   1 to 22. Defaults to \code{getOption("digits")}. Stored values are
#'   never rounded; this affects display only.
#' @param rows The number of stages shown at the start and at the end of the
#'   printed table. A whole number of at least 1; default 3. A path with no
#'   more than \code{2 * rows + 1} stages is shown in full.
#' @param row.names,optional Ignored; present for compatibility with the
#'   generic.
#' @param ... Not used.
#'
#' @return \code{print()} returns its input invisibly. \code{format()}
#'   returns a character vector, one element per printed line.
#'   \code{summary()} returns a one row \code{data.frame} with columns
#'   \code{method}, \code{estimate_term}, \code{stop_n} (\code{NA} when the
#'   rule was never met), \code{n} (the stage described), \code{estimate},
#'   \code{lower_limit}, \code{upper_limit}, \code{current_n},
#'   \code{criterion_n}, \code{is_satisfied}, \code{variance_unusable},
#'   \code{conf_level}, \code{width} (\code{NA} where the procedure has no
#'   such quantity), and \code{stages} (the number of stages checked).
#'   \code{as.data.frame()} returns a data frame (class
#'   \code{"data.frame"} only) with one row per stage and the columns
#'   \code{n}, \code{estimate}, \code{lower_limit}, \code{upper_limit},
#'   \code{current_n}, \code{criterion_n}, \code{is_satisfied}, and
#'   \code{variance_unusable}.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @references
#' Kelley, K., Darku, F. B., & Chattopadhyay, B. (2018). Accuracy in
#' parameter estimation for a general class of effect sizes: A sequential
#' approach. \emph{Psychological Methods}, \emph{23}(2), 226--243.
#' \doi{10.1037/met0000127}
#'
#' @examples
#' set.seed(113)
#' # Weekly hours of remote work reported by 90 employees, in survey order.
#' hours <- rgamma(90, shape = 9, rate = 0.6)
#' path <- seq_path(seq_aipe_cv, hours, width = 0.15)
#' print(path, digits = 4, rows = 2)
#' summary(path)
#' head(as.data.frame(path))
#' path[path$is_satisfied, c("n", "estimate")][1, ]
#'
#' @seealso \code{\link{seq_path}}, \code{\link{plot.smsd_path}},
#'   \code{\link{seq_report}}
#' @family seq_path
#' @name smsd_path_methods
#' @aliases smsd_path
NULL

#' @rdname smsd_path_methods
#' @export
format.smsd_path <- function(x, digits = getOption("digits"), rows = 3L,
                             ...) {
  .path_validate(x)
  digits <- .path_digits(digits)
  rows <- .check_count(rows, "rows", min = 1L)
  tab <- .path_strip(x)
  k <- nrow(tab)
  i_stop <- .path_stop_row(x)
  stop_n <- attr(x, "stop_n", exact = TRUE)
  by <- attr(x, "by", exact = TRUE)
  width <- attr(x, "width", exact = TRUE)
  conf_level <- attr(x, "conf_level", exact = TRUE)
  term <- attr(x, "estimate_term", exact = TRUE)
  args_pieces <- .path_args_pieces(attr(x, "arguments", exact = TRUE),
                                   omit = c("width", "conf_level"))
  num <- function(v) format(v, digits = digits)

  lines <- c(
    paste("Sequential path:", attr(x, "method", exact = TRUE)),
    "",
    sprintf("  Stages checked  : %d (n = %s)", k,
            if (k == 1L) {
              format(tab$n[1L])
            } else {
              paste0(tab$n[1L], " to ", tab$n[k],
                     if (by > 1L) paste0(" by ", by))
            }),
    sprintf("  First check at  : %s",
            switch(attr(x, "start_rule", exact = TRUE),
                   pilot = "the pilot sample size",
                   earliest = "the smallest sample size the check accepts",
                   "the sample size supplied as 'start'"))
  )
  if (!is.null(width)) {
    lines <- c(lines, sprintf(
      "  %-16s: %s",
      if (.report_width_is_bound(attr(x, "method", exact = TRUE))) {
        "Length at most"
      } else {
        "Target width"
      }, .path_arg_text(width)
    ))
  }
  if (!is.null(conf_level)) {
    lines <- c(lines, sprintf("  Confidence level: %s",
                              paste0(.path_arg_text(100 * conf_level), "%")))
  }
  if (length(args_pieces)) {
    lines <- c(lines, .path_args_lines(args_pieces))
  }
  lines <- c(lines, sprintf("  Rule first met  : %s",
                            if (is.na(stop_n)) "at no stage checked"
                            else paste("n =", stop_n)))
  est_txt <- sprintf("%s = %s", term, num(tab$estimate[i_stop]))
  # No interval beside an estimate from a stage whose variance estimate was
  # unusable (it can only be the last stage of a rule never met): the limits
  # in that row rest on no estimated variance.
  if (!(tab$variance_unusable[i_stop] %in% TRUE) &&
      is.finite(tab$lower_limit[i_stop]) &&
      is.finite(tab$upper_limit[i_stop])) {
    # The numbers follow the session's decimal mark, as everything printed
    # does. With a decimal comma the limits are separated by a semicolon,
    # so that the comma is not both decimal mark and separator.
    between <- if (identical(getOption("OutDec"), ",")) "; " else ", "
    est_txt <- paste0(est_txt, " [", num(tab$lower_limit[i_stop]), between,
                      num(tab$upper_limit[i_stop]), "]")
  }
  lines <- c(lines, sprintf("  %-16s: %s",
                            if (is.na(stop_n)) "Latest estimate"
                            else "Estimate there", est_txt))
  n_unusable <- sum(tab$variance_unusable %in% TRUE)
  if (n_unusable > 0L) {
    lines <- c(lines, sprintf("  Variance        : not usable at %d of %d %s",
                              n_unusable, k,
                              if (k == 1L) "stage" else "stages"))
  }
  # The notes of the check at the stage reported, as print() shows them for
  # the check itself: a path that stops where the check adds a caution (the
  # lower bound of the variance was used, say) must show the caution. The
  # standard note about an unusable variance estimate is left to the line
  # under the table, which says the same for every such stage by name.
  notes <- attr(attr(x, "final_check", exact = TRUE), "notes", exact = TRUE)
  if (is.character(notes) && (tab$variance_unusable[i_stop] %in% TRUE)) {
    notes <- setdiff(notes, .note_variance_unusable())
  }
  if (is.character(notes)) {
    for (note in notes) {
      lines <- c(lines, strwrap(note, width = 78,
                                initial = "  Note            : ",
                                prefix = strrep(" ", 20L)))
    }
  }

  show <- c("n", "estimate",
            if (any(is.finite(tab$lower_limit) | is.finite(tab$upper_limit)))
              c("lower_limit", "upper_limit"),
            if (!isTRUE(all(tab$current_n == tab$n))) "current_n",
            "criterion_n", "is_satisfied",
            if (n_unusable > 0L) "variance_unusable")
  lines <- c(lines, "", .path_table_lines(tab[show], digits, rows), "")
  if (n_unusable > 0L) {
    # The stages are named by their sample sizes, because the table shows
    # only the first and last few and the marked rows may be among those
    # left out.
    unusable_note <- paste0(
      "At ", .path_stages_text(tab$n, which(tab$variance_unusable %in% TRUE),
                               by),
      " the variance estimate was zero, negative, or not finite (column ",
      "variance_unusable), so the stopping rule was reported as not met ",
      "there.")
    lines <- c(lines, strwrap(unusable_note, width = 72, indent = 2,
                              exdent = 2), "")
  }

  note <- if (is.na(stop_n)) {
    paste0("The stopping rule was not met at any stage checked, so ",
           "sampling would continue.")
  } else if (i_stop < k) {
    paste0("The estimate to report is the one at n = ", stop_n, ". Later ",
           "rows show what the rule would have said had sampling ",
           "continued; they are not grounds for choosing a different ",
           "stopping point.")
  } else {
    paste0("The estimate to report is the one at n = ", stop_n, ".")
  }
  c(lines, strwrap(note, width = 72, indent = 2, exdent = 2))
}

# The stages in rows `at` of a path whose stages are at the sample sizes
# `n`, as words for a sentence: "n = 10", "n = 10 and 11", "n = 3 to 6, 10,
# and 14 to 16". A run of three or more stages in a row is given by its
# ends, followed by "by <by>" when the stages are more than one observation
# apart (as the header gives the stages checked). More than `most` pieces
# would no longer read as a sentence; the stages are then counted and the
# reader pointed to the column.
#' @noRd
.path_stages_text <- function(n, at, by, most = 8L) {
  run <- cumsum(c(1L, diff(at) != 1L))
  pieces <- unlist(lapply(split(at, run), function(rows) {
    if (length(rows) < 3L) {
      as.character(n[rows])
    } else {
      paste0(n[rows[1L]], " to ", n[rows[length(rows)]],
             if (by > 1L) paste0(" by ", by))
    }
  }), use.names = FALSE)
  k <- length(pieces)
  if (k > most) {
    return(paste0(length(at), " stages between n = ", n[at[1L]], " and n = ",
                  n[at[length(at)]], " (as.data.frame() shows which)"))
  }
  paste0("n = ", if (k == 1L) {
    pieces
  } else if (k == 2L) {
    paste(pieces, collapse = " and ")
  } else {
    paste0(paste(pieces[-k], collapse = ", "), ", and ", pieces[k])
  })
}

# The first and last `rows` stages as aligned text lines.
#' @noRd
.path_table_lines <- function(tab, digits, rows) {
  k <- nrow(tab)
  cut <- k > 2L * rows + 1L
  shown <- if (cut) c(seq_len(rows), seq.int(k - rows + 1L, k)) else
    seq_len(k)
  cols <- lapply(names(tab), function(nm) {
    v <- tab[[nm]][shown]
    txt <- if (is.logical(v)) as.character(v) else format(v, digits = digits)
    txt[is.na(v)] <- "NA"
    formatC(c(nm, txt), width = max(nchar(c(nm, txt))))
  })
  body <- paste0("  ", do.call(paste, c(cols, sep = "  ")))
  if (cut) {
    body <- c(body[seq_len(rows + 1L)], "  ...",
              body[seq.int(rows + 2L, length(body))])
  }
  body
}

#' @rdname smsd_path_methods
#' @export
print.smsd_path <- function(x, digits = getOption("digits"), rows = 3L,
                            ...) {
  cat(format(x, digits = digits, rows = rows), sep = "\n")
  invisible(x)
}

#' @rdname smsd_path_methods
#' @export
summary.smsd_path <- function(object, ...) {
  .path_validate(object)
  tab <- .path_strip(object)
  i <- .path_stop_row(object)
  conf_level <- attr(object, "conf_level", exact = TRUE)
  width <- attr(object, "width", exact = TRUE)
  data.frame(
    method = attr(object, "method", exact = TRUE),
    estimate_term = attr(object, "estimate_term", exact = TRUE),
    stop_n = attr(object, "stop_n", exact = TRUE),
    n = tab$n[i],
    estimate = tab$estimate[i],
    lower_limit = tab$lower_limit[i],
    upper_limit = tab$upper_limit[i],
    current_n = tab$current_n[i],
    criterion_n = tab$criterion_n[i],
    is_satisfied = tab$is_satisfied[i],
    variance_unusable = tab$variance_unusable[i],
    conf_level = if (is.null(conf_level)) NA_real_ else conf_level,
    width = if (is.null(width)) NA_real_ else width,
    stages = nrow(tab),
    stringsAsFactors = FALSE
  )
}

#' @rdname smsd_path_methods
#' @export
as.data.frame.smsd_path <- function(x, row.names = NULL, optional = FALSE,
                                    ...) {
  # Every column the object holds, whatever has happened to it since
  # seq_path() made it: this is the way out for an object the other methods
  # refuse.
  cols <- unclass(x)
  attributes(cols) <- list(names = names(cols))
  as.data.frame(cols, stringsAsFactors = FALSE)
}

# Subsetting. The complete path keeps its class and its attributes (the
# data frame method drops the attributes for some forms of index, such as
# path[names(path)], so they are put back); anything less is returned as a
# plain data frame (or the vector base R would return), so part of a path
# can never be printed as though it were the whole of one.
#' @export
`[.smsd_path` <- function(x, ...) {
  out <- NextMethod()
  if (!is.data.frame(out)) {
    return(out)
  }
  complete <- identical(names(out), names(x)) &&
    identical(unclass(out)$n, unclass(x)$n)
  if (complete) {
    for (a in .path_attrs) {
      attr(out, a) <- attr(x, a, exact = TRUE)
    }
    class(out) <- class(x)
    return(out)
  }
  for (a in intersect(names(attributes(out)), .path_attrs)) {
    attr(out, a) <- NULL
  }
  class(out) <- "data.frame"
  out
}

#' Plot the Stopping Path of a Sequential Procedure
#'
#' Draws the path returned by \code{\link{seq_path}} with base graphics. The
#' \code{"estimate"} panel shows the estimate, with its confidence interval
#' when the procedure has one, against the sample size, so the interval can
#' be seen narrowing toward its target. The \code{"criterion"} panel shows
#' the current sample size against the criterion sample size, so the
#' stopping rule can be seen being met. The first sample size at which the
#' rule is met is marked in both.
#'
#' In the \code{"estimate"} panel the stages after the stopping point are
#' drawn in gray. They show what the rule would have said had sampling
#' continued and are not part of the result: the estimate and interval a
#' sequential design reports are the ones at the marked sample size (see
#' \code{\link{seq_path}}).
#'
#' The \code{"criterion"} panel plots the columns \code{current_n} and
#' \code{criterion_n} of the path, which are on the scale the procedure
#' itself uses (for example the total over both groups for
#' \code{\link{seq_mr_smd}}). For most procedures the rule is met when the
#' current sample size reaches the criterion; where a procedure adds a
#' further condition, the marked sample size follows the procedure's own
#' decision. In particular, a stage at which the variance estimate was zero,
#' negative, or not finite is never the stopping point, whatever the two
#' lines show there (the column \code{variance_unusable} of the path). A
#' procedure that reports no criterion sample size has no
#' \code{"criterion"} panel.
#'
#' When both panels are drawn they are stacked in one column, and the
#' graphical parameters changed to do so are restored on exit. When a single
#' panel is requested no graphical parameter is changed, so the panel can be
#' placed in a layout of the user's own.
#'
#' @param x An object of class \code{"smsd_path"}.
#' @param which Which panels to draw: \code{"estimate"}, \code{"criterion"},
#'   or both (the default; the \code{"criterion"} panel is then left out for
#'   a procedure that reports no criterion sample size).
#' @param main The title, placed over the first panel drawn, above a line
#'   that gives the target (the target width or, for a procedure that only
#'   bounds the length of its interval, the largest length allowed) and the
#'   sample size at which the rule was first met: one string or an
#'   expression. The default is built from the method. Use
#'   \code{main = ""} for neither.
#' @param xlab The label of the horizontal axis: one string or an
#'   expression. Default \code{"Sample size (n)"}.
#' @param ylab The labels of the vertical axes: one string or expression,
#'   or one for each panel drawn (a character vector, an expression vector,
#'   or a list). The defaults are the name of the quantity estimated and
#'   \code{"Sample size"}.
#' @param legend Logical; draw a legend in each panel? Default \code{TRUE}.
#' @param ... Further graphical parameters, applied to every panel drawn:
#'   \code{col}, \code{lwd}, and \code{lty} set the estimate line (and the
#'   current sample size line); anything else (for example \code{xlim},
#'   \code{ylim}, \code{log}, \code{las}, \code{cex.axis}) is passed to
#'   \code{\link[graphics]{plot.default}}. Request a single panel with
#'   \code{which} when a parameter such as \code{ylim} should apply to that
#'   panel only. Room for the legend is left above the data unless
#'   \code{ylim} or a logarithmic vertical axis is asked for; the legend
#'   can then lie over the data, and \code{legend = FALSE} removes it.
#'
#' @return \code{x}, invisibly. Called for its side effect of drawing on the
#'   current graphics device.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @references
#' Kelley, K., Darku, F. B., & Chattopadhyay, B. (2018). Accuracy in
#' parameter estimation for a general class of effect sizes: A sequential
#' approach. \emph{Psychological Methods}, \emph{23}(2), 226--243.
#' \doi{10.1037/met0000127}
#'
#' @examples
#' set.seed(113)
#' # Satisfaction (1 to 7 scale) and weekly hours of product use for 150
#' # customers, in the order surveyed. How does the interval for their
#' # correlation narrow toward a full width of 0.30?
#' use <- rnorm(150, mean = 10, sd = 3)
#' satisfaction <- 4 + 0.15 * (use - 10) + rnorm(150, sd = 0.8)
#' path <- seq_path(seq_aipe_r, list(use, satisfaction), width = 0.30)
#' plot(path)
#'
#' # One panel, with the vertical axis and its label chosen by the user.
#' plot(path, which = "estimate", ylim = c(0, 1.2), las = 1,
#'      ylab = expression(hat(rho)))
#'
#' @seealso \code{\link{seq_path}}, \code{\link{smsd_path_methods}}
#' @family seq_path
#' @export
plot.smsd_path <- function(x, which = c("estimate", "criterion"),
                           main = NULL, xlab = NULL, ylab = NULL,
                           legend = TRUE, ...) {
  .path_validate(x)
  default_which <- missing(which)
  which <- .path_which(which)
  legend <- .check_flag(legend, "legend")
  dots <- list(...)

  tab <- .path_strip(x)
  which <- .path_plot_panels(tab, which, default_which)

  main <- if (is.null(main)) {
    paste("Sequential path:", attr(x, "method", exact = TRUE))
  } else {
    .path_plot_label(main, "main", 1L)[[1L]]
  }
  has_main <- !(is.character(main) && !nzchar(main))
  xlab <- if (is.null(xlab)) {
    "Sample size (n)"
  } else {
    .path_plot_label(xlab, "xlab", 1L)[[1L]]
  }
  ylab <- if (is.null(ylab)) {
    default_ylab <- c(
      estimate = paste0("Estimate (", attr(x, "estimate_term", exact = TRUE),
                        ")"),
      criterion = "Sample size"
    )
    as.list(unname(default_ylab[which]))
  } else {
    .path_plot_label(ylab, "ylab", length(which))
  }

  stacked <- length(which) == 2L
  if (stacked) {
    # Changing the layout also resets cex, so cex is saved and restored
    # with it. The order matters: mfrow first (it resets cex), then cex,
    # then mar (whose size in inches is worked out from cex).
    old <- graphics::par("mfrow", "cex", "mar")
    on.exit(graphics::par(old), add = TRUE)
    graphics::par(mfrow = c(2L, 1L), mar = c(2.6, 4.6, 3.6, 1.1))
  }
  for (j in seq_along(which)) {
    if (stacked && j == 2L) {
      graphics::par(mar = c(4.1, 4.6, 1.6, 1.1))
    }
    .path_panel(x, tab, which[j], main = if (j == 1L && has_main) main,
                xlab = if (stacked && j == 1L) "" else xlab,
                ylab = ylab[[j]], legend = legend, dots = dots)
  }
  invisible(x)
}

# The panels asked for: "estimate", "criterion", or both, in the order
# given. Each name goes through the shared choice check (so it may be
# abbreviated, as with match.arg()); what that check does not cover is that
# two may be given.
#' @noRd
.path_which <- function(which) {
  panels <- c("estimate", "criterion")
  one <- function(panel) {
    .path_also(.check_choice(panel, panels, "which"),
               "Give both, or leave 'which' out, to draw both panels.")
  }
  if (!is.character(which) || length(which) < 2L) {
    return(one(which))
  }
  unique(vapply(which, one, character(1), USE.NAMES = FALSE))
}

# The panels plot.smsd_path() draws: the panels asked for, without the
# "criterion" panel when it was not asked for by name and the procedure
# reports no criterion sample size.
#' @noRd
.path_plot_panels <- function(tab, which, default_which) {
  if ("criterion" %in% which && !any(is.finite(tab$criterion_n))) {
    if (!default_which) {
      .smsd_stop("This procedure reports no criterion sample size, so ",
                 "there is no \"criterion\" panel to draw. Use ",
                 "which = \"estimate\".")
    }
    which <- setdiff(which, "criterion")
  }
  if ("estimate" %in% which && !any(is.finite(tab$estimate))) {
    .smsd_stop("The estimate is not finite at any stage of this path, so ",
               "there is nothing to draw in the \"estimate\" panel.")
  }
  which
}

# A title or an axis label, as a list of `k` labels: each one string or one
# expression. `value` may be a character vector, an expression vector, a
# single call or name (as quote() gives), or a list of these, of length 1 or
# `k`.
#' @noRd
.path_plot_label <- function(value, name, k) {
  labels <- if (is.character(value) || is.expression(value)) {
    lapply(seq_along(value), function(i) value[i])
  } else if (is.call(value) || is.name(value)) {
    list(as.expression(value))
  } else if (is.list(value) && !is.data.frame(value)) {
    lapply(value, function(v) {
      if (is.call(v) || is.name(v)) as.expression(v) else v
    })
  }
  one <- function(v) {
    (is.character(v) && length(v) == 1L && !is.na(v)) ||
      (is.expression(v) && length(v) == 1L)
  }
  ok <- !is.null(labels) && length(labels) %in% c(1L, k) &&
    all(vapply(labels, one, logical(1)))
  if (!ok) {
    .smsd_stop("'", name, "' must be one string or expression",
               if (k > 1L) {
                 paste0(", or one for each panel drawn (", k, ")")
               }, "; got ", .smsd_describe(value), ".")
  }
  rep_len(labels, k)
}

# Colors of the path plot: blue for the estimate and the current sample
# size, orange for the criterion and the stopping marker, and neutral grays
# for everything else. The blue and the orange stay far apart under
# simulated protanopia, deuteranopia, and tritanopia, and the criterion is
# drawn as a step line, so the two series do not rely on hue alone. The tint
# of the interval band is opaque so that it draws on every device.
#' @noRd
.path_colors <- list(
  line = "#2A78D6",
  band = "#D4E4F7",
  criterion = "#EB6834",
  after_line = "gray60",
  after_band = "gray92",
  stop = "gray25",
  grid = "gray93",
  ink = "gray30"
)

# The line under the title: the target and where the rule was first met.
#' @noRd
.path_subtitle <- function(x, tab) {
  width <- attr(x, "width", exact = TRUE)
  conf_level <- attr(x, "conf_level", exact = TRUE)
  stop_n <- attr(x, "stop_n", exact = TRUE)
  target <- c(
    if (!is.null(width)) {
      paste(if (.report_width_is_bound(attr(x, "method", exact = TRUE))) {
        "length at most"
      } else {
        "target width"
      }, .path_arg_text(width))
    },
    if (!is.null(conf_level)) {
      paste0(.path_arg_text(100 * conf_level), "% confidence")
    }
  )
  if (is.null(target)) {
    target <- .path_args_text(attr(x, "arguments", exact = TRUE))
    if (nchar(target) > 50L) target <- ""
  }
  target <- paste(target[nzchar(target)], collapse = ", ")
  met <- if (is.na(stop_n)) {
    paste("rule not met by n =", tab$n[nrow(tab)])
  } else {
    paste("rule first met at n =", stop_n)
  }
  txt <- paste(c(target[nzchar(target)], met), collapse = "; ")
  paste0(toupper(substr(txt, 1L, 1L)), substr(txt, 2L, nchar(txt)))
}

# What one panel of the path plot shows, as data: the series drawn against
# the sample size, the stages drawn as the path proper (`before`: up to and
# including the stopping stage, or every stage when the rule is never met)
# and those drawn in gray in the estimate panel (`after`: from the stopping
# stage on, so that the gray line joins the stopping point; empty when
# nothing follows it), the stopping stage and the point that marks it
# (`marker`: its sample size and its height on the series; NULL when the
# rule is never met), and the range of the vertical axis before any room for
# the legend is added. Values that are not finite are not drawn (NA).
# .path_panel() draws from this list and from nothing else, so that what a
# panel shows can be tested without a graphics device.
#' @noRd
.path_panel_data <- function(tab, stop_n, panel) {
  finite <- function(v) {
    v[!is.finite(v)] <- NA_real_
    v
  }
  k <- nrow(tab)
  n <- tab$n
  i_stop <- if (is.na(stop_n)) NA_integer_ else match(stop_n, n)
  if (panel == "estimate") {
    y <- finite(tab$estimate)
    lower <- finite(tab$lower_limit)
    upper <- finite(tab$upper_limit)
    criterion <- NULL
  } else {
    y <- finite(tab$current_n)
    lower <- upper <- NULL
    criterion <- finite(tab$criterion_n)
  }
  y_range <- range(c(y, lower, upper, criterion), na.rm = TRUE)
  if (diff(y_range) == 0) {
    y_range <- y_range + c(-0.5, 0.5) * max(abs(y_range[1L]), 1)
  }
  list(
    n = n, y = y, lower = lower, upper = upper, criterion = criterion,
    has_band = panel == "estimate" && any(!is.na(lower) & !is.na(upper)),
    i_stop = i_stop,
    marker = if (!is.na(i_stop)) c(n = n[i_stop], y = y[i_stop]),
    before = if (is.na(i_stop)) seq_len(k) else seq_len(i_stop),
    after = if (is.na(i_stop) || i_stop == k) integer(0) else
      seq.int(i_stop, k),
    y_range = y_range
  )
}

# One panel of the path plot. `main` is NULL for a panel without a title.
#' @noRd
.path_panel <- function(x, tab, panel, main, xlab, ylab, legend, dots) {
  colors <- .path_colors
  shown <- .path_panel_data(tab, attr(x, "stop_n", exact = TRUE), panel)
  k <- length(shown$n)
  n <- shown$n
  y <- shown$y
  marker <- shown$marker
  before <- shown$before
  after <- shown$after
  line_col <- if (is.null(dots$col)) colors$line else dots$col
  line_lwd <- if (is.null(dots$lwd)) 2 else dots$lwd
  line_lty <- if (is.null(dots$lty)) 1 else dots$lty
  frame_dots <- dots[setdiff(names(dots), c("col", "lwd", "lty", "type"))]

  # Headroom above the data, so that the legend does not sit on it.
  y_range <- shown$y_range
  headroom <- legend && is.null(frame_dots$ylim) &&
    !any(grepl("y", frame_dots$log, fixed = TRUE))
  data_top <- y_range[2L]
  if (headroom) {
    y_range[2L] <- y_range[2L] + 0.28 * diff(y_range)
  }

  frame <- list(x = range(n), y = y_range, type = "n", main = "",
                xlab = xlab, ylab = ylab, las = 1, bty = "l",
                xlim = if (k == 1L) n + c(-1, 1) else range(n))
  frame[names(frame_dots)] <- frame_dots
  do.call(graphics::plot.default, frame)
  graphics::abline(h = graphics::axTicks(2L), col = colors$grid, lwd = 1)
  graphics::box(bty = frame$bty)
  if (!is.null(main)) {
    # A long method label is shrunk to fit rather than cut off. The title
    # is centered on the plot region, so the room it has is twice the
    # distance from that center to the nearer edge of the figure.
    cex_main <- frame_dots$cex.main
    if (is.null(cex_main)) {
      cex_main <- graphics::par("cex.main")
      center <- mean(graphics::par("plt")[1:2])
      room <- 0.97 * 2 * min(center, 1 - center)
      needed <- graphics::strwidth(main, units = "figure", cex = cex_main,
                                   font = graphics::par("font.main"))
      if (needed > room) cex_main <- max(0.5, cex_main * room / needed)
    }
    graphics::title(main = main, line = 2, cex.main = cex_main)
    graphics::mtext(.path_subtitle(x, tab), side = 3L, line = 0.6,
                    cex = 0.85 * graphics::par("cex"), col = colors$ink)
  }

  show_points <- k <= 30L
  if (panel == "estimate") {
    if (shown$has_band) {
      .path_band(n, shown$lower, shown$upper, after, colors$after_band)
      .path_band(n, shown$lower, shown$upper, before, colors$band)
    }
    if (length(after)) {
      graphics::lines(n[after], y[after], col = colors$after_line,
                      lwd = line_lwd, lty = line_lty,
                      type = if (show_points) "o" else "l", pch = 16,
                      cex = 0.6)
    }
    graphics::lines(n[before], y[before], col = line_col, lwd = line_lwd,
                    lty = line_lty, type = if (show_points) "o" else "l",
                    pch = 16, cex = 0.6)
  } else {
    graphics::lines(n, shown$criterion, col = colors$criterion,
                    lwd = line_lwd, type = if (k == 1L) "p" else "s",
                    pch = 16)
    graphics::lines(n, y, col = line_col, lwd = line_lwd, lty = line_lty,
                    type = if (show_points) "o" else "l", pch = 16,
                    cex = 0.6)
  }

  key_args <- NULL
  if (legend) {
    conf_level <- attr(x, "conf_level", exact = TRUE)
    if (panel == "estimate") {
      key <- list(
        label = "Estimate", col = line_col, lty = line_lty, lwd = line_lwd,
        pch = NA, pt_bg = NA, pt_cex = 1
      )
      if (shown$has_band) {
        key <- .path_key(key, if (is.null(conf_level)) "Interval" else
          paste0(.path_arg_text(100 * conf_level),
                 "% confidence interval"),
          col = colors$band, lty = NA, lwd = NA, pch = 15, pt_cex = 1.8)
      }
    } else {
      key <- list(
        label = "Current sample size", col = line_col, lty = line_lty,
        lwd = line_lwd, pch = NA, pt_bg = NA, pt_cex = 1
      )
      key <- .path_key(key, "Criterion sample size", col = colors$criterion,
                       lty = 1, lwd = line_lwd, pch = NA, pt_cex = 1)
    }
    if (!is.null(marker)) {
      key <- .path_key(key, paste0("Rule first met (n = ", marker[["n"]],
                                   ")"),
                       col = colors$stop, lty = 2, lwd = 1, pch = 21,
                       pt_bg = colors$criterion, pt_cex = 1.1)
      if (length(after) && panel == "estimate") {
        key <- .path_key(key, "After stopping (for illustration)",
                         col = colors$after_line, lty = line_lty,
                         lwd = line_lwd, pch = NA, pt_cex = 1)
      }
    }
    key_args <- list("topright", legend = key$label, col = key$col,
                     lty = key$lty, lwd = key$lwd, pch = key$pch,
                     pt.bg = key$pt_bg, pt.cex = key$pt_cex, bty = "n",
                     cex = 0.8, ncol = if (length(key$label) > 2L) 2L else
                       1L, text.col = colors$ink)
  }

  if (!is.null(marker)) {
    # The stopping marker: a dashed line up from the axis and a point on the
    # series. Positions are worked out in the coordinates of par("usr"),
    # which are logarithms on a log axis.
    y_log <- graphics::par("ylog")
    on_axis <- function(v) {
      if (y_log) log10(ifelse(is.na(v) | v <= 0, NA_real_, v)) else v
    }
    from_axis <- function(v) if (y_log) 10^v else v
    edge <- graphics::par("usr")[3:4]
    box <- if (!is.null(key_args)) {
      do.call(graphics::legend, c(key_args, list(plot = FALSE)))$rect
    }
    top <- .path_marker_top(
      limit = if (headroom) on_axis(data_top) else edge[2L],
      point = on_axis(marker[["y"]]),
      at = if (graphics::par("xlog")) log10(marker[["n"]]) else
        marker[["n"]],
      box = box
    )
    graphics::segments(marker[["n"]], from_axis(edge[1L]), marker[["n"]],
                       from_axis(top), col = colors$stop, lty = 2, lwd = 1)
    graphics::points(marker[["n"]], marker[["y"]], pch = 21,
                     bg = colors$criterion, col = "white", cex = 1.5,
                     lwd = 1.5)
  }
  if (!is.null(key_args)) {
    do.call(graphics::legend, key_args)
  }
  invisible(NULL)
}

# Where the dashed line of the stopping marker ends. It runs up to `limit`
# (the top of the data when room was left above them for the legend, the
# top of the plot otherwise), but stops below the legend's box wherever it
# would run through it (a `ylim` of the user's own leaves no room for the
# legend), though never below the point it marks. `box` is the `rect` that
# graphics::legend() returns (left, top, w, h), or NULL without a legend;
# `at` is the horizontal position of the line. All in the coordinates of
# par("usr").
#' @noRd
.path_marker_top <- function(limit, point, at, box = NULL) {
  top <- limit
  if (!is.null(box) && at >= box$left && at <= box$left + box$w) {
    top <- min(top, box$top - box$h)
  }
  max(top, point, na.rm = TRUE)
}

# Append one entry to a legend key.
#' @noRd
.path_key <- function(key, label, col, lty, lwd, pch, pt_bg = NA,
                      pt_cex = 1) {
  list(label = c(key$label, label), col = c(key$col, col),
       lty = c(key$lty, lty), lwd = c(key$lwd, lwd), pch = c(key$pch, pch),
       pt_bg = c(key$pt_bg, pt_bg), pt_cex = c(key$pt_cex, pt_cex))
}

# The interval as a filled band over each unbroken run of stages in `idx`
# (a single stage is drawn as a vertical bar).
#' @noRd
.path_band <- function(n, lo, hi, idx, col) {
  if (!length(idx)) {
    return(invisible(NULL))
  }
  ok <- !is.na(lo[idx]) & !is.na(hi[idx])
  runs <- rle(ok)
  ends <- cumsum(runs$lengths)
  starts <- ends - runs$lengths + 1L
  for (r in which(runs$values)) {
    i <- idx[starts[r]:ends[r]]
    if (length(i) == 1L) {
      graphics::segments(n[i], lo[i], n[i], hi[i], col = col, lwd = 6,
                         lend = 1)
    } else {
      graphics::polygon(c(n[i], rev(n[i])), c(lo[i], rev(hi[i])),
                        col = col, border = NA)
    }
  }
  invisible(NULL)
}

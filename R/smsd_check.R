#' Construct an smsd_check object
#'
#' Internal constructor for the tidy data frame returned by the sequential
#' check functions. Each row is one quantity (e.g., the current estimate,
#' the lower and upper confidence limits, the asymptotic variance estimate,
#' the criterion sample size).
#'
#' The function is not exported. Method authors should call it to assemble
#' the return value of a check function rather than constructing the data
#' frame by hand, so that the class, attribute, and column type contract
#' is maintained in one place.
#'
#' @param term A character vector of term names.
#' @param value A numeric vector of values, the same length as \code{term}.
#'   Each is a finite number, or \code{NA} where the quantity is not
#'   defined; \code{Inf}, \code{-Inf}, and \code{NaN} are refused with a
#'   message that names the row.
#' @param method A short character label identifying the method
#'   (e.g., \code{"standardized mean difference"}).
#' @param is_satisfied A length one logical indicating whether the
#'   sequential stopping rule has been satisfied by the supplied data.
#'   Never \code{NA}.
#' @param current_n A length one numeric giving the current per-group
#'   (or total, for univariate methods) sample size on which the check
#'   is based.
#' @param criterion_n A length one numeric giving the criterion sample
#'   size implied by the procedure given the current variance estimate.
#' @param conf_level The confidence level used for the interval (numeric
#'   in (0, 1)), or \code{NULL} for risk based procedures.
#' @param width The target full width of the confidence interval, when
#'   applicable. Numeric; pass \code{NULL} for procedures that
#'   parameterize precision differently.
#' @param extra A named list of additional attributes to attach (e.g., a
#'   correlation method label, a contrast weight vector,
#'   \code{variance_unusable}). Optional.
#' @param notes A character vector of plain-language notes for the user
#'   (for example that the variance estimate is not usable, or that the budget
#'   is exhausted). Printed under the stopping line. Optional.
#'
#' @return An object of class \code{c("smsd_check", "data.frame")} with
#'   columns \code{term} (character) and \code{value} (numeric). Metadata
#'   is stored as attributes.
#'
#' @keywords internal
#' @noRd
new_smsd_check <- function(term, value, method, is_satisfied, current_n,
                           criterion_n, conf_level = NULL, width = NULL,
                           extra = NULL, notes = NULL) {
  stopifnot(
    is.character(term),
    is.numeric(value),
    length(term) == length(value),
    is.character(method), length(method) == 1L,
    is.logical(is_satisfied), length(is_satisfied) == 1L,
    !is.na(is_satisfied),
    is.numeric(current_n), length(current_n) == 1L,
    is.numeric(criterion_n), length(criterion_n) == 1L
  )
  if (!is.null(conf_level)) {
    stopifnot(is.numeric(conf_level), length(conf_level) == 1L,
              conf_level > 0, conf_level < 1)
  }
  if (!is.null(width)) {
    stopifnot(is.numeric(width), length(width) == 1L, width > 0)
  }
  if (is.null(notes)) notes <- character(0)
  stopifnot(is.character(notes))
  # A row is a finite number, or NA where the quantity is not defined; never
  # Inf or NaN. Sums and differences of data and targets at the top of the
  # double precision range overflow (an interval from -1.7e308 to 1.7e308
  # has no representable length). That is an error that says to change the
  # units, not a result.
  bad <- which(is.infinite(value) | is.nan(value))
  if (length(bad)) {
    .smsd_stop("The '", term[bad[1L]], "' row of the result is too large to ",
               "represent (", format(value[bad[1L]]), "). Express the data ",
               "in larger units (divide them by a constant) and state the ",
               "precision target in the same units.")
  }
  out <- data.frame(term = term, value = as.numeric(value),
                    stringsAsFactors = FALSE)
  attr(out, "method") <- method
  attr(out, "is_satisfied") <- is_satisfied
  attr(out, "current_n") <- current_n
  attr(out, "criterion_n") <- criterion_n
  attr(out, "conf_level") <- conf_level
  attr(out, "width") <- width
  attr(out, "notes") <- notes
  if (!is.null(extra)) {
    for (nm in names(extra)) {
      attr(out, nm) <- extra[[nm]]
    }
  }
  class(out) <- c("smsd_check", "data.frame")
  out
}

# Names of the attributes every smsd_check carries (beyond data frame ones).
#' @noRd
.smsd_core_attrs <- c("method", "is_satisfied", "current_n", "criterion_n",
                      "conf_level", "width", "notes")

# Drop the class and every sequential attribute, leaving the plain
# term / value data frame.
#' @noRd
.smsd_strip <- function(x) {
  data.frame(term = x$term, value = x$value, stringsAsFactors = FALSE)
}

# The row holding the point estimate: the row between lower_limit and
# upper_limit for interval procedures, otherwise the first row.
#' @noRd
.smsd_estimate_row <- function(x) {
  lo <- match("lower_limit", x$term)
  hi <- match("upper_limit", x$term)
  if (!is.na(lo) && !is.na(hi) && hi == lo + 2L) lo + 1L else 1L
}

# What is wrong with an object that has class smsd_check but is no longer
# one check as a check function returned it; NULL when nothing is. Subsetting
# through `[` returns a plain data frame for anything less than the whole
# check, but other operations keep the class: rbind() of several checks
# returns their rows under the attributes of the first, so that a log of the
# stages of a study would print, summarize, and be reported as its first
# stage alone, and an attribute can be removed by hand. A check does not
# record which terms it should have, so part of a check cannot be
# recognized; what can be is a term that appears more than once (no check
# function returns one), the two columns, and the attributes every method
# reads.
#' @noRd
.smsd_check_problem <- function(x) {
  cols <- unclass(x)
  if (!is.list(cols) || !all(c("term", "value") %in% names(cols)) ||
      !is.character(cols$term) || !is.numeric(cols$value) ||
      length(cols$term) != length(cols$value)) {
    return("its columns are not 'term' (character) and 'value' (numeric)")
  }
  one <- function(v, is_type) is_type(v) && length(v) == 1L
  ok <- c(
    method = one(attr(x, "method", exact = TRUE), is.character),
    is_satisfied = one(attr(x, "is_satisfied", exact = TRUE), is.logical) &&
      !is.na(attr(x, "is_satisfied", exact = TRUE)),
    current_n = one(attr(x, "current_n", exact = TRUE), is.numeric),
    criterion_n = one(attr(x, "criterion_n", exact = TRUE), is.numeric)
  )
  if (!all(ok)) {
    bad <- names(ok)[!ok]
    return(paste0(if (length(bad) == 1L) "the attribute " else
      "the attributes ", paste(bad, collapse = ", "),
      if (length(bad) == 1L) " is" else " are",
      " missing or not as the check function set ",
      if (length(bad) == 1L) "it" else "them"))
  }
  repeated <- cols$term[duplicated(cols$term)]
  if (length(repeated)) {
    return(paste0("the term ", repeated[1L], " appears ",
                  sum(cols$term == repeated[1L]),
                  " times, so it holds more than one check"))
  }
  NULL
}

# Stop unless `x` is one check as a check function returned it.
#' @noRd
.smsd_check_validate <- function(x) {
  problem <- .smsd_check_problem(x)
  if (!is.null(problem)) {
    .smsd_stop("This object has class 'smsd_check' but is no longer a check ",
               "as the check function returned it (", problem, "). It was ",
               "changed afterward, for example by binding checks together ",
               "or by changing an attribute. Use as.data.frame() for the ",
               "plain table of terms; to keep a record of the stages of a ",
               "study, bind the one row summary() of each check, or use ",
               "seq_path().")
  }
  invisible(x)
}

# The text after "Stopping rule met:". For every procedure but one it is the
# stopping flag. The budget constrained procedure ends in two ways, and its
# flag is TRUE for both: every group has met its criterion (rule_met), or
# the budget cannot pay for another needed observation (budget_exhausted)
# while a group is still below its criterion. "TRUE" above a criterion the
# sample has not reached invites the reading that the rule was met, so the
# second ending is printed as what it is. Display only: is_satisfied and the
# other attributes are as seq_budget() set them. (With an unusable variance
# the flag is FALSE as well, and the line reads the same.)
#' @noRd
.smsd_check_status <- function(x) {
  if (isTRUE(attr(x, "budget_exhausted", exact = TRUE)) &&
      !isTRUE(attr(x, "rule_met", exact = TRUE))) {
    return("FALSE (sampling has ended: the budget is exhausted)")
  }
  if (isTRUE(attr(x, "is_satisfied", exact = TRUE))) "TRUE" else "FALSE"
}

#' Methods for Sequential Check Objects
#'
#' Methods for the \code{smsd_check} objects returned by the sequential
#' check functions (the \code{seq_aipe_*} accuracy in parameter estimation
#' family, the \code{seq_mr_*} and \code{seq_br_*} minimum risk and bounded
#' risk families, the fixed width and two stage procedures, and the budget
#' constrained method).
#'
#' \code{print()} shows the method, the current and criterion sample sizes,
#' the target width and confidence level when the procedure has them, every
#' term, whether the stopping rule is met, and any notes (for example that
#' the variance estimate is not yet usable, in which case the rule is
#' reported as not met). For a check from \code{\link{seq_budget}} at which
#' sampling has ended because the budget cannot pay for another needed
#' observation, with a group still below its criterion, the stopping line
#' reads \code{FALSE (sampling has ended: the budget is exhausted)},
#' because the rule itself was not met; the attribute \code{is_satisfied},
#' which says whether sampling has ended, is \code{TRUE}. \code{format()}
#' returns the same lines as a character vector.
#' \code{summary()} and \code{glance()} return a one row data frame
#' describing the check; \code{tidy()} returns the point estimate with its
#' confidence limits in the column layout used by the sibling package
#' \pkg{DMAR}. \code{as.data.frame()} returns the plain
#' \code{term}/\code{value} data frame without the class or the sequential
#' attributes. Subsetting with \code{[} returns a plain data frame (or
#' vector) whenever the result is no longer the complete check, so a partial
#' object is never printed as if it were a complete one.
#'
#' \strong{Width and coverage.} An AIPE check carries
#' \code{attr(x, "coverage_guarantee") = "asymptotic"}. Its display and
#' \code{\link{seq_report}} distinguish the asymptotic confidence level
#' from satisfying the width criterion. The latter does not establish
#' nominal coverage at the observed sample size. See
#' \code{vignette("coverage-and-design", package = "SMSD")} for measured
#' operating characteristics and simulation uncertainty.
#'
#' \strong{A check that cannot be evaluated.} In a small sample the estimate
#' of the variance that a stopping rule needs can be zero, negative, or not
#' finite. No variance has then been estimated: the stopping flag is
#' \code{FALSE}, \code{attr(x, "variance_unusable")} is \code{TRUE}, and a
#' note says so. This is a deliberate departure from the stopping rules as
#' printed in the source articles, under which such an estimate, or the lower
#' bound \eqn{n^{-3}} that some of them put in its place, would make the rule
#' be met at once. A positive estimate, however small, is used as the articles
#' prescribe (with that bound where an article has one). When the point
#' estimate itself is not finite, the flag is \code{FALSE} as well and a note
#' gives the reason. (The two stage \code{ts_*} functions whose total sample
#' size the first stage variance or scale fixes stop with an error on a first
#' stage with no variation instead; see \code{\link{SMSD}}.) Test
#' \code{attr(x, "is_satisfied")}, which is always \code{TRUE} or
#' \code{FALSE}, to drive a sampling loop. What the limit rows then hold
#' depends on the procedure. The \code{seq_aipe_*} functions report \code{NA}
#' for both limits. The \code{seq_fwd_*} functions still report their limits,
#' which are computed from the estimate and \code{d} alone and involve no
#' variance. \code{\link{seq_budget}} still reports an interval, computed with
#' a variance of zero for each group whose variance is unusable, which is too
#' narrow. Neither of the last two is an interval to report. \code{summary()},
#' \code{glance()}, and \code{tidy()} return the limits as stored, without the
#' notes or \code{variance_unusable}, so test
#' \code{attr(x, "variance_unusable")} before using them.
#'
#' \strong{An object that is no longer one check.} Some operations leave
#' the class on something that is not a check as the function returned it:
#' \code{rbind()} of two or more checks (which keeps the attributes of the
#' first), a data frame verb of another package that binds checks, or an
#' attribute removed by hand. \code{print()}, \code{format()},
#' \code{summary()}, \code{tidy()}, \code{glance()}, and
#' \code{\link{seq_report}} stop with a message that says so, rather than
#' describe the object as if it were one complete check.
#' \code{as.data.frame()} always works and gives the plain table. To keep a
#' record of the stages of a study, bind the one row \code{summary()} of
#' each check, or use \code{\link{seq_path}}, which records every stage.
#'
#' @param x,object An object of class \code{"smsd_check"}.
#' @param digits Number of significant digits to show for numeric values;
#'   by default the session's \code{getOption("digits")}. Stored values are
#'   never rounded; this affects display only.
#' @param row.names,optional Ignored; present for compatibility with the
#'   generic.
#' @param ... Not used.
#'
#' @return \code{print()} returns its input invisibly. \code{format()}
#'   returns a character vector, one element per printed line.
#'   \code{summary()} and \code{glance()} return a one row
#'   \code{data.frame} with columns \code{method}, \code{estimate_term},
#'   \code{estimate}, \code{lower_limit}, \code{upper_limit},
#'   \code{current_n}, \code{criterion_n}, \code{is_satisfied},
#'   \code{conf_level}, and \code{width} (\code{NA} where the procedure has
#'   no such quantity). \code{tidy()} returns a one row \code{data.frame}
#'   with columns \code{term}, \code{estimate}, \code{ci_lower},
#'   \code{ci_upper}, and \code{conf_level}. \code{as.data.frame()} returns
#'   a \code{data.frame} with columns \code{term} and \code{value}.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @examples
#' set.seed(113)
#' group_1 <- rnorm(40, mean = 0.5)
#' group_2 <- rnorm(40, mean = 0)
#' out <- seq_aipe_smd(group_1, group_2, width = 0.50)
#' out
#' summary(out)
#' generics::tidy(out)
#' as.data.frame(out)
#' out[out$term == "smd", "value"]
#'
#' @name smsd_check_methods
#' @aliases smsd_check
NULL

#' @rdname smsd_check_methods
#' @importFrom DMAR smd
#' @export
format.smsd_check <- function(x, digits = getOption("digits"), ...) {
  .smsd_check_validate(x)
  conf_level <- attr(x, "conf_level")
  width <- attr(x, "width")
  notes <- attr(x, "notes")

  lines <- c(
    sprintf("Sequential check (%s)", attr(x, "method")),
    "",
    sprintf("  Current n       : %s", format(attr(x, "current_n"))),
    sprintf("  Criterion n     : %s", format(attr(x, "criterion_n")))
  )
  if (!is.null(width)) {
    lines <- c(lines, sprintf("  Target width    : %s",
                              format(width, digits = digits)))
  }
  if (!is.null(conf_level)) {
    lines <- c(lines, sprintf("  Confidence level: %s",
                              paste0(format(100 * conf_level, trim = TRUE),
                                     "%")))
  }
  if (identical(attr(x, "coverage_guarantee", exact = TRUE), "asymptotic")) {
    lines <- c(lines, "  Coverage        : asymptotic; assess with sensitivity analysis")
  }
  lines <- c(lines, "")
  for (i in seq_along(x$term)) {
    lines <- c(lines, sprintf("  %-18s: %s", x$term[i],
                              format(x$value[i], digits = digits)))
  }
  lines <- c(lines, "",
             sprintf("  Stopping rule met: %s", .smsd_check_status(x)))
  if (length(notes)) {
    for (note in notes) {
      lines <- c(lines, "",
                 strwrap(paste("Note:", note), width = 72, indent = 2,
                         exdent = 8))
    }
  }
  lines
}

#' @rdname smsd_check_methods
#' @export
print.smsd_check <- function(x, digits = getOption("digits"), ...) {
  cat(format(x, digits = digits), sep = "\n")
  invisible(x)
}

#' @rdname smsd_check_methods
#' @export
summary.smsd_check <- function(object, ...) {
  .smsd_check_validate(object)
  i <- .smsd_estimate_row(object)
  pick <- function(nm) {
    j <- match(nm, object$term)
    if (is.na(j)) NA_real_ else object$value[j]
  }
  conf_level <- attr(object, "conf_level")
  width <- attr(object, "width")
  data.frame(
    method = attr(object, "method"),
    estimate_term = object$term[i],
    estimate = object$value[i],
    lower_limit = pick("lower_limit"),
    upper_limit = pick("upper_limit"),
    current_n = attr(object, "current_n"),
    criterion_n = attr(object, "criterion_n"),
    is_satisfied = isTRUE(attr(object, "is_satisfied")),
    conf_level = if (is.null(conf_level)) NA_real_ else conf_level,
    width = if (is.null(width)) NA_real_ else width,
    stringsAsFactors = FALSE
  )
}

#' @rdname smsd_check_methods
#' @export
as.data.frame.smsd_check <- function(x, row.names = NULL, optional = FALSE,
                                     ...) {
  .smsd_strip(x)
}

#' @rdname smsd_check_methods
#' @importFrom generics tidy
#' @exportS3Method generics::tidy
tidy.smsd_check <- function(x, ...) {
  s <- summary.smsd_check(x)
  data.frame(
    term = s$estimate_term,
    estimate = s$estimate,
    ci_lower = s$lower_limit,
    ci_upper = s$upper_limit,
    conf_level = s$conf_level,
    stringsAsFactors = FALSE
  )
}

#' @rdname smsd_check_methods
#' @importFrom generics glance
#' @exportS3Method generics::glance
glance.smsd_check <- function(x, ...) {
  summary.smsd_check(x)
}

# Subsetting. A complete check keeps its class; anything less is returned as
# a plain data frame (or the vector base R would return), so a partial
# object can never be printed as though it were the whole check.
#' @export
`[.smsd_check` <- function(x, ...) {
  out <- NextMethod()
  if (!is.data.frame(out)) {
    return(out)
  }
  complete <- identical(names(out), c("term", "value")) &&
    identical(out$term, x$term)
  if (complete) {
    return(out)
  }
  attrs <- setdiff(names(attributes(out)), c("names", "row.names", "class"))
  for (a in attrs) attr(out, a) <- NULL
  class(out) <- "data.frame"
  out
}

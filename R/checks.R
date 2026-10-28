# Internal argument checks shared by every exported function. Not exported.
#
# Each check stops with a message that names the argument, says what is
# required, shows what was received, and, where it helps, says what to do
# next. A failed check never surfaces as R's own "missing value where
# TRUE/FALSE needed" or "the condition has length > 1".

#' @noRd
.smsd_stop <- function(...) {
  stop(..., call. = FALSE)
}

# Short description of an offending value, for use in messages.
#' @noRd
.smsd_describe <- function(x) {
  if (is.null(x)) return("NULL")
  if (!is.null(dim(x))) {
    return(sprintf("a %s with dimensions %s",
                   if (is.data.frame(x)) "data frame" else "matrix or array",
                   paste(dim(x), collapse = " x ")))
  }
  if (length(x) != 1L) {
    cls <- class(x)[1L]
    return(sprintf("%s %s vector of length %d",
                   if (grepl("^[aeiou]", cls)) "an" else "a", cls,
                   length(x)))
  }
  if (is.numeric(x) || is.logical(x)) return(format(x))
  if (is.character(x)) return(sprintf("\"%s\"", x))
  sprintf("an object of class '%s'", class(x)[1L])
}

# A number that carries a class is refused by the scalar checks, with advice
# to convert it. Such a number can have comparison and arithmetic methods of
# its own: a 64 bit integer of package bit64 (what data.table::fread()
# returns for large whole numbers) compares as NA with Inf, and arithmetic
# with it truncates to whole numbers. The procedures use their scalar
# arguments in arithmetic as given, so a check that accepted such a number
# would pass it on to give a wrong answer. Data vectors are different:
# .check_data() and .check_vector() return them converted to doubles.
#' @noRd
.stop_if_classed <- function(x, name, what = "a plain number") {
  if (is.object(x)) {
    .smsd_stop("'", name, "' must be ", what, ", not an object of class '",
               class(x)[1L], "'. Pass as.numeric(", name, ").")
  }
  invisible(NULL)
}

# A single finite number, optionally bounded. Bounds are closed unless the
# matching *_open flag is TRUE. Returns the value as a double.
#' @noRd
.check_number <- function(x, name, lower = -Inf, upper = Inf,
                          lower_open = FALSE, upper_open = FALSE) {
  ok <- is.numeric(x) && is.null(dim(x)) && length(x) == 1L
  if (ok) .stop_if_classed(x, name)
  ok <- ok && is.finite(x)
  if (ok) {
    ok <- (if (lower_open) x > lower else x >= lower) &&
      (if (upper_open) x < upper else x <= upper)
  }
  if (!ok) {
    range_txt <- if (is.finite(lower) && is.finite(upper)) {
      sprintf(" %s %s and %s %s",
              if (lower_open) "greater than" else "at least", format(lower),
              if (upper_open) "less than" else "at most", format(upper))
    } else if (is.finite(lower)) {
      sprintf(" %s %s",
              if (lower_open) "greater than" else "at least", format(lower))
    } else if (is.finite(upper)) {
      sprintf(" %s %s",
              if (upper_open) "less than" else "at most", format(upper))
    } else {
      ""
    }
    .smsd_stop("'", name, "' must be a single finite number", range_txt,
               "; got ", .smsd_describe(x), ".")
  }
  as.numeric(x)
}

# A single finite number strictly greater than zero.
#' @noRd
.check_positive <- function(x, name) {
  .check_number(x, name, lower = 0, lower_open = TRUE)
}

# The target full width of a confidence interval.
#' @noRd
.check_width <- function(width, name = "width") {
  if (missing(width)) {
    .smsd_stop("'", name, "' must be supplied: the target full width of ",
               "the confidence interval.")
  }
  .check_positive(width, name)
}

# A confidence level strictly between 0 and 1.
#' @noRd
.check_conf_level <- function(conf_level, name = "conf_level") {
  .check_number(conf_level, name, lower = 0, upper = 1,
                lower_open = TRUE, upper_open = TRUE)
}

# A whole number of at least `min` (a sample size, a replication count).
# Returns an integer. A number with a class is refused (.stop_if_classed()).
#' @noRd
.check_count <- function(x, name, min = 1L) {
  ok <- is.numeric(x) && is.null(dim(x)) && length(x) == 1L
  if (ok) .stop_if_classed(x, name, "a plain whole number")
  ok <- ok && is.finite(x) && x == round(x) && x >= min &&
    x <= .Machine$integer.max
  if (!ok) {
    .smsd_stop("'", name, "' must be a single whole number of at least ",
               format(min), "; got ", .smsd_describe(x), ".")
  }
  as.integer(x)
}

# TRUE or FALSE.
#' @noRd
.check_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    .smsd_stop("'", name, "' must be TRUE or FALSE; got ",
               .smsd_describe(x), ".")
  }
  x
}

# A one dimensional numeric array is a vector that carries a dim attribute:
# it is what tapply() returns (one mean per participant, one mean or standard
# deviation per group), and mean(), sd(), and t.test() take it as a vector.
# The dim attribute and the names are dropped so that the checks below treat
# it as one. Only an array with no class is converted: a one dimensional
# table() holds counts of values, not the values, and a classed object (a
# table, a by() result, a Date with a dim) is returned as it is, to be
# refused, as are a matrix, a higher array, and a data frame.
#' @noRd
.drop_one_dim <- function(x) {
  if (length(dim(x)) == 1L && is.numeric(x) && !is.object(x)) {
    as.vector(x)
  } else {
    x
  }
}

# A numeric data vector. Returns the cleaned vector as doubles (so integer
# input cannot overflow in sums and products downstream). A one dimensional
# numeric array counts as a vector (.drop_one_dim()).
#
# na.rm: TRUE removes missing values; FALSE stops if any are present; NULL
#   is for functions with no na.rm argument and stops with advice to remove
#   or impute.
# min_n: smallest number of (non missing) observations accepted.
# positive / nonnegative: sign requirements of the method.
# advice: sentence closing the too-few-observations message; the default suits
#   a sequential check, pass NULL for a plain estimator.
#' @noRd
.check_data <- function(x, name, min_n = 2L, na.rm = NULL,
                        positive = FALSE, nonnegative = FALSE,
                        advice = paste("Collect more data, then check the",
                                       "stopping rule again.")) {
  x <- .drop_one_dim(x)
  if (!is.null(dim(x))) {
    # The advice to pass a column is given only where x[, 1] can work.
    .smsd_stop("'", name, "' must be a numeric vector, not ",
               .smsd_describe(x), ".",
               if (length(dim(x)) == 2L) {
                 paste0(" Pass a single column, for example ", name,
                        "[, 1] or ", name, "[[\"column\"]].")
               })
  }
  if (!is.numeric(x)) {
    .smsd_stop("'", name, "' must be a numeric vector; got ",
               .smsd_describe(x), ".")
  }
  x <- as.numeric(x)
  if (anyNA(x)) {
    if (isTRUE(na.rm)) {
      x <- x[!is.na(x)]
    } else if (is.null(na.rm)) {
      .smsd_stop("'", name, "' contains missing values (NA or NaN). ",
                 "Remove or impute them before calling.")
    } else {
      .smsd_stop("'", name, "' contains missing values (NA or NaN). ",
                 "Remove them, or set na.rm = TRUE.")
    }
  }
  if (any(!is.finite(x))) {
    .smsd_stop("'", name, "' contains infinite values. Remove them before ",
               "calling.")
  }
  if (length(x) < min_n) {
    .smsd_stop("'", name, "' needs at least ", min_n,
               if (min_n == 1L) " observation" else " observations",
               "; it has ", length(x), ".",
               if (!is.null(advice)) paste0(" ", advice))
  }
  if (positive && any(x <= 0)) {
    .smsd_stop("'", name, "' must contain only positive values for this ",
               "procedure; ", sum(x <= 0), " of ", length(x),
               " are zero or negative.")
  }
  if (nonnegative && any(x < 0)) {
    .smsd_stop("'", name, "' must contain no negative values for this ",
               "procedure; ", sum(x < 0), " of ", length(x),
               " are negative.")
  }
  x
}

# Round a criterion sample size up to an integer. Vectorized.
#
# Two protections. A value too large to represent (which would otherwise
# coerce to NA and leak into the stopping flag) is an error whose closing
# sentence the caller supplies through `advice`. And a criterion that is a
# whole number in exact arithmetic but lands a few units in the last place
# above it in floating point (56.00000000000001 for an exact 56) is not
# rounded up to the next integer: the ceiling is taken after shrinking the
# value by one part in 1e10, far above accumulated rounding error and far
# below any statistical meaning (it can move the result only when the
# fractional part is below 1e-10 times the value, that is, never for a
# sample size under ten million unless the excess is rounding noise).
#' @noRd
.as_sample_size <- function(x, what = "sample size",
                            advice = paste("Check the precision target and",
                                           "the other arguments.")) {
  bad <- !is.finite(x) | x > .Machine$integer.max
  if (any(bad)) {
    .smsd_stop("The implied ", what, " (", format(x[bad][1L]),
               ") is too large to represent.",
               if (!is.null(advice)) paste0(" ", advice))
  }
  as.integer(ceiling(x * (1 - 1e-10)))
}

# One of a fixed set of character choices, with partial matching. Unlike
# match.arg(), the message names the argument as the user knows it. Passing
# the full vector of choices (the unevaluated default) selects the first.
#' @noRd
.check_choice <- function(x, choices, name) {
  if (identical(x, choices)) return(choices[1L])
  if (is.character(x) && length(x) == 1L && !is.na(x)) {
    i <- pmatch(x, choices)
    if (!is.na(i)) return(choices[i])
  }
  .smsd_stop("'", name, "' must be one of ",
             paste0("\"", choices, "\"", collapse = ", "), "; got ",
             .smsd_describe(x), ".")
}

# A seed for a reproducible simulation: NULL (use the caller's stream) or a
# single whole number that set.seed() accepts. Returns NULL or an integer.
#' @noRd
.check_seed <- function(seed) {
  if (is.null(seed)) return(NULL)
  ok <- is.numeric(seed) && is.null(dim(seed)) && length(seed) == 1L
  if (ok) .stop_if_classed(seed, "seed", "NULL or a plain whole number")
  ok <- ok && is.finite(seed) && seed == round(seed) &&
    abs(seed) <= .Machine$integer.max
  if (!ok) {
    .smsd_stop("'seed' must be NULL or a single whole number that ",
               "set.seed() accepts; got ", .smsd_describe(seed), ".")
  }
  as.integer(seed)
}

# A numeric vector of finite values with no dimensions (weights, costs, true
# means, one entry per group). `what` says in the message what the entries
# are; `len`, when given, is the required length (for example the number of
# groups) and `len_what` names what it must match. Returns doubles. A one
# dimensional numeric array counts as a vector (.drop_one_dim()), so group
# means or standard deviations from tapply() are accepted.
#' @noRd
.check_vector <- function(x, name, what, len = NULL, len_what = NULL,
                          positive = FALSE) {
  x <- .drop_one_dim(x)
  ok <- is.numeric(x) && is.null(dim(x)) && length(x) >= 1L
  if (!ok) {
    .smsd_stop("'", name, "' must be a numeric vector (", what, "); got ",
               if (length(x) == 0L && is.null(dim(x)) && !is.null(x)) {
                 "a vector of length 0"
               } else {
                 .smsd_describe(x)
               }, ".")
  }
  # Converted before the values are checked, so that the checks below compare
  # doubles whatever class of number was passed.
  x <- as.numeric(x)
  bad <- which(!is.finite(x))
  if (length(bad)) {
    .smsd_stop("'", name, "' must contain only finite numbers (", what,
               "); entry ", bad[1L], " is ", format(x[bad[1L]]),
               if (length(bad) > 1L) {
                 paste0(" (", length(bad), " entries in all are missing or ",
                        "infinite)")
               }, ".")
  }
  if (!is.null(len) && length(x) != len) {
    .smsd_stop("'", name, "' must have ", len,
               if (len == 1L) " entry" else " entries",
               if (!is.null(len_what)) paste0(" (", len_what, ")"),
               "; it has ", length(x), ".")
  }
  if (positive && any(x <= 0)) {
    i <- which(x <= 0)[1L]
    .smsd_stop("'", name, "' must contain only positive numbers (", what,
               "); entry ", i, " is ", format(x[i]), ".")
  }
  x
}

# Is a sample constant up to rounding error? TRUE when its range is no more
# than 1e-10 times its largest absolute value (agreement to ten significant
# digits), which catches the noise left by arithmetic on equal values (for
# example gains computed from records near 100,000) without calling real
# variation constant. One threshold for the whole package.
#' @noRd
.smsd_constant_tol <- 1e-10

#' @noRd
.is_constant <- function(x) {
  r <- range(x)
  (r[2L] - r[1L]) <= .smsd_constant_tol * max(abs(r))
}

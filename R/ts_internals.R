# Internal helpers shared by the two stage (ts_*) and fixed width procedures.
# Not exported. The first stage ("pilot") data fix the criterion sample size;
# the second stage ("additional") data top the sample up. The check functions
# accept both and combine them. Argument checking is delegated to the shared
# checks in R/checks.R.

#' @noRd
.ts_combine <- function(pilot_data, additional_data, na.rm, min_n = 2L, ...) {
  if (missing(pilot_data)) {
    stop("'pilot_data' must be supplied.", call. = FALSE)
  }
  na.rm <- .check_flag(na.rm, "na.rm")
  pilot_data <- .check_data(pilot_data, "pilot_data", min_n = min_n,
                            na.rm = na.rm, ...)
  if (is.null(additional_data)) {
    additional_data <- numeric(0)
  } else {
    additional_data <- .check_data(additional_data, "additional_data",
                                   min_n = 0L, na.rm = na.rm, ...)
  }
  list(pilot = pilot_data,
       all = c(pilot_data, additional_data),
       m = length(pilot_data),
       n = length(pilot_data) + length(additional_data))
}

# `what` describes the role of d in the calling procedure: the half width for
# the symmetric intervals, the full width for the one-sided negative
# exponential intervals.
#' @noRd
.ts_check_width <- function(d, what = "the confidence interval half width") {
  if (missing(d)) {
    stop("'d' must be supplied (", what, ").", call. = FALSE)
  }
  ok <- is.numeric(d) && is.null(dim(d)) && length(d) == 1L
  if (ok) .stop_if_classed(d, "d")
  ok <- ok && is.finite(d) && d > 0
  if (!ok) {
    stop("'d' must be a single finite number greater than 0 (", what,
         "); got ", .smsd_describe(d), ".", call. = FALSE)
  }
  invisible(TRUE)
}

#' @noRd
.ts_alpha <- function(conf_level) {
  1 - .check_conf_level(conf_level)
}

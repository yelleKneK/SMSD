# The lower bound shared by the .var_* asymptotic variance estimators, and
# what a check function does when an estimate cannot be used.
#
# The source articles bound each variance estimate below by n^(-3). The bound
# keeps the estimate away from zero in the early, small-sample stages, and
# the articles treat a non-positive estimate as a negligible event that the
# bound absorbs.
#
# It is not negligible at small pilots: several of the estimators are
# differences of moment terms and are negative in a substantial share of
# small samples (about one pilot in six for the regression slope and the
# Pearson correlation at a target width of 0.4). Replacing a negative
# estimate by n^(-3) makes the criterion sample size tiny, so the stopping
# rule would be met at once with an interval far narrower than the data
# support. The package therefore separates two cases.
#
# * The raw estimate is positive. It is used as the articles prescribe:
#   max(raw, n^(-3)). A small positive estimate may be the truth (data on a
#   small scale, a small coefficient of variation), and the bound is then the
#   conservative value the articles intend.
# * The raw estimate is zero, negative, or not finite. No variance has been
#   estimated. The check reports the stopping rule as not met, sets
#   attr(., "variance_unusable") to TRUE, and attaches
#   .note_variance_unusable(), so the user keeps sampling. This is a
#   deliberate, documented departure from the stopping rules as printed.
#
# Every .var_* helper takes `floor = TRUE`; with `floor = FALSE` it returns
# the raw estimate, which is what a check function needs in order to tell the
# two cases apart:
#
#     raw <- .var_x(..., floor = FALSE)
#     unusable <- .var_unusable(raw)
#     xi2 <- .var_floor(raw, n)

#' @noRd
.var_floor <- function(value, n, floor = TRUE) {
  if (!floor) {
    return(value)
  }
  floor_value <- n^(-3)
  if (!is.finite(value)) {
    return(floor_value)
  }
  max(value, floor_value)
}

# Is a raw variance estimate (before the bound) unusable: zero, negative, or
# not finite? Vectorized, so a procedure with several groups can test all of
# its variances in one call.
#' @noRd
.var_unusable <- function(raw) {
  !is.finite(raw) | raw <= 0
}

# The note attached to a check whose variance estimate is unusable.
#' @noRd
.note_variance_unusable <- function() {
  paste0("The variance estimate is not yet usable at this sample size (it ",
         "is zero, negative, or not finite), so the stopping rule is ",
         "reported as not met. Collect more data and check again.")
}

# The note attached to a check whose point estimate is not finite.
#' @noRd
.note_not_estimable <- function() {
  paste0("The estimate is not finite at the current data (for example, a ",
         "sample variance or mean of zero), so the stopping rule is reported ",
         "as not met. Collect more data and check again.")
}

# Unbiased estimator of the third central moment (the h-statistic of order
# three): n sum((x - mean(x))^3) / ((n - 1)(n - 2)). Kelley, Darku, and
# Chattopadhyay (2018, Equations 75 and 78). Internal; not exported.
#
# x: numeric vector with at least 3 (non missing) values.
# na.rm: remove missing values first? Default FALSE.
#
# The estimator is exactly unbiased for any distribution with a finite third
# moment and is location invariant (the data enter only through deviations
# from the sample mean).
#
# Authors: Ken Kelley, Francis Bilson Darku, Bhargab Chattopadhyay.
#' @noRd
.mu3 <- function(x, na.rm = FALSE) {
  if (!is.numeric(x)) {
    .smsd_stop("'x' must be a numeric vector.")
  }
  .check_flag(na.rm, "na.rm")
  if (na.rm) {
    x <- x[!is.na(x)]
  }
  if (length(x) < 3) {
    .smsd_stop("'x' must contain at least 3 observations.")
  }
  n <- as.numeric(length(x))
  d <- .aipe_center(x)
  n * sum(d^3) / ((n - 1) * (n - 2))
}

# Estimator of the joint central moment mu_31 = E[(X - mu_X)^3 (Y - mu_Y)],
# the mirror image of .mu13(): mu31_hat = k31 + 3 k20 k11, Kelley, Darku,
# and Chattopadhyay (2019, Equations 73 and 74). Internal; not exported.
#
# x, y: paired numeric vectors with at least 4 pairs.
# na.rm: drop a pair when either member is missing? Default FALSE.
#
# Authors: Ken Kelley, Francis Bilson Darku, Bhargab Chattopadhyay.
#' @noRd
.mu31 <- function(x, y, na.rm = FALSE) {
  .mu13(y, x, na.rm = na.rm)
}

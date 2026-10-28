# Internal helpers shared by the minimum risk point estimation (seq_mr_*)
# procedures. Not exported. Argument checking is delegated to the shared
# checks in R/checks.R.

# The loss constant A: given directly, or derived as price / epsilon^2. A
# derived value that is not a positive finite number (price / epsilon^2
# overflows for a tiny epsilon and underflows to zero for a huge one) is
# reported in terms of the two arguments the user supplied, not as a fault
# in 'A', which the user did not pass.
#' @noRd
.mr_loss_constant <- function(A, price, epsilon) {
  if (!is.null(A)) {
    if (!is.null(price) || !is.null(epsilon)) {
      stop("Supply either 'A' or ('price', 'epsilon'), not both.",
           call. = FALSE)
    }
    return(.check_positive(A, "A"))
  }
  if (is.null(price) || is.null(epsilon)) {
    stop("Supply either 'A', or both 'price' and 'epsilon'.", call. = FALSE)
  }
  price <- .check_positive(price, "price")
  epsilon <- .check_positive(epsilon, "epsilon")
  A <- price / epsilon^2
  if (!is.finite(A) || A <= 0) {
    .smsd_stop("The loss constant implied by 'price' and 'epsilon', ",
               "price / epsilon^2, must be a positive finite number; it is ",
               format(A), " for price = ", format(price), " and epsilon = ",
               format(epsilon), ". Check both arguments and their units.")
  }
  A
}

# Returns the cleaned data vector (doubles, missing values removed when
# na.rm is TRUE; an error when they are present and na.rm is FALSE).
#' @noRd
.mr_check_data <- function(data, na.rm, min_n = 2L, arg = "x", ...) {
  if (missing(data)) {
    stop("'", arg, "' must be supplied.", call. = FALSE)
  }
  na.rm <- .check_flag(na.rm, "na.rm")
  .check_data(data, arg, min_n = min_n, na.rm = na.rm, ...)
}

#' @noRd
.mr_check_cost <- function(c) {
  .check_positive(c, "c")
  invisible(TRUE)
}

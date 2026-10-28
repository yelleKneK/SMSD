# Internal helper shared by the bounded risk (seq_br_*, ts_br_*) procedures.
# Not exported. Argument checking is delegated to the shared checks in
# R/checks.R.

#' @noRd
.br_check_Akw <- function(A, k, w) {
  if (missing(A)) stop("'A' must be supplied.", call. = FALSE)
  if (missing(k)) stop("'k' must be supplied.", call. = FALSE)
  if (missing(w)) stop("'w' must be supplied (the risk bound).", call. = FALSE)
  .check_positive(A, "A")
  if (is.numeric(k) && is.null(dim(k)) && length(k) == 1L) {
    .stop_if_classed(k, "k", "1 or 2 as a plain number")
  }
  if (!is.numeric(k) || !is.null(dim(k)) || length(k) != 1L || is.na(k) ||
      !(k %in% c(1, 2))) {
    stop("'k' must be 1 (absolute error) or 2 (squared error); got ",
         .smsd_describe(k), ".", call. = FALSE)
  }
  .check_positive(w, "w")
  invisible(TRUE)
}

# The risk constant of the bounded risk procedures for a normal mean. With
# loss A |xbar - mu|^k and xbar - mu distributed as (sigma / sqrt(n)) Z, the
# fixed sample size risk is B (sigma^2 / n)^(k / 2) with B = A E|Z|^k. In
# general E|Z|^k = 2^(k / 2) Gamma((k + 1) / 2) / sqrt(pi); for the two powers
# the procedures accept that is sqrt(2 / pi) (k = 1) and 1 (k = 2). The closed
# forms are used because the general expression evaluates to 1 + 2.2e-16 at
# k = 2, which made B a shade above A and the estimated risk a shade above
# the bound at the very sample size the stopping rule calls for.
#' @noRd
.br_constant <- function(A, k) {
  if (k == 2) A else A * sqrt(2 / pi)
}

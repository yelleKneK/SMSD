#' Construct an smsd_sensitivity object
#'
#' Internal constructor for the tidy summary returned by every
#' \code{*_sensitivity()} function: a data frame with one row per summary
#' quantity. It is an ordinary \code{term}/\code{value} data frame with one
#' more class, so that it prints readably.
#'
#' @param term A character vector of row names.
#' @param value A numeric vector of values, the same length as \code{term}.
#'
#' @return An object of class \code{c("smsd_sensitivity", "data.frame")}.
#'
#' @keywords internal
#' @noRd
new_smsd_sensitivity <- function(term, value) {
  stopifnot(is.character(term), is.numeric(value),
            length(term) == length(value))
  out <- data.frame(term = term, value = as.numeric(value),
                    stringsAsFactors = FALSE)
  class(out) <- c("smsd_sensitivity", "data.frame")
  out
}

#' Print a Sensitivity Summary
#'
#' The \code{*_sensitivity()} functions, \code{\link{seq_aipe_smd_sensitivity}}
#' among them, return a summary that these methods print. The object is a
#' data frame with columns \code{term} and \code{value} and can be used as
#' one; these methods only change how it is displayed. Each value is
#' formatted on its own, in fixed notation, so that a table mixing a
#' replication count, a coverage proportion, and a mean sample size does not
#' fall into scientific notation as the default data frame print would have
#' it.
#'
#' @param x An object of class \code{"smsd_sensitivity"}.
#' @param digits Number of significant digits to show. Defaults to
#'   \code{getOption("digits")}. Stored values are never rounded; this
#'   affects display only.
#' @param ... Not used.
#'
#' @return \code{print()} returns its input invisibly. \code{format()}
#'   returns a character vector, one element per printed line.
#'
#' @author Ken Kelley \email{kkelley@@nd.edu}
#'
#' @examples
#' # G is kept small so the example runs quickly; a study to be reported
#' # deserves G = 10000 or more.
#' out <- seq_aipe_smd_sensitivity(true_value = 0.5, width = 0.8, G = 20,
#'                                 seed = 113)
#' out
#' out$value[out$term == "coverage"]
#'
#' @name smsd_sensitivity_methods
#' @aliases smsd_sensitivity
NULL

#' @rdname smsd_sensitivity_methods
#' @export
format.smsd_sensitivity <- function(x, digits = getOption("digits"), ...) {
  value <- vapply(x$value, function(v) {
    if (is.na(v)) "NA" else format(v, digits = digits, scientific = FALSE)
  }, character(1))
  width <- max(nchar(x$term), 4L)
  c(sprintf("%-*s  %s", width, "term", "value"),
    sprintf("%-*s  %s", width, x$term, value))
}

#' @rdname smsd_sensitivity_methods
#' @export
print.smsd_sensitivity <- function(x, digits = getOption("digits"), ...) {
  cat(format(x, digits = digits), sep = "\n")
  invisible(x)
}

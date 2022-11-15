#' Calculates exact value of C
#'
#' This function replicates the function in Mukhopadhyay \& de Silva (chapter 6, section 2) for the exact sample size given the variance, significance level, and confidence interval width for a fixed width confidence interval for the population mean of normal random variables.
#'
#' @param var The population variance.
#' @param alpha The chosen Type I error rate (complement of the confidence coefficient).
#' @param d Half of the confidence interval width. That is the confidence
#' interval is 2d wide.
#'
#' @return The optimal sample size $C$, the ``optimal fixed sample size required."
#'
#' @references
#' Mukhopadhyay, N., \& de Silva, B. M. (2009). \emph{Sequential Methods and Their Applications}. New York: CRC Press.
#'
#' @export calculate_c
#'
#'
#' @author Ken Kelley \email{kkelley@nd.edu},
#' Francis Bilson Darku \email{fbilsond@nd.edu},
#' Bhargab Chattopadhyay \email{bhargab@iimv.ac.in}
#'
#' @examples
#' calculate_c(var = 1, alpha = 0.05, d = 0.75)
#'
#'
calculate_c <- function(var, alpha, d){
  # Calculate a
  a <- stats::qnorm(1-(alpha/2))
  # Calculate c
  c <- ((a^2)*var)/(d^2)
  return(c)
}

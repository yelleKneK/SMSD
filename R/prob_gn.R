#' ProbGn
#'
#' This function replicates the ProbGn function from chapter 2.
#' This program computes the probability of K(n) > c*sigma for given alpha
#' and sample size.
#'
#' @param alpha The selected significance level.
#' @param n The selected sample size.
#' @param c_min The minimum value of c to try in the calulcations, default = 0.
#' @param c_max The maximum value of c to try in the calculations,
#' default = 10.
#' @param by The increment of the sequence between c_min and c_max,
#' default = 0.04.
#' @param dec The number of decimal places to use when removing probabilities
#' of 0 and 1 from the ends of the results.
#'
#' @return A sequence of values of c and the corresponding calculated
#'  probabilities.
#'
#' @references
#' Mukhopadhyay, N., \& de Silva, B. M. (2009). \emph{Sequential Methods and Their Applications}. New York: CRC Press.
#' 
#' @author Bhargab Chattopadhyay \email{bhargab@iimv.ac.in}, Ken Kelley \email{kkelley@nd.edu}
#'
#'
#' @examples
#' alpha = 0.05
#' n = 20
#' vals <- prob_gn(alpha = alpha, n = n)
#' vals
#'
#' @export prob_gn
#'
prob_gn <- function(alpha, n, c_min = 0, c_max = 10, by = 0.04,
                    dec = 4) {
  c_seq = seq(from = c_min, to=c_max, by=by)
  p_vals <- rep(NA, length(c_seq))
  for(i in 1:length(c_seq)){
    p_vals[i] <- chi_prob(alpha, n, c_seq[i])
  }
  res_dat <- cbind.data.frame(c_seq, p_vals)
  res_mod <- res_dat[which(round(p_vals, dec) < 1 &
                             round(p_vals, dec) > 0 ),]
  return(res_mod)
}
#' chi_prob
#'
#' This function is used internally by several of the chapter 2 functions.
#'
#' @param alpha The significance level.
#' @param n The selected sample size.
#' @param c The value of c to use for this calculation.
#'
#' @return The calculated probabilities for the given parameters.
#' @examples
#' chi_prob(0.05, 5, 1)
#'
#' @export chi_prob
#'
chi_prob <- function(alpha, n, c){
  t1 = stats::qt(.5 * alpha, df=2*(n-2))
  t2 = (n * (n-1) * (c^(2)))/(4 * (t1^(2)))
  prob <- stats::pchisq(t2, df = 2 * (n-2), lower.tail = FALSE)
  return(prob)
}

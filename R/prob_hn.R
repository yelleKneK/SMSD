#' ProbHn
#'
#' prob_hn replicates the ProbHn function from chapter 2 of Mukhopadhyay and de Silva (2009).
#' 
#' This function computes the probability of K(n) > c*sigma for given alpha
#' and sample size.
#'
#' @param alpha The selected significance level.
#' @param n The sample size.
#' @param c_min The minimum value of c to try in the calulcations, default = 0.
#' @param c_max The maximum value of c to try in the calculations,
#' default = 10.
#' @param by The increment of the sequence between c_min and c_max,
#' default = 0.04.
#' @param dec The number of decimal places to use when removing probabilities
#' of 0 and 1 from the ends of the results.
#'
#' @return A sequence of values of c and the corresponding probabilities.
#'
#' @references
#' Mukhopadhyay, N., \& de Silva, B. M. (2009). \emph{Sequential Methods and Their Applications}. New York: CRC Press.
#'
#'@author Bhargab Chattopadhyay \email{Bhargab@iiitvadodara.ac.in}, Ken Kelley \email{kkelley@nd.edu}
#'
#' @examples
#' alpha = 0.05
#' n = 20
#' vals <- prob_hn(alpha = alpha, n = n)
#' vals
#'
#' @export prob_hn
#'
prob_hn <- function(alpha, n, c_min = 0, c_max = 10, by = 0.04,
                    dec = 4) {
  c_seq = seq(from = c_min, to=c_max, by=by)
  p_vals <- rep(NA, length(c_seq))
  for(i in 1:length(c_seq)){
    p_vals[i] <- chi_prob_2(alpha, n, c_seq[i])
  }
  res_dat <- cbind.data.frame(c_seq, p_vals)
  res_mod <- res_dat[which(round(p_vals, dec) < 1 &
                             round(p_vals, dec) > 0 ),]
  return(res_mod)
}
#' chi_prob_2
#'
#' This function is called internally by \code{ProbHn_2}.
#'
#' @param alpha The selected significance level.
#' @param n The sample size.
#' @param c The c value.
#'
#'
#' @return The calculated probability.
#'
#' @examples
#' chi_prob_2(0.05, 5, 1)
#'
#' @export chi_prob_2
chi_prob_2 <- function(alpha, n, c){
  t1 = stats::qf(alpha, df1=2, df2=((2 * n)-2), lower.tail = FALSE)
  t2 = (n * (n-1) * (2 *c))/t1
  prob <- stats::pchisq(t2, df = 2*(n-2), lower.tail = FALSE)
  return(prob)
}

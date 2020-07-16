#' power for two independent groups when population standard deviation known
#'
#' This function calculates the probability of Type II error when the sample size is known, and the required sample size to achieve a preassigned Type II error probability (beta) for a test of two independent groups. This function uses the $z$-distribution (note the more commonly used $t$-distribution) is used.
#
#' @param mu0 Population mean of group 0 (e.g., control).
#' @param mu1 Population mean of group 1 (e.g., experimental).
#' @param sigma The population standard deviation, assumed equal across groups.
#' @param alpha The Type I error rate.
#' @param n The number of samples, used to calculate Type II error for a given
#' sample size.
#' @param n1 The lowest number of samples to consider when generating a
#' sequence of Type II errors.
#' @param n2 The highest number of samples to consider when generating a
#' sequence of Type II errors.
#' @param n_vals The number of samples to generate Type II error for.
#' @param beta The desired Type II error rate, used when calculating the
#' necessary sample size.
#' @param method What should the function calculate, (1) "sample_size" returns
#' the calculated sample size necessary to achieve a specified Type II error
#' rate. (2) "type_2" returns the calculated Type II error rate for a given
#' sample size. (3) "type_2_seq" returns a series of length "n_vals" of
#' calculated type II error rates from n1 to n2.
#'
#' @return Depending on the method chosen the function returns either the
#' necessary sample size to achieve a certain type II error, the type II error
#' for a given sample size or a sequence of type II errors for different
#' sample sizes.
#'
#' @references
#' Mukhopadhyay, N., \& de Silva, B. M. (2009). \emph{Sequential Methods and Their Applications}. New York: CRC Press.
#'
#'
#' @examples
#' # Set parameter values
#' mu0 <- 0
#' mu1 <- 1
#' sigma <- 2
#' alpha <- 0.05
#' n <- 16
#'
#' # Calculate type II error
#' t2 <- power(mu0 = mu0, mu1 = mu1, sigma = sigma,
#'             alpha = alpha , n = n, method="type_2")
#' t2
#'
#' # Calculate sample size
#' beta <- 0.361
#' size <- power(mu0 = mu0, mu1 = mu1, sigma = sigma,
#'               alpha = alpha, beta = beta, method = "sample_size")
#' size
#'
#' # Return sequence of type II errors
#' n1 <- 5
#' n2 <- 50
#' n_vals <- 10
#' seq_vals <- power(mu0 = mu0, mu1 = mu1, sigma=sigma, alpha = alpha,
#'                   n1 = n1, n2 = n2, n_vals = n_vals,
#'                   method = "type_2_seq")
#' seq_vals
#'
#' @export power
#'
power <- function(mu0, mu1, sigma, alpha, n=NULL,
                  n1 = NULL, n2=NULL, n_vals=NULL, beta=NULL,
                  method = c("sample_size", "type_2", "type_2_seq")){
  # Warning Messages
  if(method == "sample_size" & is.null(beta)){
    print("Please provide beta value to calculate sample size")
  } else if (method == "type_2" & is.null(n)){
    print("Please provide n1 value to calculate type II error")
  } else if(method == "type_2_seq" & (is.null(n1) | is.null(n2) |
                                      is.null(n_vals))){
    print("Please provide n1, n2 (Start and end points, and n_vals
          (Number of values) to calculate the sequence of type II error")
  }
  # Sample size calculation
  if(method == "sample_size"){
    size <- sample_size(mu0, mu1, sigma, alpha, beta)
    return(size)
  }
  if(method == "type_2"){
    t_2 <- type_two_calc(mu0, mu1, sigma, alpha, n)
    return(t_2)
  }
  if(method == "type_2_seq"){
    ps <- power_seq(mu0, mu1, sigma, alpha, n1, n2, n_vals)
    return(ps)
  }
}

#' power_seq
#'
#' This function calculates the power for a sequence of sample sizes using \code{power} for two independent groups when the population standard deviation is assumed known and common across groups.
#' 
#' @param mu0 Population mean of group 0 (e.g., control).
#' @param mu1 Population mean of group 1 (e.g., experimental).
#' @param sigma The population standard deviation, assumed equal across groups.
#' @param alpha The Type I error rate.
#' @param n1 The lowest number of samples to consider when generating a
#' sequence of Type II errors.
#' @param n2 The highest number of samples to consider when generating a
#' sequence of Type II errors.
#' @param n_vals The number of samples to generate Type II error for.
#'
#' @return Returns a series of length "n_vals" of
#' calculated Type II error rates from n1 to n2.
#'
#' @examples
#' # Set parameter values
#' mu0 <- 0
#' mu1 <- 1
#' sigma <- 2
#' alpha <- 0.05
#' n1 <- 5
#' n2 <- 50
#' n_vals <- 10
#' # Calculate sequence
#' seq_vals <- power_seq(mu0 = mu0, mu1 = mu1, sigma=sigma, alpha = alpha,
#'                   n1 = n1, n2 = n2, n_vals = n_vals)
#'
#' @export power_seq
power_seq <- function(mu0, mu1, sigma, alpha, n1, n2, n_vals){
  n_seq <- seq(from=n1, to=n2, length.out = n_vals )
  t_2_err <- rep(NA, length(n_seq))
  for(i in 1:length(t_2_err)){
    t_2_err[i] <- type_two_calc(mu0, mu1, sigma, alpha, n_seq[i])
  }
  return(cbind.data.frame(n_seq, t_2_err))
}

#' type_two_calc
#'
#' This function is called internally by \code{power} and is used to
#' calculate the probability of Type II error.
#'
#' @param mu0 Mean value 0.
#' @param mu1 Mean value 1.
#' @param sigma The standard deviation.
#' @param alpha The Type I error rate.
#' @param n The number of samples, used to calculate Type II error for a given
#' sample size.
#'
#' @return The Type II error for a given samples size.
#'
#' @export type_two_calc
#'
type_two_calc <- function(mu0, mu1, sigma, alpha, n){
  z <- stats::qnorm(1-alpha)
  t2 <- (sqrt(n) * (mu1 - mu0))/sigma
  type_2_err <- stats::pnorm(z - t2, 0, 1)
  return(type_2_err)
}
#' sample_size
#'
#' This function calculates the probability of Type II error for two independent groups when the 
#' when the population standard deviation is assumed known (i.e., using the $z$-distribution)#'  in order to achieve a preassigned Type II
#' error probability (beta) of interest.
#'
#' @param mu0 Mean value 0.
#' @param mu1 Mean value 1.
#' @param sigma The standard deviation.
#' @param alpha The Type I error rate.
#' @param beta The desired Type II error rate, used when calculating the
#' necessary sample size.
#'
#' @return Returns the Type II error for a given sample size.
#'
#' @export sample_size
#'
sample_size <- function(mu0, mu1, sigma, alpha, beta){
  z_a <- stats::qnorm(1-alpha)
  z_b <- stats::qnorm(1-beta)
  size <- (((z_a + z_b)/(mu1 - mu0))^(2)) * (sigma^(2))
  return(size)
}

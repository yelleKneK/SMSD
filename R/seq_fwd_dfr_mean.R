#' seq_fwd_dfr_mean
#'
#' Purely Sequential Approach to Get Fixed Width Confidence
#' Interval for the mean of distribution free random variables. Calculates a confidence
#' interval for given data and parameters.
#'
#' @param data The data for which to calculate the confidence interval.
#' A numeric vector.
#' @param d Half of the confidence interval width, must be a non-zero positive
#'  value.
#' @param alpha The significance level. A value between 0 and 1.
#' @param gamma gamma
#' @param pilot Should a pilot sample be generated. True/False value.
#' default value is \code{FALSE}.
#' @param verbose Should the criterion be printed. Default is \code{FALSE}.
#' @param na.rm This parameter controls whether NA values are removed from
#' the data prior to calculation. Default is \code{TRUE}.
#'
#' @return The calculated confidence interval, the sample size, data mean,
#' and an indicator of if the criterion was satisfied.
#'
#' @author Bhargab Chattopadhyay \email{Bhargab@iiitvadodara.ac.in},
#' Neetu Shah \email{201451015@iiitvadodara.ac.in}, Ken Kelley \email{kkelley@nd.edu}
#'
#'
#' @references
#' Mukhopadhyay, N., \& de Silva, Basil M. (2009). \emph{Sequential Methods and Their Applications}. New York: CRC Press.
#'
#'
#' @export seq_fwd_dfr_mean
#'
#' @examples
#' pilot_ss <- seq_fwd_dfr_mean(alpha=0.1, d=0.2, gamma=1, pilot=TRUE)
#' SLS <- rnorm(pilot_ss, mean=10, sd=16)
#' seq_fwd_dfr_mean(data=SLS, d=0.2, alpha=0.1, pilot=FALSE)
seq_fwd_dfr_mean <- function(data, d, alpha, gamma, pilot=FALSE,
                             verbose=FALSE,
                             na.rm=TRUE)
{
  if(missing(d)){
    stop("You must specify \'d\'")
  }

  if(missing(alpha)){
    stop("You must specify \'alpha\'")
  }

  if(!missing(d)){
    if(d<0) stop("d should be a non-zero positive value")
  }

  if(!missing(alpha)){
    if(alpha>1 & alpha<0) stop("alpha should be between 0 and 1")
  }

  if(missing(gamma)){
    gamma <- 1
  }

  if(pilot==FALSE)
  {
    if (!is.data.frame(data) && !is.matrix(data) && !is.vector(data)){
      stop("The argument 'data' must be a data.frame
           or matrix with one column")
    }

    if (dim(data)[2] != 1 && !is.null(data) && !is.vector(data)){
      stop("The argument 'data' must have only one column,
           or be 'NULL' for pilot = TRUE")
    }


    if(is.data.frame(data))
    {
      data <- as.vector(data)
    }
    if(na.rm){
      data <-  data[!is.na( data)]
    }
    n <- length(data)

    optimal_N <- as.integer((stats::qnorm(1-
                                            alpha/2)^2*(stats::var(data)+
                                                          1/n))/d^2)

    ci <- mean(data)+c(-1,1)*d # This is the fixed width confidence interval.

    if(verbose==FALSE){
      Outcome <- list("Confidence Interval"=ci[1],ci[2],
                            N=n, mean=mean(data),
                            "Is.Satisfied?"=(n >= optimal_N))
    } else if(verbose==TRUE){
      Outcome <- list("Confidence Interval"=ci[1],ci[2],
                            N=n, mean=mean(data), Criterion=optimal_N,
                            Is.Satisfied=n >= optimal_N)
    }
  }

  if(pilot==TRUE)
  {
    Outcome <- c(Pilot.SS=max(2,
                              ceiling((stats::qnorm(1-alpha)/d)^(2/(1+
                                       gamma)))+1))
  }

  return(Outcome)
}

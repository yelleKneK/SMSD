#' ts_fwd_nexp_mean
#'
#' @description Two Stage approach to Get Fixed Width Confidence
#' Interval  for the mean of negative exponential random variables.
#'
#' @param data The data for which to calculate the confidence interval.
#' A numeric vector.
#' @param d Half of the confidence interval width, must be a non-zero positive
#' value.
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
#' @references
#' Mukhopadhyay, N., \& de Silva, Basil M. (2009). \emph{Sequential Methods and Their Applications}. New York: CRC Press.
#'
#' @export ts_fwd_nexp_mean
#'
#' @examples
#' pilot_ss <- ts_fwd_nexp_mean(alpha=0.01, d=0.5, gamma=1, pilot=TRUE)
#' \dontrun{
#' SLS <- rexp(pilot_ss, rate=1)
#' }
#' SLS <- rexp(100, rate=1)
#' ts_fwd_nexp_mean(data=SLS, d=0.5, alpha=0.01, pilot=FALSE)
#'
ts_fwd_nexp_mean <- function(data, d, alpha, gamma, pilot=FALSE,
                             verbose=FALSE, na.rm=TRUE)
{
  if(missing(d)) {
    stop("You must specify \'d\'")
  }

  if(missing(alpha)){
    stop("You must specify \'alpha\'")
  }

  if(!missing(d)){
    if(d<0){
      stop("d should be a non-zero positive value")
    }
  }

  if(!missing(alpha)){
    if(alpha>1 & alpha<0){
      stop("alpha should be between 0 and 1")
    }
  }

  if(missing(gamma)){
    gamma <- 1

  }
  if(pilot==FALSE){
    if (!is.data.frame(data) && !is.matrix(data) && !is.vector(data)){
      stop("The argument 'data' must be a data.frame
           or matrix with one column")
    }

    if (dim(data)[2] != 1 && !is.null(data) && !is.vector(data)){
      stop("The argument 'data' must have only one column,
           or be 'NULL' for pilot = TRUE")
    }


    if(is.data.frame(data)){
      data <- as.vector(data)
    }
    if(na.rm){
      data <-  data[!is.na( data)]
    }
    n <- length(data)
    U <- sum(data-min(data))/(n-1)
    optimal_N <- max(n, ceiling((stats::qf(1-alpha, 2,2*n-2)*U)/d))

    ci <- c(min(data)-d, min(data)) # This is the fixed width confidence interval.

    if(verbose==FALSE){
      Outcome <- rbind(list("Confidence Interval"=ci[1],
                            ci[2],
                            N=n, mean=mean(data)
                            ,"Is.Satisfied?"=(n >= optimal_N),
                            Stage=1))
    }
    if(verbose==TRUE){
      Outcome <- rbind(list("Confidence Interval"=ci[1],
                            ci[2],
                            N=n, mean=mean(data),
                            Criterion=optimal_N,
                            Is.Satisfied=n >= optimal_N,
                            Stage=1))
    }
    if((n >= optimal_N)==FALSE){
      print(Outcome)
    }
    while((length(data) >= optimal_N)==FALSE){
      cat(optimal_N-length(data)," more are needed ")
      obs <- as.integer(strsplit(readline(), " ")[[1]])
      data <- c(data,obs)
      if(length(data)==optimal_N){
        ci <- c(min(data)-d, min(data))
        Outcome <- rbind(list("Confidence Interval"=ci[1],
                              ci[2], N=length(data), mean=mean(data), Stage=2))
      }else if(length(data)> optimal_N){
        ci <- mean(data)+
          c(-1,1)*qt(1-alpha/2,
                     length(data)-1)*sqrt(stats::var(data)/length(data))
        Outcome <- rbind(list("Confidence Interval"=ci[1], ci[2],
                              N=length(data), mean=mean(data), Stage=2))
      }
    }
  }

  if(pilot==TRUE)
  {
    Outcome <- c(pilot_ss=max(2,ceiling((stats::qnorm(1-alpha)/d)^(2/(1+
                                                                  gamma)))))
  }

  return(Outcome)
}

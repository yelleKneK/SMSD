#' seq_blci_dfr_mean
#'
#' Purely Sequential approach to Get Bounded Length Confidence Interval for the
#' median of distribution free random variables
#'
#' @param data The data for which to calculate the confidence interval.
#' A numeric vector.
#' @param d Half of the confidence interval width, must be a non-zero positive
#'  value.
#' @param alpha The significance level. A value between 0 and 1.
#' @param pilot Should a pilot sample be generated. True/False value.
#' default value is \code{FALSE}.
#' @param verbose Should the criterion be printed. Default is \code{FALSE}.
#' @param na.rm This parameter controls whether NA values are removed from
#' the data prior to calculation. Default is \code{TRUE}.
#'
#' @return The calculated confidence interval, the sample size, data mean,
#' and an indicator of if the confidence level was satisfied.
#'
#' @author Bhargab Chattopadhyay \email{Bhargab@iiitvadodara.ac.in},
#' Neetu Shah \email{201451015@iiitvadodara.ac.in}
#'
#' @references
#' Mukhopadhyay, N., \& de Silva, B. M. (2009). \emph{Sequential Methods and Their Applications}. New York: CRC Press.
#'
#' @examples
#' pilot_ss <- seq_blci_dfr_mean(alpha=0.05, d=1, pilot=TRUE)
#' SLS <- rnorm(pilot_ss, mean=0, sd=1)
#' seq_blci_dfr_mean(data=SLS, d=1, alpha=0.5, pilot=FALSE)
#'
#' @export seq_blci_dfr_mean
#'
seq_blci_dfr_mean <- function(data, d, alpha, pilot=FALSE,
                              verbose=FALSE, na.rm=TRUE){
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

    b <- max(1, ceiling((n-stats::qnorm(1-alpha/2)*sqrt(n)-1)/2))

    a <- n-b+1

    data <- sort(data, decreasing = FALSE)

    ci <- c(data[b], data[a])

    if(verbose==FALSE){
      Outcome <- rbind(list("Confidence Interval"=ci[1],
                            ci[2], N=n, mean=mean(data),
                            "Is.Satisfied?"=((data[a]-data[b])<=2*d)))
    }
    if(verbose==TRUE){
      Outcome <- rbind(list("Confidence Interval"=ci[1],
                            ci[2], N=n, mean=mean(data),
                            2*d,
                            Is.Satisfied=(data[a]-data[b])<=2*d))
    }
  }

  if(pilot==TRUE)
  {
    Outcome <- c(Pilot.SS=max(2,ceiling((stats::qnorm(1-alpha)/d)^2)+1))
  }

  return(Outcome)
}

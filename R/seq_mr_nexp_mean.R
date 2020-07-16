#' seq_mr_nexp_mean
#'
#' @description Purely Sequential approach to Get Minimum Risk Point Estimation
#'  for negative exponential random variables
#'
#' @param data The data for which to calculate the minimum risk point.
#' @param A The loss function constant.
#' @param c The cost of unit sample.
#' @param k Loss function index 1.
#' @param t Loss function index 2.
#' @param pilot Should a pilot sample be generated. True/False value.
#' default value is \code{FALSE}.
#' @param verbose Should the criterion be printed. Default is \code{FALSE}.
#' @param na.rm This parameter controls whether NA values are removed from
#' the data prior to calculation. Default is \code{TRUE}
#'
#' @return The calculated risk, the sample size, mean, standard deviation,
#'  and
#' an indicator of if the criterion was satisfied.
#'
#' @author Bhargab Chattopadhyay \email{Bhargab@iiitvadodara.ac.in},
#' Neetu Shah \email{201451015@iiitvadodara.ac.in}, Ken Kelley \email{kkelley@nd.edu}
#'
#' @references
#' Mukhopadhyay, N., \& de Silva, B. M. (2009). \emph{Sequential Methods and Their Applications}. New York: CRC Press.
#'
#' @export seq_mr_nexp_mean
#'
#' @examples
#' pilot_ss <- seq_mr_nexp_mean(A=500, c=0.01, k=1, t=2, pilot=TRUE)
#' SLS <- rexp(pilot_ss, rate=1)
#' seq_mr_nexp_mean(data=SLS, A=500, c=0.01, k=1, t=2, pilot=FALSE)
seq_mr_nexp_mean <- function(data, A, c, k, t, pilot=FALSE, verbose=FALSE,
                             na.rm=TRUE)
{

  if(!missing(A)){
    if(A<=0) {
      stop("A should be a non-zero positive value")
    }
  }

  if(!missing(c)) {
    if(c<=0) {
      stop("c should be a non-zero positive value")
    }
  }

  if(missing(A)){
    stop("You must specify \'A\'")
  }

  if(missing(c)){
    stop("You must specify \'c\'")
  }

  if(!missing(k)){
    if(k<=0) {
      stop("k should be a non-zero positive value")
    }
  }
  if(missing(k)){
    stop("You must specify \'k\'")
  }

  if(!missing(t)){
    if(t<=0){
      stop("t should be a non-zero positive value")
    }
  }
  if(missing(t)){
    stop("You must specify \'t\'")
  }

  B <- A*k^2*gamma(k)

  if(pilot==FALSE){
    if (!is.data.frame(data) && !is.matrix(data) && !is.vector(data)){
      stop("The argument 'data' must be a data.frame or matrix with one column")
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

    optimal_N <- as.integer(((B*U)/(c*t))^(1/(t+k)))

    Rk<-c*(1+t/k)*n^t # This is the risk function.

    if(verbose==FALSE){
      Outcome <- rbind(list(Risk=Rk,
                            N=n,
                            cv=mean(data),
                            "Is.Satisfied?"=(n >= optimal_N)))
    }
    if(verbose==TRUE){
      Outcome <- rbind(list(Risk=Rk,
                            N=n,
                            cv=mean(data),
                            Criterion=optimal_N,
                            Is.Satisfied=n >= optimal_N))
    }
  }

  if(pilot==TRUE){
    Outcome <- c(Pilot.SS=max(4, ceiling((B/(c*t))^(1/(t+k)))))
  }

  return(Outcome)
}

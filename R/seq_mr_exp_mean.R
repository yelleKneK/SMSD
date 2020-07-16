#'seq_mr_exp_mean
#'
#' @description Purely Sequential approach to Get Minimum Risk
#'  Point Estimation for exponential random variables
#'
#' @param data The data for which to calculate the minimum risk point.
#' @param A The loss function constant.
#' @param c The cost of unit sample.
#' @param pilot Should a pilot sample be generated. True/False value.
#' default value is \code{FALSE}.
#' @param verbose Should the criterion be printed. Default is \code{FALSE}.
#' @param na.rm This parameter controls whether NA values are removed from
#' the data prior to calculation. Default is \code{TRUE}.
#'
#' @return The calculated risk, the sample size, mean, standard deviation, and
#' an indicator of if the criterion was satisfied.
#'
#' @author Bhargab Chattopadhyay \email{Bhargab@iiitvadodara.ac.in},
#' Neetu Shah \email{201451015@iiitvadodara.ac.in}, Ken Kelley \email{kkelley@nd.edu}
#'
#' @references
#' Mukhopadhyay, N., \& de Silva, Basil M. (2009). \emph{Sequential Methods and Their Applications}. New York: CRC Press.
#'
#' @export seq_mr_exp_mean
#'
#' @examples
#' pilot_ss <- seq_mr_exp_mean(A=100, c=0.04, pilot=TRUE)
#' SLS <- rexp(pilot_ss, rate=0.5)
#' seq_mr_exp_mean(data=SLS, A=100, c=0.04, pilot=FALSE)
seq_mr_exp_mean <- function(data, A, c, pilot=FALSE, verbose=FALSE,
                            na.rm=TRUE){

  if(!missing(A)){
    if(A<=0) stop("A should be a non-zero positive value")
  }
  if(!missing(c)) {
    if(c<=0) stop("c should be a non-zero positive value")
  }

  if(missing(A)){
    stop("You must specify \'A\'")
  }

  if(missing(c)) {
    stop("You must specify \'c\'")
  }

  if(pilot==FALSE) {
    if (!is.data.frame(data) && !is.matrix(data) && !is.vector(data)){
      stop("The argument 'data' must be a data.frame or
           matrix with one column")
    }

    if (dim(data)[2] != 1 && !is.null(data) && !is.vector(data)){
      stop("The argument 'data' must have only one column, or be 'NULL'
           for pilot = TRUE")
    }


    if(is.data.frame(data)){
      data <- as.vector(data)
    }
    if(na.rm){
      data <-  data[!is.na( data)]
    }
    n <- length(data)

    optimal_N <- as.integer(sqrt(A/c)*mean(data))

    Rk<-(A*mean(data)^2)/n +c*n

    if(verbose==FALSE) {
      Outcome <- rbind(list(Risk=Rk, N=n, cv=mean(data), "Is.Satisfied?"=(n >= optimal_N)))
    }
    if(verbose==TRUE){
      Outcome <- rbind(list(Risk=Rk,
                            N=n, cv=mean(data),
                            Criterion=optimal_N,
                            Is.Satisfied=n >= optimal_N))
    }
  }

  if(pilot==TRUE){
    Outcome <- c(Pilot.SS=max(4, ceiling((A/c)^(1/2))))
  }

  return(Outcome)
}

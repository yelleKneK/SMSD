#' seq_ubdes_norm
#'
#' @description Unbalanced Design using Purely Sequential approach
#' for two different size sets of normal random variables
#'
#'
#' @param data The data matrix for which to calculate the
#' @param C A data frame.
#' @param K K
#' @param w The risk bound,A non-negative integer greater than zero.
#' @param alpha The significance level. A value between 0 and 1.
#' @param pilot Should a pilot sample be generated. True/False value.
#' default value is \code{FALSE}.
#' @param verbose Should the criterion be printed. Default is \code{FALSE}.
#' @param na.rm This parameter controls whether NA values are removed from
#' the data prior to calculation. Default is \code{TRUE}.
#'
#' @return The calculated mean and an indicator of if the criterion is
#' satisfied.
#'
#' @author Ken Kelley \email{KKelley@nd.edu},
#' Francis Bilson Darku \email{FBilsonD@nd.edu},
#' Bhargab Chattopadhyay \email{Bhargab@iiitvadodara.ac.in},
#' Neetu Shah \email{201451015@iiitvadodara.ac.in}
#'
#' @references
#' Ken Kelley, Francis Bilson Darku, and Bhargab Chattopadhyay. (2018).
#' Accuracy in Parameter Estimation for a General Class of
#' Effect Sizes: A Sequential Approach. \emph{Psychological Methods}, \emph{23}, 226--243.
#'
#' @export seq_ubdes_norm
#'
#' @examples
#' C <- matrix(c(2,3), nrow=2)
#' pilot_ss <- seq_ubdes_norm(alpha=0.05, w=4, C=C, pilot=TRUE)
#' SLS <- matrix( rnorm(pilot_ss[1],mean=0,sd=1), rnorm(pilot_ss[2],
#' mean = 0, sd = 1), nrow=2,ncol=2)
#' seq_ubdes_norm(data=SLS, C=C, 2, w=4, alpha=0.05, pilot=FALSE)
#'
seq_ubdes_norm <- function(data, C, K, w, alpha, pilot=FALSE, verbose=FALSE,
                           na.rm=TRUE){
  if(missing(C)){
    stop("You must specify \'C\'")
  }
  if (!is.data.frame(C) && !is.matrix(C) && !is.vector(C)){
    stop("The argument 'C' must be a data.frame or matrix with one column")
  }
  if (dim(C)[2] != 1 && !is.null(C) && !is.vector(C)){
    stop("The argument 'C' must have only one column")
  }
  if(missing(w)){
    stop("You must specify \'w\'")
  }
  if(!missing(w)){
    if(w<=0) stop("w should be a non-negative integer greater than zero")
  }

  if(missing(alpha)){
    stop("You must specify \'alpha\'")
  }

  if(!missing(alpha)){
    if(alpha>1 & alpha<0) stop("alpha should be between 0 and 1")
  }

  if(pilot == FALSE){
    if(missing(K)){
      stop("You must specify \'K\'")
    }
    if(!missing(K)){
      if(K<=0) stop("K should be a non-negative integer greater than zero")
    }

    if(length(C)<K){
      stop("provide more sample cost value")
    }
    if(length(C)>K){
      stop("provided cost values are more than required")
    }
    if (!is.data.frame(data) && !is.matrix(data) && !is.vector(data)){
      stop("The argument 'data' must be a data.frame or matrix with more than one column")
    }

    if (dim(data)[2] <! 1 && !is.null(data) && !is.vector(data)){
      stop("The argument 'data' must have more than one column for pilot = TRUE")
    }
    if(na.rm){
      data <-  matrix(data[!is.na( data)], ncol=ncol(data))
    }
    sv <- c()
    mv <- c()
    optimal_N <- c()
    size_v <- c()
    temp_d <- c()

    for(i in 1:nrow(data)){
      for(j in 1:ncol(data)){
        temp_d <- c(temp_d, data[i,j])
      }
      sv <- c(sv, sqrt(stats::var(temp_d)))
      size_v <- c(size_v, length(temp_d))
      mv <- c(mv, mean(temp_d))
      temp_d <- c()
    }
    sm <- 0
    for(i in C){
      for(j in sv){
        sm <- i*j+sm
        o_N <- as.integer(((4*i*j*stats::qnorm(1-alpha/2)^2)/w^2)*sm)
      }
      optimal_N <- c(optimal_N, o_N)
    }

    if(verbose==FALSE) {
      outcome <- cbind(mean=mv, "Is.satisfied"=(size_v >= optimal_N))
    }
    if(verbose==TRUE) {
      outcome <- cbind(mean=mv,
                       optimal_N,
                       "Is.satisfied"=(size_v >= optimal_N))
    }
  }

  if(pilot==TRUE){
    Pilot.SS <- c()
    for(i in C){
      Pilot.SS <- c(Pilot.SS, max(2,(2*stats::qnorm(1-alpha/2)*sqrt(i))/w))
    }
    outcome <- (Pilot.SS=Pilot.SS)
  }
  return(outcome)
}

#' ts_br_exp_mean
#'
#' @description Two Stage approach to Get Bounded Risk Point
#'  Estimation for the mean of exponential random variables
#'
#'
#' @param data The data for which to calculate the bounded risk point.
#' @param A The loss function constant.
#' @param w The risk bound. A non-negative integer greater than zero.
#' @param pilot Should a pilot sample be generated. True/False value.
#' default value is \code{FALSE}.
#' @param verbose Should the criterion be printed. Default is \code{FALSE}.
#' @param na.rm This parameter controls whether NA values are removed from
#' the data prior to calculation. Default is \code{TRUE}.
#'
#' @return The current sample size, the mean and an indicator of if the
#' criterion is satisfied.
#'
#' @author Bhargab Chattopadhyay \email{Bhargab@iiitvadodara.ac.in},
#' Neetu Shah \email{201451015@iiitvadodara.ac.in}, Ken Kelley \email{kkelley@nd.edu}
#'
#' @references
#' Mukhopadhyay, N., \& de Silva, Basil M. (2009). \emph{Sequential Methods and Their Applications}. New York: CRC Press.
#'
#' @export ts_br_exp_mean
#'
#' @examples
#' pilot_ss <- ts_br_exp_mean(w=4, pilot=TRUE)
#' SLS <- rexp(pilot_ss, rate=2)
#' \dontrun{
#' SLS <- rexp(pilot_ss, rate=2)
#' }
#' SLS <- rexp(100, rate=2)
#' ts_br_exp_mean(SLS, A=100, w=4, pilot=FALSE)
ts_br_exp_mean <- function(data, A, w, pilot=FALSE, verbose=FALSE,
                           na.rm=TRUE)
{
  if(!missing(w)){
    if(w<=0) stop("w should be a non-zero positive value")
  }

  if(missing(w)){
    stop("You must specify \'w\'")
    }

  if(pilot==FALSE)
  {
    if(!missing(A))
    {
      if(A<=0) {
        stop("A should be a non-zero positive value")
      }
    }
    if(missing(A)){
      stop("You must specificy \'A\'")
      }

    if (!is.data.frame(data) && !is.matrix(data) && !is.vector(data)) {
      stop("The argument 'data' must be a data.frame or
           matrix with one column")
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

    B <- 0.565*(2*length(data)*(length(data)-1)*A)/((length(data)-1)
                                                    *(length(data)-2))

    optimal_N <- max(length(data),ceiling((B*mean(data)^2)/w))

    if(verbose==FALSE) {
      Outcome <- rbind(list(N=n, mean=mean(data),
                                             "Is.Satisfied?"=(n >= optimal_N),
                                             Stage=1))
    }
    if(verbose==TRUE){
      Outcome <- rbind(list(N=n, mean=mean(data),
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
        Outcome <- rbind(list(N=length(data), mean=mean(data), Stage=2))
      }
    }
  }

  if(pilot==TRUE){
    Outcome <- c(pilot_ss=max(3, ceiling(mean(stats::rchisq(3,2))/w)+1))
  }

  return(Outcome)
}

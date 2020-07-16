#' seq_mr_norm_mean
#'
#' @description Two Stage approach to Get minimum Risk Point Estimation
#'  for normal random variables.
#'
#' @param data The data for which to calculate the minimum risk point.
#' @param A The loss function constant.
#' @param c The cost of a sampling a unit.
#' @param price price
#' @param epsilon epislon
#' @param pilot Should a pilot sample be generated. TRUE/FALSE value.
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
#'
#' @import stats
#'
#' @export seq_mr_norm_mean
#'
#' @examples
#' pilot_ss <- seq_mr_norm_mean(A=100, c=0.5, pilot=TRUE)
#' SLS <- rnorm(pilot_ss)
#' seq_mr_norm_mean(data=SLS, A=100, c=0.5, pilot=FALSE)
seq_mr_norm_mean <- function(data, A, c, price, epsilon,
                             pilot=FALSE,
                             verbose=FALSE, na.rm=TRUE){

  if(!missing(A)){
    if(!missing(price)){
      stop("Because you specified \'A\' directly,
                             you should not also specify \'price\'.")
    }
    if(!missing(epsilon)){
      stop("Because you specified \'A\' directly,
                               you should not also specify \'epsilon\'.")
    }
    if(A<=0){
      stop("A should be a non-zero positive value")
    }
  }

  if(missing(A)){
    if(missing(price)) {
        stop("Because you did not specificy \'A\' directly,
                            you must specify \'price\'.")
    }
    if(missing(epsilon)){
      stop("Because you did not specificy \'A\' directly,
                              you must specify \'epsilon\'.")

    }
    A <- price/(epsilon^2)
    if(A<=0) {
      stop("A should be a non-zero positive value")
    }
  }

  if(pilot==FALSE) {
    if (!is.data.frame(data) && !is.matrix(data) && !is.vector(data)){
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

    optimal_N <- as.integer((A/c)^(1/2))*sqrt((stats::var(data)))

    Rk<-(1/n)*A*var(data)+(c*n)

    if(verbose==FALSE){
      Outcome <- rbind(list(Risk=Rk, N=n,
                            mean=mean(data),
                            "Is.Satisfied?"=(n >= optimal_N)))
    }
    if(verbose==TRUE){
      Outcome <- rbind(list(Risk=Rk, N=n,
                            mean=mean(data),
                            Criterion=optimal_N,
                            Is.Satisfied=n >= optimal_N))
    }
  }

  if(pilot==TRUE) {
    Outcome <- c(Pilot.SS=max(4, ceiling((A/c)^(1/2))))
  }

  return(Outcome)
}

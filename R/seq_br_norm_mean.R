#' seq_br_norm_mean
#'
#' @description Sequential approach to Get Bounded Risk Point Estimation
#' for normal random variables. Calculates the risk for a normal
#' random variables.
#'
#' @param data The data for which to calculate the bounded risk point.
#' @param A The loss function constant.
#' @param k Controls whether the absolute error or squared error is used.
#' K=1 uses absolute error, K=2 uses squared error.
#' @param w The risk bound.
#' @param pilot Should a pilot sample be generated. True/False value.
#' default value is \code{FALSE}.
#' @param verbose Should the criterion value be returned, True/False value,
#' default value is \code{FALSE}.
#' @param na.rm This parameter controls whether NA values are removed from
#' the data prior to calculation. Default is \code{TRUE}.
#'
#'
#' @return The calculated risk, the sample size, mean, standard deviation, and
#' an indicator of if the criterion was satisfied.
#'
#' @author Bhargab Chattopadhyay \email{Bhargab@iiitvadodara.ac.in},
#' Neetu Shah \email{201451015@iiitvadodara.ac.in}
#'
#' @references
#' Mukhopadhyay, N., \& de Silva, Basil M. (2009). \emph{Sequential Methods and Their Applications}. New York: CRC Press.
#'
#' @import stats
#'
#' @export seq_br_norm_mean
#'
#' @examples
#' pilot_ss <- seq_br_norm_mean(A=2, k=1, w=0.4, pilot=TRUE)  # k=1 absolute error,
#' # k=2 squared error
#' SLS <- rnorm(pilot_ss, mean=2, sd=3)
#' seq_br_norm_mean(data=SLS, A=2, k=1, w=0.4, pilot=FALSE)
seq_br_norm_mean <- function(data = NULL, A, k, w, pilot=FALSE, verbose=FALSE,
                             na.rm=TRUE)
{
  if(missing(w)){
    stop("You must specify  \'w\'")

  }
  if(missing(k)){
    stop("You must specify  \'k\'")
  }

  if(!missing(k)){
    if(k!=1 & k!=2) stop("k should be either 1 or 2")
  }

  B <- ((2^(k/2)*A)/sqrt(pi))*gamma((k+1)/2)

  if(pilot==FALSE)
  {
    if(!missing(A))
    {
      if(A<=0) stop("A should be a non-zero positive value")
    }

    if(missing(A)){
      stop("You must specify \'A\'")
    }

    if (!is.data.frame(data) && !is.matrix(data) && !is.vector(data)){
      stop("The argument 'data' must be a data.frame or matrix
           with one column")
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

    Criterion <- as.integer((B/w)^(2/k)*stats::var(data))

    CV <- c(mean = mean(data), sd = sqrt(stats::var(data)))

    Rk<-B*(stats::var(data)^k)

    if(verbose==FALSE){
      Outcome <- list(Risk=Rk, N=n, cv=t(CV),
                            "Is.Satisfied?"=(n >= Criterion))
    } else if(verbose==TRUE){
      Outcome <- list(Risk=Rk, N=n, cv=t(CV),
                           Criterion=Criterion,
                           "Is.Satisfied?"= (n >= Criterion))
    }
  }

  if(pilot==TRUE)
  {
    Outcome <- c(Pilot.SS=max(k, (B/w)^(2/k)))
  }

  return(Outcome)
}

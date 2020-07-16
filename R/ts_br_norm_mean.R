#' ts_br_norm_mean
#' 
#' Two Stage approach to Get Bounded Risk Point
#'  Estimation for the mean of normal random variables.
#'
#'
#' @param data The data for which to calculate the bounded risk point.
#' @param A The loss function constant.
#' @param k Controls whether the absolute error or squared error is used.
#' K=1 uses absolute error, K=2 uses squared error.
#' @param w The risk bound.
#' @param pilot Should a pilot sample be generated. True/False value.
#' default value is \code{FALSE}.
#' @param verbose Should the criterion be printed. Default is \code{FALSE}.
#' @param na.rm This parameter controls whether NA values are removed from
#' the data prior to calculation. Default is \code{TRUE}.
#'
#' @return The calculated risk, the sample size, the mean and the number of
#'  stages,
#' and an indicator of if the criterion is satisfied.
#'
#' @author Bhargab Chattopadhyay \email{Bhargab@iiitvadodara.ac.in},
#' Neetu Shah \email{201451015@iiitvadodara.ac.in}, Ken Kelley \email{kkelley@nd.edu}

#'
#' @references
#' Mukhopadhyay, N., \& de Silva, Basil M. (2009). \emph{Sequential Methods and Their Applications}. New York: CRC Press.
#'
#' @export ts_br_norm_mean
#'
#' @examples
#' pilot_ss <- ts_br_norm_mean(A=2, k=1, w=0.5, pilot=TRUE)
#' # k=1 absolute error, k=2 squared error
#' \dontrun{
#' SLS <- rnorm(pilot_ss)
#' }
#' SLS <- rnorm(100)
#' ts_br_norm_mean(data=SLS, A=2, k=1, w=0.5, pilot=FALSE)
ts_br_norm_mean <- function(data, A, k, w, pilot=FALSE, verbose=FALSE,
                            na.rm=TRUE)
{

  if(!missing(A)){
    if(A<=0) stop("A should be a non-zero positive value")
  }

  if(missing(A)){
    stop("You must specify \'A\'")
  }

  if(missing(w)){
    stop("You must specify  \'w\'")
  }

  if(missing(k)){
    stop("You must specify  \'k\'")
  }

  if(!missing(k)){
    if(k!=1 & k!=2){
      stop("k should be either 1 or 2")
    }
  }

  B <- ((2^(k/2)*A)/sqrt(pi))*gamma((k+1)/2)

  if(pilot==FALSE){
    if (!is.data.frame(data) && !is.matrix(data) && !is.vector(data)) {
      stop("The argument 'data' must be a data.frame or
           matrix with one column")
    }

    if (dim(data)[2] != 1 && !is.null(data) && !is.vector(data)) {
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

    b <- (((n-1)/2)^(k/2)*gamma((n-k-1)/2))/gamma((n-1)/2)

    optimal_N <- max(length(data),
                     ceiling((((b*B)/w)^(2/k))*var(data)))

    Rk<-B*(var(data)^k)

    if(verbose==FALSE){
      Outcome <- rbind(list(Risk=Rk, N=n, mean=mean(data),
                            "Is.Satisfied?"=(n >= optimal_N),
                            Stage=1))
    }
    if(verbose==TRUE){
      Outcome <- rbind(list(Risk=Rk, N=n, mean=mean(data),
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
        Rk<-B*(stats::var(data)^k)
        Outcome <- rbind(list(Risk=Rk, N=length(data),
                              mean=mean(data), Stage=2))
      }
    }
  }

  if(pilot==TRUE)
  {
    Outcome <- c(pilot_ss=ceiling((B/w)^(2/k)))
  }

  return(Outcome)
}

#' seq_mr_smd
#'
#' @description Purely Sequential approach to Get Minimum Risk Point Estimation
#'  for the Standardized Mean Difference
#'
#'
#' @param data1 The first data vector for which to calculate the minimum
#' risk point.
#' @param data2 The second data vector for which to calculate the minimum
#' risk point.
#' @param A The loss function constant.
#' @param c1 The cost of unit sample for the first data vector.
#' @param c2 The cost of unit sample for the second data vector.
#' @param gamma gamma
#' @param pilot Should a pilot sample be generated. TRUE/FALSE value.
#' default value is \code{FALSE}.
#' @param verbose Should the criterion be printed. Default is \code{FALSE}.
#' @param na.rm This parameter controls whether NA values are removed from
#' the data prior to calculation. Default is \code{TRUE}.
#'
#' @return The calculated minimum risk point, the sample size of each data
#' vector, the mean of each vector, and an indicator of if the
#' criterion is satisfied.
#'
#' @author Ken Kelley \email{KKelley@@nd.edu},
#' Bhargab Chattopadhyay \email{Bhargab@iiitvadodara.ac.in},
#' Neetu Shah \email{201451015@iiitvadodara.ac.in}
#'
#' @references
#' Chattopadhyay, B., & Kelley, K. (2017). Estimating the standardized mean
#' difference with minimum risk: Maximizing accuracy and minimizing
#' cost with sequential estimation. Psychological Methods, 22(1), 94-113
#'
#' @export seq_mr_smd
#'
#'@author Bhargab Chattopadhyay \email{Bhargab@iiitvadodara.ac.in}, Ken Kelley \email{kkelley@nd.edu}
#'
#' @examples
#' pilot_ss <- seq_mr_smd(gamma=1, A=100, c1=4, c2=5, pilot=TRUE)
#' SLS1 <- rnorm(pilot_ss[1], mean=0, sd=1)
#' SLS2 <- rnorm(pilot_ss[2], mean=0, sd=1)
#' seq_mr_smd(data1=SLS1, data2=SLS2, A=100, c1=4, c2=5, pilot=FALSE)
seq_mr_smd <- function(data1, data2, A, c1, c2, gamma, verbose=FALSE,
                       pilot=FALSE, na.rm=TRUE){
  if(missing(A)){
    stop("You must specify \'A\'")
  }

  if(missing(c1)){
    stop("You must specify \'c1\'")
  }

  if(missing(c1)){
    stop("You must specify \'c2\'")
  }

  if(missing(gamma)){
    gamma <- 1
  }

  if(!missing(A)){
    if(A<=0){
      stop("A should be a non-negative integer")
    }
  }

  if(!missing(c1)){
    if(c1<=0){
      stop("c1 should be a non-negative integer")
    }
  }

  if(!missing(c2)){
    if(c2<=0){
      stop("c2 should be a non-negative integer")
    }
  }

  if(pilot==FALSE){
    if (!is.data.frame(data1) && !is.matrix(data1) && !is.vector(data1)) {
      stop("The argument 'data1' must be a data.frame or
           matrix with one column")
    }

    if (dim(data1)[2] != 1 && !is.null(data1) && !is.vector(data1)) {
      stop("The argument 'data1' must have only one column,
           or be 'NULL' for pilot = TRUE")
    }


    if (!is.data.frame(data2) && !is.matrix(data2) && !is.vector(data2)) {
      stop("The argument 'data2' must be a data.frame or
           matrix with one column")
    }

    if (dim(data2)[2] != 1 && !is.null(data2) && !is.vector(data2)) {
      stop("The argument 'data2' must have only one column,
           or be 'NULL' for pilot = TRUE")
    }


    if(is.data.frame(data1)) {
      data1 <- as.vector(data1)
    }

    if(is.data.frame(data2)){
      data2 <- as.vector(data2)
    }
    if(na.rm){
      data1 <-  data1[!is.na( data1)]
      data2 <-  data2[!is.na( data2)]
    }

    delta <- (mean(data1) - mean(data2))/((((length(data1)-1)*
                                              stats::var(data1))+
                                             ((length(data2)-1)*
                                                stats::var(data2)))/
                                            (length(data1)+length(data2)-2))

    optimal_n1 <- as.integer(sqrt(A/c1)*
                               (1+delta^2/2*(1+sqrt(A/((c2-c1)*
                                                         length(data1)^2+A)))))

    optimal_n2 <- as.integer(optimal_n1/sqrt(1+((c2-c1)*optimal_n1^2)/A))

    risk <- A*((1/length(data1)+1/length(data2))+
                 delta^2/(2*(length(data1)+
                               length(data2))))+c1*length(data1)+
      c2*length(data2)

    if(verbose==FALSE){
      outcome <- rbind(list(Risk=risk,
                            n1=length(data1),
                            n2=length(data2),
                            mean1=mean(data1),
                            mean2=mean(data2),
                            "Is.satisfied"=(length(data1)>=optimal_n1)))
    }
    if(verbose==TRUE) {
      outcome <- rbind(list(Risk=risk, n1=length(data1),
                            n2=length(data2),
                            mean1=mean(data1),
                            mean2=mean(data2),
                            Criterion1=optimal_n1,
                            Criterion2=optimal_n2,
                            "Is.satisfied"=(length(data1)>=optimal_n1)))
    }
  }

  if(pilot==TRUE){
    outcome <- c(data1=(Pilot.SS=max(2,ceiling((A/c1)^(1/(2+2*gamma))))),
                 data2=(Pilot.SS=max(2,ceiling((A/c2)^(1/(2+2*gamma))))))
  }
  return(outcome)
}

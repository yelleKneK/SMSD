#' V2_cor
#'
#' This function is used internally to calculate the estimate of the asymptotic variance of the coefficient of variation.
#' 
#' @param x The first vector.
#' @param y The second vector.
#' @param method The correlation method to be used.
#' @param na.rm This parameter controls whether NA values are removed from
#' the data prior to calculation. Default is \code{FALSE}.
#'
#' @return The estimated asymptotic variance of the calculated coefficient of variation.
#'
#' @export V2_cor
#'
#' @importFrom copula F.n
#
#' @author Bhargab Chattopadhyay \email{Bhargab@iiitvadodara.ac.in}, Ken Kelley \email{kkelley@nd.edu}
#'
#' @examples
#' x <- c(1:4)
#' y <- c(2:5)
#' V2_cor(x,y, method = "pearson")
#'
#'
V2_cor <- function(x, y, method = c("pearson", "kendall", "spearman"),
                   na.rm=FALSE)
{
  method <- match.arg(method)
  if(na.rm){
    x <-  x[!is.na( x)]
    y <-  y[!is.na( y)]
  }
  n <- length(x)
  r <- stats::cor(x, y)

  if (method == "pearson")
  {

    t1 <- mu4(x)/stats::var(x)^2
    t2 <- mu4(y)/stats::var(y)^2
    t3 <- 2*mu22(x,y)/(stats::var(x)*stats::var(y))
    t4 <- 4*mu22(x,y)/stats::cov(x,y)^2
    t5 <- 4*mu31(x,y)/(stats::cov(x,y)*stats::var(x))
    t6 <- 4*mu13(x,y)/(stats::cov(x,y)*stats::var(y))
    vn2 <- r^2/4*(t1+t2+t3+t4-t5-t6)

  }else if (method == "kendall"){
    U1 <- rank(x)/(length(x)+1)
    U2 <- rank(y)/(length(y)+1)
    U <- cbind(U1, U2)
    W <- 2 * copula::F.n(U, U) - U1 - U2
    vn2 <- 16*stats::var(W)

  }else{ #spearman
    U1 <- rank(x)/(n + 1)
    U2 <- rank(y)/(n + 1)
    V <- U1*U2 + colMeans(t(outer(U1, U1, "<="))*U2 + t(outer(U2, U2, "<="))*U1)
    vn2 <- 144*stats::var(V)
  }

  vn2 <- max(vn2, n^(-3))
  return(vn2)
}

#' seq_aipe_slr_slope
#'
#' @description Sequential approach to Accuracy in Parameter
#'  Estimation for Effect Sizes
#' (AIPE): Simple Linear Regression Slope
#'
#' @param alpha The significance level., default is 0.05.
#' @param omega omega
#' @param data The data for which to calculate the slope.
#' @param x The data vector for the independent variable.
#' @param y The data vector for the response.
#' @param pilot Should a pilot sample be generated.
#' @param m0 The initial sample size.
#' @param na.rm This parameter controls whether NA values are removed from
#' the data prior to calculation. Default is \code{TRUE}.
#'
#' @return The current sample size, the calculated slope, and an indicator of
#' if the criterion is satisfied.
#'
#' @author Ken Kelley \email{KKelley@@nd.edu},
#' Francis Bilson Darku \email{FBilsonD@nd.edu},
#' Bhargab Chattopadhyay \email{Bhargab@iiitvadodara.ac.in}
#'
#' @import stats
#'
#' @export seq_aipe_slr_slope
#'
#' @examples
#' pilot_ss <- seq_aipe_slr_slope(alpha=0.05, omega=0.2, pilot=TRUE)
#' SLS <- matrix( rnorm(pilot_ss[1],mean=0,sd=1),
#'  rnorm(pilot_ss[1], mean = 0, sd = 1), nrow=20, ncol= 2)
#' seq_aipe_slr_slope(alpha=0.05, omega=0.2,data = SLS)
#'
seq_aipe_slr_slope <- function(alpha=0.05, omega, data = NULL, x=NULL, y=NULL,
                              pilot=FALSE, m0=4, na.rm=TRUE) {

  if (missing(alpha) && missing(omega)){
    stop("You must specify \'omega\' and \'alpha\'.")
  }

  # Changed to x and y
  if (pilot == FALSE && is.null(data) && (is.null(x) | is.null(y))) {
    stop("For pilot = FALSE, provide 'data'  or both 'x' and 'y'")
  }

  if (pilot == FALSE && !is.null(data) && (!is.null(x) | !is.null(y))){
    stop("For pilot = FALSE, provide 'data' only or both 'x' and 'y'")
  }

  if (!is.null(data) && !is.data.frame(data) && !is.matrix(data)){
    stop("The argument 'data' must be a data.frame or matrix with two columns")
  }

  if (!is.null(data) && dim(data)[2] != 2){
    stop("The argument 'data' must have two columns, or be 'NULL' for pilot = TRUE")
  }

  if (!is.null(data) && dim(data)[1] < 4){
    stop("The argument 'data' must have at least 4 rows, or be 'NULL' for pilot = TRUE")
  }

  if (is.null(data) & (!is.null(x) | !is.null(y))){
    if (length(x) != length(y)){
      stop("Group.1 should be of the same length as Group.2")
    }
    if (length(x) < 4){
      stop("Group.1 and Group.2 should have at least 4 observations")
    }
  }

  if (alpha <= 0 || alpha >= 1){
    stop("alpha be betweeen 0 and 1")
  }

  if (!is.null(data)){
    x <- data[, 1]
    y <- data[, 2]
  }
  if (pilot == FALSE){
    if(na.rm){
      x <-  x[!is.na( x)]
      y <-  y[!is.na( y)]
    }
    stop <- FALSE
    n <- length(x)
    beta1 <- slr_slope(x, y)
    V2n <- V2_slr_slope(x, y)
    Criterion <- ceiling((2*stats::qnorm(1 - alpha/2)/omega)^2 * (V2n + 1/n))

    if (n >= Criterion){
      Stop <- TRUE
    } else if (n < Criterion){
      Stop <- FALSE
    }
    lci <- beta1 - stats::qnorm(1 - alpha/2)*sqrt(V2n/n)
    uci <- beta1 + stats::qnorm(1 - alpha/2)*sqrt(V2n/n)
    ci <- c(lci, uci)
    if (Stop == FALSE){
      cat("The stopping rule has not yet been met;
            sample size is not large enough")
      Outcome <- list("Current.n" = n, "Current.slope" = beta1,
                      "Is.Satisfied?" = Stop)
    }

    if (Stop == TRUE){
      cat("The stopping rule has been met; sample size is large enough.")
      Outcome <- list("Current.n" = n, "Current.slope" = beta1,
                      "Is.Satisfied?" = Stop,
                      "Confidence Interval"= ci)
    }

  }
  if(pilot==TRUE){
    if(m0 < 4) stop("The value of 'm0' must be 4 or greater.")
    Outcome <- max(m0, ceiling(2*stats::qnorm(1-alpha/2)/omega))
  }
  return(Outcome)
}


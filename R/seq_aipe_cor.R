#' @title seq_aipe_cor
#'
#' @description Sequential approach to Accuracy in Parameter
#' Estimation for Effect Sizes
#' (AIPE): Correlation Coefficient
#'
#' @param alpha The significance level., default is 0.05.
#' @param omega omega
#' @param data The data set, should have two columns.
#' @param Group.1 The data for the first group.
#' @param Group.2 The data for the second group.
#' @param pilot Should a pilot sample be generated.
#' @param m0 The initial sample size.
#' @param method The correlation method to be used.
#' @param na.rm This parameter controls whether NA values are removed from
#' the data prior to calculation. Default is \code{TRUE}.
#'
#' @return The current sample size, the current correlation and an
#' indicator of if the criterion has been satisfied.
#'
#' @author Ken Kelley \email{KKelley@@nd.edu},
#' Francis Bilson Darku \email{FBilsonD@nd.edu},
#' Bhargab Chattopadhyay \email{Bhargab@iiitvadodara.ac.in}
#'
#' @references
#' Kelley, K., Darku, F. B., \& Chattopadhyay, B. (2018). Accuracy in parameter estimation for a general class of effect sizes: A sequential approach. \emph{Psychological Methods}, \emph{23}, 226–243.
#' 
#' Kelley, K., Darku, F. B., \& Chattopadhyay, B. (2019). Sequential Accuracy in parameter estimation for population correlation coefficients. \emph{Psychological Methods}, \emph{24}, 492–515.
#' 
#' @examples
#' pilot_ss <- seq_aipe_cor(alpha=0.05, omega=0.2, pilot=TRUE)
#'   SLS <- matrix(c(rexp(pilot_ss, rate=0.05),
#'   rexp(pilot_ss, rate=0.05)), ncol = 2)
#' seq_aipe_cor(alpha=0.05, omega=0.2,data = SLS)
#'
#' @export seq_aipe_cor
#'
seq_aipe_cor <- function(alpha=0.05, omega, data = NULL, Group.1=NULL,
                         Group.2=NULL,
                        pilot=FALSE, m0=4,
                        method = c("pearson", "kendall", "spearman"),
                        na.rm=TRUE)
{
  method <- match.arg(method)
  if (missing(alpha) && missing(omega)){
    stop("You must specify \'omega\' and \'alpha\'.")
  }

  if (pilot == FALSE && is.null(data) &&
      (is.null(Group.1) | is.null(Group.2))){
    stop("For pilot = FALSE, provide 'data'  or both 'Group.1' and 'Group.2'")
  }

  if (pilot == FALSE && !is.null(data) &&
      (!is.null(Group.1) | !is.null(Group.2))){
    stop("For pilot = FALSE, provide 'data' only or
         both 'Group.1' and 'Group.2'")
  }

  if (!is.data.frame(data) && !is.matrix(data) && !is.null(data)){
    stop("The argument 'data' must be a data.frame or matrix with two columns")
  }

  if (!is.null(data) && dim(data)[2] != 2){
    stop("The argument 'data' must have two columns,
         or be 'NULL' for pilot = TRUE")
  }

  if (!is.null(data) && dim(data)[1] < 4){
    stop("The argument 'data' must have at least 4 rows,
         or be 'NULL' for pilot = TRUE")
  }

  if (is.null(data) && (!is.null(Group.1) | !is.null(Group.2))){
    if (length(Group.1) != length(Group.2)){
      stop("Group.1 should be of the same length as Group.2")
    }

    if (length(Group.1) < 4){
      stop("Group.1 and Group.2 should have at least 4 observations")
    }

  }

  if (alpha <= 0 || alpha >= 1){
    stop("alpha be betweeen 0 and 1")
  }

  if (!is.null(data)){
    Group.1 <- data[, 1]
    Group.2 <- data[, 2]
  }
  if (pilot == FALSE) {
    if(na.rm){
      Group.1 <-  Group.1[!is.na( Group.1)]
      Group.2 <-  Group.2[!is.na( Group.2)]
    }
    stop <- FALSE
    n <- length(Group.1)
    r <- stats::cor(Group.1, Group.2, method = method)
    V2n <- V2_cor(Group.1, Group.2, method = method)
    Criterion <- ceiling((2*stats::qnorm(1 - alpha/2)/omega)^2 * (V2n + 1/n))

    if (n >= Criterion) {
      Stop <- TRUE
    } else if (n < Criterion) {
      Stop <- FALSE
    }
    lci <- r - stats::qnorm(1 - alpha/2)*sqrt(V2n/n)
    uci <- r + stats::qnorm(1 - alpha/2)*sqrt(V2n/n)
    ci <- c(lci, uci)
    if (Stop == FALSE){
      cat("The stopping rule has not yet been met;
            sample size is not large enough")
      Outcome <- list("Current.n" = n, "Current.cor" = r,
                      "Is.Satisfied?" = Stop)
    } else if (Stop == TRUE) {
      cat("The stopping rule has been met; sample size is large enough.")
      Outcome <- list("Current.n" = n, "Current.cor" = r,
                      "Is.Satisfied?" = Stop,
                      "Confidence Interval"= ci)
    }

  }
  if(pilot==TRUE)
  {
    if(m0 < 4) stop("The value of 'm0' must be 4 or greater.")
    Outcome <- max(m0, ceiling(2*stats::qnorm(1-alpha/2)/omega))
  }
  return(Outcome)
}

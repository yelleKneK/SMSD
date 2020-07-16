#' seq_aipe_smd
#'
#' Sequential approach to Accuracy in Parameter
#' Estimation for Effect Sizes
#' (AIPE): Standardized Mean Difference
#'
#'
#' @param alpha The significance level. default is 0.05.
#' @param omega omega
#' @param data The data set for which to calculate the standardized mean
#' difference.
#' @param Group.1 The data vector for the first group.
#' @param Group.2 The data vector for the second group.
#' @param pilot Should a pilot sample be generated.
#' @param m0 The initial sample size.
#' @param na.rm This parameter controls whether NA values are removed from
#' the data prior to calculation. Default is \code{TRUE}.
#'
#' @return The current sample size, the calculated standardized mean
#' difference, and an indicator of if the criterion has been satisfied.
#'
#' @author Ken Kelley \email{KKelley@@nd.edu},
#' Francis Bilson Darku \email{FBilsonD@nd.edu},
#' Bhargab Chattopadhyay \email{Bhargab@iiitvadodara.ac.in}
#'
#' @references
#' Kelley, K., Darku, F. B., \& Chattopadhyay, B. (2018). Accuracy in parameter estimation for a general class of effect sizes: A sequential approach. \emph{Psychological Methods}, \emph{23}, 226–243.
#'
#' @examples
#' pilot_ss <- seq_aipe_smd(alpha=0.05, omega=0.2, pilot=TRUE)
#' SLS <- matrix( rnorm(pilot_ss[1],mean=0,sd=1),
#' rnorm(pilot_ss[1], mean = 0, sd = 1), nrow=20, ncol= 2)
#' seq_aipe_smd(alpha=0.05, omega=0.2,data = SLS)
#'
#' @export seq_aipe_smd
#'
seq_aipe_smd <- function(alpha=0.05, omega, data = NULL, Group.1=NULL,
                        Group.2=NULL, pilot=FALSE, m0=4,
                        na.rm=TRUE)
{

  if (missing(alpha) && missing(omega)){
    stop("You must specify \'omega\' and \'alpha\'.")
  }

  if (pilot == FALSE && is.null(data) && (is.null(Group.1) | is.null(Group.2)))
    {
    stop("For pilot = FALSE, provide 'data'  or both 'Group.1' and 'Group.2'")
    }
  if (pilot == FALSE && !is.null(data) && (!is.null(Group.1) | !is.null(Group.2))){
    stop("For pilot = FALSE, provide 'data' only or both 'Group.1' and 'Group.2'")
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

  if (is.null(data) & (!is.null(Group.1) | !is.null(Group.2))) {
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
    SMD <- V2_smd(Group.1, Group.2)
    V2n <- V2_smd(Group.1, Group.2)
    Criterion <- ceiling((2*stats::qnorm(1 - alpha/2)/omega)^2 * (V2n + 1/n))

    if (n >= Criterion) Stop <- TRUE
    if (n < Criterion) Stop <- FALSE
    lci <- SMD - stats::qnorm(1 - alpha/2)*sqrt(V2n/n)
    uci <- SMD + stats::qnorm(1 - alpha/2)*sqrt(V2n/n)
    ci <- c(lci, uci)
    if (Stop == FALSE)
    {
      cat("The stopping rule has not yet been met;
            sample size is not large enough")
      Outcome <- list("Current.n" = n, "Current.smd" = SMD,
                      "Is.Satisfied?" = Stop)
    }

    if (Stop == TRUE){
      cat("The stopping rule has been met; sample size is large enough.")
      Outcome <- list("Current.n" = n, "Current.smd" = SMD,
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

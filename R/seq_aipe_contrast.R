#' seq_aipe_contrast
#'
#' @description Sequential approach to Accuracy in Parameter Estimation
#' for a Contrast
#'
#' @param alpha The significance level., default is 0.05
#' @param omega omega,
#' @param data The dataset, default is NULL
#' @param coef The coefficients, default is NULL
#' @param balanced Is the data balanced, default is TRUE
#' @param Group.1 The first group
#' @param Group.2 The second group
#' @param pilot Is pilot of interest? Default is FALSE
#' @param m0 The initial sample size.
#'
#' @return The current sample size, the current contrast, and an indicator of
#' if the criterion is satisfied.
#'
#' @author Ken Kelley \email{KKelley@@nd.edu},
#' Francis Bilson Darku \email{FBilsonD@nd.edu},
#' Bhargab Chattopadhyay \email{Bhargab@iiitvadodara.ac.in}
#'
#' @export seq_aipe_contrast
#'
#' @examples
#' pilot_ss <- seq_aipe_contrast(alpha=0.05, omega=0.2, pilot=TRUE)
#' SLS <- matrix(c(rexp(pilot_ss, rate=0.05),
#' rexp(pilot_ss, rate=0.05)), ncol = 2)
#' seq_aipe_contrast(alpha=0.05, omega=0.2,data = SLS)
#'
seq_aipe_contrast <- function(alpha=0.05, omega, data = NULL, coef=NULL,
                             balanced=TRUE, Group.1 =NULL, Group.2 = NULL,
                             pilot=FALSE, m0=4){
  if (missing(alpha) && missing(omega)){
    stop("You must specify \'omega\' and \'alpha\'.")
  }

  if (pilot == FALSE && is.null(data) && (is.null(Group.1) |
                                          is.null(Group.2))){
    stop("For pilot = FALSE, provide 'data'  or both 'Group.1' and 'Group.2'")
  }

  if (pilot == FALSE && !is.null(data) &&
      (!is.null(Group.1) | !is.null(Group.2))){
    stop("For pilot = FALSE, provide 'data' only or both
         'Group.1' and 'Group.2'")
  }

  if (!is.data.frame(data) && !is.matrix(data) && !is.null(data)){
    stop("The argument 'data' must be a data.frame or matrix with two columns")
  }

  if (!is.null(data) && dim(data)[2] != 2){
    stop("The argument 'data' must have two columns, or be 'NULL' for pilot = TRUE")
  }

  if (!is.null(data) && dim(data)[1] < 4){
    stop("The argument 'data' must have at least 4 rows, or be 'NULL' for pilot = TRUE")
  }

  if (alpha <= 0 || alpha >= 1){
    stop("alpha must be betweeen 0 and 1")
  }


  if (pilot == FALSE) {
    stop <- FALSE
    if (is.null(colnames(data))){
      colnames(data) <- paste("x", 1:ncol(data), sep="")
    }
    CONT <- sum(coef*colMeans(data, na.rm = TRUE))

    if(balanced) {
      n <- dim(data)[1]
      V2n <- V2_contrast(data, coef)
      Criterion <- ceiling((2*stats::qnorm(1 - alpha/2)/omega)^2 * (V2n + 1/n))

      if (n >= Criterion) {
        Stop <- TRUE
      } else if (n < Criterion){
        Stop <- FALSE
      }
      lci <- CONT - stats::qnorm(1 - alpha/2)*sqrt(V2n/n)
      uci <- CONT + stats::qnorm(1 - alpha/2)*sqrt(V2n/n)
      ci <- c(lci, uci)
      total.n <- n*dim(data)[2]

    }else{
      n.k <- colSums(!is.na(data))
      S.k <- matrixStats::colSds(data, na.rm = TRUE)
      V2.k <- coef*S.k*sum(coef*S.k)
      Criterion <- ceiling((2*stats::qnorm(1 -
                                             alpha/2)/omega)^2 * (V2.k + 1/n))
      if (all(n.k >= Criterion)) {
        Stop <- TRUE
      } else if (any(n.k < Criterion)){
        Stop <- FALSE
      }
      lci <- CONT - stats::qnorm(1 - alpha/2)*sqrt(sum((coef*S.k)^2/n.k))
      uci <- CONT + stats::qnorm(1 - alpha/2)*sqrt(sum((coef*S.k)^2/n.k))
      ci <- c(lci, uci)
      total.n <- sum(n.k)
    }


    if (Stop == FALSE) {
      cat("The stopping rule has not yet been met; sample size is not large enough")
      if(balanced) {
        Outcome <- list("Current.n" = n, "Current.contrast" = CONT,
                        "Is.Satisfied?" = Stop)
      }else{
        Outcome <- list("Current.n" = n.k, "Current.contrast" = CONT,
                        "Is.Satisfied?" = Stop)
      }

    } else if (Stop == TRUE) {
      cat("The stopping rule has been met; sample size is large enough.")
      Outcome <- list("Current.n" = n, "Current.contrast" = CONT,
                      "Is.Satisfied?" = Stop,
                      "Confidence Interval"= ci)
    }

  }
  if(pilot == TRUE){
    if(m0 < 4) {
      stop("The value of 'm0' must be 4 or greater.")
    }
    Outcome <- max(m0, ceiling(2*stats::qnorm(1-alpha/2)/omega))
  }
  return(Outcome)
}

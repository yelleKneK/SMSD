% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/seq_aipe_std_mean.R
\name{seq_aipe_std_mean}
\alias{seq_aipe_std_mean}
\title{seq_aipe_std_mean}
\usage{
seq_aipe_std_mean(
  alpha = 0.05,
  omega,
  data = NULL,
  pilot = FALSE,
  m0 = 4,
  na.rm = TRUE
)
}
\arguments{
\item{alpha}{The significance level. default is 0.05.}

\item{omega}{omega}

\item{data}{The data for which to calculate the standardized.}

\item{pilot}{Should a pilot sample be generated.}

\item{m0}{The initial sample size.}

\item{na.rm}{This parameter controls whether NA values are removed from
the data prior to calculation. Default is \code{FALSE}.}
}
\value{
The current sample size, the standardized mean, an indicator of if
the criterion is satisfied and the confidence interval.
}
\description{
Sequential approach to Accuracy in Parameter Estimation for Effect Sizes
(AIPE): Standardized Mean
}
\examples{
pilot_ss <- seq_aipe_std_mean(alpha=0.05, omega=0.2, pilot=TRUE)
SLS <- rexp(pilot_ss, rate=0.05)
seq_aipe_std_mean(alpha=0.05, omega=0.2,data = SLS)

}
\references{
Kelley, K., Darku, F. B., & Chattopadhyay, B. (2018). Accuracy in parameter estimation for a general class of effect sizes: A sequential approach. \emph{Psychological Methods}, \emph{23}, 226–243.
}
\author{
Ken Kelley \email{KKelley@nd.edu},
Francis Bilson Darku \email{FBilsonD@nd.edu},
Bhargab Chattopadhyay \email{Bhargab@iiitvadodara.ac.in}
}

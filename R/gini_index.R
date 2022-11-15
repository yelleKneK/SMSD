#' Gini Index
#'
#' This function calculates the Gini Index for a vector \code{x}.
#'
#' @param x A vector containing the data for which the Gini index will be
#' calculated.
#' @param na.rm This parameter controls whether NA values are removed from
#' the vector (or not). Default is \code{TRUE}.
#' 
#' @return The calculated value of the Gini index. 
#'
#' @author Ken Kelley \email{kkelley@nd.edu},
#' Francis Bilson Darku \email{fbilsond@nd.edu},
#' Bhargab Chattopadhyay \email{bhargab@iimv.ac.in}
#'
#' @examples
#' # Household money from income brackets (2010) from Table C of
#' # https://en.wikipedia.org/wiki/Gini_coefficient. 
#' x <- c(13.7, 12.0, 10.9, 13.9, 17.7, 11.4, 12.1, 4.5, 3.9) 
#' gini_index(x)
#'
#' @export gini_index
#' 
gini_index <- function(x, na.rm=TRUE){
  if(na.rm){
    x <- x[!is.na(x)]
  }
  y <- sort(x)
  n <- length(x)
  i <- 1:n
  g <- (2*sum(as.numeric(i*y))/sum(y) - (n + 1))/(n - 1)
  return(g)
}

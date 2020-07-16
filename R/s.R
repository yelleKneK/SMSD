#' s
#'
#' This function is used internally to raise two vectors to the specified powers
#' and return the sum of their element-wise product. 
#'
#' @param dat A matrix with two columns
#' @param i power 1
#' @param j power 2
#'
#' @return The sum of \code{dat[,1]^i} and \code{dat[,2]^j}
#'
#' @examples
#' dat <- cbind(c(1:4),c(2:5))
#' s(dat,3,4)
#'
#' @export s
#'
s <- function(dat, i,j){
  return(sum(dat[,1]^i * dat[,2]^j))
}

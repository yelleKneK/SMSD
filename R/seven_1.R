#' seven_1
#'
#' This function replicates program 1 in seq 7 from the textbook.
#'
#' @param d Half the confidence interval width.
#' @param var The variance.
#' @param alpha The significance level.
#'
#' @return The calculate value of c.
#'
#' @export seven_1
#'
#' @references
#' Mukhopadhyay, N., \& de Silva, Basil M. (2009). \emph{Sequential Methods and Their Applications}. New York: CRC Press.
#'
#' @examples
#' seven_1(d = 0.2, var = 10, alpha=0.01)
#'
seven_1 <- function(d, var, alpha){
  a = log(1/alpha)
  c = (a * var)/d
  return(c)
}

################################################################################
#' Returns a function to retrieve the \eqn{\alpha}-quantile from a vector. To be used in
#' conjunction with apply.
#'
#' @name   alphaq
#' @export
#' @param  a  the quantile level \eqn{\alpha}
#'
#' @return Returns a function that takes an argument \code{x} and calculates the
#'         \eqn{\alpha}-quantile from the vector \code{x}.
#' @examples
#'
#'##############################################################################

alphaq <- function(a){
  aq <- function(x){
    quantile(x,prob = a)
  }
  return(aq)
}




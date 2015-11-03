################################################################################
#' Returns a function to retrieve the \eqn{\alpha}-quantile from a vector.
#' 
#' To be used in
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

################################################################################
#' Returns the argument names of a function.
#' 
#' To check if an estimate or simulate function is valid.
#'
#' @name   arg.names
#' @export
#' @param  fun a function from which the arguments should be retrieved from
#'
#' @return a vector of characters containing the names of the arguments of the function \code{fun} 
#'
#'##############################################################################

arg.names <- function(fun){
 if(class(fun) == "function"){
   return(names(as.list(args(fun)))[-length(as.list(args(fun)))])
 }else{
   stop("Argument is not a function")
 }
}

################################################################################
#' Compare arguments with character vector
#' 
#' To check if an estimate or simulate function is valid.
#'
#' @name   compare.arg.names
#' @export
#' @param  fun a function from which the arguments should be retrieved from
#' @param names names the arguments should be compared to
#'
#' @return Logical if the names of the arguments are identical to \code{names}
#'
#'##############################################################################

compare.arg.names <- function(fun,names){
  if(class(fun) == "function"){
    return(isTRUE(all.equal(arg.names(fun),names)))
  }else{
    stop("First argument is not a function")
  }
}





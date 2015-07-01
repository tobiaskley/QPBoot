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
################################################################################
#' Test if bootstrap holds the level on each frequence
#'
#' @name test_level
#'
#' @export
#'
#' @param model the model to perform the test on
#' @param alpha the quantile level \eqn(\alpha)
#'
#'
################################################################################
test_level <- function(model,param,alpha,n = 256,rep = 100){
for(r in 1:rep){
    data = model$simulate(param,n)
  levels = c(.25,.5,.75)
  ln = length(levels)
  erg = main(data,alpha = alpha,levels = levels,model = model)
  res = array(rep(0,length(getFrequencies(erg[[1]]))*ln^2),dim = c(length(getFrequencies(erg[[1]])),ln,ln))
  for(ln1 in 1:ln){
    for(ln2 in 1:ln){
      if (ln1 >= ln2){
        res[,ln1,ln2] = Re(getValues(erg[[1]],frequencies = getFrequencies(erg[[1]]))[,ln1,ln2,1])
      }
      else{
        res[,ln1,ln2] = Im(getValues(erg[[1]],frequencies = getFrequencies(erg[[1]]))[,ln2,ln1,1])
      }
    }
  }
  low = erg[[2]]
  up = erg[[3]]
}
  return(apply((res < low) | (res > up),c(2,3),sum))
}





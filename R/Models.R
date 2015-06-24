################################################################################
#' Autoregressive time series of order p
#'
#' @name models
#' @aliases ARp
#' @export
#'
#' @param p numeric sets the order of the AR(p) process
#'
#' @examples
#' model = getAR(2)
################################################################################
getAR <- function(p){
p = p
# Yule-Walker estimator
estimator <- function(x){
  ar(x,aic = FALSE,order.max = p)$ar
}
# for simulation
simulate <- function(param,n){
  arima.sim(list(ar = param),n,rand.gen = rnorm)
}

return(list(estimator = estimator, simulate = simulate) )
}



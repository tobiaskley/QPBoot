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

################################################################################
#' Moving-Averages time series of order q
#'
#' @name models
#' @aliases MAq
#' @export
#'
#' @param q numeric sets the order of the MA(q) process
#'
#' @examples
#' model = getMA(2)
################################################################################
getMA <- function(q){

  # ML estimator
  estimator <- function(x){
    coef(arima(x,order = c(0,0,q),include.mean = FALSE))
  }
  # for simulation
  simulate <- function(param,n){
    arima.sim(list(ma = param),n,rand.gen = rnorm)
  }

  return(list(estimator = estimator, simulate = simulate) )
}

################################################################################
#' GARCH process of order p,q
#'
#' @name models
#' @aliases garch
#' @export
#'
#' @param q numeric sets the order of the MA(q) process
#'
#' @examples
#' model = getMA(2)
################################################################################
getGARCH <- function(p,q){

  # ML estimator
  estimator <- function(x){
    coef(garchFit(formula = ~garch(p,q), include.mean = FALSE, data = x))
  }
  # for simulation
  simulate <- function(param,n){
    spec = garchSpec(model = param)
    as.ts(garchSim(spec = spec,n = n))
  }

  return(list(estimator = estimator, simulate = simulate) )
}


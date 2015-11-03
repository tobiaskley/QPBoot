NULL
#' @include tsModel-class.R
#' 
################################################################################
#' @title Predefined Time-Series Models
#' 
#' @description  Creates a \link{tsModel-class} object representing a time-series model
#' 
#'   
#' @param spec a list specifying the structure of the parameters of the model
#' 
#' @name Models
#' @aliases ARMA
#' @export
#'
#'
################################################################################
getARMA <- function(spec = list(ar.order=NA,ma.order=NA)){
  
  arma = new("tsModel",spec = spec,name = paste("ARMA","(",spec$ar.order,",",spec$ma.order,")",sep = ""))
  
  
  # Yule-Walker estimator
  estimate <- function(data,spec){
    est = arima(data,order = c((spec$ar.order),0,(spec$ma.order)),include.mean = FALSE)
    #return(list(as.numeric(coef(est))))
    return(list(ar = as.numeric(coef(est))[1:spec$ar.order],ma = as.numeric(coef(est))[(spec$ar.order+1):(spec$ar.order+spec$ma.order)]))
  }
  # for simulation
  Simulate <- function(n,spec,par){
   arima.sim(list(order = c(spec$ar.order,0,spec$ma.order),ar = par$ar,ma=par$ma),n = n,rand.gen = rnorm)}
  
  setEstimate(arma,estimate)
  setSimulate(arma,Simulate)
 return(arma)
}


################################################################################
#' @name Models
#' @aliases AR
#' @export
#'
################################################################################
getAR <- function(spec = list(ar.order = 1)){
 ar = getARMA(spec = list(ar.order = spec$ar.order, ma.order = 0))
 ar@name = paste("AR","(",spec$ar.order,")",sep = "")
 return(ar)
}

################################################################################
#' @name Models
#' @aliases MA
#' @export
#'
################################################################################
getMA <- function(spec = list(ma.order = 1)){
  ma = getARMA(spec = list(ma.order = spec$ma.order, ar.order = 0))
  ma@name = paste("MA","(",spec$ma.order,")",sep = "")
  return(ar)
}
################################################################################
#' @name Models
#' @aliases GARCH
#' @export
#'
################################################################################
getGARCH <- function(spec = list(alpha = 1,beta = 1)){
  
  garch = new("tsModel",name = paste("GARCH","(",spec$alpha,",",spec$beta,")",sep = ""),spec = spec) 
  
  # QML estimator
  estimate <- function(data,spec){
    
    gspec = ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(spec$alpha,spec$beta)), mean.model = 
                                              list(armaOrder = c(0,0),include.mean = FALSE))
    fit = ugarchfit(gspec,data)
    par = as.numeric(coef(fit))
    return(list(mu = par[1],alpha = par[2],beta = par[3]))
    }
  
  # for simulation
  Simulate <- function(n,spec,par){
    gspec = ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(spec$alpha,spec$beta)), mean.model = 
                         list(armaOrder = c(0,0),include.mean = FALSE))
    setfixed(gspec) = list(omega = par$mu,alpha1 = par$alpha,beta1 = par$beta)
    gpath = ugarchpath(spec = gspec,n.sim = n)
    return(as.numeric(rugarch::fitted(gpath)))
  }
  setEstimate(garch,estimate)
  setSimulate(garch,Simulate)
  return(garch)
}

################################################################################
#' @name Models
#' @aliases EGARCH
#' @export
################################################################################
getEGARCH <- function(spec = list(alpha = 1,beta = 1)){
  
  
  garch = new("tsModel",name = paste("EGARCH","(",spec$alpha,",",spec$beta,")",sep = ""),spec = spec) 
  
  # QML estimator
  estimate <- function(data,spec){
    
    gspec = ugarchspec(variance.model = list(model = "eGARCH", garchOrder = c(spec$alpha,spec$beta)), mean.model = 
                         list(armaOrder = c(0,0),include.mean = FALSE))
    fit = ugarchfit(gspec,data)
    par = as.numeric(coef(fit))
    return(list(mu = par[1],alpha = par[2],beta = par[3],gamma = par[4]))
  }
  
  # for simulation
  Simulate <- function(n,spec,par){
    gspec = ugarchspec(variance.model = list(model = "eGARCH", garchOrder = c(spec$alpha,spec$beta)), mean.model = 
                         list(armaOrder = c(0,0),include.mean = FALSE))
    setfixed(gspec) = list(omega = par$mu,alpha1 = par$alpha,beta1 = par$beta,gamma1 = par$gamma)
    gpath = ugarchpath(spec = gspec,n.sim = n)
    return(as.numeric(fitted(gpath)))
  }
  setEstimate(garch,estimate)
  setSimulate(garch,Simulate)
  return(garch)
}

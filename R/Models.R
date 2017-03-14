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
    param = coef(est)
    #ar = as.numeric(est$coef[grepl("ar",names(est$coef))])    #get the ar Parameters
    #ma = as.numeric(est$coef[grepl("ma",names(est$coef))])    #get the ma Parameters
   
    return(as.numeric(param))
  }
  # for simulation
  Simulate <- function(n,spec,par){
   if(spec$ar.order > 0){ar = par[1:spec$ar.order]}else{ar=c()}
   if(spec$ma.order > 0){ma = par[(spec$ar.order + 1):(spec$ar.order + spec$ma.order)]}else{ma=c()}  
   param = list(ar = ar,ma = ma)
   arima.sim(param,n = n,rand.gen = rnorm)
    }
  
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
#' @aliases Noise
#' @export
#'
################################################################################
getNoise <- function(){
  noise = new("tsModel",name = "noise",spec = list()) 
estimate <- function(data,spec=list()){
  return(c(mean(data),var(data)))
}
simulate <- function(n,spec = list(),par){
  return(sqrt(par[2])*rnorm(n) + par[1])
  
}
setSimulate(noise,simulate)
setEstimate(noise,estimate)
return(noise)
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
    
    gspec = rugarch::ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(spec$alpha,spec$beta)), mean.model = 
                                              list(armaOrder = c(0,0),include.mean = FALSE))
    fit = rugarch::ugarchfit(gspec,data)
    par = coef(fit)
    return(par)
    }
  
  # for simulation
  Simulate <- function(n,spec,par){
    fixed.pars = par
    gspec = rugarch::ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(spec$alpha,spec$beta)), mean.model = 
                         list(armaOrder = c(0,0),include.mean = FALSE),fixed.par = fixed.pars)
    # setfixed(gspec) = par
    gpath = rugarch::ugarchpath(spec = gspec,n.sim = n)
    return(as.numeric(rugarch::fitted(gpath)))
  }
  setEstimate(garch,estimate)
  setSimulate(garch,Simulate)
  return(garch)
}
################################################################################
#' @name Models
#' @aliases fGARCH
#' @export
#'
################################################################################
getARCH <- function(spec = list(alpha = 1)){
  
  garch = new("tsModel",name = paste("ARCH","(",spec$alpha,")",sep = ""),spec = spec) 
  
  # QML estimator
  estimate <- function(data,spec){
    
    gspec = rugarch::ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(spec$alpha,0)), mean.model = 
                                  list(armaOrder = c(0,0),include.mean = FALSE))
    fit = rugarch::ugarchfit(gspec,data)
    par = coef(fit)
    return(par)
  }
  
  # for simulation
  Simulate <- function(n,spec,par){
    fixed.pars = par
    gspec = rugarch::ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(spec$alpha,0)), mean.model = 
                                  list(armaOrder = c(0,0),include.mean = FALSE),fixed.par = fixed.pars)
    # setfixed(gspec) = par
    gpath = rugarch::ugarchpath(spec = gspec,n.sim = n)
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
    par = coef(fit)
    return(par)
  }
  
  # for simulation
  Simulate <- function(n,spec,par){
    gspec = ugarchspec(variance.model = list(model = "eGARCH", garchOrder = c(spec$alpha,spec$beta)), mean.model = 
                         list(armaOrder = c(0,0),include.mean = FALSE))
    setfixed(gspec) = par
    gpath = ugarchpath(spec = gspec,n.sim = n)
    return(as.numeric(fitted(gpath)))
  }
  setEstimate(garch,estimate)
  setSimulate(garch,Simulate)
  return(garch)
}

################################################################################
#' @name Models
#' @aliases ARMA_GARCH
#' @export
#'
################################################################################
getARMA_GARCH <- function(spec = list(ar = 1, ma = 1,alpha = 1,beta = 1)){
  
  garch = new("tsModel",name = paste("ARMA_GARCH","(",spec$ar,",",spec$ma,",",spec$alpha,",",spec$beta,")",sep = ""),spec = spec) 
  
  # QML estimator
  estimate <- function(data,spec){
    gspec = ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(spec$alpha,spec$beta)), mean.model = 
                         list(armaOrder = c(spec$ar,spec$ma),include.mean = FALSE))
    fit = ugarchfit(gspec,data)
    par = coef(fit)
    return(par)
  }
  
  # for simulation
  simulate <- function(n,spec,par){
    gspec = ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(spec$alpha,spec$beta)), mean.model = 
                         list(armaOrder = c(spec$ar,spec$ma),include.mean = FALSE))
    setfixed(gspec) = par
    gpath = ugarchpath(spec = gspec,n.sim = n)
    return(as.numeric(rugarch::fitted(gpath)))
  }
  setEstimate(garch,estimate)
  setSimulate(garch,simulate)
  return(garch)
}

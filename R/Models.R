NULL

################################################################################
#' ARMA time series of order p,q
#' 
#' Returns two functions \code{estimate(x)} und \code{simulate(param,n)} and a string
#' that contains the name of the model and its parameters.
#' 
#' \enumerate{
#' \item \code{estimate(x)} estimates the parameters of the \code{ARMA(p,q)} Model by using
#' the \code{\link[stats]{arima}} function
#' }
#' @name getARMA
#' @aliases ARMApq
#' @export
#'
#' @param p numeric sets the order of the AR part
#' @param q numeric sets the order of the MA part
#'
################################################################################
getARMA <- function(p,q){
  # Yule-Walker estimator
  estimator <- function(x){
    est = arima(x,order = c(p,0,q),include.mean = FALSE)
    return(coef(est))
  }
  # for simulation
  simulate <- function(param,n){
   arima.sim(list(ar = param[1:p],ma= param[(p+1):(p+q)]),n = n,rand.gen = rnorm)}
  
  return(list(estimator = estimator, simulate = simulate, name = paste("ARMA(",p,q,")",sep = "")))
}


################################################################################
#' Autoregressive time series of order p
#'
#' @name getAR
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

return(list(estimator = estimator, simulate = simulate, name = paste("AR(",p,")",sep = "")) )
}

################################################################################
#' Moving-Averages time series of order q
#'
#' @name getMA
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

  return(list(estimator = estimator, simulate = simulate, name = paste("MA(",q,")",sep = "")) )
}

################################################################################
#' GARCH process of order (1,1)
#' 
#' @name getGARCH
#' @aliases GARCH
#' @export
#'
#'
################################################################################
getGARCH <- function(){

  # ML estimator
  estimator <- function(x){
    param = coef(garchFit(formula = ~garch(1,1), include.mean = FALSE,trace = FALSE, data = as.numeric(x)))
    if(param[2] + param[3] < 1){
      return(param)
    }else{
      while(param[2] + param[3] >= 1){
      warning(paste("Estimated parametes do not suffice the stationarity condition. Using a stationary approximation"),sep = "")
      p2 = param[2]
      p3 = param[3]
      param[2] = p2/(p2+p3) - 10^(-3)
      param[3] = p3/(p2+p3) - 10^(-3)
      }
      return(param)
    }
  }
  # for simulation
  simulate <- function(param,n){
    spec = garchSpec(model = list(omega = param[1],alpha = param[2],beta = param[3]))
    as.ts(garchSim(spec = spec,n = n))
  }

  return(list(estimator = estimator, simulate = simulate, name = paste("GARCH(1,1)",sep = "")) )
}

################################################################################
#' EGARCH process of order (1,1)
#' 
#' @name getEGARCH
#' @aliases EGARCH
#' @export
#'
#'
################################################################################
getEGARCH <- function(){
  
  # ML estimator
  estimator <- function(x){
    j = median(which(egarch_data == x,arr.ind = TRUE)[,2])
    return(egarch_param[,j])
    }
  # for simulation
  simulate <- function(param,n){
    Y <- rep(0,n+500)
    sig <- rep(0,n+500)
    Z = rnorm(n+500)
    Y[1] = Z[1]
    sig[1] = 0
    for(t in 2:(n+500)){
      sig[t] = param[1] + param[2]*abs(Z[t-1]) +  param[4]*Z[t-1] + param[3]*sig[t-1]
      Y[t] = exp(sig[t]/2)*Z[t]
    }
    return(Y[501:(n+500)])
     
  }
  
  return(list(estimator = estimator, simulate = simulate, name = paste("EGARCH(1,1)",sep = "")) )
}
################################################################################
#' AGARCH process of order (1,1)
#' 
#' @name getAGARCH
#' @aliases AGARCH
#' @export
#'
#'
################################################################################
getAGARCH <- function(){
  
  # ML estimator
  estimator <- function(x){
    
  }
  # for simulation
  simulate <- function(param,n){
    Y <- rep(0,n+500)
    sig <- rep(0,n+500)
    Z = rnorm(n+500)
    Y[1] = Z[1]
    sig[1] = 0
    for(t in 2:(n+500)){
      sig[t] = param[1] + param[2]*(abs(Y[t-1]) -  param[4]*Y[t-1])^2 + param[3]*sig[t-1]
      Y[t] = sqrt(sig[t])*Z[t]
    }
    return(Y[501:(n+500)])
    
  }
  
  return(list(estimator = estimator, simulate = simulate, name = paste("AGARCH(1,1)",sep = "")) )
}
################################################################################
#' GJR process of order (1,1)
#' 
#' @name getGJR
#' @aliases GJR
#' @export
#'
#'
################################################################################
getGJR <- function(){
  
  # ML estimator
  estimator <- function(x){
    
  }
  # for simulation
  simulate <- function(param,n){
    Y <- rep(0,n+500)
    sig <- rep(0,n+500)
    Z = rnorm(n+500)
    Y[1] = Z[1]
    sig[1] = 0
    for(t in 2:(n+500)){
      sig[t] = param[1] + param[2]*(Y[t-1])^2 +  param[4]*(Y[t-1] < 0)*(Y[t-1])^2 + param[3]*sig[t-1]
      Y[t] = sqrt(sig[t])*Z[t]
    }
    return(Y[501:(n+500)])
    
  }
  
  return(list(estimator = estimator, simulate = simulate, name = paste("GJR(1,1)",sep = "")) )
}
################################################################################
#' TGARCH process of order (1,1)
#' 
#' @name getTGARCH
#' @aliases TGARCH
#' @export
#'
#'
################################################################################
getTGARCH <- function(){
  
  # ML estimator
  estimator <- function(x){
    
  }
  # for simulation
  simulate <- function(param,n){
    Y <- rep(0,n+500)
    sig <- rep(0,n+500)
    Z = rnorm(n+500)
    Y[1] = Z[1]
    sig[1] = 0
    for(t in 2:(n+500)){
      sig[t] = param[1] + param[2]*max(Y[t-1],0) +  param[3]*max(-Y[t-1],0) + param[4]*sig[t-1]
      Y[t] = sig[t]*Z[t]
    }
    return(Y[501:(n+500)])
    
  }
  
  return(list(estimator = estimator, simulate = simulate, name = paste("TGARCH(1,1,1)",sep = "")) )
}
################################################################################
#' MGARCH process of order (1,1)
#' 
#' @name getMGARCH
#' @aliases MGARCH
#' @export
#'
#'
################################################################################
getMGARCH <- function(){
  
  # ML estimator
  estimator <- function(x){
    
  }
  # for simulation
  simulate <- function(param,n){
    Y <- rep(0,n+500)
    sig <- rep(0,n+500)
    Z = rnorm(n+500)
    Y[1] = Z[1]
    sig[1] = 0
    for(t in 2:(n+500)){
      sig[t] = param[1] + param[2]*Y[t-1]^2 + param[3]*sig[t-1]
      Y[t] = param[4]*sqrt(sig[t]) + sqrt(sig[t])*Z[t]
    }
    return(Y[501:(n+500)])
    
  }
  
  return(list(estimator = estimator, simulate = simulate, name = paste("TGARCH(1,1,1)",sep = "")) )
}

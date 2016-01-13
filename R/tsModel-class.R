#' @include Generics.R
#' @include aux-functions.R
################################################################################
#' @title Class for a Parametric Time-Series Model
#'
#' @description \code{tsModel} is a class to contain parametric time-series models (like ARMA
#'  oder GARCH) so that they can be used as arguments for \link{qpBoot}. There
#'  are some premade \link{Models} 
#'
#'
#' @name tsModel-class
#' @aliases tsModel
#' @exportClass tsModel
#'
#' @keywords S4-classes
#'
#' @slot name   the name of the model (e.g. "GARCH")
#' @slot spec   a list containing additional specification of the model
#' @slot env An environment to allow for slots which need to be
#'           accessable in a call-by-reference manner:
#'       \describe{
#'         \item{\code{est_function}}{a function implementing an estimator for the
#'         parameters of the model. It has the argumens \code{object} and \code{data} and 
#'         returns the estimated parameter. Also it sets \code{par} to the estimated value.}
#'         \item{\code{sim_function}}{a function implementing a way to simulate from the
#'         the model. It has the argumens \code{object} and \code{n}. Note that \code{par}
#'         has to be set in order to simulate.}
#'         \item{\code{par}}{a numeric vector that contains the parameters of the model. Can be
#'         empty at the beginning.}
#'          }
################################################################################
setClass(
  Class = "tsModel",
  representation=representation(
      name = "character",
      spec = "list",
      env = "environment"
      )
)



setMethod(
  f = "initialize",
  signature = "tsModel",
  definition = function(.Object,spec = list(),sim_function = function(){},est_function= function(){},par = NA,name=""){
    .Object@env <- new.env(parent=emptyenv())
    
    
    .Object@env$sim_function= sim_function
    .Object@env$est_function= est_function
    .Object@env$par = par
    
    .Object@spec = spec
    .Object@name = name
    return(.Object)
  }
)

################################################################################
#' @title
#' Estimates the parameter of a \link{tsModel-class}
#' 
#' @description 
#' The methods estimates the parameters from the time-series in \code{data}, using the method
#' defined in the \link{tsModel-class}.
#'
#' @return 
#' Returns the estimated parameters \code{par} and sets the \code{par} slot of the \link{tsModel-class}.
#' 
#' @name Estimate-tsModel
#' @aliases Estimate,tsModel-method
#' @export
#'
#' 
#'
#' @param object the \link{tsModel-class} for that the parameters shall be estimated
#' @param data a univariate numeric vector containing the time-series data
#'
#'
################################################################################
setMethod(
  f = "Estimate",
  signature = signature("tsModel"),
  definition = function(object,data){
    if(class(object@env$est_function)=="function"){
    object@env$par = object@env$est_function(data,object@spec)
    return(object@env$par)
    }
    else
    {
      stop("No proper estimation function defined, use ?setEstimate to learn more")
    }
  }
)

################################################################################
#' @title
#' Sets the estimation Method of a \link{tsModel-class}
#' 
#' @description 
#' Defines the method that should be used when calling \link{Estimate} on the \link{tsModel-class}.
#' The \code{Estimate} function must have exactly the following two arguments: \code{data} (numeric vector) and 
#' \code{spec} (a list) and return a list containing the estimated parameters.
#'
#' @return 
#' Sets the \code{est_function} slot of the \link{tsModel-class}.
#' 
#' @name setEstimate-tsModel
#' @aliases setEstimate,tsModel-method
#' @export
#'
#' 
#'
#' @param object the \link{tsModel-class} for that the estimation method shall be defined
#' @param estimate a function that has exactly two argumens: data (numeric vector) and spec (a list)
#'                 and returns a list containing the estimated parameters.
#'
#'
################################################################################


setMethod(
  f = "setEstimate",
  signature = signature("tsModel"),
  definition = function(object,estimate){
      if(!(class(estimate)=="function")){
        stop("estimator needs to be a function")
      } 
      if(!(compare.arg.names(estimate,c("data","spec")))){
        stop("estimator needs to be a function with exactly two arguments: data and spec")
      }
      
      object@env$est_function = estimate
  }
  
)

################################################################################
#' @title
#' Simulates from a \link{tsModel-class}
#' 
#' @description 
#' \code{Simulate} produces a numeric vector of length \code{n} with data simulated according to
#' the \code{sim_function} defined in the \link{tsModel-class}. To use this function the \code{par} slot
#' has to be set, either by calling \link{Estimate} or directly via \link{setParameter}.
#'
#' @return 
#' Returns a numeric vector of length \code{n}, that contains data, simulated according to the slot 
#' \code{sim_function} and \code{par} in the \link{tsModel-class} object.
#' 
#' @name Simulate-tsModel
#' @aliases Simulate,tsModel-method
#' @export
#'
#' 
#'
#' @param object the simulation will be based on this \link{tsModel-class}
#' @param n length of the simulated data, will be rounded down if it is not a whole number
#'
#'
################################################################################
setMethod(
  f = "Simulate",
  signature = signature("tsModel"),
  definition = function(object,n){
    if(class(object@env$sim_function)=="function"){
      if((n >= 1)){
        if(is.na(object@env$par[1])){stop("Parameter not set, use setParameter to do so")}
      return(object@env$sim_function(floor(n),spec = object@spec,par = object@env$par))
      }
      else
      {
        stop("n needs to be greater then 1")
      }
    }
    else
    {
      stop("No proper simulation function defined, use ?setSimulate to learn more")
    }
  }
)


################################################################################
#' @title
#' Sets the simulation Method of a \link{tsModel-class}
#' 
#' @description 
#' Defines the method that will be used when calling \link{Simulate} on the \link{tsModel-class}.
#' The passed \code{Simulate} function must have exactly the following three arguments: \code{n} (numeric),
#' \code{spec} (a list) and \code{par} (another list). It returns a numeric vector with the 
#' simulated data.
#'
#' @return 
#' Nothing, it sets the \code{sim_function} slot of the \link{tsModel-class}.
#' 
#' @name setSimulate-tsModel
#' @aliases setSimulate,tsModel-method
#' @export
#'
#' 
#'
#' @param object the \link{tsModel-class} for that the estimation method shall be defined
#' @param Simulate a function that has exactly three argumens: \code{n} (numeric),
#'                  \code{spec} (a list) and \code{par} (another list).
#'
#'
################################################################################


setMethod(
  f = "setSimulate",
  signature = signature("tsModel"),
  definition = function(object,Simulate){
    if(!(class(Simulate)=="function")){
      stop("Simulate needs to be a function")
    } 
    if(!(compare.arg.names(Simulate,c("n","spec","par")))){
      stop("Simulate needs to be a function with exactly three arguments: n, spec and par")
    }
    object@env$sim_function = Simulate
  }
  
)

################################################################################
#' @title
#' Sets the Parameter of a \link{tsModel-class} manually
#' 
#' @description 
#' This can be used to set the \code{par} slot of a \link{tsModel-class} by hand, in contrast
#' to call \link{Estimate} on a given data set. The parameter should be a list with the
#' name an the values of the parameters to be set. See the example below.
#'
#' @return 
#' Nothing, it sets the \code{par} slot of the \link{tsModel-class}.
#' 
#' @name setParameter-tsModel
#' @aliases setParameter,tsModel-method
#' @export
#'
#' 
#'
#' @param object the \link{tsModel-class} for that the Parameter will be set
#' @param par a list that contains the names and values of the parameters
#' 
#' 
#'
################################################################################

setMethod(
  f = "setParameter",
  signature = signature("tsModel"),
  definition = function(object,par){
    if(!(class(par)=="numeric")){
      par = as.list(par)
    } 
    object@env$par = par
  }
  
)

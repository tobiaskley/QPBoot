################################################################################
#' @title Generic functions for accessing attributes of objects
#'
#' These generic functions are needed to access or set the objects' attributes.
#'
#' @name generics-accessors
#'
#' @param object object from which to get the value
#' @param ... optional parameters; for documentation see the documentation of
#'             the methods to each of the generic.
#'


## Class-tsModel

#' @name generics-accessors
#' @aliases Estimate
#' @export
setGeneric("Estimate", function(object, ...){standardGeneric("Estimate")})

#' @name generics-accessors
#' @aliases Simulate
#' @export
setGeneric("Simulate", function(object, ...){standardGeneric("Simulate")})

#' @name generics-accessors
#' @aliases setEstimate
#' @export
setGeneric("setEstimate", function(object, ...){standardGeneric("setEstimate")})

#' @name generics-accessors
#' @aliases setSimulate
#' @export
setGeneric("setSimulate", function(object, ...){standardGeneric("setSimulate")})

#' @name generics-accessors
#' @aliases setParameter
#' @export
setGeneric("setParameter", function(object, ...){standardGeneric("setParameter")})
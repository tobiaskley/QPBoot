#' @include Models.R
NULL

################################################################################
#' Main function, which does 3 things
#' \enumerate{
#'  \item Estimates a parametric \code{model} from a given set of \code{data}
#'  \item Simulates from that \code{model} and computes the smoothed Quantile
#'        Periodogram (\link[quantspec]{smoothedPG}) for each simulated time
#'        series and the given \code{data}
#'  \item
#' }
#'
#'
#' @name main
#' @aliases main
#' @export
#'
#' @param data
#' @param alpha
#' @param model
#' @param levels
#' @param weight
#' @param SimNum
#'
#'
#' @examples
#' model = getAR(2)
################################################################################

main <- function(data,alpha = 0.05,model = getAR(2),levels = c(.25,.5,.75),weight = kernelWeight(bw = 0.025),SimNum = 100){
# Set
ln = length(levels)
n = length(data)

# Compute smoothed Periodogram of the data
sPG = smoothedPG(data,levels.1 = levels,weight = weight)

# Estimate the parametric model and Simulate from there
param = model$estimator(data)
fhat = list()
S = array(rep(0,SimNum*length(data)),dim = c(SimNum,length(data)))
res = array(rep(0,SimNum*(n/2+1)*ln^2),dim = c(SimNum,n/2+1,ln,ln))
for (i in 1:SimNum){
  S[i,] = model$simulate(param,length(data))
  fhat[[i]] = smoothedPG(S[i,],levels.1 = levels,weight = weight)
  res[i,,,] = getValues(fhat[[i]],frequencies = getFrequencies(fhat[[i]]))
}
 return(res)
}

################################################################################
#' Computes the \eqn{\alpha/2} and \eqn{1-\alpha/2} quantile for each frequence
#' from a list of \link[quantspec]{SmoothedPG}.
#'
#' @name getQuantiles-Main
#' @aliases getQuantiles,Main
#'
#' @param sPG
#' @param alpha
################################################################################
getQuantiles <- function(sPG,alpha){
  freq = getFrequencies(sPG[[1]])

  }

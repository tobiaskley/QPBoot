#' @include Models.R
#' @include generics.R
NULL
################################################################################
#' Class for a Parametric Bootstrap based on Quantile Spectral Analysis
#'
#' \code{QPBoot} is a class to compute und contain the results of a parametric
#' bootstrap based on Quantile Spectral Analysis.
#'
#'
#' @name QPBoot-class
#' @aliases QPBoot
#' @exportClass QPBoot
#'
#' @keywords S4-classes
#'
#' @slot data   the original data where the parametric bootstrap is based upon
#' @slot sPG    a smoothed quantile Periodogram from the \code{quantspec}-package
#'              calculatet from \code{data}
#' @slot model  parametric model for the bootstrap, this is a function that returns a
#'              simulate function and an estimate function see \link{models}
#' @slot param  parameter estimated for the \code{model} from the \code{data}
#' @slot sPGsim smoothed quantile periodograms of the simulated time series
################################################################################
setClass(
  Class = "QPBoot",
  representation=representation(
    sPG = "SmoothedPG",
    data = "numeric",
    model = "list",
    param = "numeric",
    sPGsim = "list"
  )
)

setMethod(
  f = "initialize",
  signature = "QPBoot",
  definition = function(.Object, data, sPG, model, sPGsim, param) {
    .Object@sPG = sPG
    .Object@data = data
    .Object@model = model
    .Object@sPGsim = sPGsim
    .Object@param = param
    
    return(.Object)
}
)
################################################################################
#' @title
#' Pointwise Confidence Intervalls
#' 
#' @description 
#' Depending on \code{method} this calculates pointwise confidence intervalls for
#' a smoothed periodgram that belongs to a time series defined by \code{model} and
#' \code{param}. If \code{(method = "quantiles")} it computes the \eqn{\alpha/2}
#' and \eqn{1-\alpha/2} quantiles from the Values of the simulated smoothed periodograms and
#' returns those. If \code{(method = "norm")} it uses the asymptotic normality
#' of the smoothed periodograms by estimating \code{mean} and 
#' \code{standard deviation} for each frequency and computing the \eqn{\alpha/2}
#' and \eqn{1-\alpha/2} quantiles from a normal distribution with the estimated
#' parameters.
#'
#' @return 
#' Returns a \code{list} with four elements
#' \item{\code{q_up}}{}
#' \item{\code{q_low}}{}
#' \item{\code{mean}}{}
#' \item{\code{sd}}{}
#' @name computeCIs-QPBoot
#' @aliases computeCIs,QPBoot-method
#' @export
#'
#' 
#'
#' @param object
#' @param alpha
#' @param method
#'
#'
################################################################################
computeCIs <- function(object,alpha = 0.05,method = c("quantiles","norm"),levels = object@sPG@levels[[1]]){
            #Extract Values
            ln = length(levels)
            freq = getFrequencies(object@sPGsim[[1]])
            fhat = object@sPGsim
            SimNum = length(fhat)
            n = length(freq)
            #Define Array
            SimValues = array(rep(0,SimNum*n*ln^2),dim = c(SimNum,n,ln,ln))
            
            for(i in 1:length(fhat)){
              fhatValues = getValues(fhat[[i]],levels.1 = levels,levels.2 = levels,frequencies = getFrequencies(fhat[[i]]))
            for(ln1 in 1:ln){
              for(ln2 in 1:ln){
                if (ln1 >= ln2){
                  SimValues[i,,ln1,ln2] = Re(fhatValues)[,ln1,ln2,1]
                }
                else{
                  SimValues[i,,ln1,ln2] = Im(fhatValues)[,ln2,ln1,1]
                }
              }
            }
            }

            method <- match.arg(method)[1]
            switch(method,
                   "quantiles" = {
                     #Compute Quantiles
                     q_low = rep(0,length(freq))
                     q_up = rep(0,length(freq))
                     
                     q_low = apply(SimValues,c(2,3,4),alphaq(alpha/2))
                     q_up = apply(SimValues,c(2,3,4),alphaq(1-alpha/2))
                     },
                   "norm" = {
                     #Estimate mean and variance
                     mean = apply(SimValues,c(2,3,4),mean)
                     sd  = apply(SimValues,c(2,3,4),sd)
                     q_low = qnorm(alpha/2,mean = mean, sd = sd)
                     q_up = qnorm(1-alpha/2,mean = mean, sd = sd)
                   }
            )
            
            #Return:
            return(list(q_low = q_low,q_up = q_up,mean = mean, sd = sd))
          }
################################################################################
#' Create an instance of the \code{QPBoot} class by doing 3 things
#' \enumerate{
#'  \item Estimates a parametric \code{model} from a given set of \code{data},
#'        this estimate can be overwritten by using the parameter \code{trueparam}
#'  \item Simulates from that \code{model} and computes the smoothed Quantile
#'        Periodogram (\link[quantspec]{smoothedPG}) for each simulated time
#'        series and the given \code{data}
#'  \item Returns an object of the class \link{QPBoot} with the calculated
#'        smoothed Periodograms
#' }
#'
#'
#' @name QPBoot-constructor
#' @aliases qpBoot
#' @export
#'
#' @keywords Constructors
#'
#' @param data
#' @param model
#' @param levels
#' @param weight
#' @param SimNum
#' @param trueparam
#'
#'
################################################################################

qpBoot <- function(data,model = getAR(2),levels = c(.1,.5,.9),weight = kernelWeight(bw = 0.1),SimNum = 1000,trueparam = NULL){
# Set
ln = length(levels)
n = length(data)

# Compute smoothed Periodogram of the data
sPG = smoothedPG(data,levels.1 = levels,weight = weight)

# Estimate the parametric model and Simulate from there
if (is.null(trueparam)){
param = model$estimator(data)
}else{
  param = as.numeric(trueparam)
}
#param =c(.1,.8,0)  #use true parameter

sPGsim = list()
S = array(rep(0,SimNum*length(data)),dim = c(SimNum,length(data)))
#Simulate SimNum times
pb = txtProgressBar()
for (i in 1:SimNum){
  setTxtProgressBar(pb,i/SimNum)
  S[i,] = model$simulate(param,length(data))
  sPGsim[[i]] = smoothedPG(S[i,],levels.1 = levels,weight = weight)
}
close(pb)
obj <- new("QPBoot",as.numeric(data),sPG,model,sPGsim,param)
return(obj)
}


################################################################################
#' Plot the values of a \code{\link{QPBoot}}.
#'
#' Creates a \code{K} x \code{K} plot depicting a smoothed quantile periodogram.
#' Optionally pointwise confidence intervals from the parametric bootstrap can
#' be displayed.
#' In each of the subplots either the real part (on and below the diagonal;
#' i. e., \eqn{\tau_1 \leq \tau_2}{tau1 <= tau2}) or the imaginary parts
#' (above the diagonal; i. e., \eqn{\tau_1 > \tau_2}{tau1 > tau2}) of
#' \itemize{
#'   \item the smoothed quantile periodogram (blue line),
#'   \item pointwise confidence intervals from the parametric bootstrap (light gray area),
#' }
#' for the combination of levels \eqn{\tau_1}{tau1} and \eqn{\tau_2}{tau2}
#' denoted on the left and bottom margin of the plot are displayed.
#'
#' @name plot-QPBoot
#' @aliases plot,QPBoot,ANY-method
#' @export
#'
#' @importFrom abind abind
#'
#' @param x  The \code{\link{SmoothedPG}} object to plot
#' @param ptw.CIs the confidence level for the conspec = garchSpec(model = param)fidence intervals to be
#'                 displayed; must be a number from [0,1]; if null, then no
#'                 confidence intervals will be plotted.
#' @param ratio quotient of width over height of the subplots; use this
#'               parameter to produce landscape or portrait shaped plots.
#' @param widthlab width for the labels (left and bottom); default is
#'                  \code{lcm(1)}, cf. \code{\link[graphics]{layout}}.
#' @param xlab label that will be shown on the bottom of the plots; can be
#'              an expression (for formulas), characters or \code{NULL} to
#'              force omission (to save space).
#' @param ylab label that will be shown on the left side of the plots;
#'              can be an expression (for formulas), characters or
#'              \code{NULL} to force omission (to save space).
#' @param type.scaling a method for scaling of the subplots; currently there
#'                      are three options: \code{"individual"} will scale each of the
#'                      \code{K^2} subplots to minimum and maximum of the values
#'                      in that plot, \code{"real-imaginary"} will scale each of the
#'                      subplots displaying real parts and each of the subplots
#'                      displaying imaginary parts to the minimum and maximum of
#'                      the values display in these subportion of plots. The
#'                      option \code{"all"} will scale the subplots to the minimum and
#'                      maximum in all of the subplots.
#' @param frequencies a set of frequencies for which the values are to be
#'                    plotted.
#' @param levels a set of levels for which the values are to be plotted.
#'
#' @return Returns the plot described in the Description section.
################################################################################

setMethod(f = "plot",
          signature = signature("QPBoot"),
          definition = function(x, ptw.CIs = 0.1,method = "quantiles", ratio = 3/2, widthlab = lcm(1), xlab = expression(omega/2*pi), ylab = NULL,
                                type.scaling = c("individual", "real-imaginary", "all"),
                                frequencies=x@sPG@frequencies,
                                levels=intersect(x@sPG@levels[[1]], x@sPG@levels[[2]])) {
            
            def.par <- par(no.readonly = TRUE) # save default, for resetting...
            
            # workaround: default values don't seem to work for generic functions?
            if (!hasArg(frequencies)) {
              frequencies <- x@sPG@frequencies
            }
            if (!hasArg(levels)) {
              levels <- intersect(x@sPG@levels[[1]], x@sPG@levels[[2]])
            }
            if (!hasArg(plotPG)) {
              plotPG <- FALSE
            }
            if (!hasArg(ptw.CIs)) {
              ptw.CIs <- 0.1
            }
            if (!hasArg(type.CIs)) {
              type.CIs <- "naive.sd"
            }
            if (!hasArg(ratio)) {
              ratio <- 3/2
            }
            if (!hasArg(widthlab)) {
              widthlab <- lcm(1)
            }
            if (!hasArg(xlab)) {
              xlab <- expression(omega/2*pi)
            }
            if (!hasArg(ylab)) {
              ylab <- NULL
            }
            if (!hasArg(type.scaling)) {
              type.scaling <- c("individual", "real-imaginary", "all")
            }
            # end: workaround
            
            if (length(levels) == 0) {
              stop("There has to be at least one level to plot.")
            }
            
            tryCatch({
              
              K <- length(levels)
              values <- getValues(x@sPG, frequencies = frequencies,
                                  levels.1=levels, levels.2=levels)
              if (ptw.CIs > 0) {
                CI <- computeCIs(x,alpha = ptw.CIs,method = method,levels = levels)
                lowerCIs  <- CI$q_low
                upperCIs  <- CI$q_up
                #text.headline <- (paste(text.headline, ", includes ",1-ptw.CIs,"-CI (ptw. of type '",type.CIs,"')",sep=""))
              }
              
              X <- frequencies/(2*pi)
              
              allVals <- array(values[,,,1], dim=c(length(X), K, K))
              if (plotPG)  {
                allVals <- abind(allVals, array(PG[,,,1], dim=c(length(X), K, K)), along=1)
              }
              if (hasArg(qsd)) {
                allVals <- abind(allVals, csd, along=1)
              }
              if (ptw.CIs > 0) {
                allVals <- abind(allVals, lowerCIs, upperCIs, along=1)
              }
              type.scaling <- match.arg(type.scaling)[1]
              
              p <- K
              M <- matrix(1:p^2, ncol=p)
              M.heights <- rep(1,p)
              M.widths  <- rep(ratio,p)
              
              # Add places for tau labels
              M <- cbind((p^2+1):(p^2+p),M)
              M.widths <- c(widthlab,M.widths)
              M <- rbind(M,c(0,(p^2+p+1):(p^2+2*p)))
              M.heights <- c(M.heights, widthlab)
              
              i <- (p^2+2*p+1)
              # Add places for labels
              if (length(xlab)>0) {
                M.heights <- c(M.heights, widthlab)
                M <- rbind(M,c(rep(0,length(M.widths)-p),rep(i,p)))
                i <- i + 1
              }
              
              if (length(ylab)>0) {
                M <- cbind(c(rep(i,p),rep(0,length(M.heights)-p)),M)
                M.widths <- c(widthlab,M.widths)
              }
              
              nf <- layout(M, M.widths, M.heights, TRUE)
              
              par(mar=c(2,2,1,1))
              
              for (i1 in 1:K) {
                for (i2 in 1:K) {
                  if (i1 <= i2) {
                    switch(type.scaling,
                           "individual" = {
                             y.min <- min(c(Re(allVals[,i1,i2]),lowerCIs[,i2,i1]))
                             y.max <- max(c(upperCIs[,i2,i1],Re(allVals[,i1,i2])))},
                           "real-imaginary" = {
                             y.min <- min(Re(allVals))
                             y.max <- max(Re(allVals))},
                           "all" = {
                             y.min <- min(Re(allVals),Im(allVals))
                             y.max <- max(Re(allVals),Im(allVals))}
                    )
                    
                    plot(x=0,y=0, type="n", xlab="", ylab="", #xlab=xl, ylab=yl,
                         xlim=c(min(X), max(X)), ylim=c(y.min, y.max))
                    if (ptw.CIs > 0) {
                      polygon(x=c(X,rev(X)), y=c(Re(lowerCIs[,i2,i1]),rev(Re(upperCIs[,i2,i1]))),
                              col="lightgray", border=NA)
                    }
                    if (plotPG) {
                      lines(x=X, y=Re(PG[,i1,i2,1]), col=gray(0.5))
                    }
                    if (hasArg(qsd)) {
                      lines(x=freq.csd/(2*pi), y=Re(csd[,i1,i2]), col="red")
                    }
                    lines(x=X, y=Re(values[,i1,i2,1]),
                          ylim=c(min(Re(allVals)), max(Re(allVals))),
                          type="l", col="blue")
                  } else {
                    switch(type.scaling,
                           "individual" = {
                             y.min <- min(c(Im(allVals[,i1,i2]),lowerCIs[,i2,i1]))
                             y.max <- max(c(upperCIs[,i2,i1],Im(allVals[,i1,i2])))},
                           "real-imaginary" = {
                             y.min <- min(Im(allVals))
                             y.max <- max(Im(allVals))},
                           "all" = {
                             y.min <- min(Re(allVals),Im(allVals))
                             y.max <- max(Re(allVals),Im(allVals))}
                    )
                    plot(x=0,y=0, type="n", xlab="", ylab="", #xlab=xl, ylab=yl,
                         xlim=c(min(X), max(X)), ylim=c(y.min, y.max))
                    if (ptw.CIs > 0) {
                      polygon(x=c(X,rev(X)), y=c((lowerCIs[,i2,i1]),rev((upperCIs[,i2,i1]))),
                              col="lightgray", border=NA)
                    }
                    if (plotPG) {
                      lines(x=X, y=Im(PG[,i1,i2,1]), col=gray(0.5))
                    }
                    if (hasArg(qsd)) {
                      lines(x=freq.csd/(2*pi), y=Im(csd[,i1,i2]), col="red")
                    }
                    lines(x=X, y=Im(values[,i1,i2,1]),
                          ylim=c(min(Im(allVals)), max(Im(allVals))),
                          type="l", col="blue")
                  }
                }
              }
              
              par(mar=c(0,0,0,0))
              for (i in 1:p) {
                plot.new()
                text(0.5,0.5,substitute(paste(tau[1],"=",k),list(k=levels[i])), srt=90)
              }
              
              for (i in 1:p) {
                plot.new()
                text(0.5,0.5,substitute(paste(tau[2],"=",k),list(k=levels[i])))
              }
              if (length(xlab)>0) {
                plot.new()
                text(0.5, 0.5, xlab)
              }
              if (length(ylab)>0) {
                plot.new()
                text(0.5, 0.5, ylab, srt=90)
              }
              
            },  error = function(e) e,
            warning = function(w) w,
            finally = {
              par(def.par)  #- reset to default
            })
          }
)







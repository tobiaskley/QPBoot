#' @include generic-accessors.R
NULL
################################################################################
#' Class for a Parametric Bootstrap based on Quantile Spectral Analysis
#'
#' \code{QPBoot} is a class to compute und contain the results of a parametric
#' bootstrap based on Quantile Spectral Analysis.
#'
#'
#' @name QPBoot-class
#' @aliases QPBoot-class
#' @exportClass QPBoot
#'
#' @keywords S4-classes
#'
#' @slot data   the original data where the parametric bootstrap is based upon
#' @slot sPG    a smoothed quantile Periodogram from the \code{quantspec}-package
#'              calculatet from \code{data}
#' @slot model  parametric model for the bootstrap from \link{tsModel-class}, some Examples
#'              can be found in  \link{Models}
#' @slot param  parameter estimated for the \code{model} from the \code{data}
#' @slot sPGsim smoothed quantile periodograms of the simulated time series
#' 
#' @importFrom methods setMethod
################################################################################
setClass(
  Class = "QPBoot",
  representation=representation(
    sPG = "SmoothedPG",
    data = "numeric",
    model = "tsModel",
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
#' If \code{(method = "uniform-wrt-taus")} then [TODO: explain the procedure].
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
#' @param object   the \link{QPBoot} object that will be plotted
#' @param alpha    the significiant level of the confidence intervalls, defaults
#' 								 to \code{0.05}
#' @param method   either "quantile" or "norm" or "uniform-wrt-taus",
#' 								 determines how the confidence intervalls are calculated.
#'                 see description for details
#' @param freq a vector of frequencies for which to compute the CIs
#' @param levels   numeric vector containing values between 0 and 1 for which the 
#'                 \link[quantspec]{smoothedPG}. Will be estimated. These are the
#'                 quantiles levels that are used for the validation
#'
#'
################################################################################
computeCIs <- function(object, alpha = 0.05,
                       method = c("quantiles", "norm", "uniform-wrt-taus"),
                       freq = getFrequencies(object@sPGsim[[1]]),
                       levels = object@sPG@levels[[1]]){
                     
            #Extract Values
            ln = length(levels)
            fhat = object@sPGsim
            SimNum = length(fhat)
            n = length(freq)
            #Define Array
            SimValues = array(rep(0,SimNum*n*ln^2),dim = c(SimNum,n,ln,ln))
            
            for(i in 1:length(fhat)){
              fhatValues = getValues(fhat[[i]],levels.1 = levels,levels.2 = levels,frequencies = freq)
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
                     aux <- apply(SimValues,c(2,3,4),quantile,
                              probs=c(1-alpha/2,alpha/2))
                     q_low <- aux[2,,,]
                     q_up  <- aux[1,,,]
                     },
                   "norm" = {
                     #Estimate mean and variance
                     mean = apply(SimValues,c(2,3,4),mean)
                     sd  = apply(SimValues,c(2,3,4),sd)
                     q_low = qnorm(alpha/2,mean = mean, sd = sd)
                     q_up = qnorm(1-alpha/2,mean = mean, sd = sd)
                   },
                   "multi-norm" = {
                     
                   },
                   "uniform-wrt-taus" = {
                     # compute mean over bootstrap replications
                     fbar <- apply(SimValues, c(2, 3, 4), mean)
                     # replicate SimNum times
                     fbar <- array(fbar, dim = c(dim(fbar), SimNum) )
                     fbar <- aperm(fbar, c(4, 1, 2, 3) )
                     
                     dif <- abs(SimValues - fbar)
                     
                     # define an array to scale the difference with
                     var_fak <- array(dim = c(ln, ln))
                     for(ln1 in 1:ln){
                       for(ln2 in 1:ln){
                         var_fak[ln1, ln2] <- min(levels[ln1], levels[ln2]) -
                                              levels[ln1] * levels[ln2]
                       }
                     }
                     var_fak <- array(var_fak, c(dim(var_fak), SimNum, n))
                     var_fak <- aperm(var_fak, c(3, 4, 1, 2))
                     
                     # two separate versions with only the Re or Im part
                     Re_dif <- dif
                     Im_dif <- dif
                     for (i in 2:ln) {
                       # block out the other part
                       Re_dif[ , , i-1, i:ln] <- 0
                       Im_dif[ , , i:ln, i] <- 0
                     }
                     Im_dif[ , , 1:ln, 1] <- 0
                     
                     max_Re_dif <- apply(Re_dif / var_fak, c(1, 2), max)
                     max_Im_dif <- apply(Im_dif / var_fak, c(1, 2), max)
                     
                     max_Re_dif <- apply(max_Re_dif, c(2), quantile,
                                               prob = 1-alpha/2)
                     max_Im_dif <- apply(max_Im_dif, c(2), quantile,
                                             prob = 1-alpha/2)
                     max_dif <- array(dim = c(n, ln, ln))
                     for (ln1 in 1:ln) {
                       for (ln2 in 1:ln) {
                         if (ln1 >= ln2) {
                           max_dif[ , ln1, ln2] <- max_Re_dif
                         } else {
                           max_dif[ , ln1, ln2] <- max_Im_dif
                         }
                         
                       }
                     }
                     
                     q_low = rep(0, length(freq))
                     q_up = rep(0, length(freq))
                     
                     q_low = fbar[1,,,,drop=TRUE] -
                             var_fak[1,,,,drop=TRUE] * max_dif
                     q_up = fbar[1,,,,drop=TRUE] +
                         var_fak[1,,,,drop=TRUE] * max_dif
                     
                   }
            )
            
            #Return:
            return(list(q_low = q_low,q_up = q_up,mean = mean, sd = sd))
          }
####################################################################################################################
#' qpBoot
#' 
#' Create an instance of the \code{QPBoot} class by doing 3 things
#' \enumerate{
#'  \item Estimates a parametric \code{model} from a given set of \code{data},
#'        this estimate can be overwritten by using the parameter \code{fix.param}
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
#' @param data         numeric vector, containing the time-series data
#' @param model        an object from the class \link{tsModel-class}. 
#' @param levels       numeric vector containing values between 0 and 1 for which the 
#'                     \link[quantspec]{smoothedPG}. Will be estimated. These are the
#'                     quantiles levels that are used for the validation
#' @param frequencies  a vector containing frequencies at which to determine the smoothed periodogram.
#' @param weight       an object of the class \link[quantspec]{KernelWeight} that is used to
#'                     in the estimation of the \link[quantspec]{smoothedPG}.
#' @param SimNum       number of bootstrap 
#' @param fix.param    defaults to \code{NULL}. In this case the parameters for the simulations are
#'                     estimated via the methode defined in the argument \code{model}. If this is not
#'                     \code{NULL}, it has to contain a list that can be used to set the parameters
#'                     in the \link{tsModel-class}. All simulations are then done with these fixed parameters.
#'
#'
####################################################################################################################

qpBoot <- function(data,model = getARMA(list(ar.order = 2,ma.order =0)),levels = c(.1,.5,.9),frequencies = 2 * pi/length(data) * 0:(length(data) -
                                                                                                         1), weight = kernelWeight(bw = 0.1),SimNum = 1000,fix.param = NULL){
# Set
ln = length(levels)
n = length(data)

# Compute smoothed Periodogram of the data
sPG = smoothedPG(data,levels.1 = levels,weight = weight,frequencies = frequencies)

# Estimate the parametric model and simulate from there
if (is.null(fix.param)){
param = Estimate(model,data)
print("Estimated Parameter:")
print(param)
}else{
  param = fix.param
  setParameter(model,param)
}

sPGsim = list()
S = array(rep(0,SimNum*length(data)),dim = c(SimNum,length(data)))

#Simulate SimNum times
pb = txtProgressBar()
for (i in 1:SimNum){
  setTxtProgressBar(pb,i/SimNum)
  S[i,] = Simulate(model,length(data))
  sPGsim[[i]] = smoothedPG(S[i,],levels.1 = levels,weight = weight,frequencies = frequencies)
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
#' \describe{
#'   \item{}{the smoothed quantile periodogram (blue line)}
#'   \item{}{pointwise confidence intervals from the parametric bootstrap (light gray area)}
#' }
#' for the combination of levels \eqn{\tau_1}{tau1} and \eqn{\tau_2}{tau2}
#' denoted on the left and bottom margin of the plot are displayed.
#' The \code{method} argument determines how the confidence intervalls are calculated.
#'  \describe{
#'   \item{quantile}{calculates the (1-\eqn{\alpha/2}) and \eqn{\alpha/2} quantiles from the bootstrap} 
#'   \item{norm}{asymptotic normality of the smoothed Periodograms is used, mean and standard deviation are
#'                 estimated from the bootstrap}
#' }
#' 
#'
#' @name plot-QPBoot
#' @aliases plot,QPBoot,ANY-method
#' @export
#'
#' @importFrom abind abind
#' 
#'
#' @param x  The \code{\link{SmoothedPG}} object to plot
#' @param ptw.CIs the confidence level for the conspec = garchSpec(model = param)fidence intervals to be
#'                 displayed; must be a number from [0,1]; if null, then no
#'                 confidence intervals will be plotted.
#' @param method either "quantile" or "norm", determines how the confidence intervalls are calculated.
#'               see description for details
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







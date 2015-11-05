##############################################################################################
#
#
#             QPBoot Demo analyzing the log-returns of the DAX 2000 - 2010 
#
#
##############################################################################################

# Plot the closing Values of the dataset
plot(dax,type = "l")

# Compute log-returns and plot those
rtn = log(dax[2:length(dax)]) - log(dax[1:(length(dax)-1)])
plot(rtn, type = "l")

# Try to fit a GARCH(1,1) Model to the data (this may take some minutes)
garch = getGARCH(spec = list(alpha = 1,beta = 1))
levels = c(.1,.5,.9)
bw = 0.1
qpB <- qpBoot(rtn, model = garch, levels = levels, weight = kernelWeight(bw = bw), SimNum = 100)

# Plot confidence intervalls using quantiles from the parametric bootstrap
plot(qpB, ptw.CIs = 0.1, method = "quantiles")

# And again using the asymptotic normality of the smoothed periodogram estimator
# Differences between the two methods are hard to spot
plot(qpB, ptw.CIs = 0.1, method = "norm")

# 






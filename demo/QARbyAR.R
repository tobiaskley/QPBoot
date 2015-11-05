##############################################################################################
#
#
#             QPBoot Demo analyzing an AR(2) and an QAR(1) process
#
#
##############################################################################################


# Simulate and plot one path of the QAR(1) process
Y <- ts1(2^12)
plot(Y, type = "l")

# Try to fit an AR(2) Model to the data (this may take some minutes)
arma = getARMA(spec = list(ar.order=2, ma.order=0))
levels = c(.25, .5, .75)
bw = 0.1
qpB <- qpBoot(Y, model = arma, levels = levels, weight = kernelWeight(bw = bw), SimNum = 100)

# Plot confidence intervalls using quantiles from the parametric bootstrap
plot(qpB, ptw.CIs = 0.1, method = "quantiles")

## Repeat with an AR(2) process

# Simulate and plot one path of the AR(2) process
Y <- ts2(2^12)
plot(Y, type = "l")

# Try to fit an AR(2) Model to the data (this may take some minutes)
arma = getARMA(spec = list(ar.order=2, ma.order=1))
levels = c(.25, .5, .75)
bw = 0.1
qpB <- qpBoot(Y, model = arma, levels = levels, weight = kernelWeight(bw = bw), SimNum = 100)

# Plot confidence intervalls using quantiles from the parametric bootstrap
plot(qpB, ptw.CIs = 0.1, method = "quantiles")

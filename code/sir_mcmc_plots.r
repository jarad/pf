# Load mcmc data
load("../data/sir_mcmc_test.rdata")

# Diagnostics
require(coda)
samps = mcmc.list(mcmc(out[[1]]$theta), mcmc(out[[2]]$theta), mcmc(out[[3]]$theta))
gelman.diag(samps)
summary(samps)

# Traceplots
n.samp = dim(samps[[1]])[1]
iter = 1001:100000
par(mfrow=c(2,1))
minbeta = min(sapply(samps, function(x) min(x[iter,1])))
maxbeta = max(sapply(samps, function(x) max(x[iter,1])))
mingamma = min(sapply(samps, function(x) min(x[iter,2])))
maxgamma = max(sapply(samps, function(x) max(x[iter,2])))
plot(iter,samps[[1]][iter,1],type="l",ylim=c(minbeta,maxbeta),xlab="",ylab=expression(beta))
lines(iter,samps[[2]][iter,1],col=2)
lines(iter,samps[[3]][iter,1],col=4)
plot(iter,samps[[1]][iter,2],type="l",ylim=c(mingamma,maxgamma),xlab="Iteration",ylab=expression(gamma))
lines(iter,samps[[2]][iter,2],col=2)
lines(iter,samps[[3]][iter,2],col=4)
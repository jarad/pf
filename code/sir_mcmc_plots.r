# Load mcmc data
load("../data/sir_mcmc_test.rdata")

# Diagnostics
require(coda)
samps = mcmc(list(chain1=out$theta))
gelman.diag(list(samps))

# Traceplots
iter = 1001:100000
par(mfrow=c(2,1))
plot(iter,out$theta[iter,1],type="l",xlab="",ylab=expression(beta))
plot(iter,out$theta[iter,2],type="l",xlab="Iteration",ylab=expression(gamma))
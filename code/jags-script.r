require(rjags)

# Load simulated data
load("../data/sim-orig.rdata")

# Set data values
y = t(mysim$sim$y[,1:60])
N = dim(y)[2]
L = length(mysim$true.params$b)
b = mysim$true.params$b
varsigma = mysim$true.params$varsigma
sigma = mysim$true.params$sigma
eta = mysim$true.params$eta
P = mysim$true.params$P
nu = mysim$true.params$theta[3]

# Test JAGS model
d = list(y=y,N=N,L=L,b=b,varsigma=varsigma,sigma=sigma,eta=eta,P=P,nu=nu)
inits = list(list(beta=.2399, gamma=.1066),list(beta=.15, gamma=.13),list(beta=.45, gamma=.095))
mod = jags.model("jags-model.txt", data=d, n.chains=3, n.adapt=1e3, inits=inits)
samps = coda.samples(mod, c("beta","gamma"), n.iter=1e5)
save(samps, file="../data/jags-script.rdata")

# Diagnostics
gelman.diag(samps)
summary(samps)

# Trace plots
n.samp = dim(samps[[1]])[1]
par(mfrow=c(2,1))
minbeta = min(sapply(samps, function(x) min(x[,1])))
maxbeta = max(sapply(samps, function(x) max(x[,1])))
mingamma = min(sapply(samps, function(x) min(x[,2])))
maxgamma = max(sapply(samps, function(x) max(x[,2])))
plot(1:n.samp,samps[[1]][,1],type="l",ylim=c(minbeta,maxbeta),xlab="",ylab=expression(beta))
lines(1:n.samp,samps[[2]][,1],col=2)
lines(1:n.samp,samps[[3]][,1],col=4)
plot(1:n.samp,samps[[1]][,2],type="l",ylim=c(mingamma,maxgamma),xlab="Iteration",ylab=expression(gamma))
lines(1:n.samp,samps[[2]][,2],col=2)
lines(1:n.samp,samps[[3]][,2],col=4)
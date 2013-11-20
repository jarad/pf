source("sir_functions.r")
source("pf_functions.r")

# Set data path
dpath = "../data/"

# Load simulated data and redefine data / known parameter values
load(paste(dpath,"sim-orig.rdata",sep=""))
y = mysim$sim$y
P = mysim$true.params$P
b = mysim$true.params$b
varsigma = mysim$true.params$varsigma
sigma = mysim$true.params$sigma
eta = mysim$true.params$eta
nu = mysim$true.params$theta[3]

# Define function to draw prior values of parameters and transform to original scale
theta.mean = c(-1.3296, -2.1764)
theta.sd = c(.3248, .1183)
ftheta = function(theta,param=1) exp(theta)
rtheta = function() rnorm(2, theta.mean, theta.sd)
n = 1000

# Run kernel density particle filter
dllik_kd = function(y, x, theta) dllik(y, x, b, varsigma, sigma, eta)
pstate_kd = function(x, theta) revo(x, c(ftheta(theta,1:2),nu), P, FALSE)
revo_kd = function(x, theta) revo(x, c(ftheta(theta,1:2),nu), P)
rprior_kd = function() rprior(rtheta)
source("kd_pf.r")
set.seed(62)
out_kd = kd_pf(y, dllik_kd, pstate_kd, revo_kd, rprior_kd, n, log=F, method="residual", nonuniformity = "ess", threshold = 0.8)
save(out_kd, file=paste(dpath,"sir_test_kd-1000.rdata",sep=""))

# Run resample-move particle particle
source("sir_mcmc_functions.r")
dllik_rm = function(y, x, theta) dllik(y, x, b, varsigma, sigma, eta)
revo_rm = function(x, theta) revo(x, c(theta[1:2],nu), P)
rtheta = function() exp(rnorm(2, theta.mean, theta.sd))
rprior_rm = function() rprior(rtheta)
rmove_rm = function(y, x, theta) rmove(y, x, c(theta[1:2], nu), list(b=b, varsigma=varsigma, sigma=sigma, eta=eta, P=P)) 
source("rm_pf.r")
set.seed(62)
out_rm = rm_pf(y, dllik_rm, revo_rm, rprior_rm, rmove_rm, n, log=F, method="residual", nonuniformity = "ess", threshold = 0.8)
save(out_rm, file=paste(dpath,"sir_test_rm-1000.rdata",sep=""))

# Overwrite old trajectories
require(plyr)
out.t = laply(out_rm$state, function(x) x[,,dim(x)[3]])
out_rm$state = aperm(out.t, c(2,3,1))

# Compute quantiles of states and plot
state.quant.kd = pf.quantile(out_kd$state, out_kd$weight, function(x, p=1) x, c(.025, .975))
state.quant.rm = pf.quantile(out_rm$state, out_rm$weight, function(x, p=1) x, c(.025, .975))
nt = dim(state.quant.kd)[1] - 1
plot(0:nt, state.quant.kd[,1,1], type="l", ylim=c(0,1), col = 3, xlab="Time",ylab="% Population",main="Estimated Epidemic Curves")
mtext(paste(n," particles",sep=""),side=3)
lines(0:nt, mysim$sim$x[1,])
lines(0:nt, mysim$sim$x[2,])
lines(0:nt, state.quant.kd[,1,2], col = 3)
lines(0:nt, state.quant.kd[,2,1], col = 3)
lines(0:nt, state.quant.kd[,2,2], col = 3)
lines(0:nt, state.quant.rm[,1,1], col = 6)
lines(0:nt, state.quant.rm[,1,2], col = 6)
lines(0:nt, state.quant.rm[,2,1], col = 6)
lines(0:nt, state.quant.rm[,2,2], col = 6)
legend("topright",c("Truth","KD","RM"),lty=c(1,1,1),col=c(1,3,6))

# Compute quantiles of parameters and plot
theta.quant.kd = pf.quantile(out_kd$theta, out_kd$weight, ftheta, c(.025, .975))
theta.quant.rm = pf.quantile(out_rm$theta, out_rm$weight, function(x,p=1) x, c(.025, .975))
windows(width = 10, height = 5)
par(mfrow=c(1,2))
plot(0:nt, theta.quant.kd[,1,1], type = "l", col = 3, main=paste(n, " particles",sep=""), xlab="Time", ylab=expression(beta), ylim = c(min(theta.quant.kd[,1,],theta.quant.rm[,1,]),max(theta.quant.rm[,1,],theta.quant.rm[,1,])))
lines(0:nt, theta.quant.kd[,1,2], col = 3)
lines(0:nt, theta.quant.rm[,1,1], col = 6)
lines(0:nt, theta.quant.rm[,1,2], col = 6)
abline(h=mysim$true.params$theta[1])
plot(0:nt, theta.quant.kd[,2,1], type = "l", col = 3, main=paste(n, " particles",sep=""), xlab="Time", ylab=expression(gamma), ylim = c(min(theta.quant.kd[,2,],theta.quant.rm[,2,]),max(theta.quant.kd[,2,],theta.quant.rm[,2,])))
lines(0:nt, theta.quant.kd[,2,2], col = 3)
lines(0:nt, theta.quant.rm[,2,1], col = 6)
lines(0:nt, theta.quant.rm[,2,2], col = 6)
abline(h=mysim$true.params$theta[2])


# Set data and graphics path
dpath = "../data/"
gpath = "../graphs/"

# Load mcmc data and simulated data
load(paste(dpath,"sir_mcmc_test.rdata",sep=""))
load(paste(dpath,"sim-orig.rdata",sep=""))

# Diagnostics on unknown parameters
require(coda)
samps = mcmc.list(mcmc(out[[1]]$theta), mcmc(out[[2]]$theta), mcmc(out[[3]]$theta))
gelman.diag(samps)
summary(samps)

# Traceplots on unknown parameters
n.samp = dim(samps[[1]])[1]
n.params = dim(samps[[1]])[2]
n.chains = 3
iter = 1:n.samp
ylabs = expression(beta,gamma,nu)
par(mfrow=c(n.params,1))
mins = apply(sapply(samps, function(x) apply(x, 2, min)), 1, min)
maxs = apply(sapply(samps, function(x) apply(x, 2, max)), 1, max)
for(i in 1:n.params)
{
  plot(iter,samps[[1]][iter,i],type="l",ylim=c(mins[i],maxs[i]),xlab="",ylab=ylabs[i])
  abline(h=mysim$true.params$theta[i])
  if(n.chains > 1)
  {
    for(j in 2:n.chains) lines(iter,samps[[j]][iter,i],col=2*(j-1))
  }
}

# Diagnostics on states
mystates = floor(seq(1, 126, len=126))
samps.states = mcmc.list(mcmc(cbind(out[[1]]$x[,1,mystates],out[[1]]$x[,2,mystates])), mcmc(cbind(out[[2]]$x[,1,mystates],out[[2]]$x[,2,mystates])), mcmc(cbind(out[[3]]$x[,1,mystates],out[[3]]$x[,2,mystates])))
gelman.diag(samps.states)

# Traceplots on (some) states
n.states = 6
mystates = floor(seq(1, 126, len=n.states))
samps.states = mcmc.list(mcmc(cbind(out[[1]]$x[,1,mystates],out[[1]]$x[,2,mystates])), mcmc(cbind(out[[2]]$x[,1,mystates],out[[2]]$x[,2,mystates])), mcmc(cbind(out[[3]]$x[,1,mystates],out[[3]]$x[,2,mystates])))
n.samp = dim(samps.states[[1]])[1]
n.params = dim(samps.states[[1]])[2]
n.chains = 3
iter = 1:n.samp
ylabs = paste(c(rep("s", n.states), rep("i", n.states)), rep(mystates,2), sep=" ")
windows(width=7.5,height=10)
par(mfrow=c(4,3))
mins = apply(sapply(samps.states, function(x) apply(x, 2, min)), 1, min)
maxs = apply(sapply(samps.states, function(x) apply(x, 2, max)), 1, max)
for(i in 1:n.params)
{
  plot(iter,samps.states[[1]][iter,i],type="l",ylim=c(mins[i],maxs[i]),xlab="",ylab=ylabs[i])
  abline(h=mysim$sim[[1]]$x[(i > n.states) + 1, mystates[6*(!(i %% n.states)) + (i %% n.states)]])
  if(n.chains > 1)
  {
    for(j in 2:n.chains) lines(iter,samps.states[[j]][iter,i],col=2*(j-1))
  }
}

# Acceptance rates
overall.chain = sapply(out, function(x) mean(cbind(x$accept.x,x$accept.theta)))
overall = mean(overall.chain)
states.chain = sapply(out, function(x) apply(x$accept.x, 2, mean))
states.overall.chain = apply(states.chain, 2, mean)
states.overall = mean(states.overall.chain)
theta.chain = sapply(out, function(x) apply(x$accept.theta, 2, mean))
theta.overall = apply(theta.chain, 1, mean)

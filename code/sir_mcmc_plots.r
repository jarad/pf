# Set data and graphics path
dpath = "../data/1-10-14/"
gpath = "../graphs/"

# Load simulated data
load(paste(dpath,"sim-orig.rdata",sep=""))

sir_mcmc_plots <- function(n.chains, x, beta, gamma, nu)
{
  # Load MCMC data
  params = which(as.logical(c(beta,gamma,nu)))
  states = as.logical(x)
  this.out = list()
  for(i in 1:n.chains)
  {
    load(paste(dpath,"sir_mcmc_test-",paste(i,x,beta,gamma,nu,sep="-"),".rdata",sep=""))
    this.out[[i]] = list()
    n.sims = out$mcmc.details$n.sims
    n.burn = out$mcmc.details$n.burn
    n.thin = out$mcmc.details$n.thin
    n.keep = (n.sims - n.burn) %/% n.thin
    if(states)
    {
      this.out[[i]]$x = out$x
      this.out[[i]]$accept.x = out$accept.x
    }
    if(length(params) > 0)
    {
      this.out[[i]]$theta = as.matrix(out$theta[,params])
      this.out[[i]]$accept.theta = out$accept.theta[params]
    }
  }
  out = this.out

  # Traceplots on unknown parameters
  if(length(params) > 0)
  {
    iter = (1:n.keep)*n.thin
    ylabs = expression(beta,gamma,nu)[params]
    pdf(file=paste(gpath,"sir_mcmc_test-",paste(x,beta,gamma,nu,sep="-"),"-traceplots-params.pdf",sep=""))
    par(mfrow=c(length(params),1), mar=c(5,6,4,2)+.1)
    mins = apply(as.matrix(sapply(out, function(x) apply(x$theta, 2, min))), 1, min)
    maxs = apply(as.matrix(sapply(out, function(x) apply(x$theta, 2, max))), 1, max)
    for(i in 1:length(params))
    {
      plot(iter,out[[1]]$theta[,i],type="l",ylim=c(mins[i],maxs[i]),xlab="Iteration",ylab=ylabs[i],cex.lab=1.5)
      if(n.chains > 1)
      {
        for(j in 2:n.chains) lines(iter,out[[j]]$theta[,i],col=2*(j-1))
      }
      abline(h=mysim$true.params$theta[params[i]])
    }
    dev.off()

    # Plot acceptance rates for parameters
    sink(file = paste(gpath,"sir_mcmc_test-",paste(x,beta,gamma,nu,sep="-"),"-acceptance-params.txt",sep=""))
    print("Acceptance rates of unknown parameters for each chain")
    theta.chain = sapply(out, function(x) x$accept.theta / n.sims)
    print(theta.chain)
    print("Acceptance rates of unknown parameters over all chains")
    theta.overall = mean(theta.chain)
    print(theta.overall)
    sink()
  }

  # Joint traceplots on (some) states
  if(states){
    nt = dim(out[[1]]$x)[3] - 1
    n.states = 9
    mystates = floor(seq(0, nt, len=n.states))
    xlabs = paste(rep("s", n.states), mystates, sep=" ")
    ylabs = paste(rep("i", n.states), mystates, sep=" ")
    mins.s = apply(sapply(out, function(x) apply(x$x[,1,mystates+1], 2, min)), 1, min)
    maxs.s = apply(sapply(out, function(x) apply(x$x[,1,mystates+1], 2, max)), 1, max)
    mins.i = apply(sapply(out, function(x) apply(x$x[,2,mystates+1], 2, min)), 1, min)
    maxs.i = apply(sapply(out, function(x) apply(x$x[,2,mystates+1], 2, max)), 1, max)
    pdf(file=paste(gpath,"sir_mcmc_test-",paste(x,beta,gamma,nu,sep="-"),"-traceplots-states.pdf",sep=""),width=4*sqrt(n.states),height=4*sqrt(n.states))
    par(mfrow=c(sqrt(n.states),sqrt(n.states)),mar=c(5,6,4,2)+.1)
    for(i in 1:n.states)
    {
      plot(out[[1]]$x[,1,mystates[i]+1],out[[1]]$x[,2,mystates[i]+1],type="l",xlim=c(mins.s[i],maxs.s[i]),ylim=c(mins.i[i],maxs.i[i]),xlab=xlabs[i],ylab=ylabs[i],cex.lab=2,cex.axis=1.5)
      mtext(paste("s = ",round(mysim$sim[[1]]$x[1,mystates[i]+1],3),", i = ",round(mysim$sim[[1]]$x[2,mystates[i]+1],3),sep=""),side=3,cex=.8)
      if(n.chains > 1) for(j in 2:n.chains) lines(out[[j]]$x[,1,mystates[i]+1],out[[j]]$x[,2,mystates[i]+1],col=2*(j-1))
      points(mysim$sim[[1]]$x[1,mystates[i]+1],mysim$sim[[1]]$x[2,mystates[i]+1],pch=20,cex=2,col="gray")
    }
    dev.off()

    # Marginal traceplots on (some) states
    iter = (1:n.keep)*n.thin
    xlabs = c("Iteration",rep("",n.states-1))
    ylabs.s = paste(rep("s", n.states), mystates, sep=" ")
    ylabs.i = paste(rep("i", n.states), mystates, sep=" ")
    mins.s = apply(sapply(out, function(x) apply(x$x[,1,mystates+1], 2, min)), 1, min)
    maxs.s = apply(sapply(out, function(x) apply(x$x[,1,mystates+1], 2, max)), 1, max)
    mins.i = apply(sapply(out, function(x) apply(x$x[,2,mystates+1], 2, min)), 1, min)
    maxs.i = apply(sapply(out, function(x) apply(x$x[,2,mystates+1], 2, max)), 1, max)
    pdf(file=paste(gpath,"sir_mcmc_test-",paste(x,beta,gamma,nu,sep="-"),"-traceplots-states-s.pdf",sep=""),width=4*sqrt(n.states),height=4*sqrt(n.states))
    par(mfrow=c(sqrt(n.states),sqrt(n.states)),mar=c(5,6,4,2)+.1)
    for(i in 1:n.states)
    {
      plot(iter,out[[1]]$x[,1,mystates[i]+1],type="l",ylim=c(mins.s[i],maxs.s[i]),xlab=xlabs[i],ylab=ylabs.s[i],cex.lab=2,cex.axis=1.5)
      mtext(paste("s = ", round(mysim$sim[[1]]$x[1,mystates[i]+1],3),sep=""),side=3)
      if(n.chains > 1) for(j in 2:n.chains) lines(iter,out[[j]]$x[,1,mystates[i]+1],col=2*(j-1))
      abline(h=mysim$sim[[1]]$x[1,mystates[i]+1])
    }
    dev.off()
    pdf(file=paste(gpath,"sir_mcmc_test-",paste(x,beta,gamma,nu,sep="-"),"-traceplots-states-i.pdf",sep=""),width=4*sqrt(n.states),height=4*sqrt(n.states))
    par(mfrow=c(sqrt(n.states),sqrt(n.states)),mar=c(5,6,4,2)+.1)
    for(i in 1:n.states)
    {
      plot(iter,out[[1]]$x[,2,mystates[i]+1],type="l",ylim=c(mins.i[i],maxs.i[i]),xlab=xlabs[i],ylab=ylabs.i[i],cex.lab=2,cex.axis=1.5)
      mtext(paste("i = ", round(mysim$sim[[1]]$x[2,mystates[i]+1],3),sep=""),side=3)
      if(n.chains > 1) for(j in 2:n.chains) lines(iter,out[[j]]$x[,2,mystates[i]+1],col=2*(j-1))
      abline(h=mysim$sim[[1]]$x[1,mystates[i]+1])
    }
    dev.off()

    # Print acceptance rates for states
    sink(file = paste(gpath,"sir_mcmc_test-",paste(x,beta,gamma,nu,sep="-"),"-acceptance-states.txt",sep=""))
    print("Acceptance rates of each state for each chain")
    states.chain = sapply(out, function(x) x$accept.x / n.sims)
    print(states.chain)
    print("Acceptance rates over all states for each chain")
    states.overall.chain = apply(states.chain, 2, mean)
    print(states.overall.chain)
    print("Acceptance rates over all states over all chains")
    states.overall = mean(states.overall.chain)
    print(states.overall)
    sink()
  }
}

mydata = data.frame(n.chains=3,x=1,beta=1,gamma=1,nu=1)
require(plyr)
m_ply(.data = mydata, .fun = sir_mcmc_plots)

# Compute 95% credible intervals and compare with KD particle filter
sir_mcmc_quant <- function(n.chains, x, beta, gamma, nu)
{
  for(i in 1:n.chains)
  {
    load(paste(dpath,"sir_mcmc_test-",paste(i,x,beta,gamma,nu,sep="-"),".rdata",sep=""))
    if(i == 1)
    {
      if(as.logical(x)) out.x = array(NA, dim = c(n.chains*dim(out$x)[1], dim(out$x)[2], dim(out$x)[3]))
      if(as.logical(sum(c(beta,gamma,nu)))) out.theta = matrix(NA, nr = n.chains*dim(out$theta)[1], nc = dim(out$theta)[2])
    }
    if(as.logical(x)) out.x[((i-1)*dim(out$x)[1] + 1):(i*dim(out$x)[1]),,] = out$x
    if(as.logical(sum(c(beta,gamma,nu)))) out.theta[((i-1)*dim(out$theta)[1] + 1):(i*dim(out$theta)[1]),] = out$theta
  }
  out.quant = list()
  if(as.logical(x)) out.quant$x = apply(out.x, 2:3, function(x) quantile(x, probs = c(0.025, 0.975)))
  if(as.logical(sum(c(beta,gamma,nu)))) out.quant$theta = apply(out.theta, 2, function(x) quantile(x, probs = c(0.025, 0.975)))

  return(out.quant)
}

# Calculate 95% CI for mcmc estimates of states and unknown parameters
mcmc.quant = sir_mcmc_quant(3,1,1,1,1)
colnames(mcmc.quant$theta) = c("beta","gamma","nu")

# Calculate 95% CI for kd_pf estimates of states and unknown parameters
load(paste(dpath,"PF-quant-KD-lognormal-stratified-20000-orig-0.99-61-1.rdata",sep=""))
kd.quant = list()
kd.quant$theta = t(pf.quant.out$theta.quant[126,,4:5])
dimnames(kd.quant$theta)[[1]] = c("2.5%","97.5%")
dimnames(kd.quant$theta)[[2]] = c("beta","gamma","nu")

# Plot 95% credible intervals for states over time
ymins.mcmc = apply(mcmc.quant$x[1,,], 1, min)
ymaxs.mcmc = apply(mcmc.quant$x[2,,], 1, max)
ymins.kd = apply(pf.quant.out$state.quant[,1:2,4], 2, min)
ymaxs.kd = apply(pf.quant.out$state.quant[,1:2,5], 2, max)
ymins = c(min(ymins.mcmc[1], ymins.kd[1], mysim$sim[[1]]$x[1,]), min(ymins.mcmc[2], ymins.kd[2], mysim$sim[[1]]$x[2,]))
ymaxs = c(max(ymaxs.mcmc[1], ymaxs.kd[1], mysim$sim[[1]]$x[1,]), max(ymaxs.mcmc[2], ymaxs.kd[2], mysim$sim[[1]]$x[2,]))
nt = dim(mcmc.quant$x)[3] - 1
windows(width = 10, height = 5)
par(mfrow = c(1,2))
plot(0:nt, mysim$sim[[1]]$x[1,], type = "l", xlab = expression(t), ylab = expression(s), main = "% Susceptible Over Time")
lines(0:nt, mcmc.quant$x[1,1,], col = 2)
lines(0:nt, mcmc.quant$x[2,1,], col = 2)
lines(0:nt, pf.quant.out$state.quant[,1,4], col = 4)
lines(0:nt, pf.quant.out$state.quant[,1,5], col = 4)
legend("topright", c("MCMC", "KD"), lty = c(1, 1), col = c(2, 4))
plot(0:nt, mysim$sim[[1]]$x[2,], type = "l", , xlab = expression(t), ylab = expression(i), main = "% Infected Over Time")
lines(0:nt, mcmc.quant$x[1,2,], col = 2)
lines(0:nt, mcmc.quant$x[2,2,], col = 2)
lines(0:nt, pf.quant.out$state.quant[,2,4], col = 4)
lines(0:nt, pf.quant.out$state.quant[,2,5], col = 4)

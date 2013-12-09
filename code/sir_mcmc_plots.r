# Set data and graphics path
dpath = "/storage/sheinson_research/"
gpath = "../graphs/"

# Load simulated data
load(paste(dpath,"sim-orig.rdata",sep=""))

sir_mcmc_plots <- function(n.chains, n.thin, tune.type, x, beta, gamma, nu)
{
  # Load MCMC data
  params = which(as.logical(c(beta,gamma,nu)))
  states = as.logical(x)
  this.out = list()
  for(i in 1:n.chains)
  {
    load(paste(dpath,"sir_mcmc_test-",paste(i,tune.type,x,beta,gamma,nu,sep="-"),".rdata",sep=""))
    this.out[[i]] = list()
    if(states)
    {
      this.out[[i]]$x = out$x
      this.out[[i]]$accept.x = out$accept.x
      this.out[[i]]$tuning$states = out$tuning$states
    }
    if(length(params) > 0)
    {
      this.out[[i]]$theta = as.matrix(out$theta[,params])
      this.out[[i]]$accept.theta = as.matrix(out$accept.theta[,params])
      this.out[[i]]$tuning$params = as.matrix(out$tuning$params[,params])
    }
  }
  out = this.out

  # Traceplots on unknown parameters
  if(length(params) > 0)
  {
    n.keep = dim(out[[1]]$theta)[1]
    iter = (1:n.keep)*n.thin
    ylabs = expression(beta,gamma,nu)[params]
    pdf(file=paste(gpath,"sir_mcmc_test-",paste(tune.type,x,beta,gamma,nu,sep="-"),"-traceplots-params.pdf",sep=""))
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
    sink(file = paste(gpath,"sir_mcmc_test-",paste(tune.type,x,beta,gamma,nu,sep="-"),"-acceptance-params.txt",sep=""))
    print("Acceptance rates of unknown parameters for each chain")
    theta.chain = sapply(out, function(x) apply(x$accept.theta, 2, mean))
    print(theta.chain)
    print("Acceptance rates of unknown parameters over all chains")
    theta.overall = apply(as.matrix(theta.chain), 1, mean)
    print(theta.overall)
    sink()
  }

  # Joint traceplots on (some) states
  if(states){
    nt = dim(mysim$sim[[1]]$y)[2]
    n.states = 9
    mystates = floor(seq(0, nt, len=n.states))
    n.keep = dim(out[[1]]$x)[1]
    iter = 1:n.keep
    xlabs = paste(rep("s", n.states), mystates, sep=" ")
    ylabs = paste(rep("i", n.states), mystates, sep=" ")
    mins.s = apply(sapply(out, function(x) apply(x$x[,1,mystates+1], 2, min)), 1, min)
    maxs.s = apply(sapply(out, function(x) apply(x$x[,1,mystates+1], 2, max)), 1, max)
    mins.i = apply(sapply(out, function(x) apply(x$x[,2,mystates+1], 2, min)), 1, min)
    maxs.i = apply(sapply(out, function(x) apply(x$x[,2,mystates+1], 2, max)), 1, max)
    pdf(file=paste(gpath,"sir_mcmc_test-",paste(tune.type,x,beta,gamma,nu,sep="-"),"-traceplots-states.pdf",sep=""),width=4*sqrt(n.states),height=4*sqrt(n.states))
    par(mfrow=c(sqrt(n.states),sqrt(n.states)),mar=c(5,6,4,2)+.1)
    for(i in 1:n.states)
    {
      plot(out[[1]]$x[iter,1,mystates[i]+1],out[[1]]$x[iter,2,mystates[i]+1],type="l",xlim=c(mins.s[i],maxs.s[i]),ylim=c(mins.i[i],maxs.i[i]),xlab=xlabs[i],ylab=ylabs[i],cex.lab=2,cex.axis=1.5)
      mtext(paste("s = ",round(mysim$sim[[1]]$x[1,mystates[i]+1],3),", i = ",round(mysim$sim[[1]]$x[2,mystates[i]+1],3),sep=""),side=3,cex=.8)
      if(n.chains > 1) for(j in 2:n.chains) lines(out[[j]]$x[iter,1,mystates[i]+1],out[[j]]$x[iter,2,mystates[i]+1],col=2*(j-1))
      points(out[[1]]$x[iter[1],1,mystates[i]+1],out[[1]]$x[iter[1],2,mystates[i]+1],pch="S",cex=2,col="gray")
      points(out[[1]]$x[iter[length(iter)],1,mystates[i]+1],out[[1]]$x[iter[length(iter)],2,mystates[i]+1],pch="F",cex=2,col="gray")
      if(n.chains > 1)
      {
        for(j in 2:n.chains)
        {
          points(out[[j]]$x[iter[1],1,mystates[i]+1],out[[j]]$x[iter[1],2,mystates[i]+1],pch="S",cex=2,col="gray")
          points(out[[j]]$x[iter[length(iter)],1,mystates[i]+1],out[[j]]$x[iter[length(iter)],2,mystates[i]+1],pch="F",cex=2,col="gray")
        }
      }
      points(mysim$sim[[1]]$x[1,mystates[i]+1],mysim$sim[[1]]$x[2,mystates[i]+1],pch=20,cex=2,col="gray")
    }
    dev.off()

    # Print acceptance rates for states
    sink(file = paste(gpath,"sir_mcmc_test-",paste(tune.type,x,beta,gamma,nu,sep="-"),"-acceptance-states.txt",sep=""))
    print("Acceptance rates of each state for each chain")
    states.chain = sapply(out, function(x) apply(x$accept.x, 2, mean))
    print(states.chain)
    print("Acceptance rates over all states for each chain")
    states.overall.chain = apply(states.chain, 2, mean)
    print(states.overall.chain)
    print("Acceptance rates over all states over all chains")
    states.overall = mean(states.overall.chain)
    print(states.overall)
    sink()
  }

  # Plot tuning parameters over time for parameters
  if(length(params) > 0)
  {
    n.sims = dim(out[[1]]$tuning$params)[1]
    iter = (1:n.sims)[1:n.sims %% n.thin == 0]
    ylabs = expression(tau[beta],tau[gamma],tau[nu])[params]
    pdf(file=paste(gpath,"sir_mcmc_test-",paste(tune.type,x,beta,gamma,nu,sep="-"),"-tunings-params.pdf",sep=""))
    par(mfrow=c(length(params),1),mar=c(5,6,4,2)+.1)
    mins = apply(as.matrix(sapply(out, function(x) apply(x$tuning$params, 2, min))),1, min)
    maxs = apply(as.matrix(sapply(out, function(x) apply(x$tuning$params, 2, max))),1, max)
    for(i in 1:length(params))
    {
      plot(iter,out[[1]]$tuning$params[iter,i],type="l",ylim=c(mins[i],maxs[i]),xlab="Iteration",ylab=ylabs[i],cex.lab = 1.5, main = "Tuning parameters for fixed parameters")
      if(n.chains > 1)
      {
        for(j in 2:n.chains) lines(iter,out[[j]]$tuning$params[iter,i],col=2*(j-1))
      }
    }
    dev.off()

    # Print final tuning parameter values
    sink(file = paste(gpath,"sir_mcmc_test-",paste(tune.type,x,beta,gamma,nu,sep="-"),"-finaltunings-params.txt",sep=""))
    print("Final tuning parameters of fixed parameters")
    print(sapply(out, function(x) x$tuning$params[dim(x$tuning$params)[1],]))
    sink()
  }

  # Plot tuning parameters over time for (some) states
  if(states)
  {
    nt = dim(mysim$sim[[1]]$y)[2]
    n.states = 9
    mystates = floor(seq(0, nt, len=n.states))
    n.sims = dim(out[[1]]$tuning$states)[3]
    iter = (1:n.sims)[1:n.sims %% n.thin == 0]
    ylabs = paste("t = ",mystates,sep="")
    pdf(file=paste(gpath,"sir_mcmc_test-",paste(tune.type,x,beta,gamma,nu,sep="-"),"-tunings-states-s.pdf",sep=""),width=4*sqrt(n.states),height=4*sqrt(n.states))
    par(mfrow=c(sqrt(n.states),sqrt(n.states)),mar=c(5,6,4,2)+.1)
    mins = apply(as.matrix(sapply(out, function(x) apply(x$tuning$states[1, mystates+1, ], 1, min))),1, min)
    maxs = apply(as.matrix(sapply(out, function(x) apply(x$tuning$states[1, mystates+1, ], 1, max))),1, max)
    for(i in 1:n.states)
    {
      plot(iter,out[[1]]$tuning$states[1,mystates[i]+1,iter],type="l",ylim=c(mins[i],maxs[i]),xlab="Iteration",ylab=ylabs[i],cex.lab = 2, cex.axis = 1.5, main="Tuning parameters for states")
      if(n.chains > 1)
      {
        for(j in 2:n.chains) lines(iter,out[[j]]$tuning$states[1,mystates[i]+1,iter],col=2*(j-1))
      }
    }
    dev.off()

  # Print final tuning parameter values
  sink(file = paste(gpath,"sir_mcmc_test-",paste(tune.type,x,beta,gamma,nu,sep="-"),"-finaltunings-states.txt",sep=""))
  print("Final tuning parameters of states")
  print(sapply(out, function(x) x$tuning$states[1,,dim(x$tuning$states)[3]]))
  sink()
  }
}

mydata = matrix(nr=0,nc=0)
#data = data.frame(x=c(0,0,0,1,1),beta=c(1,0,0,0,1),gamma=c(0,1,0,0,1),nu=c(0,0,1,0,1))
data = data.frame(x=c(1,1),beta=c(0,1),gamma=c(0,1),nu=c(0,1))
for(k in 1:dim(data)[1])
{
  for(i in 1:2)
  {
    mydata = rbind(mydata, data.frame(n.chains = 3, tune.type=i, n.thin = 15, data[k,]))
  }
}
rownames(mydata) = 1:dim(mydata)[1]

require(plyr)
m_ply(.data = mydata, .fun = sir_mcmc_plots)
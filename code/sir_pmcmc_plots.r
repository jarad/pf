require(pomp)

# Set data and graphics path
dpath = "../data/"
gpath = "../graphs/"

# Load simulated data and redefine data / known parameter values
load(paste(dpath,"sim-orig.rdata",sep=""))
y = mysims[[1]]$sim$y
P = mysims[[1]]$true.params$P
b = mysims[[1]]$true.params$b
varsigma = mysims[[1]]$true.params$varsigma
sigma = mysims[[1]]$true.params$sigma
eta = mysims[[1]]$true.params$eta

# Function to create traceplots and tables with 95% credible intervals
sir_pmcmc_plots <- function(n.chains, niter, np, y.max, nburn = 0, nthin = 1)
{
  # Load pmcmc objects
  this.out = list(); length(this.out) = n.chains
  for(i in 1:n.chains)
  {
    load(paste(dpath,"sir_pmcmc_test-",paste(i,niter,np,y.max,sep="-"),".rdata",sep=""))
    this.out[[i]] = out
    rm(out)
  }
  out = this.out
  rm(this.out)
  
  # Create traceplots for fixed parameters
  iter = seq(nburn+nthin+1,niter+1,nthin)
  xlabs = c("","","Iteration")
  ylabs = expression(beta,gamma,nu)
  pdf(paste(gpath,"sir_pmcmc_test-",paste(i,niter,np,y.max,sep="-"),"-traceplots.pdf",sep=""))
  par(mfrow=c(3,1),mar=c(5,7,4,2)+0.1)
  for(i in 1:3)
  {
    ymin = min(sapply(out, function(x) conv.rec(x)[iter,3+i]))
    ymax = max(sapply(out, function(x) conv.rec(x)[iter,3+i]))
    plot(iter,conv.rec(out[[1]])[iter,3+i],type = "l",ylim=c(ymin,ymax),xlab=xlabs[i],ylab=ylabs[i],cex.lab=1.75)
    if(n.chains > 1) 
    {
      for(j in 2:n.chains) lines(iter, conv.rec(out[[j]])[iter,3+i], col = 2*(j-1))
    }
  }
  dev.off()
  
  # Plot filtered means of states
  pdf(paste(gpath,"sir_pmcmc_test-",paste(i,niter,np,y.max,sep="-"),"-meanFilteredStates.pdf",sep=""))
  par(mar=c(5,7,4,2)+0.1)
  for(i in 1:3)
  {
    filt = filter.mean(out[[i]])
    if(i == 1) plot(1:dim(filt)[2],filt[1,],type="l",ylim=c(0,1),xlab="Time",ylab="Filtered Mean",main="Filtered Means of s and i",cex.lab=1.75) else lines(1:dim(filt)[2],filt[1,], col = 2*(i-1))
    if(i == 1) lines(1:dim(filt)[2],filt[2,]) else lines(1:dim(filt)[2],filt[2,],col=2*(i-1))
  }
  legend("topright",legend=1:n.chains,lty=c(1,1,1),col=c(1,2*(1:(n.chains-1))),cex=1.75)
  dev.off()
  
  # Calculate 95% credible intervals
  cred.int = matrix(NA, nr=2, nc = 3)
  colnames(cred.int) = c("beta","gamma","nu")
  for(i in 1:3) cred.int[,i] = quantile(sapply(out, function(x) conv.rec(x)[iter,3+i]),c(0.025,0.975))
  print(cred.int)
}

# Process pmcmc objects
require(plyr)
mydata = expand.grid(n.chains = 3, niter = 20000, np = 100, nburn = 0, nthin = 1, y.max = c(30,60,90,125), stringsAsFactors = FALSE)
m_ply(mydata, sir_pmcmc_plots)

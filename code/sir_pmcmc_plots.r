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
  rownames(cred.int) = c("2.5%","97.5%")
  for(i in 1:3) cred.int[,i] = quantile(sapply(out, function(x) conv.rec(x)[iter,3+i]),c(0.025,0.975))
  return(cred.int)
}

# Process pmcmc objects
require(plyr)
y.max = c(30,seq(60,100,5),125)
mydata = expand.grid(n.chains = 3, niter = 20000, np = 100, nburn = 0, nthin = 1, y.max = y.max, stringsAsFactors = FALSE)
cred.int.pmcmc = maply(mydata, sir_pmcmc_plots)

# Load KDPF credible intervals
load(paste(dpath,"PF-quant-1-20000-KD-stratified-orig-log-0.99-61.rdata",sep=""))
kd_quant <- function(y.max)
{
  cred.int = t(pf.quant.out$theta.quant[y.max+1,,4:5])
  rownames(cred.int) = c("2.5%","97.5%")
  colnames(cred.int) = c("beta","gamma","nu")
  return(cred.int)
}
cred.int.kdpf = maply(.data = data.frame(y.max = y.max), .fun = kd_quant)

# Plot 95% CI for pmcmc versus kdpf
params=expression(beta,gamma,nu)
xlabs = c("Time","","")
ylabs = c("Parameter Value","","")
nt = dim(mysims[[1]]$sim$y)[2]
pdf(file=paste(gpath,"sir-pmcmc-kdpf.pdf",sep=""),width=30,height=10)
par(mfrow=c(1,3),mar=c(9,11,7,1)+.1,mgp=c(7,2,0))
for(i in 1:3)
{
  ymin = min(cred.int.pmcmc[,1,i],cred.int.kdpf[,1,i],mysims[[1]]$true.params$theta[i])
  ymax = max(cred.int.pmcmc[,2,i],cred.int.kdpf[,2,i],mysims[[1]]$true.params$theta[i])
  plot(y.max,cred.int.pmcmc[,1,i],ylim=c(ymin,ymax),type="b",lwd=4,col=6,xlab=xlabs[i],ylab=ylabs[i],main=params[i],cex.lab=6,cex.main=7,cex.axis=4,cex=4)
  lines(y.max,cred.int.pmcmc[,2,i],type="b",col=6,lwd=4,cex=4)
  lines(y.max,cred.int.kdpf[,1,i],type="b",col=3,lwd=4,cex=4)
  lines(y.max,cred.int.kdpf[,2,i],type="b",col=3,lwd=4,cex=4)
  abline(h=mysims[[1]]$true.params$theta[i],lwd=6,col="gray47")
  if(i == 1) legend("bottomright",legend=c("KDPF","PMCMC"),lty=c(1,1),pch=c(1,1),col=c(3,6),cex=4,pt.cex=c(4,4),lwd=c(4,4),bg="white",pt.bg="white")
}
dev.off()
require(pomp)
require(mcmcse)
require(smcUtils)
require(coda)

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
sir_pmcmc_plots <- function(chains, niter, np, y.max, nburn = 0, nthin = 1)
{
  # Load pmcmc objects
  this.out = list(); 
  length(this.out) = length(chains)
  for(i in 1:length(chains))
  {
    load(paste(dpath,"sir_pmcmc_test-",paste(chains[i],niter,np,y.max,sep="-"),".rdata",sep=""))
    this.out[[i]] = out
    rm(out)
  }
  out = this.out
  rm(this.out)
  
  # Create traceplots for fixed parameters and calculate effective sample size
  iter = seq(nburn+nthin+1,niter+1,nthin)
  xlabs = c("","","Iteration")
  ylabs = expression(beta,gamma,nu)
  ess.mcmcse = rep(NA,3)
  ess.coda = rep(NA,3)
  pdf(paste(gpath,"sir_pmcmc_test-",paste(chains,collapse="-"),"-",paste(length(chains),niter,np,y.max,sep="-"),"-traceplots.pdf",sep=""))
  par(mfrow=c(3,1),mar=c(5,7,4,2)+0.1)
  for(i in 1:3)
  {
    ymin = min(sapply(out, function(x) conv.rec(x)[iter,3+i]))
    ymax = max(sapply(out, function(x) conv.rec(x)[iter,3+i]))
    param = conv.rec(out[[1]])[iter,3+i]
    param.all = param
    plot(iter,param,type = "l",ylim=c(ymin,ymax),xlab=xlabs[i],ylab=ylabs[i],cex.lab=1.75)
    if(length(chains) > 1) 
    {
      for(j in 2:length(chains))
      {
        param = conv.rec(out[[j]])[iter,3+i]
        param.all = c(param.all,param)
        lines(iter, param, col = 2*(j-1))
      }
    }
    ess.mcmcse[i] = ess(param.all)
    ess.coda[i] = effectiveSize(param.all)
    mtext(paste("ESS:",round(ess.coda[i],2)),side=3,cex=0.85)
    if(i == 1) title("PMCMC",cex.main=1.75)
  }
  dev.off()
  ess.all = rbind(ess.mcmcse,ess.coda)
  rownames(ess.all) = c("mcmcse","coda")
  colnames(ess.all) = c("beta","gamma","nu")
  
  # Plot filtered means of states
  pdf(paste(gpath,"sir_pmcmc_test-",paste(chains,collapse="-"),"-",paste(length(chains),niter,np,y.max,sep="-"),"-meanFilteredStates.pdf",sep=""))
  par(mar=c(5,7,4,2)+0.1)
  for(i in 1:length(chains))
  {
    filt = filter.mean(out[[i]])
    if(i == 1) plot(1:dim(filt)[2],filt[1,],type="l",ylim=c(0,1),xlab="Time",ylab="Filtered Mean",main="Filtered Means of s and i",cex.lab=1.75) else lines(1:dim(filt)[2],filt[1,], col = 2*(i-1))
    if(i == 1) lines(1:dim(filt)[2],filt[2,]) else lines(1:dim(filt)[2],filt[2,],col=2*(i-1))
  }
  legend("topright",legend=chains,lty=rep(1,length(chains)),col=c(1,2*(1:(length(chains)-1))),cex=1.75)
  dev.off()
  
  # Calculate 95% credible intervals
  cred.int = matrix(NA, nr=2, nc = 3)
  colnames(cred.int) = c("beta","gamma","nu")
  rownames(cred.int) = c("2.5%","97.5%")
  for(i in 1:3) cred.int[,i] = quantile(sapply(out, function(x) conv.rec(x)[iter,3+i]),c(0.025,0.975))
  return(rbind(cred.int,ess.all))
}

# Process pmcmc objects
require(plyr)
y.max.pmcmc = seq(5,125,5)
n.iter = 30000
mydata = expand.grid(chains = 1:3, niter = n.iter, np = 100, y.max = y.max.pmcmc, nburn = 0, nthin = 1, stringsAsFactors = FALSE)
cred.int.pmcmc = maply(mydata, sir_pmcmc_plots)
n.pmcmc = dim(cred.int.pmcmc)[1]
# mydata = expand.grid(niter = n.iter, np = 100, y.max = y.max.pmcmc, nburn = 0, nthin = 1, stringsAsFactors = FALSE)
# cred.int.pmcmc = maply(mydata, function(niter, np, y.max, nburn, nthin) sir_pmcmc_plots(1:3, niter, np, y.max, nburn, nthin))
# cred.int.pmcmc = array(cred.int.pmcmc, c(1,dim(cred.int.pmcmc)))
# n.pmcmc = 1

# Calculate ESS for KDPF runs
y.max = 0:125
n.sims = 20
np = 20000
ess_kdpf <- function(y.max,n.sim,np)
{
  load(paste(dpath,"PF-1-",np,"-KD-stratified-orig-log-0.99-",60+n.sim,".rdata",sep=""))
  ess.kdpf = rep(NA,length(y.max))
  for(j in 1:length(y.max)) ess.kdpf[j] = ess.weights(pf.out$out$weight[,y.max[j]+1])
  return(ess.kdpf)
}
require(plyr)
kdpf.ess <- maply(data.frame(n.sim=1:n.sims,np=np), function(n.sim,np) ess_kdpf(y.max,n.sim,np))

# Plot ESS for KDPF over time
pdf(file=paste(gpath,"sir-kdpf-ess-",n.sims,"-",np,".pdf",sep=""),width=10,height=10)
par(mar = c(5,7,4,2)+0.1)
ymax = max(kdpf.ess); ymin = min(kdpf.ess)
plot(y.max,kdpf.ess[1,],ylim=c(ymin,ymax),xlim=c(0,max(y.max,y.max.pmcmc)),type="l",main="KDPF ESS for all parameters",xlab="Time",ylab="ESS",cex.lab=1.5,cex.main=2,cex.axis=1.25)
if(n.sims > 1) for(i in 2:n.sims) lines(y.max,kdpf.ess[i,],col=i)
dev.off()

# Plot ESS for PMCMC over time
pdf(file=paste(gpath,"sir-pmcmc-ess-",n.pmcmc,"-",n.iter,".pdf",sep=""),width=9,height=3)
par(mfrow=c(1,3),mar=c(5,7,4,2)+0.1)
ymax = apply(cred.int.pmcmc,4,max)
ymin = apply(cred.int.pmcmc,4,min)
params = expression(beta,gamma,nu)
xlabs = c("Time","","")
ylabs = c("ESS","","")
for(i in 1:3)
{
  plot(y.max.pmcmc,cred.int.pmcmc[1,,3,i],ylim=c(ymin[i],ymax[i]),col=4,xlab=xlabs[i],ylab=ylabs[i],main=params[i],cex.lab=1.75,cex.main=2,cex.axis=1.25)
  if(n.pmcmc > 1) for(j in 2:n.pmcmc) points(y.max.pmcmc,cred.int.pmcmc[j,,3,i],col=4)
  for(j in 1:dim(cred.int.pmcmc)[1]) points(y.max.pmcmc,cred.int.pmcmc[j,,4,i],col=2)
  if(i == 1) legend("topleft",c("mcmcse","coda"),pch=c(1,1),col=c(4,2))
}
dev.off()

## Create plot comparing KDPF credible intervals with PMCMC
kd_quant <- function(y.max)
{
  cred.int = t(pf.quant.out$theta.quant[y.max+1,,4:5])
  rownames(cred.int) = c("2.5%","97.5%")
  colnames(cred.int) = c("beta","gamma","nu")
  return(cred.int)
}

# Find maximum and minimum kdpf quantiles
ymax.kdpf = rep(-Inf,3)
ymin.kdpf = rep(Inf,3)
for(i in 1:n.sims)
{
  # Load KDPF credible intervals
  load(paste(dpath,"PF-quant-1-",np,"-KD-stratified-orig-log-0.99-",60+i,".rdata",sep=""))
  cred.int.kdpf = maply(.data = data.frame(y.max = y.max), .fun = kd_quant)
  ymins = apply(cred.int.kdpf[,1,],2,min)
  ymaxs = apply(cred.int.kdpf[,2,],2,max)
  ymax.kdpf = c(max(ymax.kdpf[1],ymaxs[1]),max(ymax.kdpf[2],ymaxs[2]),max(ymax.kdpf[3],ymaxs[3]))
  ymin.kdpf = c(min(ymin.kdpf[1],ymins[1]),min(ymin.kdpf[2],ymins[2]),min(ymin.kdpf[3],ymins[3]))
}

# Plot 95% CI for pmcmc versus kdpf
params=expression(beta,gamma,nu)
xlabs = c("Time","","")
ylabs = c("Parameter Value","","")
nt = dim(mysims[[1]]$sim$y)[2]
pdf(file=paste(gpath,"sir-pmcmc-kdpf-",n.pmcmc,"-",n.iter,"-",n.sims,"-",np,".pdf",sep=""),width=30,height=10)
par(mfrow=c(1,3),mar=c(9,11,7,1)+.1,mgp=c(7,2,0))
for(i in 1:3)
{
  ymin = min(cred.int.pmcmc[,,1,i],ymin.kdpf[i],mysims[[1]]$true.params$theta[i])
  ymax = max(cred.int.pmcmc[,,2,i],ymax.kdpf[i],mysims[[1]]$true.params$theta[i])
  load(paste(dpath,"PF-quant-1-",np,"-KD-stratified-orig-log-0.99-",61,".rdata",sep=""))
  cred.int.kdpf = maply(.data = data.frame(y.max = y.max), .fun = kd_quant)
  plot(y.max,cred.int.kdpf[,1,i],type="l",ylim=c(ymin,ymax),lwd=2,col=3,xlab=xlabs[i],ylab=ylabs[i],main=params[i],cex.lab=6,cex.main=7,cex.axis=4,cex=2)
  lines(y.max,cred.int.kdpf[,2,i],lwd=2,col=3,cex=2)
  abline(h=mysims[[1]]$true.params$theta[i],lwd=6,col="gray47")
  if(n.sims > 1)
  {
    for(j in 1:n.sims)
    {
      load(paste(dpath,"PF-quant-1-",np,"-KD-stratified-orig-log-0.99-",60+j,".rdata",sep=""))
      cred.int.kdpf = maply(.data = data.frame(y.max = y.max), .fun = kd_quant)
      lines(y.max,cred.int.kdpf[,1,i],col=3,lwd=2,cex=2)
      lines(y.max,cred.int.kdpf[,2,i],col=3,lwd=2,cex=2)
    }
  }
  for(k in 1:n.pmcmc)
  {
    points(y.max.pmcmc,cred.int.pmcmc[k,,1,i],lwd=4,cex=4)
    points(y.max.pmcmc,cred.int.pmcmc[k,,2,i],lwd=4,cex=4)
    if(i == 1) legend("bottomright",legend=c("KDPF","PMCMC"),lty=c(1,NA),pch=c(NA,1),col=c(3,1),cex=4,pt.cex=c(NA,4),lwd=c(2,4),bg="white",pt.bg="white")
  }
}
dev.off()
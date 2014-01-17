source("pf_functions.r")

# Set data and graphics path
dpath = "../data/"
gpath = "../graphs/"

## Compare particle filters over different # particles using systematic resampling, uniform priors
n = c(100, 1000, 10000, 20000)
filt = c("BF", "APF", "KD")
cols = c(2, 4, 3)
probs = c(4, 5)
n.sim = 1:20
load.label <- function(filt, n, n.sim) paste(dpath,"PF-quant-",n.sim,"-",n,"-",filt,"-systematic-unif-logit-0.99-61.rdata",sep="")
states = c(TRUE, FALSE)
states.label <- c("states", "params")
for(i in n.sim)
{
  for(j in 1:2)
  {
    if(j == 1)
    {
      params <- expression(s,i,r)
      burn = 0
    } else {
      params <- expression(beta,gamma,nu)
      burn = c(15, 0, 0)
    }
    create.label <- paste(gpath,"PF-",i,"-systematic-unif-logit-0.99-61-",states.label[j],".pdf",sep="")
    pf_plot(n, params, filt, i, probs, cols, create.label, load.label, states[j], burn = c(15, 0, 0))
  }
}

## Compare resampling schemes over different # particles using KD pf and lognormal priors
n = c(100, 1000, 10000, 20000)
filt = c("multinomial", "residual", "stratified", "systematic")
cols = rainbow(length(filt))
probs = c(4, 5)
n.sim = 1:20
load.label <- function(filt, n, n.sim) paste(dpath,"PF-quant-",n.sim,"-",n,"-KD-",filt,"-orig-log-0.99-61.rdata",sep="")
states = c(TRUE, FALSE)
states.label <- c("states", "params")
for(i in n.sim)
{
  for(j in 1:2)
  {
    if(j == 1)
    {
      params <- expression(s,i,r)
      burn = 0
      create.label <- paste(gpath,"PF-",i,"-KD-resamp-orig-log-0.99-61-",states.label[j],".pdf",sep="")
      pf_plot(n, params, filt, i, probs, cols, create.label, load.label, states[j], burn = c(15, 0, 0))
    } else {
      params <- expression(beta,gamma,nu)
      burn = c(15, 0, 0)
      
      pf_average <- function(n.param)
      {
        load(paste(dpath,"sim-orig.rdata",sep=""))
        tt = dim(mysims[[i]]$sim$x)[2]
        avg.quant = matrix(0, nr=tt, nc = 2)
        for(k in 1:length(filt))
        {
          load(load.label(filt[k], n[length(n)], i))
          avg.quant[,1] = avg.quant[,1] + pf.quant.out$theta.quant[,n.param,probs[1]]
          avg.quant[,2] = avg.quant[,2] + pf.quant.out$theta.quant[,n.param,probs[2]]
        }
        avg.quant = avg.quant / length(filt)
        return(avg.quant)
      }  
      
      require(plyr)
      mydata = expand.grid(n.param = 1:length(params))
      out.avg = maply(mydata, pf_average)
      create.label <- paste(gpath,"PF-",i,"-KD-resamp-orig-log-0.99-61-",states.label[j],".pdf",sep="")
      pf_plot(n, params, filt, i, probs, cols, create.label, load.label, states[j], out.avg = out.avg, burn = c(15, 0, 0))
    }
  }
}

## Compare delta parameter over different # particles using lognormal priors and stratified resampling
n = c(100, 1000, 10000, 20000)
filt = c(0.9, 0.95, 0.96, 0.97, 0.98, 0.99)
cols = rainbow(length(filt))
probs = c(4, 5)
n.sim = 1:20
load.label <- function(filt, n, n.sim) paste(dpath,"PF-quant-",n.sim,"-",n,"-KD-stratified-orig-log-",filt,"-61.rdata",sep="")
states = c(TRUE, FALSE)
states.label <- c("states", "params")
for(i in n.sim)
{
  for(j in 1:2)
  {
    if(j == 1)
    {
      params <- expression(s,i,r)
      burn = 0
      create.label <- paste(gpath,"PF-",i,"-KD-stratified-orig-log-delta-61-",states.label[j],".pdf",sep="")
      pf_plot(n, params, filt, i, probs, cols, create.label, load.label, states[j], burn = c(15, 0, 0))
    } else {
      params <- expression(beta,gamma,nu)
      burn = c(15, 0, 0)
      
      pf_average <- function(n.param)
      {
        load(paste(dpath,"sim-orig.rdata",sep=""))
        tt = dim(mysims[[i]]$sim$x)[2]
        avg.quant = matrix(0, nr=tt, nc = 2)
        for(k in 1:length(filt))
        {
          load(load.label(filt[k], n[length(n)], i))
          avg.quant[,1] = avg.quant[,1] + pf.quant.out$theta.quant[,n.param,probs[1]]
          avg.quant[,2] = avg.quant[,2] + pf.quant.out$theta.quant[,n.param,probs[2]]
        }
        avg.quant = avg.quant / length(filt)
        return(avg.quant)
      }  
      
      require(plyr)
      mydata = expand.grid(n.param = 1:length(params))
      out.avg = maply(mydata, pf_average)
      
      create.label <- paste(gpath,"PF-",i,"-KD-stratified-orig-log-delta-61-",states.label[j],".pdf",sep="")
      pf_plot(n, params, filt, i, probs, cols, create.label, load.label, states[j], out.avg = out.avg, burn = c(15, 0, 0))      
    }
  }
}

## Compare particle filters over different # particles using stratified resampling, original priors
n = c(100, 1000, 10000, 20000)
filt = c("KD", "RM")
cols = c(2, 4, 3, 6)
probs = c(4, 5)
n.sim = 1
load.label <- function(filt, n, n.sim)
{
  if(filt != "RM") return(paste(dpath,"PF-quant-",n.sim,"-",n,"-",filt,"-stratified-orig-log-0.99-61.rdata",sep=""))
  if(filt == "RM") return(paste(dpath,"PF-quant-",n.sim,"-",n,"-",filt,"-stratified-orig-none-0.99-61.rdata",sep=""))
}
states = c(TRUE, FALSE)
states.label <- c("states", "params")
for(i in n.sim)
{
  for(j in 1:2)
  {
    if(j == 1)
    {
      params <- expression(s,i,r)
      burn = 0
    } else {
      params <- expression(beta,gamma,nu)
      burn = c(15, 0, 0)
    }
    create.label <- paste(gpath,"PF-",n.sim,"-stratified-orig-0.99-61-",states.label[j],".pdf",sep="")
    pf_plot(n, params, filt, i, probs, cols, create.label, load.label, states[j], burn = c(15, 0, 0))
  }
}

# Plot coverage probabilities for different particle filters with uniform priors, systematic resampling
quantiles <- c(0.5, 0.25, 0.75, 0.025, 0.975, 0.05, 0.95)
probs <- c(2, 3)
my_pf_coverage <- function(n, filt, states)
{
  load.label <- function(filt, n, n.sim)
  {
    if(filt != "RM") return(paste(dpath,"PF-quant-",n.sim,"-",n,"-",filt,"-systematic-unif-logit-0.99-61.rdata",sep=""))
    if(filt == "RM") return(paste(dpath,"PF-quant-",n.sim,"-",n,"-",filt,"-systematic-unif-none-0.99-61.rdata",sep=""))
  }
  pf_coverage(20, n, filt, probs, load.label, states) 
}
mydata = expand.grid(n = c(100, 1000, 10000, 20000), filt = c("BF", "APF", "KD"), states = c(TRUE, FALSE), stringsAsFactors = FALSE)
require(plyr)
coverage <- maply(mydata, my_pf_coverage)
alpha = quantiles[probs[2]]-quantiles[probs[1]]
params = expression(beta, gamma, nu)
states = expression(s, i, r)
cols = c(2,4,3)
create.label <- paste(gpath,"PF-coverage-",alpha,"-systematic-unif-filt-params.pdf",sep="")
pf_coverage_plot(coverage[,,1,,], alpha, 20, params, cols, create.label)
create.label <- paste(gpath,"PF-coverage-",alpha,"-systematic-unif-filt-states.pdf",sep="")
pf_coverage_plot(coverage[,,2,,], alpha, 20, states, cols, create.label)

# Plot coverage probabilities for resampling schemes with KD pf, original priors
quantiles <- c(0.5, 0.25, 0.75, 0.025, 0.975, 0.05, 0.95)
probs <- c(2, 3)
my_pf_coverage <- function(n, filt, states)
{
  load.label <- function(filt, n, n.sim) paste(dpath,"PF-quant-",n.sim,"-",n,"-KD-",filt,"-orig-log-0.99-61.rdata",sep="")
  pf_coverage(20, n, filt, probs, load.label, states) 
}
mydata = expand.grid(n = c(100, 1000, 10000, 20000), filt = c("multinomial", "residual", "stratified", "systematic"), states = c(TRUE, FALSE), stringsAsFactors = FALSE)
require(plyr)
coverage <- maply(mydata, my_pf_coverage)
alpha = quantiles[probs[2]]-quantiles[probs[1]]
params = expression(beta, gamma, nu)
states = expression(s, i, r)
cols = rainbow(4)
create.label <- paste(gpath,"PF-coverage-",alpha,"-KD-orig-resamp-params.pdf",sep="")
pf_coverage_plot(coverage[,,1,,], alpha, 20, params, cols, create.label)
create.label <- paste(gpath,"PF-coverage-",alpha,"-KD-orig-resamp-states.pdf",sep="")
pf_coverage_plot(coverage[,,2,,], alpha, 20, states, cols, create.label)

# Plot coverage probabilities for delta values with lognormal priors, stratified resampling (KD filter only)
quantiles <- c(0.5, 0.25, 0.75, 0.025, 0.975, 0.05, 0.95)
probs <- c(2, 3)
my_pf_coverage <- function(n, filt, states)
{
  load.label <- function(filt, n, n.sim) paste(dpath,"PF-quant-",n.sim,"-",n,"-KD-stratified-orig-log-",filt,"-61.rdata",sep="")
  pf_coverage(20, n, filt, probs, load.label, states) 
}
mydata = expand.grid(n = c(100, 1000, 10000, 20000), filt = c(0.9,0.95,0.96,0.97,0.98,0.99), states = c(TRUE, FALSE), stringsAsFactors = FALSE)
require(plyr)
coverage <- maply(mydata, my_pf_coverage)
alpha = quantiles[probs[2]]-quantiles[probs[1]]
params = expression(beta, gamma, nu)
states = expression(s, i, r)
cols = rainbow(6)
create.label <- paste(gpath,"PF-coverage-",alpha,"-KD-stratified-orig-delta-params.pdf",sep="")
pf_coverage_plot(coverage[,,1,,], alpha, 20, params, cols, create.label)
create.label <- paste(gpath,"PF-coverage-",alpha,"-KD-stratified-orig-delta-states.pdf",sep="")
pf_coverage_plot(coverage[,,2,,], alpha, 20, states, cols, create.label)

# Plot coverage probabilities for original versus disperse priors, stratified resampling, delta = .99 (KD pf)
quantiles <- c(0.5, 0.25, 0.75, 0.025, 0.975, 0.05, 0.95)
probs <- c(2, 3)
my_pf_coverage <- function(n, filt, states)
{
  load.label <- function(filt, n, n.sim) paste(dpath,"PF-quant-",n.sim,"-",n,"-KD-stratified-",filt,"-log-0.99-61.rdata",sep="")
  pf_coverage(20, n, filt, probs, load.label, states) 
}
mydata = expand.grid(n = c(100, 1000, 10000, 20000), filt = c("orig","disp"), states = c(TRUE, FALSE), stringsAsFactors = FALSE)
require(plyr)
coverage <- maply(mydata, my_pf_coverage)
alpha = quantiles[probs[2]]-quantiles[probs[1]]
params = expression(beta, gamma, nu)
states = expression(s, i, r)
cols = rainbow(6)
create.label <- paste(gpath,"PF-coverage-",alpha,"-KD-stratified-priors-params.pdf",sep="")
pf_coverage_plot(coverage[,,1,,], alpha, 20, params, cols, create.label)
create.label <- paste(gpath,"PF-coverage-",alpha,"-KD-stratified-priors-states.pdf",sep="")
pf_coverage_plot(coverage[,,2,,], alpha, 20, states, cols, create.label)

## Create scatterplots of beta v gamma over time
source("sir_functions.r")

# Set graphical parameters
xlim=c(.135,.33); ylim=c(.08,.16)
borderx=c(.14,.5); bordery=c(.09,.143)
msize=5; labsize=5; axsize=3; ptsize=3
ptsty=20; ptcol="gray75"; rline = -2.3; rsize = 1.8

# Load simulated data to get true values of beta and gamma
load(paste(dpath,"sim-orig.rdata",sep=""))

# Panel 1 - uniform prior draws, logit transformation
# Panel 2 - uniform prior draws, log transformation
trans = c("logit", "log")
for(k in 1:length(trans))
{
  # Load particle filtered data
  load(paste(dpath,"PF-1-10000-KD-systematic-unif-",trans[k],"-0.99-61.rdata",sep=""))
  
  # Resample particles at cutoff points to have equal weights
  cutoff = seq(16, 61, len=4)
  betas = pf.out$ftheta(pf.out$out$theta[1,,],1)
  gammas = pf.out$ftheta(pf.out$out$theta[2,,],2)
  myout = array(NA,dim=c(2,dim(betas)[1],dim(betas)[2]))
  myout[1,,] = betas
  myout[2,,] = gammas
  myscat = pf.scat(myout,pf.out$out$weight,cutoff)
  
  # Scatterplots over time
  file = paste("../graphs/PF-betaGammaScat-1-10000-KD-systematic-unif-",trans[k],"-0.99-61.pdf",sep="")
  pdf(file,width=20,height=5)
  par(mfrow=c(1,4),mar=c(7,10,5,1)+.1,mgp=c(6,1.55,0))
  for(i in 1:length(cutoff))
  {
    if(k == 1)
    {
      if(i == 1)
      {
        plot(myscat$xrw1[,i],myscat$xrw2[,i],col="white",xlim=xlim,ylim=ylim,xlab=expression(beta),ylab=expression(gamma),main=paste("t = ",cutoff[i]-1,sep=""),cex.main=msize,cex.lab=labsize,cex.axis=axsize)
      } else {
        plot(myscat$xrw1[,i],myscat$xrw2[,i],col="white",xlim=xlim,ylim=ylim,xlab="",ylab="",main=paste("t = ",cutoff[i]-1,sep=""),cex.main=msize,axes=FALSE)
        box()
      }
      abline(v=borderx,lty=2)
      abline(h=bordery,lty=2)
      points(myscat$xrw1[,i],myscat$xrw2[,i],col=ptcol,pch=ptsty,cex=ptsize)
      points(mysims[[1]]$true.params$theta[1],mysims[[1]]$true.params$theta[2],col=2,pch=3,lwd=5,cex=1.5*ptsize)
    } else {
      plot(myscat$xrw1[,i],myscat$xrw2[,i],col="white",xlim=xlim,ylim=ylim,xlab="",ylab="",axes=FALSE)
      box()
      abline(v=borderx,lty=2)
      abline(h=bordery,lty=2)
      points(myscat$xrw1[,i],myscat$xrw2[,i],col=ptcol,pch=ptsty,cex=ptsize)
      points(mysims[[1]]$true.params$theta[1],mysims[[1]]$true.params$theta[2],col=2,pch=3,lwd=5,cex=1.5*ptsize)
    }
  }
  dev.off()
}

## Figure 4 - Extended model: KD PF for n particles with original priors and stratified resampling

# Load data
n = 20000
n.sim = 1
filt = "KD"
resamp = "stratified"
prior = "orig"
transform = "log"
delta = 0.99
seed = 61
probs = c(4, 5)

load(paste(dpath,"sim-ext.rdata",sep=""))
load(paste(dpath,"PF-ext-",n.sim,"-",n,"-",filt,"-",resamp,"-",prior,"-",transform,"-",delta,"-",seed,".rdata",sep=""))

# Get points where resampling done for extended model
tt = dim(mysim$sim$x)[2]; nt = tt - 1
dpts = which(!is.na(mysim$sim$y[1,]))
resampled = rep(0,nt)
parents = pf.out$out$parent
for(i in 1:nt) resampled[i] = !all(parents[,i+1] == 1:n)
spts.ext = which(as.logical(resampled))
resampled = rep(0,nt)

# Get points where resampling done for original model
load(paste(dpath,"PF-ext.orig-",n.sim,"-",n,"-",filt,"-",resamp,"-",prior,"-",transform,"-",delta,"-",seed,".rdata",sep=""))
parents = pf.out$out$parent
for(i in 1:nt) resampled[i] = !all(parents[,i+1] == 1:n)
spts.org = which(as.logical(resampled))

# Set graphical parameters
params = expression(beta,gamma,nu,b,varsigma,sigma,eta)
ymins = c(0.10,0.075,0.85,0.05,.75,0.0004,1.99)
ymaxs = c(0.40,0.15,1.55,.45,1.2,.0016,2.01)
cex.lab = 6
cex.main = 7
cex.axis = 4
cex.leg = 4

# Construct plot
pdf(paste(gpath,"PF-ext-",n.sim,"-",n,"-",filt,"-",resamp,"-",prior,"-",transform,"-",delta,"-",seed,".pdf",sep=""),width=30,height=30)
par(mfrow=c(3,3),mar=c(9,11,7,1)+.1,mgp=c(7,2,0))
for(k in 1:length(params))
{
  load(paste(dpath,"PF-ext-quant-",n.sim,"-",n,"-",filt,"-",resamp,"-",prior,"-",transform,"-",delta,"-",seed,".rdata",sep=""))
  out = pf.quant.out$theta.quant
  tt = dim(out)[1]; nt = tt - 1
  if(k == 1) # label y axis, title, legend
  {
     plot(1:nt,out[-1,k,probs[1]],type="l",ylim=c(ymins[k],ymaxs[k]),col=4,xlab="",ylab=paste("J = ",n,sep=""),main=params[k],cex.lab=cex.lab,cex.main=cex.main,cex.axis=cex.axis)
     lines(1:nt,out[-1,k,probs[2]],col=4)
     points(spts.ext,rep(ymins[k],length(spts.ext)),pch="|",cex=2,col=4)
     points(spts.org,rep(ymins[k]+.03*(ymaxs[k]-ymins[k]),length(spts.org)),pch="|",cex=2,col=2)
     points(dpts,rep(ymins[k]+.06*(ymaxs[k]-ymins[k]),length(dpts)),pch="|",cex=2,col="gray47")
     legend("topright",legend=c("Truth","Original","Extended"),col=c("gray47",2,4),lty=c(1,1,1),cex=cex.leg)
  } else { # label title only
     plot(1:nt,out[-1,k,probs[1]],type="l",ylim=c(ymins[k],ymaxs[k]),col=4,xlab="",ylab="",main=params[k],cex.main=cex.main,cex.axis=cex.axis)
     lines(1:nt,out[-1,k,probs[2]],col=4)
     points(spts.ext,rep(ymins[k],length(spts.ext)),pch="|",cex=2,col=4)
     points(spts.org,rep(ymins[k]+.03*(ymaxs[k]-ymins[k]),length(spts.org)),pch="|",cex=2,col=2)
     points(dpts,rep(ymins[k]+.06*(ymaxs[k]-ymins[k]),length(dpts)),pch="|",cex=2,col="gray47")
  }
  load(paste(dpath,"sim-ext.rdata",sep=""))
  theta = c(mysim$true.params$theta,mysim$true.params$b,mysim$true.params$varsigma,mysim$true.params$sigma,mysim$true.params$eta)
  abline(h=theta[k],col="gray47")
  if(k %in% 1:3)
  {
    load(paste(dpath,"PF-ext.orig-quant-",n.sim,"-",n,"-",filt,"-",resamp,"-",prior,"-",transform,"-",delta,"-",seed,".rdata",sep=""))
    out = pf.quant.out$theta.quant
    tt = dim(out)[1]; nt = tt - 1
    lines(1:nt,out[-1,k,probs[1]],col=2)
    lines(1:nt,out[-1,k,probs[2]],col=2)
  }
}

params = expression(s,i,r)
ymin = 0; ymax = 1
for(k in 1:2)
{
  load(paste(dpath,"PF-ext-quant-",n.sim,"-",n,"-",filt,"-",resamp,"-",prior,"-",transform,"-",delta,"-",seed,".rdata",sep=""))
  out = pf.quant.out$state.quant
  tt = dim(out)[1]; nt = tt - 1
  if(k == 2) # label x axis, title
  {
     plot(1:nt,out[-1,k,probs[1]],type="l",ylim=c(ymin,1),col=4,xlab="Time (days)",ylab="",main=params[k],cex.lab=cex.lab,cex.main=cex.main,cex.axis=cex.axis)
     lines(1:nt,out[-1,k,probs[2]],col=4)
     lines(1:nt,mysim$sim$x[k,-1],col="gray47")
     points(spts.ext,rep(0,length(spts.ext)),pch="|",cex=2,col=4)
     points(spts.org,rep(.03,length(spts.org)),pch="|",cex=2,col=2)
     points(dpts,rep(.06,length(dpts)),pch="|",cex=2,col="gray47")
  } else { # label title only
     plot(1:nt,out[-1,k,probs[1]],type="l",ylim=c(0,ymax),col=4,xlab="Time (days)",ylab="",main=params[k],cex.lab=cex.lab,cex.main=cex.main,cex.axis=cex.axis)
     lines(1:nt,out[-1,k,probs[2]],col=4)
     lines(1:nt,mysim$sim$x[k,-1],col="gray47")
     points(spts.ext,rep(0,length(spts.ext)),pch="|",cex=2,col=4)
     points(spts.org,rep(.03,length(spts.org)),pch="|",cex=2,col=2)
     points(dpts,rep(.06,length(dpts)),pch="|",cex=2,col="gray47")
  }
  load(paste(dpath,"PF-ext.orig-quant-",n.sim,"-",n,"-",filt,"-",resamp,"-",prior,"-",transform,"-",delta,"-",seed,".rdata",sep=""))
  out = pf.quant.out$state.quant
  tt = dim(out)[1]; nt = tt - 1
  lines(1:nt,out[-1,k,probs[1]],col=2)
  lines(1:nt,out[-1,k,probs[2]],col=2)
}
#k = 3
#load(paste(dpath,"PF-quant-ext-ext-KD-lognormal-lognormal-stratified-",n,"-ess-80.rdata",sep=""))
#out = pf.quant.out$state.quant
#tt = dim(out)[1]; nt = tt - 1
#truex = 1 - apply(mysim$sim$x[,-1],2,sum)
#plot(1:nt,out[-1,k,probs[1]],type="l",ylim=c(0,ymax),col=4,xlab="Time (days)",ylab="",main=params[k],cex.lab=cex.lab,cex.main=cex.main,cex.axis=cex.axis)
#lines(1:nt,out[-1,probs[2],3],col=4)
#lines(1:nt,truex,col="gray47")
##points(spts.ext,rep(0,length(spts.ext)),pch="|",cex=2,col=4)
##points(spts.org,rep(.03,length(spts.org)),pch="|",cex=2,col=2)
##points(dpts,rep(.06,length(dpts)),pch="|",cex=2,col="gray47")
#load(paste(dpath,"PF-quant-ext-orig-KD-lognormal-lognormal-stratified-",n,"-ess-80.rdata",sep=""))
#out = pf.quant.out$state.quant
#tt = dim(out)[1]; nt = tt - 1
#lines(1:nt,out[-1,k,2],col=2)
#lines(1:nt,out[-1,k,3],col=2)
dev.off()

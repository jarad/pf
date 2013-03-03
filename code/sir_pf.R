# Set data path
dpath = "../data/"

# Load data
load(paste(dpath,"sim-xy.rdata",sep=""))

# Which parameters unknown?
p = 1:3
s = rep(0,3); s[p] = 1
sinv = rep(1,3); sinv[p] = 0

# pf - function to run particle filter given n number of particles, filt = "BF", "APF", or "KD" for which filter to run, resamp = "multinomial", "residual", "stratified", or "systematic" for which resampling method to use, and prior = "normal" or "uniform" for which prior to use on unknown parameters
# Returns nothing; saves .rdata data file
pf <- function(n, filt, resamp, prior, ...)
{
  # Create function to sample from prior distribution of theta and map theta to original scale
  if(prior == "uniform")
  {
    thetal = c(0.1400, 0.0900, 0.9500)
    thetau = c(0.5000, 0.1430, 1.3000)
    rtheta = function(){ u2theta(runif(3,thetal,thetau),thetal,thetau)}
    ftheta = function(theta,param=1) theta2u(theta,thetal[param],thetau[param])
  } else {
    theta.mean = c(-1.3296, -2.1764, 0.1055)
    theta.sd = sqrt(c(0.1055, 0.0140, 0.0064))
    rtheta = function(){ rnorm(3,theta.mean,theta.sd)}
    ftheta = function(theta,param=1) exp(theta)
  }

  # Get index of first non-empty observation
  empty = TRUE; ind = 1
  while(empty) if(all(is.na(sim$y[,ind]))) ind = ind + 1 else empty = FALSE 

  # Run one of the particle filters
  if(filt == "BF")
  {
    # Run bootstrap filter
    dllik_bf = function(y, x){ dllik(y, x, b, varsigma, sigma, dpower)}
    revo_bf = function(x){ revo(x, P, d, s*ftheta(x[p+2],p)+sinv*theta)}
    rprior_bf = function()
    { 
      myprior = rprior(sim$y[,ind], rtheta, b, varsigma, sigma, dpower)
      return(c(myprior$x,myprior$theta))
    }
    source("bf.r")
    out = bf(sim$y, dllik_bf, revo_bf, rprior_bf, n, progress=FALSE, method=resamp, log=F, ...)
  } else if(filt == "APF"){
    # Run auxiliary particle filter
    dllik_apf = function(y, x){ dllik(y, x, b, varsigma, sigma, dpower)}
    pstate_apf = function(x) { revo(x, P, d, s*ftheta(x[p+2],p)+sinv*theta, FALSE)}
    revo_apf = function(x){ revo(x, P, d, s*ftheta(x[p+2],p)+sinv*theta)}
    rprior_apf = function()
    { 
      myprior = rprior(sim$y[,ind], rtheta, b, varsigma, sigma, dpower)
      return(c(myprior$x,myprior$theta))
    }
    source("apf.r")
    out = apf(sim$y, dllik_apf, pstate_apf, revo_apf, rprior_apf, n, progress=FALSE, method=resamp, log=F, ...)
  } else {
    # Run kernel density particle filter
    dllik_kd = function(y, x, theta=NULL){ dllik(y, x, b, varsigma, sigma, dpower)}
    pstate_kd = function(x, mytheta) { revo(x, P, d, s*ftheta(mytheta,p)+sinv*theta, FALSE)}
    revo_kd = function(x, mytheta){ revo(x, P, d, s*ftheta(mytheta,p)+sinv*theta)}
    rprior_kd = function(){ rprior(sim$y[,ind], rtheta, b, varsigma, sigma, dpower)}
    source("kd_pf.r")
    out = kd_pf(sim$y, dllik_kd, pstate_kd, revo_kd, rprior_kd, n, progress=FALSE, method=resamp, log=F, ...)
  }

  # Save output
  pf.out = list(out=out,ftheta=ftheta,p=p)
  save(pf.out, file=paste(dpath,"PF-",filt,"-",prior,"-",resamp,"-",n,".rdata",sep=""))
}

# Apply pf to combination of pfs
require(plyr)
mydata = expand.grid(n = c(100, 1000, 10000), filt = c("BF","APF","KD"), resamp = c("multinomial","residual","stratified","systematic"), prior = c("normal","uniform"), nonuniformity="ess", threshold=0.8, stringsAsFactors=FALSE)
save(mydata, file=paste(dpath,"filt_instruc.rdata",sep=""))
m_ply(mydata,pf)

# Clear objects
rm(list=ls(all=TRUE))
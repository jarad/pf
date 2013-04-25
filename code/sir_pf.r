# Set data path
dpath = "../data/"

# Load data
load(paste(dpath,"sim-xy.rdata",sep=""))

# pf - function to run particle filter given n number of particles, filt = "BF", "APF", or "KD" for which filter to run, resamp = "multinomial", "residual", "stratified", or "systematic" for which resampling method to use, and prior = "normal" or "uniform" for which prior to use on unknown parameters
# Returns nothing; saves .rdata data file
pf <- function(n, filt, resamp, prior, nonunif = "ess", thresh = 0.8, progress)
{
  # Create function to sample from prior distribution of theta and map theta to original scale
  if(prior == "uniform")
  {
    thetal = c(0.1400, 0.0900, 0.9500)
    thetau = c(0.5000, 0.1430, 1.3000)
    rtheta = function(){ u2theta(runif(3,thetal,thetau),thetal,thetau)}
    ftheta = function(theta,param=1) theta2u(theta,thetal[param],thetau[param])
  } else if(prior == "semi-uniform") {
    thetal = c(0.1400, 0.0900, 0.9500)
    thetau = c(0.5000, 0.1430, 1.3000)
    rtheta = function(){ runif(3,thetal,thetau)}
    ftheta = function(theta,param=1) theta
  } else if(prior == "semi-normal") {
    theta.mean = c(-1.3296, -2.1764, 0.1055)
    theta.sd = sqrt(c(0.1055, 0.0140, 0.0064))
    rtheta = function(){ exp(rnorm(3,theta.mean,theta.sd))}
    ftheta = function(theta,param=1) theta
  } else {
    theta.mean = c(-1.3296, -2.1764, 0.1055)
    theta.sd = sqrt(c(0.1055, 0.0140, 0.0064))
    rtheta = function(){ rnorm(3,theta.mean,theta.sd)}
    ftheta = function(theta,param=1) exp(theta)
  }

  # Run one of the particle filters
  if(filt == "BF")
  {
    # Run bootstrap filter
    dllik_bf = function(y, x){ dllik(y, x, b, varsigma, sigma, dpower)}
    revo_bf = function(x){ revo(x, P, d, ftheta(x[3:5],1:3))}
    rprior_bf = function()
    { 
      myprior = rprior(rtheta)
      return(c(myprior$x,myprior$theta))
    }
    source("bf.r")
    out = bf(sim$y, dllik_bf, revo_bf, rprior_bf, n, progress=progress, method=resamp, log=F, nonuniformity = nonunif, threshold = thresh)
  } else if(filt == "APF"){
    # Run auxiliary particle filter
    dllik_apf = function(y, x){ dllik(y, x, b, varsigma, sigma, dpower)}
    pstate_apf = function(x) { revo(x, P, d, ftheta(x[3:5],1:3), FALSE)}
    revo_apf = function(x){ revo(x, P, d, ftheta(x[3:5],1:3))}
    rprior_apf = function()
    { 
      myprior = rprior(rtheta)
      return(c(myprior$x,myprior$theta))
    }
    source("apf.r")
    out = apf(sim$y, dllik_apf, pstate_apf, revo_apf, rprior_apf, n, progress=progress, method=resamp, log=F, nonuniformity = nonunif, threshold = thresh)
  } else {
    # Run kernel density particle filter
    dllik_kd = function(y, x, theta=NULL){ dllik(y, x, b, varsigma, sigma, dpower)}
    pstate_kd = function(x, mytheta) { revo(x, P, d, ftheta(mytheta,1:3), FALSE)}
    revo_kd = function(x, mytheta){ revo(x, P, d, ftheta(mytheta,1:3))}
    rprior_kd = function(){ rprior(rtheta)}
    source("kd_pf.r")
    out = kd_pf(sim$y, dllik_kd, pstate_kd, revo_kd, rprior_kd, n, progress=progress, method=resamp, log=F, nonuniformity = nonunif, threshold = thresh)
  }

  # Save output
  pf.out = list(out=out,ftheta=ftheta)
  save(pf.out, file=paste(dpath,"PF-",filt,"-",prior,"-",resamp,"-",n,"-",nonunif,"-",100*thresh,".rdata",sep=""))
}

# Apply pf to combination of pfs
#require(plyr)
#mydata = expand.grid(n = c(100, 1000, 10000), filt = "KD", resamp = "systematic", prior = "semi-normal", progress=TRUE, stringsAsFactors=FALSE)
#m_ply(mydata,pf)

require(plyr)
mydata = expand.grid(n = 10000, filt = "KD", resamp = "stratified", prior = "normal", progress=FALSE, nonunif=c("ess","cov","entropy"), thresh=seq(.05,.95,.05), stringsAsFactors=FALSE)
m_ply(mydata,pf)

# Clear objects
rm(list=ls(all=TRUE))
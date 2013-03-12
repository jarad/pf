# Set data path
dpath = "../data/"

# Load data
load(paste(dpath,"sim-xy-ext.rdata",sep=""))

# pf - function to run particle filter given n number of particles, filt = "BF", "APF", or "KD" for which filter to run, resamp = "multinomial", "residual", "stratified", or "systematic" for which resampling method to use, and prior = "normal" or "uniform" for which prior to use on unknown parameters
# Returns nothing; saves .rdata data file
pf <- function(n, filt, resamp, prior, progress, ...)
{
  # Create function to sample from prior distribution of theta and map theta to original scale
  if(prior == "uniform")
  {
    thetal = c(0.1400, 0.0900, 0.9500, 0.0001, 0.85, 0.0001)
    thetau = c(0.5000, 0.1430, 1.3000, 0.5000, 1.15, 0.0100)
    rtheta = function(){ u2theta(runif(3,thetal,thetau),thetal[1:3],thetau[1:3])}
    rthetaPlus = function()
    {
      .current.seed = .Random.seed
      b = u2theta(runif(1,thetal[4],thetau[4]),thetal[4],thetau[4])
      .Random.seed = .current.seed
      varsigma = u2theta(runif(1,thetal[5],thetau[5]),thetal[5],thetau[5])
      .Random.seed = .current.seed
      sigma = u2theta(runif(1,thetal[6],thetau[6]),thetal[6],thetau[6])
      return(list(b=b,varsigma=varsigma,sigma=sigma))
    }
    ftheta = function(theta,param=1) theta2u(theta,thetal[param],thetau[param])
  } else {
    theta.mean = c(-1.3296, -2.1764, 0.1055, -4.9517, 0, -6.9078)
    theta.sd = sqrt(c(0.1055, 0.0140, 0.0064, 4.7209, 0.0059, 1.3801))
    rtheta = function(){ rnorm(3,theta.mean,theta.sd)}
    rthetaPlus = function()
    {
      .current.seed = .Random.seed
      b = rnorm(1,theta.mean[4],theta.sd[4])
      .Random.seed = .current.seed
      varsigma = rnorm(1,theta.mean[5],theta.sd[5])
      .Random.seed = .current.seed
      sigma = rnorm(1,theta.mean[6],theta.sd[6])
      return(list(b=b,varsigma=varsigma,sigma=sigma))
    }
    ftheta = function(theta,param=1) exp(theta)
  }

  # Run one of the particle filters
  if(filt == "BF")
  {
    # Run bootstrap filter
    dllik_bf = function(y, x){ dllik(y, x[1:2], ftheta(x[6],4), ftheta(x[7],5), ftheta(x[8],6), dpower)}
    revo_bf = function(x){ revo(x, P, d, ftheta(x[3:5],4:6))}
    rprior_bf = function()
    { 
      myprior = rprior(rtheta, obsparam=TRUE, rthetaPlus=rthetaPlus)
      return(c(myprior$x,myprior$theta))
    }
    source("bf.r")
    out = bf(sim$y, dllik_bf, revo_bf, rprior_bf, n, progress=progress, method=resamp, log=F, ...)
  } else if(filt == "APF"){
    # Run auxiliary particle filter
    dllik_apf = function(y, x){ dllik(y, x[1:2], ftheta(x[6],4), ftheta(x[7],5), ftheta(x[8],6), dpower)}
    pstate_apf = function(x) { revo(x, P, d, ftheta(x[3:5],4:6), FALSE)}
    revo_apf = function(x){ revo(x, P, d, ftheta(x[3:5],4:6))}
    rprior_apf = function()
    { 
      myprior = rprior(rtheta, obsparam=TRUE, rthetaPlus=rthetaPlus)
      return(c(myprior$x,myprior$theta))
    }
    source("apf.r")
    out = apf(sim$y, dllik_apf, pstate_apf, revo_apf, rprior_apf, n, progress=progress, method=resamp, log=F, ...)
  } else {
    # Run kernel density particle filter
    dllik_kd = function(y, x, mytheta){ dllik(y, x, ftheta(mytheta[4],4), ftheta(mytheta[5],5), ftheta(mytheta[6],6), dpower)}
    pstate_kd = function(x, mytheta) { revo(x, P, d, ftheta(mytheta[1:3],1:3), FALSE)}
    revo_kd = function(x, mytheta){ revo(x, P, d, ftheta(mytheta[1:3],1:3))}
    rprior_kd = function() rprior(rtheta, obsparam=TRUE, rthetaPlus=rthetaPlus)
    source("kd_pf.r")
    out = kd_pf(sim$y, dllik_kd, pstate_kd, revo_kd, rprior_kd, n, progress=progress, method=resamp, log=F, ...)
  }

  # Save output
  pf.out = list(out=out,ftheta=ftheta)
  save(pf.out, file=paste(dpath,"PF-ext-",filt,"-",prior,"-",resamp,"-",n,".rdata",sep=""))
}

# Apply pf to combination of pfs
require(plyr)
mydata = expand.grid(n = c(100, 1000, 10000, 20000), filt = c("BF","APF","KD"), resamp = c("multinomial","residual","stratified","systematic"), prior = "normal", progress=FALSE, nonuniformity="ess", threshold=0.8, stringsAsFactors=FALSE)
m_ply(mydata,pf)

# Clear objects
rm(list=ls(all=TRUE))
# Set data path
dpath = "../data/"

# Initialize .Random.seed
rnorm(1)

# prior_samp - function to sample prior draws for a given number of particles, choice of prior, and number of unknown parameters (3 or 6)
prior_samp <- function(n, prior, mod)
{
  # Sample from prior distribution of theta and create function to map theta to original scale
  if(mod == "orig")
  {
    # Load data
    load(paste(dpath,"sim-orig.rdata",sep=""))
    theta.samp = matrix(NA,nr=n,nc=3)
    if(prior == "uniform")
    {
      thetal = c(0.1400, 0.0900, 0.9500)
      thetau = c(0.5000, 0.1430, 1.3000)
      for(i in 1:n) theta.samp[i,] = runif(3,thetal,thetau)
    } else if(prior == "lognormal") {
      theta.mean = c(-1.3296, -2.1764, 0.1055)
      theta.sd = c(.3248, .1183, .0800)
      for(i in 1:n) theta.samp[i,] = exp(rnorm(3,theta.mean,theta.sd))
    } else { stop("prior must be either 'uniform' or 'lognormal'")}
  } else if(mod == "ext") {
    # Load data
    load(paste(dpath,"sim-ext.rdata",sep=""))
    theta.samp = matrix(NA,nr=n,nc=7)
    if(prior == "uniform")
    {
      thetal = c(0.1400, 0.0900, 0.9500, 0.0001, 0.85, 0.0001, 0)
      thetau = c(0.5000, 0.1430, 1.3000, 0.5000, 1.15, 0.0100, 5)
      for(i in 1:n) theta.samp[i,] = runif(7,thetal,thetau)
    } else if(prior == "lognormal") {
      theta.mean = c(-1.3296, -2.1764, 0.1055, -1.609, -0.0114, -7.0516, 2.5)
      theta.sd = c(.3248, .1183, .0800, .3536, .0771, .2803, 1)
      for(i in 1:n)
      {
        prior.draw = rnorm(7,theta.mean,theta.sd)
        theta.samp[i,] = c(exp(prior.draw[1:6]),prior.draw[7])
      }
    } else { stop("prior must be either 'uniform' or 'lognormal'")}
  } else { stop("mod must be either 'orig' or 'ext'")}

  # Save output
  save(theta.samp,file=paste(dpath,n,"-",prior,"-",mod,"-prior-samp.rdata",sep=""))
}

# Sample prior draws for every combination of n, prior, and mod
mydata = expand.grid(n=c(100,1000,10000,20000,40000),prior=c("uniform","lognormal"),mod=c("orig","ext"),stringsAsFactors=FALSE)
require(plyr)
m_ply(mydata,prior_samp)

# pf - function to run particle filter given n number of particles, filt = "BF", "APF", or "KD" for which filter to run, resamp = "multinomial", "residual", "stratified", or "systematic" for which resampling method to use, and prior = "normal" or "uniform" for which prior to use on unknown parameters
# Returns nothing; saves .rdata data file
pf <- function(n, filt, resamp, prior, prior.samp, mod.data, mod.fit, nonunif = "ess", thresh = 0.8, progress = TRUE)
{
  # Load simulated data
  load(paste(dpath,"sim-",mod.data,".rdata",sep=""))

  # Create functions to load prior sample for given n, prior, and mod (on transformed scale) and to map theta to original scale
  load(paste(dpath,n,"-",prior.samp,"-",mod.fit,"-prior-samp.rdata",sep=""))
  if(mod.fit == "orig")
  {
    if(prior == "uniform")
    { 
      thetal = c(0.1400, 0.0900, 0.9500)
      thetau = c(0.5000, 0.1430, 1.3000)
      rtheta = function(j) u2theta(theta.samp[j,],thetal,thetau)
      ftheta = function(theta,param=1) theta2u(theta,thetal[param],thetau[param])
    } else if(prior == "semi-uniform") {
      thetal = c(0.1400, 0.0900, 0.9500)
      thetau = c(0.5000, 0.1430, 1.3000)
      rtheta = function(j) theta.samp[j,]
      ftheta = function(theta,param=1) theta
    } else if(prior == "lognormal") {
      rtheta = function(j) log(theta.samp[j,])
      ftheta = function(theta,param=1) exp(theta)
    } else if(prior == "semi-lognormal") {
      rtheta = function(j) theta.samp[j,]
      ftheta = function(theta,param=1) theta
    } else { stop("prior must be either 'uniform', 'semi-uniform', 'lognormal', or 'semi-lognormal'")}
  } else if(mod.fit == "ext") {
    if(prior == "uniform")
    { 
      thetal = c(0.1400, 0.0900, 0.9500, 0.0001, 0.85, 0.0001, 0)
      thetau = c(0.5000, 0.1430, 1.3000, 0.5000, 1.15, 0.0100, 5)
      rtheta = function(j) u2theta(theta.samp[j,],thetal,thetau)
      ftheta = function(theta,param=1) theta2u(theta,thetal[param],thetau[param])
    } else if(prior == "semi-uniform") {
      thetal = c(0.1400, 0.0900, 0.9500, 0.0001, 0.85, 0.0001, 0)
      thetau = c(0.5000, 0.1430, 1.3000, 0.5000, 1.15, 0.0100, 5)
      rtheta = function(j) theta.samp[j,]
      ftheta = function(theta,param=1) theta
    } else if(prior == "lognormal") {
      rtheta = function(j) c(log(theta.samp[j,1:6]),theta.samp[j,7])
      ftheta = function(theta,param=1) if(param == 7) theta else exp(theta)
    } else if(prior == "semi-lognormal") {
      rtheta = function(j) theta.samp[j,]
      ftheta = function(theta,param=1) theta
    } else { stop("prior must be either 'uniform', 'semi-uniform', 'lognormal', or 'semi-lognormal'")}
  } else { stop("mod must be either 'orig' or 'ext'")}

  # Run one of the particle filters
  if(mod.fit == "orig")
  {
    if(filt == "BF")
    {
      # Run bootstrap filter
      dllik_bf = function(z, x){ dllik(z, x, b, varsigma, sigma, mu)}
      revo_bf = function(x){ revo(x, P, d, ftheta(x[3:5],1:3))}
      rprior_bf = function(j)
      { 
        myrtheta = function() rtheta(j)
        myprior = rprior(myrtheta)
        return(c(myprior$x,myprior$theta))
      }
      source("bf.r")
      out = bf(sim$y, dllik_bf, revo_bf, rprior_bf, n, progress=progress, method=resamp, log=F, nonuniformity = nonunif, threshold = thresh)
    } else if(filt == "APF"){
      # Run auxiliary particle filter
      dllik_apf = function(z, x){ dllik(z, x, b, varsigma, sigma, mu)}
      pstate_apf = function(x) { revo(x, P, d, ftheta(x[3:5],1:3), FALSE)}
      revo_apf = function(x){ revo(x, P, d, ftheta(x[3:5],1:3))}
      rprior_apf = function(j)
      { 
        myrtheta = function() rtheta(j)
        myprior = rprior(myrtheta)
        return(c(myprior$x,myprior$theta))
      }
      source("apf.r")
      out = apf(sim$y, dllik_apf, pstate_apf, revo_apf, rprior_apf, n, progress=progress, method=resamp, log=F, nonuniformity = nonunif, threshold = thresh)
    } else {
      # Run kernel density particle filter
      dllik_kd = function(z, x, mytheta=NULL){ dllik(z, x, b, varsigma, sigma, mu)}
      pstate_kd = function(x, mytheta) { revo(x, P, d, ftheta(mytheta,1:3), FALSE)}
      revo_kd = function(x, mytheta){ revo(x, P, d, ftheta(mytheta,1:3))}
      rprior_kd = function(j)
      { 
        myrtheta = function() rtheta(j)
        rprior(myrtheta)
      }
      source("kd_pf.r")
      out = kd_pf(sim$y, dllik_kd, pstate_kd, revo_kd, rprior_kd, n, progress=progress, method=resamp, log=F, nonuniformity = nonunif, threshold = thresh)
    }
  } else if(mod.fit == "ext") {
    if(filt == "BF")
    {
      # Run bootstrap filter
      dllik_bf = function(z, x){ dllik(z, x[1:2], ftheta(x[6],4), ftheta(x[7],5), ftheta(x[8],6), ftheta(x[9],7))}
      revo_bf = function(x){ revo(x, P, d, c(ftheta(x[3],4),ftheta(x[4],5),ftheta(x[5],6)))}
      rprior_bf = function(j)
      { 
        myrtheta = function() rtheta(j)
        myprior = rprior(myrtheta)
        return(c(myprior$x,myprior$theta))
      }
      source("bf.r")
      out = bf(sim$y, dllik_bf, revo_bf, rprior_bf, n, progress=progress, method=resamp, log=F, nonuniformity = nonunif, threshold = thresh)
    } else if(filt == "APF"){
      # Run auxiliary particle filter
      dllik_apf = function(z, x){ dllik(z, x[1:2], ftheta(x[6],4), ftheta(x[7],5), ftheta(x[8],6), ftheta(x[9],7))}
      pstate_apf = function(x) { revo(x, P, d, c(ftheta(x[3],4),ftheta(x[4],5),ftheta(x[5],6)), FALSE)}
      revo_apf = function(x){ revo(x, P, d, c(ftheta(x[3],4),ftheta(x[4],5),ftheta(x[5],6)))}
      rprior_apf = function(j)
      {
        myrtheta = function() rtheta(j) 
        myprior = rprior(myrtheta)
        return(c(myprior$x,myprior$theta))
      }
      source("apf.r")
      out = apf(sim$y, dllik_apf, pstate_apf, revo_apf, rprior_apf, n, progress=progress, method=resamp, log=F, nonuniformity = nonunif, threshold = thresh)
    } else {
      # Run kernel density particle filter
      dllik_kd = function(z, x, mytheta){ dllik(z, x, ftheta(mytheta[4],4), ftheta(mytheta[5],5), ftheta(mytheta[6],6), ftheta(mytheta[7],7))}
      pstate_kd = function(x, mytheta) { revo(x, P, d, c(ftheta(mytheta[1],1),ftheta(mytheta[2],2),ftheta(mytheta[3],3)), FALSE)}
      revo_kd = function(x, mytheta){ revo(x, P, d, c(ftheta(mytheta[1],1),ftheta(mytheta[2],2),ftheta(mytheta[3],3)))}
      rprior_kd = function(j)
      { 
        myrtheta = function() rtheta(j)
        rprior(myrtheta)
      }
      source("kd_pf.r")
      out = kd_pf(sim$y, dllik_kd, pstate_kd, revo_kd, rprior_kd, n, progress=progress, method=resamp, log=F, nonuniformity = nonunif, threshold = thresh)
    }
  } else { stop("mod must be either 'orig' or 'ext'")}

  # Save output
  pf.out = list(out=out,ftheta=ftheta)
  save(pf.out, file=paste(dpath,"PF-",mod.data,"-",mod.fit,"-",filt,"-",prior.samp,"-",prior,"-",resamp,"-",n,"-",nonunif,"-",100*thresh,".rdata",sep=""))
}

# Apply pf to combination of pfs
require(plyr)
data1 = expand.grid(mod.dat = "orig", mod.fit = "orig", n = c(100, 1000, 10000, 20000, 40000), filt = c("BF","APF","KD"), resamp = "systematic", prior = "uniform", prior.samp = "uniform", progress=FALSE, stringsAsFactors=FALSE)
data2 = data.frame(mod.dat = "orig", mod.fit = "orig", n = c(10000, 20000), filt = "KD", resamp = "systematic", prior = "semi-uniform", prior.samp = "uniform", progress=FALSE, stringsAsFactors=FALSE)
data3 = expand.grid(mod.dat = "orig", mod.fit = "orig", n = c(100, 1000, 10000, 20000, 40000), filt = "KD", resamp = c("multinomial","residual","stratified","systematic"), prior = "lognormal", prior.samp = "lognormal", progress=FALSE, stringsAsFactors=FALSE)
data4 = expand.grid(mod.dat = "ext", mod.fit = c("ext","orig"), n = c(100, 1000, 10000, 20000, 40000), filt = "KD", resamp = "stratified", prior = "lognormal", prior.samp = "lognormal", progress=FALSE, stringsAsFactors=FALSE)
mydata = data.frame(rbind(data1,data2,data3,data4))
m_ply(mydata,pf)

# Clear objects
rm(list=ls(all=TRUE))
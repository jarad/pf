source("sir_functions.r")

# Set data path
dpath = "/storage/sheinson_research/"

# pf - function to run particle filter given n number of particles, n.sim-th data set, resamp = "multinomial", "residual", "stratified", or "systematic", prior = "orig" or "disp", and delta amount of jitter to particles
# KD particle filter, lognormal priors on beta, gamma, and nu
# Returns nothing; saves .rdata data file
pf <- function(n, filt, resamp, prior, transform, delta, seed, progress = TRUE)
{
  # Load simulated data
  load(paste(dpath,"sim-ext.rdata",sep=""))
  y = mysim$sim$y
  P = mysim$true.params$P  
  
  # Define function to transform theta to original scale
  if(transform == "log") {
    ftheta = function(theta,param=1:7)
    {
      if(length(param) > 1)
      {
        out.theta = rep(NA, length(param))
        for(k in 1:length(param))
        {
          if(param[k] < 7) out.theta[k] = exp(theta[param[k]]) else out.theta[k] = theta[param[k]]
        }
      } else {
        if(param < 7) return(exp(theta)) else return(theta)
      }
      return(out.theta)
    }
  } else if(transform == "none"){
    ftheta = function(theta, param=NULL) theta 
  } else { stop("transform 'log', or 'none'")}
  
  # Define functions to sample prior draws of fixed parameters
  if(prior == "orig")
  {
    if(transform == "log")
    {
      rtheta <- function()
      {
        theta <- rep(NA, 7)
        log.params <- find.mu.sigma(c(1.5, .09, .95, .1, .85, .0005), c(3, .143, 1.3, .4, 1.15, .0015))
        theta[2:7] <- exp(rnorm(6, c(log.params[[1]][2:6], 2.5), c(log.params[[2]][2:6], 1)))
        theta[1] <- theta[2]*exp(rnorm(1, log.params[[1]][1], log.params[[2]][1]))
        return(log(theta))
      }
    } else if(transform == "none") { 
      rtheta <- function()
      {
        theta <- rep(NA, 7)
        log.params <- find.mu.sigma(c(1.5, .09, .95, .1, .85, .0005), c(3, .143, 1.3, .4, 1.15, .0015))
        theta[2:7] <- exp(rnorm(6, c(log.params[[1]][2:6], 2.5), c(log.params[[2]][2:6], 1)))
        theta[1] <- theta[2]*exp(rnorm(1, log.params[[1]][1], log.params[[2]][1]))
        return(theta)
      }
    } else { stop("Must use log or no transformation with original prior") }
  } else { stop("prior must be 'orig'") }
  
  # Set seed
  set.seed(seed)
  
  if(filt == "BF" | filt == "APF")
  {
    mydllik = function(y, x) dllik(y, x, ftheta(x[6], 4), ftheta(x[7], 5), ftheta(x[8], 6), ftheta(x[9], 7))
    myrevo = function(x) revo(x, ftheta(x[3:5], 1:3), P)
    myrprior = function()
    {
      myprior = rprior(rtheta)
      return(c(myprior$x,myprior$theta))
    }
    if(filt == "BF")
    {
      source("bf.r")
      out = bf(y, mydllik, myrevo, myrprior, n, progress=progress, method="systematic", nonuniformity = "ess", threshold = 0.8, log=F)
    } else {
      pstate = function(x) revo(x, ftheta(x[3:5], 1:3), P, FALSE)
      source("apf.r")
      out = apf(y, mydllik, pstate, myrevo, myrprior, n, progress=progress, method="systematic", nonuniformity = "ess", threshold = 0.8, log=F)
    }
  } else if(filt == "KD"){
    mydllik = function(y, x, theta) dllik(y, x, ftheta(theta[4], 4), ftheta(theta[5], 5), ftheta(theta[6], 6), ftheta(theta[7], 7))
    myrevo = function(x, theta) revo(x, ftheta(theta), P)
    pstate = function(x, theta) revo(x, ftheta(theta), P, FALSE)
    myrprior = function() rprior(rtheta)
    source("kd_pf.r")
    out = kd_pf(y, mydllik, pstate, myrevo, myrprior, n, delta, progress, method=resamp, nonuniformity = "ess", threshold = 0.8, log=F)
} else { stop("filt must be one of 'BF', 'APF', or 'KD'") }
  
  # Save output
  pf.out = list(out=out,ftheta=ftheta)
  file = paste(dpath,"PF-ext-1-",n,"-",filt,"-",resamp,"-",prior,"-",transform,"-",delta,"-",seed,".rdata",sep="")
  print(file)
  save(pf.out, file=file)
}

# Create data frame and use plyr to run particle filters in parallel
mydata = expand.grid(n = c(100, 1000, 10000, 20000), filt = c("KD"), resamp = "stratified", prior = "orig", transform = "log", delta = 0.99, seed = 61, progress = FALSE, stringsAsFactors=FALSE)

require(plyr)
require(doMC)
registerDoMC()
m_ply(.data = mydata, .fun = pf, .parallel = TRUE)
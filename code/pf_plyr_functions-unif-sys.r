source("sir_functions.r")

# Set data path
dpath = "/storage/sheinson_research/"

# pf - function to run particle filter given n number of particles, n.sim-th data set, filt = "BF", "APF", or "KD", transform = "logit" or "log" (used with KD only)
# Systematic resampling, uniform priors on beta, gamma, and nu
# Returns nothing; saves .rdata data file
pf <- function(n.sim, n, filt, transform, seed, progress = TRUE)
{
  # Load simulated data
  load(paste(dpath,"sim-orig.rdata",sep=""))
  y = mysim$sim[[n.sim]]$y
  P = mysim$true.params$P
  b = mysim$true.params$b
  varsigma = mysim$true.params$varsigma
  sigma = mysim$true.params$sigma
  eta = mysim$true.params$eta

  # Define uniform bounds on unknown parameters
  thetal = c(0.1400, 0.0900, 0.9500)
  thetau = c(0.5000, 0.1430, 1.3000)

  # Set seed
  set.seed(seed)

  if(filt == "BF" | filt == "APF")
  {
    ftheta = function(x, param=1) x # for computing quantiles with pf.quantile() in pf_functions.r
    mydllik = function(y, x) dllik(y, x, b, varsigma, sigma, eta)
    myrevo = function(x) revo(x, x[3:5], P)
    myrprior = function()
    {
      myprior = rprior(function() runif(3,thetal,thetau))
      return(c(myprior$x,myprior$theta))
    }
    if(filt == "BF")
    {
      source("bf.r")
      out = bf(y, mydllik, myrevo, myrprior, n, progress=progress, method="systematic", nonuniformity = "ess", threshold = 0.8, log=F)
    } else {
      pstate = function(x) revo(x, x[3:5], P, FALSE)
      source("apf.r")
      out = apf(y, mydllik, pstate, myrevo, myrprior, n, progress=progress, method="systematic", nonuniformity = "ess", threshold = 0.8, 
      log=F)
    }
  } else if(filt == "KD"){
    mydllik = function(y, x, theta) dllik(y, x, b, varsigma, sigma, eta)
    if(transform == "logit"){
      ftheta = function(theta,param=1:3) theta2u(theta,thetal[param],thetau[param])
      rtheta = function() u2theta(runif(3,thetal,thetau),thetal,thetau)
    } else if(transform == "log") {
      ftheta = function(theta,param=1) exp(theta)
      rtheta = function() log(runif(3,thetal,thetau))
    } else { stop("transform must be either 'logit' or 'log'") }
    pstate = function(x, theta) revo(x, ftheta(theta), P, FALSE)
    myrevo = function(x, theta) revo(x, ftheta(theta), P)
    myrprior = function() rprior(rtheta)
    source("kd_pf.r")
    out = kd_pf(y, mydllik, pstate, myrevo, myrprior, n, progress=progress, method="systematic", nonuniformity = "ess", threshold = 0.8, 
    log=F)
  } else { stop("filt must be one of 'BF', 'APF', or 'KD'") }

  # Save output
  pf.out = list(out=out,ftheta=ftheta)
  save(pf.out, file=paste(dpath,"PF-uniform-systematic-",filt,"-",n,"-",transform,"-",seed,"-",n.sim,".rdata",sep=""))
}

# Create data frame and use plyr to run particle filters in parallel
mydata = expand.grid(n.sim = 1:20, n = c(100,1000,10000,20000,40000), filt=c("BF","APF","KD"), transform = "logit", seed = 61, progress=FALSE, stringsAsFactors=FALSE)
mydata = rbind(mydata, data.frame(n.sim = 1, n = 10000, filt = "KD", transform = "log", seed = 61, progress=FALSE, stringsAsFactors=FALSE))

require(plyr)
require(doMC)
registerDoMC()
m_ply(.data = mydata, .fun = pf, .parallel = TRUE)

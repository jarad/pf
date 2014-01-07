source("sir_functions.r")

# Set data path
dpath = "/storage/sheinson_research/"

# pf - function to run particle filter given n number of particles, n.sim-th data set, resamp = "multinomial", "residual", "stratified", or "systematic", prior = "orig" or "disp", and delta amount of jitter to particles
# KD particle filter, lognormal priors on beta, gamma, and nu
# Returns nothing; saves .rdata data file
pf <- function(n.sim, n, resamp, prior, delta, seed, progress = TRUE)
{
  # Load simulated data
  load(paste(dpath,"sim-orig.rdata",sep=""))
  y = mysim$sim[[n.sim]]$y
  P = mysim$true.params$P
  b = mysim$true.params$b
  varsigma = mysim$true.params$varsigma
  sigma = mysim$true.params$sigma
  eta = mysim$true.params$eta

  # Define prior mean and sd of unknown parameters on log scale
  if(prior == "orig")
  {
    log.params <- find.mu.sigma(c(.14, .09, .95), c(.5, .143, 1.3))
    rtheta <- function() rnorm(3, log.params[[1]], log.params[[2]])
    ftheta <- function(x, param=1) exp(x)
  } else if(prior == "disp") {
    log.params <- find.mu.sigma(c(.05, .01, .75), c(2, 0.25, 1.75))
    rtheta <- function() rnorm(3, log.params[[1]], log.params[[2]])
    ftheta <- function(x, param=1) exp(x)
  } else if(prior == "unit"){
    thetal = c(0, 0)
    thetau = c(1, 1)
    log.nu <- find.mu.sigma(.75, 1.75)    
    rtheta <- function() c(u2theta(runif(2,thetal,thetau),thetal,thetau), rnorm(1,log.nu[[1]],log.nu[[2]]))
    ftheta <- function(theta, param=1:3)
    {
      if(length(param) < 1) stop("param must have at least 1 element")
      u = rep(NA, length(param))
      for(d in 1:length(param))
      {
        if(param[d] < 3)
        {
          u[d] = theta2u(theta[param[d]],thetal[param[d]],thetau[param[d]])
        } else { u[d] = exp(theta[param[d]]) }
      }
      return(u)
    }
  } else { stop("prior must be 'orig', 'disp', or 'unit'") }

  # Set seed
  set.seed(seed)

  mydllik = function(y, x, theta) dllik(y, x, b, varsigma, sigma, eta)
  pstate = function(x, theta) revo(x, ftheta(theta), P, FALSE)
  myrevo = function(x, theta) revo(x, ftheta(theta), P)
  myrprior = function() rprior(rtheta)
  source("kd_pf.r")
  out = kd_pf(y, mydllik, pstate, myrevo, myrprior, n, delta, progress, method=resamp, nonuniformity = "ess", threshold = 0.8, log=F)

  # Save output
  pf.out = list(out=out,ftheta=ftheta)
  file = paste(dpath,"PF-KD-lognormal-",resamp,"-",n,"-",prior,"-",delta,"-",seed,"-",n.sim,".rdata",sep="")
  print(file)
  save(pf.out, file=file)
}

# Create data frame and use plyr to run particle filters in parallel
data1 = expand.grid(n.sim = 1:20, n = c(100, 1000, 10000, 20000, 40000), resamp = c("multinomial","residual","stratified","systematic"), prior="orig", delta = .99, seed = 61, progress=FALSE, stringsAsFactors=FALSE)
data2 = expand.grid(n.sim = 1:20, n = c(100, 1000, 10000, 20000, 40000), resamp = "stratified", prior="orig", delta = c(0.9,0.95,0.96,0.97,0.98,0.99), seed = 61, progress = FALSE, stringsAsFactors=FALSE)
mydata = rbind(data1, data2)

require(plyr)
require(doMC)
registerDoMC()
m_ply(.data = mydata, .fun = pf, .parallel = TRUE)

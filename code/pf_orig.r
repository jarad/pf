source("sir_functions.r")
source("sir_mcmc_functions.r")

# Set data path
dpath = "/storage/sheinson_research/"

# pf - function to run particle filter given n number of particles, n.sim-th data set, resamp = "multinomial", "residual", "stratified", or "systematic", prior = "orig" or "disp", and delta amount of jitter to particles
# KD particle filter, lognormal priors on beta, gamma, and nu
# Returns nothing; saves .rdata data file
pf <- function(n.sim, n, filt, resamp, prior, transform, delta, seed, progress = TRUE)
{
  # Load simulated data
  load(paste(dpath,"sim-orig.rdata",sep=""))
  y = mysims[[n.sim]]$sim$y
  P = mysims[[n.sim]]$true.params$P
  b = mysims[[n.sim]]$true.params$b
  varsigma = mysims[[n.sim]]$true.params$varsigma
  sigma = mysims[[n.sim]]$true.params$sigma
  eta = mysims[[n.sim]]$true.params$eta  
  
  # Define functions to transform theta to original scale
  if(transform == "logit"){
    ftheta = function(theta,param=1:3) theta2u(theta,c(0.1400, 0.0900, 0.9500)[param], c(0.5, 0.143, 1.3)[param])
  } else if(transform == "log") {
    ftheta = function(theta,param=1) exp(theta)
  } else if(transform == "none"){
    ftheta = function(theta, param=1) theta 
  } else { stop("transform must be one of 'logit', 'log', or 'none'")}
  
  # Define functions to sample prior draws of fixed parameters
  if(prior == "orig")
  {
    if(transform == "log")
    {
      rtheta <- function()
      {
        theta <- rep(NA, 3)
        log.params <- find.mu.sigma(c(.09, .95), c(.143, 1.3))
        theta[2:3] <- exp(rnorm(2, log.params[[1]], log.params[[2]]))
        theta[1] <- theta[2]*runif(1, 1.2, 3)
        return(log(theta))
      }
    } else if(transform == "none") { 
      rtheta <- function()
      {
        theta <- rep(NA, 3)
        log.params <- find.mu.sigma(c(.09, .95), c(.143, 1.3))
        theta[2:3] <- exp(rnorm(2, log.params[[1]], log.params[[2]]))
        theta[1] <- theta[2]*runif(1, 1.2, 3)
        return(theta)
      }
    } else { stop("Must use log or no transformation with original prior") }
  } else if(prior == "disp") {
    if(transform == "log") {
      rtheta <- function()
      {
        theta <- rep(NA, 3)
        log.params <- find.mu.sigma(c(.05, .75), c(.25, 1.5))
        theta[2:3] <- exp(rnorm(2, log.params[[1]], log.params[[2]]))
        theta[1] <- theta[2]*runif(1, 1.05, 3)
        return(log(theta))
      }
    } else if(transform == "none") {
      rtheta <- function()
      {
        theta <- rep(NA, 3)
        log.params <- find.mu.sigma(c(.05, .75), c(.25, 1.5))
        theta[2:3] <- exp(rnorm(2, log.params[[1]], log.params[[2]]))
        theta[1] <- theta[2]*runif(1, 1.05, 3)
        return(theta)
      }
    } else {stop("Must use log or no transformation with dispersed prior")}
  } else if(prior == "unif") {
    if(transform == "logit") {
      rtheta <- function()
      {
        theta = runif(3, c(0.1400, 0.0900, 0.9500), c(0.5000, 0.1430, 1.3000))
        return(u2theta(theta, c(0.1400, 0.0900, 0.9500), c(0.5000, 0.1430, 1.3000)))
      }
    } else if(transform == "none") {
      rtheta <- function()
      {
        theta = runif(3, c(0.1400, 0.0900, 0.9500), c(0.5000, 0.1430, 1.3000))
        return(theta)
      }
    } else {stop("Must use logit or no transformation with uniform prior")}
  } else { stop("prior must be 'orig', 'disp', or 'unif'") }
  
  # Set seed
  set.seed(seed)
  
  if(filt == "BF" | filt == "APF")
  {
    mydllik = function(y, x) dllik(y, x, b, varsigma, sigma, eta)
    myrevo = function(x) revo(x, ftheta(x[3:5]), P)
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
      pstate = function(x) revo(x, ftheta(x[3:5]), P, FALSE)
      source("apf.r")
      out = apf(y, mydllik, pstate, myrevo, myrprior, n, progress=progress, method="systematic", nonuniformity = "ess", threshold = 0.8, log=F)
    }
  } else if(filt == "KD"){
    mydllik = function(y, x, theta) dllik(y, x, b, varsigma, sigma, eta)
    myrevo = function(x, theta) revo(x, ftheta(theta), P)
    pstate = function(x, theta) revo(x, ftheta(theta), P, FALSE)
    myrprior = function() rprior(rtheta)
    source("kd_pf.r")
    out = kd_pf(y, mydllik, pstate, myrevo, myrprior, n, delta, progress, method=resamp, nonuniformity = "ess", threshold = 0.8, log=F)
  } else if(filt == "RM"){
    mydllik = function(y, x, theta) dllik(y, x, b, varsigma, sigma, eta)
    myrevo = function(x, theta) revo(x, ftheta(theta), P)
    myrprior = function() rprior(rtheta)
    myrmove = function(y, x, theta) rmove(y, x, ftheta(theta), list(b=b, varsigma=varsigma, sigma=sigma, eta=eta, P=P))
    source("rm_pf.r")
    out = rm_pf(y, mydllik, myrevo, myrprior, myrmove, n, progress=progress, method=resamp, nonuniformity = "ess", threshold = 0.8, log=F)
  } else { stop("filt must be one of 'BF', 'APF', or 'KD'") }
  
  # Save output
  pf.out = list(out=out,ftheta=ftheta)
  file = paste(dpath,"PF-",n.sim,"-",n,"-",filt,"-",resamp,"-",prior,"-",transform,"-",delta,"-",seed,".rdata",sep="")
  print(file)
  save(pf.out, file=file)
}

# Create data frame and use plyr to run particle filters in parallel
#data1 = expand.grid(n.sim = 1:20, n = c(100, 1000, 10000, 20000), filt = c("BF","APF","KD"), resamp = "systematic", prior = "unif", transform = "logit", delta = 0.99, seed = 61, progress = FALSE, stringsAsFactors=FALSE)
#data2 = expand.grid(n.sim = 1:20, n = c(100, 1000, 10000, 20000), filt = "KD", resamp = c("multinomial","residual","stratified","systematic"), prior="orig", transform="log", delta = .99, seed = 61, progress = FALSE, stringsAsFactors = FALSE)
#data3 = expand.grid(n.sim = 1:20, n = c(100, 1000, 10000, 20000), filt = "KD", resamp = "stratified", prior="orig", transform="log", delta = c(0.9,0.95,0.96,0.97,0.98), seed = 61, progress = FALSE, stringsAsFactors = FALSE)
#data4 = expand.grid(n.sim = 1:20, n = c(100, 1000, 10000, 20000), filt = "KD", resamp = "stratified", prior = "disp", transform = "log", delta = 0.99, seed = 61, progress = FALSE, stringsAsFactors=FALSE)
data5 = data.frame(n.sim = 1, n = 10000, filt = "KD", resamp = "systematic", prior = "unif", transform = "none", delta = 0.99, seed = 61, progress = FALSE, stringsAsFactors = FALSE)
data6 = data.frame(n.sim = 1, n = c(100, 1000, 10000, 20000), filt = "RM", resamp = "stratified", prior = "orig", transform = "none", delta = 0.99, seed = 61, progress = FALSE, stringsAsFactors = FALSE)
#mydata = rbind(data1, data2, data3, data4)
mydata = rbind(data5, data6)

require(plyr)
require(doMC)
registerDoMC()
m_ply(.data = mydata, .fun = pf, .parallel = TRUE)
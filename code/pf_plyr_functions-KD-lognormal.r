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
    log.params <- find.mu.sigma(c(.05, .05, .75), c(0.75, 0.25, 1.5))
    rtheta <- function()
    {
      log.theta = rnorm(3, log.params[[1]], log.params[[2]])
      while(!all(exp(log.theta) > 0) | !all(exp(log.theta[1:2]) < 1)) log.theta = rnorm(3, log.params[[1]], log.params[[2]])
      return(log.theta)
    }
    ftheta <- function(x, param=1) exp(x)
  } else { stop("prior must be 'orig' or 'disp'") }

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
data1 = expand.grid(n.sim = 1:20, n = c(100, 1000, 10000, 20000, 40000), resamp = c("multinomial","residual","stratified","systematic"), prior="orig", delta = .99, seed = 61, progress = FALSE, stringsAsFactors = FALSE)
data2 = expand.grid(n.sim = 1:20, n = c(100, 1000, 10000, 20000, 40000), resamp = "stratified", prior="orig", delta = c(0.9,0.95,0.96,0.97,0.98), seed = 61, progress = FALSE, stringsAsFactors = FALSE)
data3 = expand.grid(n.sim = 1:20, n = c(100, 1000, 10000, 20000, 40000), resamp = "stratified", prior = "disp", delta = 0.99, seed = 61, progress = FALSE, stringsAsFactors=FALSE)
data4 = expand.grid(n.sim = 1, n = c(60000, 80000), resamp = c("multinomial","residual","stratified","systematic"), prior = "orig", delta = 0.99, seed = 61, progress = FALSE, stringsAsFactors = FALSE)
data5 = expand.grid(n.sim = 1, n = c(60000, 80000), resamp = "stratified", prior = "orig", delta = c(0.9,0.95,0.96,0.97,0.98), seed = 61, progress = FALSE, stringsAsFactors = FALSE)
data6 = expand.grid(n.sim = 1, n = c(60000, 80000), resamp = "stratified", prior = "disp", delta = 0.99, seed = 61, progress = FALSE, stringsAsFactors = FALSE)
mydata = rbind(data1, data2, data3, data4, data5, data6)

require(plyr)
require(doMC)
registerDoMC()
m_ply(.data = mydata, .fun = pf, .parallel = TRUE)

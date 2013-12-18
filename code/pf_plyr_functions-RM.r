source("sir_functions.r")
source("pf_functions.r")
source("sir_mcmc_functions.r")

# Set data path
dpath = "/storage/sheinson_research/"

# pf_rm - function to run resample-move particle filter given n number of particles, n.sim-th data set, resamp = "stratified", prior = "orig", and delta = 0.99
# lognormal priors on beta, gamma, and nu
# Returns nothing; saves .rdata data file
pf <- function(n.sim, n, seed, progress = TRUE)
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
  log.params <- find.mu.sigma(c(.14, .09, .95), c(.5, .143, 1.3))
  rtheta <- function() exp(rnorm(3, log.params[[1]], log.params[[2]]))
  ftheta = function(x, param=1) x

  # Set seed
  set.seed(seed)

  mydllik = function(y, x, theta) dllik(y, x, b, varsigma, sigma, eta)
  myrevo = function(x, theta) revo(x, theta, P)
  myrprior = function() rprior(rtheta)
  myrmove = function(y, x, theta) rmove(y, x, theta, list(b=b, varsigma=varsigma, sigma=sigma, eta=eta, P=P))
  source("rm_pf.r")
  out = rm_pf(y, mydllik, myrevo, myrprior, myrmove, n, progress=progress, method="stratified", nonuniformity = "ess", threshold = 0.8, log=F)

  # Save output
  pf.out = list(out=out,ftheta=ftheta)
  file = paste(dpath,"PF-RM-lognormal-stratified-",n,"-orig-0.99-",seed,"-",n.sim,".rdata",sep="")
  print(file)
  save(pf.out, file=file)
}

pf.quant.rm = function(n.sim, n, seed)
{
  # Load data
  load(paste(dpath,"PF-RM-lognormal-stratified-",n,"-orig-0.99-",seed,"-",n.sim,".rdata",sep=""))
  out = pf.out$out
  ftheta = pf.out$ftheta

  # Calculate 2.5%, 50%, and 97.5% quantiles of states over time
  probs = c(.5,.025,.975)
  tt = dim(out$state)[3]
  states = array(NA,dim=c(3,n,tt))
  states[1:2,,] = out$state[1:2,,]
  for(i in 1:tt)
  {
    states[3,,i] = 1 - apply(out$state[1:2,,i],2,sum)
  }
  state.quant = pf.quantile(states, out$weight, function(x,param=1) x, probs)

  # Calculate 2.5%, 50%, and 97.5% quantiles of parameters over time
  if("theta" %in% names(out))
  {
    theta.quant = pf.quantile(out$theta, out$weight, ftheta, probs)
  } else if(dim(out$state)[1] > 2){
    theta.quant = pf.quantile(out$state[3:(dim(out$state)[1]),,], out$weight, ftheta, probs)
  } else { theta.quant = NULL}

  # Save data
  pf.quant.out = list(state.quant=state.quant,theta.quant=theta.quant,probs=probs)
  file = paste(dpath,"PF-quant-RM-lognormal-stratified-",n,"-orig-0.99-",seed,"-",n.sim,".rdata",sep="")
  print(file)
  save(pf.quant.out, file=file)
}

# Create data frame and use plyr to run particle filters in parallel
mydata = expand.grid(n.sim = 1:20, n = c(100, 1000, 10000, 20000, 40000, 60000, 80000), seed=61, progress=FALSE, stringsAsFactors=FALSE)

require(plyr)
require(doMC)
registerDoMC()
m_ply(.data = mydata, .fun = pf_rm, .parallel = TRUE)

# Calculate quantiles
mydata.quant = mydata[,-dim(mydata)[2]]
m_ply(.data = mydata.quant, .fun = pf.quant.rm, .parallel = TRUE)

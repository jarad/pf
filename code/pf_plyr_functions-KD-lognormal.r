source("sir_functions.r")
source("pf_functions.r")

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
        if(d < 3)
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
  save(pf.out, file=paste(dpath,"PF-KD-lognormal-",resamp,"-",n,"-",prior,"-",delta,"-",seed,"-",n.sim,".rdata",sep=""))
}

pf.quant = function(n.sim, n, resamp, prior, delta, seed)
{
  # Load data
  load(paste(dpath,"PF-KD-lognormal-",resamp,"-",n,"-",prior,"-",delta,"-",seed,"-",n.sim,".rdata",sep=""))
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
  save(pf.quant.out, file=paste(dpath,"PF-quant-KD-lognormal-",resamp,"-",n,"-",prior,"-",delta,"-",seed,"-",n.sim,".rdata",sep=""))
}

# Create data frame and use plyr to run particle filters in parallel
data1 = expand.grid(n.sim = 1:20, n = c(100,1000,10000,20000,40000), resamp = c("multinomial","residual","stratified","systematic"), prior="orig", delta = .99, seed = 61, progress=FALSE, stringsAsFactors=FALSE)
data2 = expand.grid(n.sim = 1:20, n = c(100,1000,10000,20000,40000), resamp = "stratified", prior="orig", delta = c(.9,.95,.96,.97,.98), seed = 61, progress=FALSE, stringsAsFactors=FALSE)
data3 = expand.grid(n.sim = 1:20, n = c(100,1000,10000,20000,40000), resamp = "stratified", prior = "disp", delta = 0.99, seed = 61, progress=FALSE, stringsAsFactors=FALSE)
data4 = expand.grid(n.sim = 1:20, n = c(100,1000,10000,20000,40000), resamp = "stratified", prior = "unit", delta = 0.99, seed = 61, progress=FALSE, stringsAsFactors=FALSE)
mydata = rbind(data1, data2, data3, data4)
require(plyr)
require(doMC)
registerDoMC()
m_ply(.data = mydata, .fun = pf, .parallel = TRUE)

# Calculate quantiles
mydata.quant = mydata[,-dim(mydata)[2]]
m_ply(.data = mydata.quant, .fun = pf.quant, .parallel = TRUE)

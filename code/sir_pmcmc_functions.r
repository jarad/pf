source("sir_functions.r")
source("sir_mcmc_functions.r")

## Functions for creating a 'pomp' object for sir model ##
## See 'pomp' method in pomp package for additional information on these functions

# rprocess generates a particle trajectory given an initial state 'xstart', vector of observation times 'times', and fixed parameters 'params'
rprocess <- function(xstart, times, params, P, log.par = FALSE, ...)
{
  ns <- dim(xstart)[1]
  np <- dim(xstart)[2]
  tt <- length(times)
  if(log.par) params = exp(params)
  
  x <- array(NA, dim = c(ns, np, tt))
  x[,,1] = xstart
  for(i in 2:tt)
  {
    for(j in 1:np)
    {
      x[,j,i] = revo(x[,j,i-1], params, P)
    }
  }
  rownames(x) = c("s","i")
  return(x)
}

# rmeasure generates an observation given the current state x, time t, and fixed parameters params
rmeasure <- function(x, t, params, b, varsigma, sigma, eta, ...)
{
  y <- robs(x,b,varsigma,sigma,eta)
  names(y) = 1:length(y)
  return(y)
}

# dmeasure evaluates the probability density (or log of pdf if log = TRUE) of an observation y given the current state x, time t, and fixed parameters params
dmeasure <- function(y, x, t, params, log = TRUE, b, varsigma, sigma, eta, ...)
{
  log.pdf = dllik(y, x, b, varsigma, sigma, eta)
  if(log) return(log.pdf) else return(exp(log.pdf))
}

# skeleton provides the deterministic one-step map from state x at time t given fixed parameters params
skeleton <- function(x, t, params, P, log.par = FALSE, ...)
{
  if(log.par) params = exp(params)
  x.map = revo(x, params, P, random = FALSE)
  names(x.map) = c("s","i")
  return(x.map)
}

# initializer generates initial value of the state given fixed parameters params and initial time t0
initializer <- function(params, t0, ...)
{
  i0 = rnorm(1,0.002,0.0005)
  while(i0 < 0 | i0 > 1) i0 = rnorm(1,0.002,0.0005)
  s0 = 1 - i0
  x.init = c(s0, i0)
  names(x.init) = c("s","i")
  return(x.init)
}

# rprior.pomp generates fixed parameter values from their joint prior distribution
rprior.pomp <- function(params, log=FALSE, ...)
{
  theta <- rep(NA, 3)
  log.params <- find.mu.sigma(c(1.5, .09, .95), c(3, .143, 1.3))
  theta[2:3] <- exp(rnorm(2, log.params[[1]][2:3], log.params[[2]][2:3]))
  theta[1] <- theta[2]*exp(rnorm(1, log.params[[1]][1], log.params[[2]][1]))
  names(theta) = c("beta","gamma","nu")
  if(!log) return(theta) else return(log(theta))
}

# dprior evaluatues the prior density (or log density if log=TRUE) of the fixed parameters params
dprior <- function(params, log=FALSE, log.par = FALSE, ...)
{
  if(log.par) params = exp(params)
  log.prior <- dlbeta0gamma0(params[1], params[2]) + dlnu0(params[3])
  if(log) return(log.prior) else return(exp(log.prior))
}

# par.trans and par.inv.trans transform the fixed parameters to the log scale and inverse-log scale, respectively
par.trans = function(params, ...)
{
  log.params = log(params)
  names(log.params) = c('beta','gamma','nu')
  return(log.params)
}

par.inv.trans = function(params, ...)
{
  exp.params = exp(params)
  names(exp.params) = c('beta','gamma','nu')
  return(exp.params)
}
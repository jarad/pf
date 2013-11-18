source("sir_functions.r")

rmove <- function(y, x, theta, psi)
{
  nt = dim(y)[2]
  stopifnot(nt == dim(x)[2] - 1)

  # Move states
  states = sample.x(y, x, theta, psi)
  
  # Move parameters
  params = sample.theta(y, x, theta, psi)
  
  return(list(state=states, theta=params[1:2]))
}

sample.x <- function(y, x, theta, psi)
{
  nt = dim(y)[2]
  stopifnot(nt == dim(x)[2] - 1)  
  x.prop = rx(x, theta, psi)
  logMH = dljoint(y, x.prop, theta, psi) + dlx(x, x.prop, theta, psi) - dljoint(y, x, theta, psi) - dlx(x.prop, x, theta, psi)
  if(log(runif(1)) < logMH) x = x.prop
  return(x)
}

sample.theta <- function(y, x, theta, psi, tuning = c(0.003, .001))
{
  # Sample beta
  nt = dim(y)[2]
  stopifnot(nt == dim(x)[2] - 1)
  beta.prop = rnorm(1, theta[1], tuning[1])
  while(!(beta.prop > 0)) beta.prop = rnorm(1, theta[1], tuning[1])
  logMH = dljoint(y, x, c(beta.prop, theta[2:3]), psi) + dnorm(theta[1], beta.prop, tuning[1], log=TRUE) - dljoint(y, x, theta, psi) -  dnorm(beta.prop, theta[1], tuning[1])
  if(log(runif(1)) < logMH) theta[1] = beta.prop
  
  # Sample gamma
  nt = dim(y)[2]
  stopifnot(nt == dim(x)[2] - 1)
  gamma.prop = rnorm(1, theta[2], tuning[2])
  while(!(gamma.prop > 0)) gamma.prop = rnorm(1, theta[2], tuning[2])
  logMH = dljoint(y, x, c(theta[1], gamma.prop, theta[3]), psi) + dnorm(theta[2], gamma.prop, tuning[2], log=TRUE) - dljoint(y, x, theta, psi) -  dnorm(gamma.prop, theta[2], tuning[2])
  if(log(runif(1)) < logMH) theta[2] = gamma.prop
  
  return(theta)
}

#### Utility functions ####

rx <- function(x, theta, psi, sdfac = 1)
{
  s <- rnorm(dim(x)[2], x[1,], sdfac*sqrt(theta[1] / psi$P^2))
  i <- rnorm(dim(x)[2], x[2,] + -1*(s - x[1,]), sdfac*sqrt(theta[2] / psi$P^2))
  while(!all(apply(cbind(s,i), 1, in.Omega)))
  {
    s <- rnorm(dim(x)[2], x[1,], sdfac*sqrt(theta[1] / psi$P^2))
    i <- rnorm(dim(x)[2], x[2,] + -1*(s - x[1,]), sdfac*sqrt(theta[2] / psi$P^2))
  }
  return(rbind(s,i))
}

dlx <- function(x, x.curr, theta, psi, sdfac=1) dnorm(x[1], x.curr[1], sdfac*sqrt(theta[1] / psi$P^2), log=TRUE) + dnorm(x[2], x.curr[2] + -1*(x[1] - x.curr[1]), sdfac*sqrt(theta[2] / psi$P^2), log=TRUE)

dljoint <- function(y, x, theta, psi)
{
  nt = dim(y)[2]
  stopifnot(nt == dim(x)[2] - 1)
  log.density = dlprior(x[,1], theta)
  for(i in 1:nt) log.density = log.density + dllik(y[,i], x[,i+1], psi$b, psi$varsigma, psi$sigma, psi$eta) + dlevo(x[,i+1], x[,i], theta, 
  psi$P)
  return(log.density)
}

# dlevo - function to evaluate the logarithm of the state transition density
# Arguments:
# Arguments:
# x - 2-element vector of state (s,i) at which to evaluate the density
# x.curr - 2-element vector, current state (s,i)
# theta - 3-element vector, current parameter (beta,gamma,nu)
# P - scalar, population size
dlevo <- function(x,x.curr,theta,P,sdfac=1)
{
  if(all(x >= 0) & sum(x) <= 1)
  {
    dnorm(x[1], x.curr[1] - theta[1]*x.curr[2]*x.curr[1]^theta[3], sdfac*sqrt(theta[1]/P^2), log=TRUE) + dnorm(x[2], x.curr[2]*(1 - theta[2]) + x.curr[1] - x[1], sdfac*sqrt(theta[2] / P^2), log=TRUE)
  } else { return(-Inf) }
}

# dlprior - function that evaluates the logarithm of the prior density of the initial states and unknown parameters
# Arguments
# x - vector, state at which to evaluate density (s,i)
# theta - vector, unknown parameter values at which to evaluate the density (beta, gamma, nu)
dlprior <- function(x, theta)
{
  if(all(x >= 0) & sum(x) <= 1 & all(theta > 0))
  {
    dnorm(x[2], .002, .00025, log=TRUE) + sum(dlnorm(theta[1:2], c(-1.3296, -2.1764, 0.1055)[1:2], c(.3248, .1183, .0800)[1:2], log=TRUE))
  } else { return(-Inf) }
}

in.Omega <- function(x) all(x >= 0) & sum(x) <= 1
